"""Python counterpart of write_v2_lut.pro: write the ORAC V2 NetCDF LUT.

IDL:  pro write_v2_lut, V2_LUT_Filename, lutstr, inststr, srfstrarr, Vavg, BextOut, BextRatOut,
         SSAOUT, GOUT, TD, TfD, RD, RfD, RBD=RBD, RfBD=RfBD, TfBD=TfBd, TB=TB, EM=EM,
         include_pressure=include_pressure

The dimensions, variables and attributes are defined in the same order as the
IDL ncdf_dimdef / ncdf_vardef / ncdf_attput sequence and then written with the
tested NetCDF primitive oraclut.io.v2.write_v2_lut.  IDL arrays are
column-major, so an IDL variable defined on [chn_dim, efr_dim, opd_dim, saz_dim]
appears in the file with dimensions (satellite_zenith, optical_depth,
effective_radius, channels); the Python arrays are indexed [l, r, a, s] as
the IDL and are transposed (.T reverses all axes) on output.

Numerical conventions retained from the validated Python:
   * the IDL writes TD/100 etc.; here the operators are multiplied by
     float32(0.01), which is what the validated products contain;
   * E_md is written unscaled because the calling routine accumulates UU/BBE
     (the IDL accumulates 100*UU/BBE and divides by 100 here).
"""

from pathlib import Path

import numpy as np

from ..io.v2 import write_v2_lut as ncdf_write_v2_lut


_FILL = np.float32(9.96921e36)          # NetCDF default float fill, as left by the IDL partial writes


def _chars(text, length):
   """A fixed-length character array (IDL /char variable)."""

   encoded = text.encode("ascii", errors="replace")[:length]
   values = np.full(length, b" ", dtype="S1")
   values[:len(encoded)] = np.frombuffer(encoded, dtype="S1")
   return values


def _channel_optical_property(values, inststr):
   """Shape an optical-property vector like the legacy dual-view writer.

   The IDL forward-view replication copies the radiative-transfer operators,
   but it does not populate the forward-view entries of Bext, BextRat, SSA, or
   G.  NetCDF therefore retains its default float fill value in those slots.
   """

   array = np.asarray(values, dtype=np.float32)
   nchan = int(inststr.number_of_channels)
   nadir = int(getattr(inststr, "number_of_nadir_channels", nchan))
   if nchan <= nadir or array.shape[0] == nchan:
      return array
   if array.shape[0] != nadir:
      raise ValueError(
         "optical-property channel dimension does not match the instrument: "
         f"{array.shape[0]} values for {nadir} nadir/{nchan} total channels"
      )
   padded = np.full((nchan, *array.shape[1:]), _FILL, dtype=np.float32)
   padded[:nadir] = array
   return padded


def _solar_uncertainty_variables(inststr, solar_index):
   """Return the legacy solar uncertainty variable definitions."""

   if int(getattr(inststr, "view", 0)) > 0:
      return (
         ("rua", np.asarray(inststr.rua[solar_index], dtype=np.float32),
          "radiance uncertainty coefficient a", "dimensionless"),
         ("rub", np.asarray(inststr.rub[solar_index], dtype=np.float32),
          "radiance uncertainty coefficient b", "[W/(m^2 sr um)]^(1/2)"),
         ("ruc", np.asarray(inststr.ruc[solar_index], dtype=np.float32),
          "radiance uncertainty coefficient c", "[W/(m^2 sr um)]"),
      )
   return (("snr", np.asarray(inststr.snr[solar_index], dtype=np.float32),
            "signal-to-noise ratio", "dimensionless"),)


def write_v2_lut(v2_lut_filename, lutstr, inststr, srfstrarr, vavg, bextout, bextratout, ssaout, gout,
                 td, tfd, rd, rfd, rbd=None, rfbd=None, tfbd=None, tb=None, em=None, include_pressure=False,
                 global_attributes=None, em_valid_range=None):
   """Write the V2 LUT file ``v2_lut_filename``.

   ``global_attributes`` (none for the legacy product) are written as NetCDF
   global attributes; ``em_valid_range`` replaces the legacy 0-1 valid range of
   E_md (V25 adiabatic clouds, whose E_md is relative to the cloud-top Planck
   radiance and may exceed 1).
   """

   numberofsolarchannels = int(np.sum(inststr.solar_channel_flag))
   solar_channels_exist = numberofsolarchannels > 0
   solar_index = np.flatnonzero(inststr.solar_channel_flag)

   numberofthermalchannels = int(np.sum(inststr.thermal_channel_flag))
   thermal_channels_exist = numberofthermalchannels > 0
   thermal_index = np.flatnonzero(inststr.thermal_channel_flag)

   numberofmixedchannels = int(np.sum(inststr.mixed_channel_flag))
   mixed_channels_exist = numberofmixedchannels > 0
   mixed_index = np.flatnonzero(inststr.mixed_channel_flag)

   xmax = np.float32(np.finfo(np.float32).max)          # IDL (machar()).xmax
   nchan = inststr.number_of_channels

   # Define the dimensions to be used
   dimensions = {
      "optical_depth": lutstr.opd_n,
      "effective_radius": lutstr.efr_n,
      "satellite_zenith": lutstr.saz_n,
      "solar_zenith": lutstr.soz_n,
      "relative_azimuth": lutstr.raa_n,
   }
   if include_pressure:
      dimensions["surface_pressure"] = lutstr.prs_n
   dimensions["channels"] = nchan
   dimensions["length"] = max(len(name) for name in inststr.srf_file)
   dimensions["st1"] = len(inststr.instrument_filename)
   dimensions["st2"] = len(inststr.platform)
   dimensions["st3"] = len(inststr.instrument)
   dimensions["st4"] = len(inststr.instrument_version)
   if solar_channels_exist:
      dimensions["solar_channels"] = numberofsolarchannels
   if thermal_channels_exist:
      dimensions["thermal_channels"] = numberofthermalchannels
   if mixed_channels_exist:
      dimensions["mixed_channels"] = numberofmixedchannels

   variables = {}
   dims = {}
   attrs = {}

   def define(name, value, dimnames, **attributes):
      variables[name] = value
      dims[name] = tuple(dimnames)
      attrs[name] = attributes

   channel_range = np.asarray([0, nchan], dtype=np.int32)

   # instrument properties
   define("instrument_filename", _chars(inststr.instrument_filename, dimensions["st1"]), ("st1",))
   define("platform", _chars(inststr.platform, dimensions["st2"]), ("st2",))
   define("instrument", _chars(inststr.instrument, dimensions["st3"]), ("st3",))
   define("instrument_version", _chars(inststr.instrument_version, dimensions["st4"]), ("st4",))
   define("max_sat_zenith", np.float32(inststr.max_sat_zenith), (),
          units="degrees", valid_range=np.asarray([0.0, 90.0], dtype=np.float32))
   define("number_of_channels", np.int32(nchan), (),
          units="dimensionless", valid_range=channel_range)

   # channel specific instrument properties
   define("SRF_file", np.stack([_chars(name, dimensions["length"]) for name in inststr.srf_file]),
          ("channels", "length"), long_name="file containing the spectral response for the channel")
   define("channel_id", np.asarray(inststr.channelid, dtype=np.int16), ("channels",),
          long_name="Instrument channel identifier", units="dimensionless", valid_range=channel_range)
   define("central_wavelength", np.asarray([s.wvl_centre for s in srfstrarr], dtype=np.float32), ("channels",),
          long_name="effective central wavelength for the channel", units="microns",
          valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("central_wavenumber", np.asarray([s.wvn_centre for s in srfstrarr], dtype=np.float32), ("channels",),
          long_name="effective central wavenumber for the channel", units="cm^{-1}",
          valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("solar_channel_flag", np.asarray(inststr.solar_channel_flag, dtype=np.int16), ("channels",),
          long_name="Flag set to 1 if channel measures reflected solar radiation otherwise 0",
          units="dimensionless", valid_range=np.asarray([0, 1], dtype=np.int16))
   define("mixed_channel_flag", np.asarray(inststr.mixed_channel_flag, dtype=np.int16), ("channels",),
          long_name="Flag set to 1 if channel measures reflected solar and emitted infrared radiation otherwise 0",
          units="dimensionless", valid_range=np.asarray([0, 1], dtype=np.int16))
   define("thermal_channel_flag", np.asarray(inststr.thermal_channel_flag, dtype=np.int16), ("channels",),
          long_name="Flag set to 1 if channel measures emitted infrared radiation otherwise 0",
          units="dimensionless", valid_range=np.asarray([0, 1], dtype=np.int16))

   if solar_channels_exist:
      define("solar_channel_id", np.asarray(inststr.channelid[solar_index], dtype=np.int16), ("solar_channels",),
             long_name="Instrument channel identifier for solar channels", units="dimensionless",
             valid_range=channel_range)
      define("oldf0", np.asarray(inststr.oldf0[solar_index], dtype=np.float32), ("solar_channels",),
             long_name="SAD file F0", units="W/(m^2 um)", valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      define("oldf1", np.asarray(inststr.oldf1[solar_index], dtype=np.float32), ("solar_channels",),
             long_name="SAD file F1", units="W/(m^2 um)", valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      define("F0", np.asarray([srfstrarr[i].f0 for i in solar_index], dtype=np.float32), ("solar_channels",),
             long_name="in-band solar radiance", units="W/(m^2 um)",
             valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      for name, values, long_name, units in _solar_uncertainty_variables(inststr, solar_index):
         attributes = {"long_name": long_name, "units": units}
         if name == "snr":
            attributes["valid_range"] = np.asarray([0.0, xmax], dtype=np.float32)
         define(name, values, ("solar_channels",), **attributes)
      # ***** TEMPORARY FOR BACK COMPATIBILITY ******
      # IDL defines oldnefr on [chn_dim] but writes only the solar channels; the
      # remaining elements hold the NetCDF fill value.
      oldnefr = np.full(nchan, _FILL, dtype=np.float32)
      oldnefr[:numberofsolarchannels] = inststr.oldnefr[solar_index]
      define("oldnefr", oldnefr, ("channels",),
             long_name="SAD file NeFr (measurement uncertainty)", units="W/(m^2 sr um)",
             valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      # *********************************************

   if mixed_channels_exist:
      define("mixed_channel_id", np.asarray(inststr.channelid[mixed_index], dtype=np.int16), ("mixed_channels",),
             long_name="Instrument channel identifier for mixed channels", units="dimensionless",
             valid_range=channel_range)

   if thermal_channels_exist:
      define("thermal_channel_id", np.asarray(inststr.channelid[thermal_index], dtype=np.int16), ("thermal_channels",),
             long_name="Instrument channel identifier for thermal channels", units="dimensionless",
             valid_range=channel_range)
      define("refbt", np.asarray(inststr.refbt[thermal_index], dtype=np.float32), ("thermal_channels",),
             long_name="reference brightness temperate at neDT has been calculated", units="K",
             valid_range=np.asarray([100, 500], dtype=np.int16))
      define("nedt", np.asarray(inststr.nedt[thermal_index], dtype=np.float32), ("thermal_channels",),
             long_name="noise equivalent delta temperature", units="K",
             valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      define("B1", np.asarray([srfstrarr[i].b1 for i in thermal_index], dtype=np.float32), ("thermal_channels",))
      define("B2", np.asarray([srfstrarr[i].b2 for i in thermal_index], dtype=np.float32), ("thermal_channels",))
      define("T1", np.asarray([srfstrarr[i].t1 for i in thermal_index], dtype=np.float32), ("thermal_channels",))
      define("T2", np.asarray([srfstrarr[i].t2 for i in thermal_index], dtype=np.float32), ("thermal_channels",))
      # ***** TEMPORARY FOR BACK COMPATIBILITY ****** (defined on [chn_dim], thermal channels written)

      def thermal_vector(values):
         vector = np.full(nchan, _FILL, dtype=np.float32)
         vector[:numberofthermalchannels] = values[thermal_index]
         return vector

      define("oldwvn", thermal_vector(inststr.oldwvn), ("channels",),
             long_name="SAD file Wvn", units="1/cm", valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      define("oldb1", thermal_vector(inststr.oldb1), ("channels",), long_name="SAD file B1")
      define("oldb2", thermal_vector(inststr.oldb2), ("channels",), long_name="SAD file B2")
      define("oldt1", thermal_vector(inststr.oldt1), ("channels",), long_name="SAD file T1")
      define("oldt2", thermal_vector(inststr.oldt2), ("channels",), long_name="SAD file T2")
      define("oldnebt", thermal_vector(inststr.oldnebt), ("channels",),
             long_name="SAD file NeBT", units="K", valid_range=np.asarray([0.0, xmax], dtype=np.float32))
      # *********************************************

   # MICROPHYSICAL properties (IDL [chn_dim, efr_dim] -> file (effective_radius, channels))
   define("average_volume_per_particle", np.asarray(vavg, dtype=np.float32), ("effective_radius",),
          long_name="average volume per particle", units="to be investigated",
          valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("extinction_coefficient", _channel_optical_property(bextout, inststr).T, ("effective_radius", "channels"),
          long_name="volume extinction coefficient", units="to be investigated",
          valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("extinction_coefficient_ratio", _channel_optical_property(bextratout, inststr).T, ("effective_radius", "channels"),
          long_name="ratio of volume extinction coefficient to the volume extinction coefficient at 550 nm",
          units="dimensionless", valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("single_scatter_albedo", _channel_optical_property(ssaout, inststr).T, ("effective_radius", "channels"),
          long_name="single scatter albedo", units="dimensionless",
          valid_range=np.asarray([0.0, 1.0], dtype=np.float32))
   define("asymmetry_parameter", _channel_optical_property(gout, inststr).T, ("effective_radius", "channels"),
          long_name="asymmetry parameter", units="dimensionless",
          valid_range=np.asarray([-1.0, 1.0], dtype=np.float32))

   # LUT axes
   define("optical_depth", np.asarray(lutstr.opd, dtype=np.float32), ("optical_depth",),
          long_name="optical depth", spacing=lutstr.opd_spacing, units="dimensionless",
          valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("effective_radius", np.asarray(lutstr.efr, dtype=np.float32), ("effective_radius",),
          long_name="particle effective radius", spacing=lutstr.efr_spacing, units="microns",
          valid_range=np.asarray([0.0, xmax], dtype=np.float32))
   define("satellite_zenith", np.asarray(lutstr.saz, dtype=np.float32), ("satellite_zenith",),
          long_name="satellite zenith angle", spacing=lutstr.saz_spacing, units="degrees",
          valid_range=np.asarray([0.0, 180.0], dtype=np.float32))
   define("solar_zenith", np.asarray(lutstr.soz, dtype=np.float32), ("solar_zenith",),
          long_name="solar zenith angle", spacing=lutstr.soz_spacing, units="degrees",
          valid_range=np.asarray([0.0, 180.0], dtype=np.float32))
   define("relative_azimuth", np.asarray(lutstr.raa, dtype=np.float32), ("relative_azimuth",),
          long_name="satellite azimuth relative to the Sun", spacing=lutstr.raa_spacing, units="degrees",
          valid_range=np.asarray([0.0, 180.0], dtype=np.float32))
   if include_pressure:
      define("surface_pressure", np.asarray(lutstr.prs, dtype=np.float32), ("surface_pressure",),
             long_name="surface pressure", spacing=lutstr.prs_spacing, units="hPa",
             valid_range=np.asarray([900.0, 1100.0], dtype=np.float32))

   # LUT operators.  IDL: [chn_dim, (prs_dim,) efr_dim, opd_dim, saz_dim] etc.; the
   # Python arrays are [l, (k,) r, a, s] and .T gives the file order.
   prs = ("surface_pressure",) if include_pressure else ()
   unit_range = np.asarray([0.0, 1.0], dtype=np.float32)
   scale = np.float32(0.01)          # IDL: /100 (see module docstring)
   define("T_dv", scale * np.asarray(td, dtype=np.float32).T,
          ("satellite_zenith", "optical_depth", "effective_radius", *prs, "channels"),
          long_name="diffuse transmission of direct light", units="dimensionless", valid_range=unit_range)
   define("T_dd", scale * np.asarray(tfd, dtype=np.float32).T,
          ("optical_depth", "effective_radius", *prs, "channels"),
          long_name="diffuse transmission", units="dimensionless", valid_range=unit_range)
   define("R_dv", scale * np.asarray(rd, dtype=np.float32).T,
          ("satellite_zenith", "optical_depth", "effective_radius", *prs, "channels"),
          long_name="direct reflection of diffuse light", units="dimensionless", valid_range=unit_range)
   define("R_dd", scale * np.asarray(rfd, dtype=np.float32).T,
          ("optical_depth", "effective_radius", *prs, "channels"),
          long_name="diffuse reflection of diffuse light", units="dimensionless", valid_range=unit_range)
   if solar_channels_exist:
      define("R_0v", scale * np.asarray(rbd, dtype=np.float32)[solar_index].T,
             ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth", "effective_radius", *prs, "solar_channels"),
             long_name="bi-directional reflectance", units="dimensionless", valid_range=unit_range)
      define("R_0d", scale * np.asarray(rfbd, dtype=np.float32)[solar_index].T,
             ("solar_zenith", "optical_depth", "effective_radius", *prs, "solar_channels"),
             long_name="diffuse reflectance of direct beam", units="dimensionless", valid_range=unit_range)
      define("T_0d", scale * np.asarray(tfbd, dtype=np.float32)[solar_index].T,
             ("solar_zenith", "optical_depth", "effective_radius", *prs, "solar_channels"),
             long_name="diffuse transmission of diffuse light", units="dimensionless", valid_range=unit_range)
      define("T_00", scale * np.asarray(tb, dtype=np.float32)[solar_index].T,
             ("solar_zenith", "optical_depth", "effective_radius", *prs, "solar_channels"),
             long_name="direct transmission", units="dimensionless", valid_range=unit_range)
   if thermal_channels_exist:
      # em already holds the emissivity as a fraction (UU/BBE): no scale here
      define("E_md", np.asarray(em, dtype=np.float32)[thermal_index].T,
             ("satellite_zenith", "optical_depth", "effective_radius", *prs, "thermal_channels"),
             long_name="diffuse emissivity", units="dimensionless",
             valid_range=unit_range if em_valid_range is None else np.asarray(em_valid_range, dtype=np.float32))

   Path(v2_lut_filename).parent.mkdir(parents=True, exist_ok=True)
   ncdf_write_v2_lut(v2_lut_filename, lut_level=2, revision=0, dimensions=dimensions, variables=variables,
                     variable_dimensions=dims, variable_attributes=attrs, global_attributes=global_attributes)
