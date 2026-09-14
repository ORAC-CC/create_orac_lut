# Current V2 NetCDF LUT format

This document describes the selected current V2 product from the authorised
archive, based on `write_v2_lut.pro` and direct read-only inspection of
`meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc`. It does not describe older
SAD/V11/V12 formats.

## Container and dimensions

The file is NetCDF4 with no global attributes. The selected reference has:

| dimension | size |
|---|---:|
| `optical_depth` | 17 |
| `effective_radius` | 20 |
| `satellite_zenith` | 10 |
| `solar_zenith` | 10 |
| `relative_azimuth` | 11 |
| `channels` | 11 |
| `solar_channels` | 4 |
| `thermal_channels` | 8 |
| `mixed_channels` | 1 |
| `length` | 32 |
| `st1`/`st2`/`st3`/`st4` | 26/11/6/7 |

The current aerosol writer can insert `surface_pressure` into operator
dimensions, but it is absent from this cloud reference. NetCDF variable
dimension order is preserved exactly by the Python reader; it is not converted
to a channel-first convention.

## Variable groups

### Instrument and channel metadata

`instrument_filename`, `platform`, `instrument`, and `instrument_version` are
NetCDF character arrays. `max_sat_zenith` and `number_of_channels` are scalar
values. `SRF_file` is `[channels, length]`; character arrays should be decoded
only as a presentation operation, not silently reshaped.

Per-channel variables include `channel_id`, `central_wavelength` (microns),
`central_wavenumber` (cm^-1), and the solar/mixed/thermal flags. Solar-channel
metadata includes `solar_channel_id`, `oldf0`, `oldf1`, `F0`, and `snr`.
Thermal-channel metadata includes `thermal_channel_id`, `refbt`, `nedt`, `B1`,
`B2`, `T1`, and `T2`. Compatibility metadata includes `oldnefr`, `oldwvn`,
`oldb1`, `oldb2`, `oldt1`, `oldt2`, and `oldnebt`.

### Particle optical properties

The writer defines:

| variable | dimensions in reference | meaning from `write_v2_lut.pro` |
|---|---|---|
| `average_volume_per_particle` | `[effective_radius]` | average particle volume; units are left as “to be investigated” by the legacy writer |
| `extinction_coefficient` | `[effective_radius, channels]` | volume extinction coefficient |
| `extinction_coefficient_ratio` | `[effective_radius, channels]` | volume extinction relative to the 550 nm coefficient |
| `single_scatter_albedo` | `[effective_radius, channels]` | single-scattering albedo |
| `asymmetry_parameter` | `[effective_radius, channels]` | asymmetry parameter |

These quantities are produced before the RT loop by the legacy particle-optics
generator and are not recomputed by the reader.

### Coordinates

The reference coordinate variables are all `float32` except `channel_id`, which
is `int16`:

- `optical_depth`: dimensionless, uneven logarithmic; values begin at
  `1.000000013351432e-10`, then `0.0078125` through `256`;
- `effective_radius`: microns, linear, `1, 3, ..., 39`;
- `satellite_zenith`: degrees, linear, `0, 10, ..., 90`;
- `solar_zenith`: degrees, linear, `0, 10, ..., 90`;
- `relative_azimuth`: degrees, linear, `0, 18, ..., 180`;
- `channel_id`: integers `1` through `11`.

The complete coordinate byte checksums and representative values are in
`validation/reference_cases/meteosat10_seviri_liquid_water_stg_v21/fingerprint.json`.

### Radiative-transfer operators

The reference operator variables and exact dimensions are:

| variable | dimensions, in NetCDF order | writer meaning |
|---|---|---|
| `T_dv` | `[satellite_zenith, optical_depth, effective_radius, channels]` | diffuse transmission of direct light |
| `T_dd` | `[optical_depth, effective_radius, channels]` | diffuse transmission |
| `R_dv` | `[satellite_zenith, optical_depth, effective_radius, channels]` | direct reflection of diffuse light |
| `R_dd` | `[optical_depth, effective_radius, channels]` | diffuse reflection of diffuse light |
| `R_0v` | `[relative_azimuth, satellite_zenith, solar_zenith, optical_depth, effective_radius, solar_channels]` | bidirectional reflectance |
| `R_0d` | `[solar_zenith, optical_depth, effective_radius, solar_channels]` | diffuse reflectance of direct beam |
| `T_0d` | `[solar_zenith, optical_depth, effective_radius, solar_channels]` | diffuse transmission of diffuse light, direct-beam case |
| `T_00` | `[solar_zenith, optical_depth, effective_radius, solar_channels]` | direct transmission |
| `E_md` | `[satellite_zenith, optical_depth, effective_radius, thermal_channels]` | diffuse emissivity |

`write_v2_lut.pro` inserts `[surface_pressure]` between `effective_radius` and
`channels` for the pressure-aware aerosol path. The current Python reader exposes the arrays
exactly in the order stored here and retains each variable’s attributes and
dtype.

## Fill values, compression, and observed statistics

The inspected variables do not define a NetCDF `_FillValue`. Character and
numeric variables are stored contiguously and have no compression filters in
the reference file. The fingerprint records NaN/Inf counts and deterministic
statistics for each scientific optical-property and operator variable; all
observed NaN and Inf counts are zero.

Some operator values lie marginally below zero or above one despite the writer’s
declared valid ranges. These are retained as evidence of current output and
must not be silently clipped by an input reader.

## Reader contract

`src/oraclut/io/lut.py` provides `read_lut(path)`, and
`src/oraclut/io/v2.py` provides `write_v2_lut(...)`. The reader:

- requires an existing file and reports a clear `FileNotFoundError` otherwise;
- reads through the established `netCDF4` library;
- copies arrays into memory without transposition;
- exposes dimensions, coordinates, variable dimensions, dtypes, variable
  attributes, and global attributes;
- validates the required current V2 dimensions and scientific variables;
- supports exploratory non-validation reads with `validate_v21=False`.

It does not write files, alter values, decode character arrays into a different
shape, or add scientific defaults.

The writer accepts `lut_level=2` and an independent integer `revision` (for
example `revision=21`). It validates supplied dimensions and array shapes and
does not create scientific values that have not been supplied by an upstream
stage.

## Pressure-aware (aerosol formulation) products

When the LUT definition carries a sixth, surface-pressure grid the legacy
writer (`write_v2_lut.pro`, `/include_pressure`) adds:

- a `surface_pressure` dimension declared **between `relative_azimuth` and
  `channels`** — the full declaration order is `optical_depth`,
  `effective_radius`, `satellite_zenith`, `solar_zenith`, `relative_azimuth`,
  `surface_pressure`, `channels`, `length`, `st1`–`st4`, then the channel-class
  dimensions (`solar_channels`, `thermal_channels`, `mixed_channels`) that exist;
- a float32 coordinate `surface_pressure(surface_pressure)` with no fill value
  and the attributes `long_name = "surface pressure"`, `spacing` (the LUT
  header word, e.g. `uneven_linear`), `units = "hPa"`,
  `valid_range = [900., 1100.]` (float32);
- `surface_pressure` inserted as the second dimension (IDL order) of every
  radiative-transfer operator, i.e. immediately before `channels` /
  `solar_channels` / `thermal_channels` in the NetCDF (C) order used by the
  Python reader.

The Python writer reproduces this layout: `oraclut.io.v2.V2_DIMENSION_ORDER`
fixes the declaration order for every product (cloud products, which have no
pressure grid, are unaffected), and `V2_COORDINATE_DEFAULTS` supplies the
literal `long_name`/`units`/`valid_range` for `surface_pressure` when the
caller passes none. `spacing` is never inferred by the writer; the LUT reader
retains the pressure-block header word as `LutGrid.surface_pressure_spacing`
and the pipeline forwards it into the coordinate metadata, so the written
attribute is the LUT definition's own keyword (`uneven_linear` for
`aerosol_test.lut`).
Verified against the captured legacy aerosol references (2026-09-14).
