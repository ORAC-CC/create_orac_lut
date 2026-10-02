"""Routine-by-routine tests of the IDL-structured path against the validated Python path.

Each IDL-mirroring routine in src/oraclut/idl_mirror is compared with the
reader or kernel of the current validated implementation (src/oraclut/config,
pipeline, optics) on the validated inputs, requiring identical values.  The
liquid-water radius limit and the Legendre expansion of Mie classes are
deliberate departures (2026-10): generate_scattering_properties is compared
with the validated path with the IDL radius limits imposed and, for the
moments, where the fixed expansion is converged (tests/test_legendre_expansion.py
and tests/test_size_integration_limits.py cover the new behaviour).
"""

import importlib
from pathlib import Path

import numpy as np
import pytest

import create_orac_luts
from oraclut.config import read_atmosphere, read_instrument, read_lut_grid, read_microphysics, read_refractive_index
from oraclut.idl_mirror import (
   call_disort, create_bwgp, create_range, generate_scattering_properties, int_tabulated, interpol, legpexp, load_atmstr, load_gasstr,
   load_inststr, load_lutstr, load_mmdat, load_srfstrarr, mie_size_dist_new, quadrature, setup_disort,
   read_dubovik_filenames, read_dubovik_grid, read_dubovik_kernel_kext,
   read_dubovik_kernel_scatt_matrix, segment,
)
from oraclut.optics.legacy import _mie_distribution, legacy_stg_optics, legendre_moments
from oraclut.pipeline import _atmosphere_path, _gas_levels, _interpolate_profile, read_channels
from oraclut.radiative_transfer import call_disort as old_call_disort, getmom


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
INST = INPUTS / "inst" / "meteosat-10_seviri_v1.inst"
RUNS = ROOT / "runs"
CONCRETE_RUNS = sorted(RUNS.glob("*.run"))


def identical(a, b):
   a = np.asarray(a)
   b = np.asarray(b)
   return a.shape == b.shape and bool(np.array_equal(a, b, equal_nan=(a.dtype.kind == "f")))


# ---------------------------------------------------------------------------
# run file (makerunfile_v2.pro + driver file)
# ---------------------------------------------------------------------------

def test_template_and_validation_run_files_parse():
   for path in CONCRETE_RUNS:
      settings = create_orac_luts.read_runfile(path)
      assert settings["forward_model"] in ("cloud", "aerosol")
      assert settings["srf_quad"] == 1 and settings["nstreams"] == 60
      # nmom is the deprecated fixed expansion (legacy run files, or Baum / T-matrix)
      assert settings["nmom"] in (None, 1000)
   assert create_orac_luts.read_runfile(RUNS / "template.run")["nmom"] is None
   settings = create_orac_luts.read_runfile(RUNS / "meteosat-10_seviri_aerosol_a79_test_ch01_ch09.run")
   assert settings["channelid"] == [1, 9] and settings["gas"] == 1 and settings["lutfile"] == "aerosol_test.lut"


def test_run_file_missing_required_setting_is_an_error(tmp_path):
   text = (RUNS / "template.run").read_text().replace("nstreams      = 60", "")
   bad = tmp_path / "bad.run"
   bad.write_text(text)
   with pytest.raises(ValueError, match="required settings missing: nstreams"):
      create_orac_luts.read_runfile(bad)


@pytest.mark.parametrize("channel_line", ["channelid = []", "# channelid = [1]"])
def test_run_file_requires_explicit_nonempty_channels(tmp_path, channel_line):
   text = (RUNS / "template.run").read_text()
   text = text.replace("channelid     = [1]", channel_line)
   bad = tmp_path / "bad_channels.run"
   bad.write_text(text)
   with pytest.raises(ValueError, match="channelid must be specified explicitly"):
      create_orac_luts.read_runfile(bad)


@pytest.mark.parametrize("channels", ([1], [1, 9], [1, 4, 9], list(range(1, 12))))
def test_run_file_accepts_explicit_channel_selections(tmp_path, channels):
   text = (RUNS / "template.run").read_text().replace("channelid     = [1]", f"channelid = {channels}")
   path = tmp_path / "channels.run"
   path.write_text(text)
   assert create_orac_luts.read_runfile(path)["channelid"] == channels


def test_run_file_rejects_duplicate_channels(tmp_path):
   text = (RUNS / "template.run").read_text().replace("channelid     = [1]", "channelid = [1, 1]")
   path = tmp_path / "duplicate.run"
   path.write_text(text)
   with pytest.raises(ValueError, match="duplicate channel IDs"):
      create_orac_luts.read_runfile(path)


def test_run_file_reports_requested_and_valid_channels(tmp_path):
   text = (RUNS / "template.run").read_text().replace("channelid     = [1]", "channelid = [1, 99]")
   path = tmp_path / "invalid.run"
   path.write_text(text)
   with pytest.raises(ValueError, match=r"requested channels \[1, 99\].*valid channels are"):
      create_orac_luts.run(path)


def test_every_shipped_run_file_has_channels_defined_by_its_instrument():
   for path in CONCRETE_RUNS:
      settings = create_orac_luts.read_runfile(path)
      instrument = load_inststr(INPUTS / "inst" / settings["instfile"])
      assert set(settings["channelid"]).issubset(set(map(int, instrument.channelid))), path


def test_earthcare_a70_run_is_an_explicit_read_only_preflight_case():
   path = RUNS / "earthcare_msi_aerosol_a70.run"
   settings = create_orac_luts.read_runfile(path)
   assert settings["channelid"] == [1, 2, 3, 4, 5, 6, 7]
   assert settings["mmfile"] == "aerosol_a70.mm"
   assert settings["lutfile"] == "aerosol.lut"
   assert settings["gas"] == 1 and settings["no_rayleigh"] == 0
   assert settings["srf_quad"] == 1 and settings["nstreams"] == 60
   assert settings["tmatrix_path"] == "/network/aopp/matin/eodg/shared/dubovik_tmatrix"
   assert settings["nmom"] == 1000 and settings["version"] == 22
   assert (INPUTS / "inst" / settings["instfile"]).is_file()
   assert (INPUTS / "microphysics" / settings["mmfile"]).is_file()
   assert (INPUTS / "lut" / settings["lutfile"]).is_file()
   instrument = load_inststr(INPUTS / "inst" / settings["instfile"])
   assert settings["channelid"] == instrument.channelid.tolist()
   assert settings["out_path"] == "luts"
   assert "earthcare_msi_m_aerosol_a2_pa70_v22.nc" == (
      f"{settings['platform']}_{settings['instrument']}_m_aerosol_a2_pa70_v22.nc"
   )


def test_create_range_single_component_and_shared_mode():
   radii = np.asarray([0.01, 0.1, 1.0, 10.0])
   nrat, nmr = create_range([1.0], [0.07], [1.7], radii)
   np.testing.assert_allclose(nrat, np.ones((1, 4)))
   np.testing.assert_allclose(nmr[0], radii * 0.07 / (0.07 * np.exp(2.5 * np.log(1.7) ** 2)))

   nrat, nmr = create_range([0.8, 0.2], [0.07, 0.07], [1.7, 1.7], radii)
   np.testing.assert_allclose(nrat, [[0.8] * 4, [0.2] * 4])
   np.testing.assert_allclose(nmr[0], nmr[1])


def test_create_range_a70_matches_idl_gdl_reference_nodes():
   # Values printed by validation/diagnostics/create_range_a70.pro using GDL
   # and the authoritative create_orac_lut/create_range.pro.
   radii = np.asarray([0.010000, 0.06158483, 0.3792691, 2.335721, 10.00000])
   expected_nrat = np.asarray([
      [0.874999940395, 0.874999940395, 0.874100327492, 0.0, 0.0],
      [0.125000000000, 0.125000000000, 0.124871484935, 0.0, 0.0],
      [0.0, 0.0, 0.001028120518, 1.0, 1.0],
   ])
   expected_nmr = np.asarray([
      [0.004946444649, 0.030462596565, 0.070000000298, 0.070000000298, 0.070000000298],
      [0.004946444649, 0.030462596565, 0.070000000298, 0.070000000298, 0.070000000298],
      [0.787999987602, 0.787999987602, 0.787999987602, 0.949818968773, 4.066491603851],
   ])
   nrat, nmr = create_range([0.866250, 0.123750, 0.010000], [0.070000, 0.070000, 0.788000],
                            [1.700000, 1.700000, 1.822000], radii)
   np.testing.assert_allclose(nrat, expected_nrat, rtol=3e-6, atol=3e-7)
   np.testing.assert_allclose(nmr, expected_nmr, rtol=3e-6, atol=3e-7)


@pytest.mark.parametrize(
   "args",
   [([], [], [], [1.0]), ([1.0], [0.0], [1.7], [1.0]), ([1.0], [0.07], [0.0], [1.0]),
    ([1.0, 1.0, 1.0, 1.0, 1.0], [0.01, 0.02, 0.03, 0.04, 0.05], [1.2] * 5, [1.0])],
)
def test_create_range_rejects_invalid_mode_inputs(args):
   with pytest.raises(ValueError, match="create_range"):
      create_range(*args)


def test_a70_component_classification_and_aspect_ratio_data():
   mmstr = load_mmdat(INPUTS / "microphysics" / "aerosol_a70.mm", INPUTS)
   assert mmstr.substance == "aerosol" and mmstr.shortname == "a70"
   assert mmstr.compname == ["waf", "saf", "mdc"]
   assert [component.code for component in mmstr.comp] == ["mie", "mie", "tmatrix"]
   assert mmstr.distname == ["log_normal"] * 3
   assert not hasattr(mmstr.comp[0], "eps") and not hasattr(mmstr.comp[1], "eps")
   assert mmstr.comp[2].eps.shape == (25,)
   assert mmstr.comp[2].neps.shape == (25,)
   np.testing.assert_allclose(mmstr.comp[2].eps[0], 0.3349, rtol=1e-6)
   np.testing.assert_allclose(mmstr.comp[2].eps[-1], 2.9860, rtol=1e-6)
   np.testing.assert_allclose(mmstr.comp[2].neps.sum(), 1.0, rtol=3e-6)


def test_port_status_inventory_records_all_unavailable_a70_science():
   status = (ROOT / "docs" / "idl_python_port_status.md").read_text()
   for marker in (
      "Ported and validated", "Partially ported", "Not ported", "Intentionally unsupported",
      "External dependency", "dubovik_lognormal_multiple_eps", "ROUTINE/name.dat",
      "KERNEL_n22_181/", "tmatrix_path",
   ):
      assert marker in status


def test_run_file_unknown_setting_is_an_error(tmp_path):
   bad = tmp_path / "bad.run"
   bad.write_text((RUNS / "template.run").read_text() + "\nphase_order = 1000\n")
   with pytest.raises(ValueError, match="unknown setting 'phase_order'"):
      create_orac_luts.read_runfile(bad)


# ---------------------------------------------------------------------------
# load_inststr / load_lutstr / load_srfstrarr / load_mmdat / load_atmstr / load_gasstr
# ---------------------------------------------------------------------------

def test_load_inststr_matches_validated_instrument_reader():
   old = read_instrument(INST)
   inststr = load_inststr(INST)
   assert inststr.platform == "meteosat-10" and inststr.instrument == "seviri" and inststr.view == 0
   assert identical(inststr.channelid, old.available_channels)
   assert identical(inststr.solar_channel_flag, [int(c in old.solar_channels) for c in old.available_channels])
   assert identical(inststr.thermal_channel_flag, [int(c in old.thermal_channels) for c in old.available_channels])
   assert identical(inststr.mixed_channel_flag, [int(c in old.solar_channels and c in old.thermal_channels) for c in old.available_channels])
   assert inststr.srf_file == [old.srf_files[c] for c in old.available_channels]
   for tag in ("oldf0", "oldf1", "oldnefr", "oldwvn", "oldb1", "oldb2", "oldt1", "oldt2", "oldnebt", "snr", "refbt", "nedt"):
      expected = np.asarray([getattr(old, tag).get(c, 0.0) for c in old.available_channels], dtype=np.float32)
      assert identical(getattr(inststr, tag), expected), tag
   subset = load_inststr(INST, requestedchannelid=[1, 9])
   assert subset.channelid.tolist() == [1, 9] and subset.number_of_nadir_channels == 2
   assert subset.solar_channel_flag.tolist() == [1, 0] and subset.thermal_channel_flag.tolist() == [0, 1]
   with pytest.raises(ValueError, match="channel 99 not found"):
      load_inststr(INST, requestedchannelid=[99])


@pytest.mark.parametrize("lutfile", ["liquid-water-cloud_test.lut", "aerosol_test.lut", "liquid-water-cloud.lut", "aerosol.lut"])
def test_load_lutstr_matches_validated_grid_reader(lutfile):
   old = read_lut_grid(INPUTS / "lut" / lutfile)
   pressure = old.surface_pressure is not None
   lutstr = load_lutstr(INPUTS / "lut" / lutfile, np.float32(90.0), include_pressure=pressure)
   assert identical(lutstr.opd, old.optical_depth) and identical(lutstr.efr, old.effective_radius)
   assert identical(lutstr.soz, old.solar_zenith) and identical(lutstr.saz, old.satellite_zenith)
   assert identical(lutstr.raa, old.relative_azimuth)
   assert (lutstr.opd_spacing, lutstr.efr_spacing, lutstr.soz_spacing, lutstr.saz_spacing, lutstr.raa_spacing) == old.spacings
   if pressure:
      assert identical(lutstr.prs, old.surface_pressure) and lutstr.prs_spacing == old.surface_pressure_spacing


def test_load_lutstr_uses_the_instrument_maximum_satellite_zenith_like_idl(tmp_path):
   # ash-plume_test.lut lists 0 89; IDL lut_quadrature replaces the end point by max_sat_zenith.
   lutstr = load_lutstr(INPUTS / "lut" / "ash-plume_test.lut", np.float32(90.0))
   assert lutstr.saz.tolist() == [0.0, 90.0]
   # Relative azimuth must be linear (IDL stop), and every count must match its values.
   text = (INPUTS / "lut" / "liquid-water-cloud_test.lut").read_text()
   bad = tmp_path / "bad_raa.lut"
   bad.write_text(text.replace("2 linear # 11 Relative azimuth angles", "2 uneven_linear # azimuth"))
   with pytest.raises(ValueError, match="Relative azimuth spacing must be linear"):
      load_lutstr(bad, 90.0)
   bad.write_text(text.replace("3 linear # 20 Effective radius", "4 uneven_linear # radius"))
   with pytest.raises(ValueError, match="Effective radius N does not match"):
      load_lutstr(bad, 90.0)


@pytest.mark.parametrize("channels", [(1,), (9,), (1, 9)])
def test_load_srfstrarr_matches_validated_channel_reader(channels):
   old = read_channels(ROOT, read_instrument(INST), channels, srf_quad=1)
   inststr = load_inststr(INST, requestedchannelid=list(channels))
   srfstrarr, nwvl_max = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
   assert nwvl_max == 1
   for srf, channel in zip(srfstrarr, old):
      assert srf.wvl_centre == channel.wavelength_microns and srf.wvn_centre == channel.wavenumber_cm_inverse
      assert srf.f0 == channel.f0
      assert (srf.b1, srf.b2, srf.t1, srf.t2) == (channel.b1, channel.b2, channel.t1, channel.t2)
      assert srf.nwvl == 1 and srf.wvl[0] == srf.wvl_centre and srf.val[0] == 1.0


def test_load_srfstrarr_rejects_unported_qm0():
   inststr = load_inststr(INST, requestedchannelid=[1])
   with pytest.raises(NotImplementedError):
      load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 0, INPUTS)


def test_int_tabulated_matches_the_idl_five_point_newton_cotes_example():
   x = np.asarray([0.0, .12, .22, .32, .36, .40, .44, .54, .64, .70, .80], dtype=np.float32)
   y = np.asarray([.200000, 1.30973, 1.30524, 1.74339, 2.07490, 2.45600,
                   2.84299, 3.50730, 3.18194, 2.36302, .231964], dtype=np.float32)
   np.testing.assert_allclose(int_tabulated(x, y), np.float32(1.6232316), rtol=2e-6, atol=2e-7)


def test_segment_preserves_descending_srf_order_and_minimum_point_count():
   raw_reader = __import__("oraclut.idl_mirror.load_srfstrarr", fromlist=["read_srfstr"])
   source = raw_reader.read_srfstr(INPUTS / "srf" / "rtcoef_msg_3_seviri_srf_ch01.txt")
   x, y = segment(source.wvl, source.srf, np.float32(0.001), minn=12)
   assert x.size == 20
   assert x.size >= 12 and y.size == x.size
   assert np.all(np.diff(x) < 0.0)
   assert x[0] == source.wvl[0] and x[-1] == source.wvl[-1]
   assert np.all(y >= 0.0)


def test_qm2_seviri_channels_have_idl_segmented_point_counts_and_weights():
   inststr = load_inststr(INST, requestedchannelid=[1, 10])
   segmented, nwvl_max = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 2, INPUTS)
   assert nwvl_max == 23
   assert [srf.nwvl for srf in segmented] == [20, 23]
   for srf in segmented:
      assert srf.nwvl >= 12
      assert np.all(np.diff(srf.wvl[:srf.nwvl]) < 0.0)
      assert np.all(srf.val[:srf.nwvl] >= 0.0)
      assert np.isfinite(srf.f0)
      assert np.float32(0.0) < np.sum(srf.val[:srf.nwvl], dtype=np.float32)


def test_qm2_solar_and_thermal_channel_metadata_is_retained():
   inststr = load_inststr(INST, requestedchannelid=[1, 10])
   segmented, _ = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 2, INPUTS)
   assert inststr.solar_channel_flag.tolist() == [1, 0]
   assert inststr.thermal_channel_flag.tolist() == [0, 1]
   assert segmented[0].f0 > 100.0
   assert 0.0 < segmented[1].f0 < 10.0


@pytest.mark.parametrize("mmfile", ["liquid-water_stg.mm", "aerosol_a79.mm"])
def test_load_mmdat_matches_validated_microphysics_reader(mmfile):
   old = read_microphysics(INPUTS / "microphysics" / mmfile)
   mmstr = load_mmdat(INPUTS / "microphysics" / mmfile, INPUTS)
   assert mmstr.substance == old.substance and mmstr.shortname == old.shortname
   assert mmstr.ncomp == len(old.components) and mmstr.nlayer == old.profile_height_km.size
   assert identical(mmstr.height, old.profile_height_km) and identical(mmstr.rext, old.profile_relative_amount)
   assert identical(mmstr.mrat, [c.mixing_ratio for c in old.components])
   assert identical(mmstr.rm, [c.size_parameters[0] for c in old.components])
   assert identical(mmstr.s, [c.size_parameters[1] for c in old.components])
   assert mmstr.distname == [c.size_distribution for c in old.components]
   for comp, old_comp in zip(mmstr.comp, old.components):
      wl, n, k = read_refractive_index(INPUTS / "ri" / old_comp.refractive_index_file)
      assert comp.code == old_comp.scattering_code
      assert identical(comp.wl, wl) and identical(np.real(comp.cm), n) and identical(np.imag(comp.cm), -k)


def test_load_atmstr_matches_validated_atmosphere_reader():
   old = read_atmosphere(_atmosphere_path(INPUTS, 2), 2)
   atmstr = load_atmstr(INPUTS / "atm" / "mls.atm", 2)
   assert atmstr.nlevels == 46
   assert identical(atmstr.height, old.height_km) and identical(atmstr.pressure, old.pressure_hpa)
   assert identical(atmstr.temperature, old.temperature_k)
   assert atmstr.height[0] == 100.0 and atmstr.height[-1] == 0.0      # top of atmosphere first


def test_load_gasstr_matches_validated_gas_reader_and_the_atmosphere_grid():
   atmstr = load_atmstr(INPUTS / "atm" / "mls.atm", 2)
   gasstr = load_gasstr(2, "meteosat-10", "seviri", [1, 9], INPUTS / "gas")
   for onegasstr, channel in zip(gasstr, (1, 9)):
      assert int(onegasstr.channelid) == channel
      assert identical(onegasstr.tau_gas, _gas_levels(ROOT, read_instrument(INST), 2, channel))
      assert identical(onegasstr.height, atmstr.height)              # the IDL ARRAY_EQUAL check


def test_interpol_matches_the_validated_idl_compatible_profile_interpolation():
   height = np.array([5.5, 4.5, 3.5, 2.5, 1.5], dtype=np.float32)
   amount = np.array([0.0, 0.0, 0.0, 472.37, 778.80], dtype=np.float32)
   query = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 5.5, 6.5, 26.25, 97.5], dtype=np.float32)
   assert identical(interpol(amount, height, query), _interpolate_profile(height, amount, query))
   # IDL: INTERPOL([1,2,4,8], [1,2,3,4], [0, 5]) -> 0.0, 12.0 (linear extrapolation both ends)
   np.testing.assert_allclose(interpol([1, 2, 4, 8], [1, 2, 3, 4], [0.0, 5.0]), [0.0, 12.0], rtol=1e-6)


# ---------------------------------------------------------------------------
# scattering: quadrature / mie_size_dist_new / legpexp / create_bwgp / generate_scattering_properties
# ---------------------------------------------------------------------------

def test_quadrature_gauss_points_match_the_validated_nodes():
   abscissa, weight = quadrature("g", 64)
   expected_abscissa, expected_weight = np.polynomial.legendre.leggauss(64)
   assert identical(abscissa, expected_abscissa) and identical(weight, expected_weight)
   abscissa, weight = quadrature("T", 5)
   assert identical(abscissa, np.linspace(-1, 1, 5)) and identical(weight, [0.25, 0.5, 0.5, 0.5, 0.25])


@pytest.mark.parametrize("distname", ["modified_gamma", "log_normal"])
def test_mie_size_dist_new_and_legpexp_match_the_validated_kernels(distname):
   abscissa, weight = quadrature("g", 64)
   qv = -abscissa
   refractive_index = complex(np.complex64(1.331 - 1e-8j))
   wavelength = 0.6381768
   if distname == "modified_gamma":
      params = [12.0, 0.1111111, 0.001, 100.0]
      old = _mie_distribution(params[0], params[1], wavelength, refractive_index, qv, xres=0.4)
   else:
      params = [0.091, 1.7, 0.001, 100.0]
      old = _mie_distribution(params[0], 0.0, wavelength, refractive_index, qv, xres=0.4,
                              distribution="log_normal", size_parameters=(params[0], params[1]))
   bext, bsca, w, g, spm, vavg = mie_size_dist_new(distname, 1.0, params, 1.0 / wavelength, refractive_index, qv, xres=0.4)
   assert (bext, w, g, vavg) == (old[0], old[1], old[2], old[4])
   assert identical(spm[0], old[3])
   inlc, lc = legpexp(64, qv, weight, spm[0])
   assert identical(lc, legendre_moments(qv, weight, old[3]))
   assert 0 < inlc < 64


def test_create_bwgp_runs_one_mie_calculation_per_wavelength():
   abscissa, _ = quadrature("g", 32)
   qv = -abscissa
   bext, w, g, phi, vavg = create_bwgp("modified_gamma", 5.5, 0.1111111, [1.33 - 1e-8j, 1.32 - 1e-3j], [0.55, 0.64], qv)
   assert bext.shape == (2,) and phi.shape == (32, 2) and vavg > 0
   assert 0 < w[1] < w[0] <= 1.0
   with pytest.raises(FileNotFoundError, match="tmatrix_path"):
      create_bwgp("log_normal", 0.07, 1.7, [1.5 - 0.04j], [0.55], qv, scode="tmatrix",
                  eps=[1.0], neps=[1.0])


def test_dubovik_metadata_readers_follow_the_idl_record_order(tmp_path):
   names = tmp_path / "name.dat"
   names.write_text("1\nK11\nK12\nK22\nK33\nK34\nK44\nKEXT\n2\nignored\n")
   result = read_dubovik_filenames(names)
   assert result.kext == ("KEXT",) and result.k11 == ("K11",) and result.k44 == ("K44",)
   grid = tmp_path / "grid1.dat"
   grid.write_text("3 0.340\n0.1 0.2 0.4\n")
   assert read_dubovik_grid(grid).radius.tolist() == [0.1, 0.2, 0.4]


def test_dubovik_kernel_readers_parse_small_wrapped_records(tmp_path):
   (tmp_path / "grid1.dat").write_text("2 0.340\n0.1 0.2\n")
   kext = tmp_path / "Rkext"
   kext.write_text("0.1 1.0 1.0\n-2\nreal\nimag\n1 -1\n0 1\n0.340 1.3 -0.01\nEXTINCTION\n1 2\nABSORPTION\n0.1 0.2\n")
   matrix = tmp_path / "Rkernel.11"
   matrix.write_text("0.1 1.0 1.0\n-2\n3\n0 90 180\nreal\nimag\n1 -1\n11 1\n0.340 1.3 -0.01\n1 2 3\n4 5 6\n")
   result = read_dubovik_kernel_kext(kext)
   scattering = read_dubovik_kernel_scatt_matrix(matrix)
   assert result.kext.shape == (1, 1, 2) and result.kext[0, 0].tolist() == [1.0, 2.0]
   assert result.kabs[0, 0].tolist() == [0.1, 0.2]
   assert scattering.k11.shape == (1, 1, 2, 3)
   assert scattering.k11[0, 0, 1].tolist() == [4.0, 5.0, 6.0]


def test_dubovik_branch_reports_a_missing_database(tmp_path):
   with pytest.raises(FileNotFoundError):
      create_bwgp("log_normal", 0.07, 1.7, [1.5 - 0.04j], [0.55], [-0.5], scode="tmatrix",
                  tmatrix_path=tmp_path, eps=[1.0], neps=[1.0])


@pytest.mark.slow
def test_generate_scattering_properties_matches_the_validated_cloud_optics(monkeypatch):
   inststr = load_inststr(INST, requestedchannelid=[1, 9])
   lutstr = load_lutstr(INPUTS / "lut" / "liquid-water-cloud_test.lut", inststr.max_sat_zenith)
   srfstrarr, nwvl_max = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
   mmstr = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS)
   with pytest.raises(ValueError, match="nmom is obsolete for Mie"):
      generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, 1000)
   # The validated path integrates liquid water over the IDL's fixed 0.001-100 um
   # and expands it in a fixed number of moments.  Port check of the bulk
   # optics (bitwise) and of the moments (where the fixed expansion is
   # converged, radii 1-10 um) with the IDL radius limits imposed:
   gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
   with monkeypatch.context() as patch:
      patch.setattr(gsp, "radius_upper_factor", lambda mmstr, c: None)
      (lmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom) = \
         generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, None)
   wavelengths = np.asarray([0.55] + [float(s.wvl_centre) for s in srfstrarr])
   ri_wavelength, ri_real, ri_imaginary = read_refractive_index(INPUTS / "ri" / "H2O_Segelstein_1981.ri")
   real = np.interp(wavelengths, ri_wavelength, ri_real).astype(np.float32)
   imaginary = np.interp(wavelengths, ri_wavelength, ri_imaginary).astype(np.float32)
   optics = legacy_stg_optics(lutstr.efr, wavelengths, real.astype(np.complex64) - 1j * imaginary.astype(np.complex64), phase_order=1000)
   assert bext.shape == (1, 2, 3) and lmom.shape == (1, 2, 3) and amom.shape == (int(lmom.max()), 1, 2, 3)
   assert identical(bext550, optics.extinction_coefficient[:, 0].astype(np.float32))
   assert identical(bext[0].T, optics.extinction_coefficient[:, 1:].astype(np.float32))
   assert identical(bextrat[0].T, (optics.extinction_coefficient[:, 1:] / optics.reference_extinction_coefficient[:, None]).astype(np.float32))
   assert identical(w[0].T, optics.single_scatter_albedo[:, 1:].astype(np.float32))
   assert identical(g[0].T, optics.asymmetry_parameter[:, 1:].astype(np.float32))
   assert identical(vavg, optics.average_volume_per_particle.astype(np.float32))
   fixed = np.transpose(optics.phase_moments[:, :, 1:], (0, 2, 1))       # (1000, channel, radius)
   for l in range(2):
      for r in range(3):
         # With the 100 um integration the far tail leaves coefficients near the
         # rounding-noise level up to the degree bound; when that noise exceeds
         # King's 1e-9 the expansion conservatively keeps them (L <= D + 1).
         length = int(lmom[0, l, r])
         common = min(length, 1000)
         assert length >= 10
         assert np.allclose(amom[:common, 0, l, r], fixed[:common, l, r], rtol=0.0, atol=1e-7)
         assert np.all(amom[length:, 0, l, r] == 0.0)
         if length < 1000:
            assert np.max(np.abs(fixed[length:, l, r] * (2.0 * np.arange(length, 1000) + 1.0))) < 1e-6
         else:
            assert np.max(np.abs(amom[1000:length, 0, l, r] * (2.0 * np.arange(1000, length) + 1.0))) < 1e-6
   # Production radius limit (3.5 x effective radius): only the removed tail changes the bulk optics.
   adopted = generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, None)
   assert np.allclose(adopted[7], bext, rtol=1e-4, atol=0.0) and np.allclose(adopted[1], bext550, rtol=1e-4, atol=0.0)
   assert np.allclose(adopted[8], w, rtol=0.0, atol=2e-5) and np.allclose(adopted[9], g, rtol=0.0, atol=2e-5)


# ---------------------------------------------------------------------------
# DISORT: setup_disort / call_disort
# ---------------------------------------------------------------------------

def test_setup_disort_mirrors_the_idl_common_block():
   disort_vars = setup_disort(60, 45, 2, 2, 1000)
   assert (disort_vars.nlyr, disort_vars.nmom, disort_vars.ntau, disort_vars.numu, disort_vars.nphi, disort_vars.nstr) == (45, 999, 2, 4, 2, 60)
   assert disort_vars.usrtau == 1 and disort_vars.usrang == 1 and disort_vars.lamber == 1 and disort_vars.albedo == 0.0
   with pytest.raises(NotImplementedError):
      setup_disort(60, 45, 2, 2, 1000, alb=0.1)


def test_call_disort_compression_and_kernel_match_the_validated_wrapper():
   nmom = 61
   rayleigh = getmom(2, 0.0, nmom - 1)
   particle = np.linspace(1.0, 0.0, nmom, dtype=np.float32)
   dtau = np.asarray([0.01, 0.02, 0.5, 0.03], dtype=np.float32)
   ssalb = np.asarray([1.0, 1.0, 0.9, 1.0], dtype=np.float32)
   pmom = np.asfortranarray(np.stack([rayleigh, rayleigh, particle, rayleigh], axis=1))
   utau = np.asarray([0.0, np.sum(dtau, dtype=np.float32)], dtype=np.float32)
   umu = np.asarray([-0.99999, -0.5, 0.5, 0.99999], dtype=np.float32)
   phi = np.asarray([0.0, 180.0], dtype=np.float32)
   disort_vars = setup_disort(16, 4, 2, 2, nmom)
   with pytest.raises(ValueError, match="DTAUC dimension missmatch"):
      call_disort(disort_vars, dtau[:3], ssalb[:3], pmom[:, :3], utau, umu, phi, 0.0, 0.6, 100.0)
   new = call_disort(disort_vars, dtau, ssalb, pmom, utau, umu, phi, 0.0, 0.6, 100.0)
   old = old_call_disort(dtau, ssalb, pmom, utau, umu, phi, 0.0, 0.6, 100.0, nstreams=16)
   for value, name in zip(new, ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")):
      assert identical(value, old[name]), name
   # thermal call (no compression)
   layers = [2]
   new = call_disort(disort_vars, dtau[layers], ssalb[layers], np.asfortranarray(pmom[:, layers]),
                     np.asarray([0.0, dtau[2]], dtype=np.float32), umu, phi, 0.0, 0.6, 0.0,
                     plank=True, wnlo=920.0, wnhi=930.0, temp=np.full(2, 250.0, dtype=np.float32), nlayer=1)
   old = old_call_disort(dtau[layers], ssalb[layers], np.asfortranarray(pmom[:, layers]),
                         np.asarray([0.0, dtau[2]], dtype=np.float32), umu, phi, 0.0, 0.6, 0.0, nstreams=16,
                         plank=True, wavenumber_low=920.0, wavenumber_high=930.0,
                         temperature=np.full(2, 250.0, dtype=np.float32))
   assert identical(new[5], old["uu"])


# ---------------------------------------------------------------------------
# end to end: which Legendre-expansion definition a LUT may use
# ---------------------------------------------------------------------------
# The V21 compact validation products that this file previously required
# bitwise were made with the IDL's fixed nmom = 1000 expansion and 0.001-100 um
# liquid integration; they are reproduced with the source revision that made
# them (9d663e9).  Current code refuses nmom for Mie classes; Baum / T-matrix
# classes still need it.

def test_mie_luts_refuse_the_fixed_nmom_expansion(tmp_path):
   with pytest.raises(ValueError, match="nmom is obsolete for Mie"):
      create_orac_luts.create_orac_cloud_lut(
         INPUTS, "meteosat-10_seviri_v1.inst", "liquid-water_stg.mm", "liquid-water-cloud_test.lut", tmp_path, 2,
         channelid=[1], srf_quad=1, version=21, nstreams=60, nmom=1000)
   with pytest.raises(ValueError, match="nmom is obsolete for Mie"):
      create_orac_luts.create_orac_aerosol_lut(
         INPUTS, "meteosat-10_seviri_v1.inst", "aerosol_a79.mm", "aerosol_test.lut", tmp_path, 2,
         channelid=[1], gas=1, srf_quad=1, version=21, nstreams=60, nmom=1000)


def test_tabulated_baum_luts_still_require_nmom(tmp_path):
   with pytest.raises(ValueError, match="nmom is required for Baum"):
      create_orac_luts.create_orac_cloud_lut(
         INPUTS, "earthcare_msi_v1.inst", "water-ice_agg.mm", "ice-cloud_test.lut", tmp_path, 2,
         channelid=[1], srf_quad=1, version=None, nstreams=60, nmom=None)
