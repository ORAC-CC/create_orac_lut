"""V25 vertically inhomogeneous cloud profile (src/oraclut/cloud_temperature.py) and its use by create_orac_luts.py.

The V25 formulation: the supplied cirrostratus profile set
(references/data/ocalut_cloudprofile_Cirrostratus.dat) gives, for 13 total
optical depths, the distribution of optical depth with depth below the cloud
top and the temperature departure dT = 8 K/km x depth; the cloud depth H(COT)
is interpolated linearly in COT; T = 240 K + dT at the boundaries of
equal-optical-depth emission layers.  No representative extinction
coefficient, saturated adiabat, depth cap or absolute altitude is involved.
"""

import inspect
from pathlib import Path

import numpy as np
import pytest

import create_orac_luts
from oraclut import cloud_temperature as ct
from oraclut.idl_mirror.load_lutstr import load_lutstr

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
DATA = ROOT / "references" / "data" / "ocalut_cloudprofile_Cirrostratus.dat"
NATIVE_COT = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]
# Cloud-base depth (km) and temperature departure (K) stated for the supplied data
SUPPLIED_BASE = {0.0625: (2.0000, 16.0000), 0.125: (2.0774, 16.6194), 0.25: (2.2323, 17.8581), 0.5: (2.5419, 20.3355),
                 1: (3.1613, 25.2903), 2: (4.4000, 35.2000), 4: (5.2000, 41.6000), 8: (7.0000, 56.0000),
                 16: (9.7200, 77.7600), 32: (10.8067, 86.4533), 64: (10.9100, 87.2800), 128: (11.0000, 88.0000),
                 256: (11.0000, 88.0000)}


@pytest.fixture(scope="module")
def profiles():
   return ct.load_cloud_profiles(DATA)


@pytest.fixture(scope="module")
def model():
   return ct.cloud_vertical_profile_model("cirrostratus", 256.0)


# ---------------------------------------------------------------------------
# The supplied file
# ---------------------------------------------------------------------------

def test_file_provenance_and_shape(profiles):
   assert DATA.is_file() and profiles.md5 == "7d7b472f3ff9681901c03afe85c3202b"
   assert profiles.cloud_type == "Cirrostratus"
   assert profiles.origin == "/tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf"
   assert profiles.cot.tolist() == NATIVE_COT                      # all 13 COTs load
   for name in ("dcot", "cumulative", "depth_km", "delta_t_k"):
      assert getattr(profiles, name).shape == (13, 99)             # all 99 rows of every profile
   assert profiles.extinction_percent.shape == (99,)


def test_layer_optical_depths_sum_to_the_stated_total(profiles):
   totals = profiles.dcot.sum(axis=1)
   assert np.allclose(totals, profiles.cot, rtol=5e-5)
   assert np.all(profiles.dcot >= 0.0)
   assert np.allclose(profiles.cumulative[:, -1], profiles.cot, rtol=5e-5)   # cumulative ends at the total
   assert np.all(np.diff(profiles.cumulative, axis=1) > 0.0)


def test_extinction_shape_is_common_and_normalised(profiles):
   assert profiles.extinction_percent.sum() == pytest.approx(100.0, abs=1e-3)
   # dCOT_i = COT f_i / 100 for every profile (to the file's six-decimal rounding)
   for k, cot in enumerate(profiles.cot):
      assert np.allclose(profiles.native_layer_optical_depths(cot), profiles.dcot[k], atol=6e-7 * max(1.0, cot) + 1e-6)


def test_depth_and_temperature_structure(profiles):
   # rows equally spaced from the cloud top to the base; dT = 8 K/km x z in every row
   for k in range(13):
      z, dt = profiles.depth_km[k], profiles.delta_t_k[k]
      assert z[0] == 0.0 and dt[0] == 0.0
      assert np.allclose(z, z[-1] * np.arange(99) / 98.0, atol=1e-5)
      assert np.allclose(dt, ct.LAPSE_K_PER_KM * z, atol=1e-4)
      assert np.all(np.diff(z) > 0.0) and np.all(np.diff(dt) > 0.0)
   assert ct.LAPSE_K_PER_KM == 8.0


def test_supplied_cloud_base_values_are_reproduced(profiles, model):
   for cot, (depth, dt) in SUPPLIED_BASE.items():
      assert profiles.depth_at(cot) == pytest.approx(depth, abs=5e-5)
      assert model.cloud_base_temperature_k(cot) - 240.0 == pytest.approx(dt, abs=5e-4)
   assert profiles.depth_scale_km[-1] == 11.0                       # the cloud never deepens beyond 11 km


def test_native_profiles_are_reproduced_exactly_at_their_own_cot(profiles, model):
   # the row temperatures at the rows' (midpoint) cumulative optical-depth fractions
   fractions = profiles.optical_fraction_midpoints
   for k, cot in enumerate(profiles.cot):
      dt = profiles.delta_t_at(cot, fractions)
      assert np.allclose(dt, profiles.delta_t_k[k], atol=1e-4)
      assert np.allclose(model.boundary_temperatures_k(cot, fractions), 240.0 + profiles.delta_t_k[k], atol=1e-4)
   assert profiles.depth_fraction(0.0) == 0.0 and profiles.depth_fraction(1.0) == 1.0


def test_parser_rejects_a_corrupted_file(tmp_path):
   text = DATA.read_text().splitlines()
   bad = tmp_path / "bad.dat"
   bad.write_text("\n".join(text[:200]) + "\n")                     # truncated: fewer than 13 profiles
   with pytest.raises(ValueError):
      ct.load_cloud_profiles(bad)
   lines = list(text)
   lines[10] = lines[10].replace(lines[10].split()[1], "9.0")       # a layer dCOT that breaks the total
   bad.write_text("\n".join(lines) + "\n")
   with pytest.raises(ValueError):
      ct.load_cloud_profiles(bad)


# ---------------------------------------------------------------------------
# Interpolation to the LUT grid and the emission layers
# ---------------------------------------------------------------------------

def test_interpolation_preserves_the_structure_on_grid_b(model):
   lutstr = load_lutstr(INPUTS / "lut" / "ice-cloud-grid-b.lut", 75.0)
   assert float(np.max(lutstr.opd)) == 256.0                        # no extrapolation above the supplied maximum
   for tau in lutstr.opd:
      tau = float(tau)
      f, dt = model.profiles.layer_boundaries(tau, 46)
      assert f[0] == 0.0 and f[-1] == 1.0 and dt[0] == 0.0           # cloud-top departure zero
      assert np.all(np.diff(dt) > 0.0)                               # monotonic temperature
      assert np.all(np.diff(f) > 0.0)
      assert 2.0 <= model.profiles.depth_at(tau) <= 11.0
   assert model.profiles.depth_at(0.03125) == 2.0                    # held below the smallest supplied COT
   with pytest.raises(ValueError):
      model.profiles.depth_at(300.0)


def test_linear_in_cot_interpolation_between_native_profiles(profiles):
   # H is exactly linear in COT in the supplied data from 0.0625 to 1
   h = profiles.depth_scale_km
   slope = np.diff(h[:5]) / np.diff(profiles.cot[:5])
   assert np.allclose(slope, slope[0], rtol=1e-3)
   assert profiles.depth_at(0.75) == pytest.approx(h[3] + 0.5 * (h[4] - h[3]))
   assert profiles.depth_at(2.0) == 4.4


def test_zero_optical_depth_is_safe(model):
   with np.errstate(all="raise"):
      t = model.boundary_temperatures_k(0.0, np.linspace(0.0, 1.0, 13))
   assert t[0] == 240.0 and np.all(np.isfinite(t)) and np.all(np.diff(t) >= 0.0)
   small = ct.cloud_vertical_profile_model("cirrostratus", 0.0)
   assert small.cloud_base_temperature_k(0.0) == pytest.approx(240.0 + 16.0)   # the 2 km profile, whose emission vanishes with tau


def test_emission_layers_conserve_optical_depth_and_properties(model):
   dtau = np.asarray([2.0, 1.0], np.float32)
   ssalb = np.asarray([0.5, 0.6], np.float32)
   pmom = np.asarray([[1.0, 1.0], [0.8, 0.7]], np.float32)
   weights = np.asarray([0.75, 0.25], np.float32)
   emtau, emssa, empmo, temper = ct.emission_layers(dtau, ssalb, pmom, weights, 64.0, model, nlayers=4)
   assert emtau.shape == (8,) and emssa.shape == (8,) and empmo.shape == (2, 8) and temper.shape == (9,)
   assert np.allclose(emtau[:4].sum(), 2.0) and np.allclose(emtau[4:].sum(), 1.0)    # total optical depth preserved
   assert np.all(emtau >= 0.0)
   assert np.all(emssa[:4] == 0.5) and np.all(emssa[4:] == 0.6)
   assert np.all(empmo[1, :4] == 0.8) and np.all(empmo[1, 4:] == 0.7)
   assert temper[0] == 240.0                                                          # cloud top
   assert np.all(np.diff(temper) > 0.0)                                               # monotonic
   assert temper[4] == pytest.approx(model.boundary_temperatures_k(64.0, [0.75])[0], abs=1e-3)
   assert temper[-1] == pytest.approx(240.0 + 87.28, abs=1e-3)
   assert ct.EMISSION_LAYERS <= ct.DISORT_MAX_LAYERS
   with pytest.raises(ValueError, match="MXCLY"):
      ct.emission_layers(dtau, ssalb, pmom, weights, 64.0, model, nlayers=24)


def test_all_channels_share_one_profile(model):
   pmom_a = np.ones((5, 1), np.float32)
   pmom_b = np.full((5, 1), 0.5, np.float32)
   _, _, _, ta = ct.emission_layers([3.0], [0.5], pmom_a, [1.0], 16.0, model)
   _, _, _, tb = ct.emission_layers([0.7], [0.99], pmom_b, [1.0], 16.0, model)
   assert np.array_equal(ta, tb) and ta.size == ct.EMISSION_LAYERS + 1
   assert ta[-1] == pytest.approx(240.0 + 77.76, abs=1e-3)


def test_no_beta_adiabat_cap_or_altitude_in_the_production_path():
   source = inspect.getsource(ct).lower()
   for word in ("beta_ext", "lapse_rate(", "adiabatic_profile", "hydrostatic(", "h_max", "p_top_hpa", "tau / beta",
                "saturation_vapour", "latent_heat", "load_atmstr", "atmstr", "altitude_km"):
      assert word not in source, word
   assert set(inspect.signature(ct.cloud_vertical_profile_model).parameters) == {"name", "max_tau_055", "t_top_k", "path"}
   generator = inspect.getsource(create_orac_luts)
   assert "adiabat" not in generator.lower() and "beta_ext" not in generator


# ---------------------------------------------------------------------------
# Run-file setting and production run files
# ---------------------------------------------------------------------------

def write_run(tmp_path, extra, mm="water-ice_sph.mm", lut="ice-cloud_test.lut"):
   text = "\n".join([
      "platform = 'meteosat-10'", "instrument = 'seviri'", "forward_model = 'cloud'",
      "in_path = 'create_orac_lut/input_files'", "instfile = 'meteosat-10_seviri_v1.inst'",
      f"mmfile = '{mm}'", f"lutfile = '{lut}'", "atmospheres = 2",
      "channelid = [9]", "gas = 0", "no_rayleigh = 0", "srf_quad = 1", "nstreams = 60", "version = 25",
      "out_path = 'validation/tmp/v25/runfile'"] + extra) + "\n"
   run = tmp_path / "case.run"
   run.write_text(text)
   return run


def test_run_file_defaults_to_the_legacy_isothermal_cloud(tmp_path):
   assert create_orac_luts.read_runfile(write_run(tmp_path, []))["cloud_vertical_profile"] == "isothermal"
   assert create_orac_luts.read_runfile(write_run(tmp_path, ["cloud_vertical_profile = 'cirrostratus'"]))[
      "cloud_vertical_profile"] == "cirrostratus"


def test_run_file_rejects_unknown_profiles_and_aerosol_use(tmp_path):
   with pytest.raises(ValueError, match="cloud_vertical_profile"):
      create_orac_luts.read_runfile(write_run(tmp_path, ["cloud_vertical_profile = 'adiabatic'"]))
   run = write_run(tmp_path, ["cloud_vertical_profile = 'cirrostratus'"])
   run.write_text(run.read_text().replace("forward_model = 'cloud'", "forward_model = 'aerosol'"))
   with pytest.raises(ValueError, match="cloud forward model only"):
      create_orac_luts.read_runfile(run)


def test_generator_refuses_the_cirrostratus_profile_for_liquid_water(tmp_path):
   with pytest.raises(ValueError, match="ice-cloud profile"):
      create_orac_luts.create_orac_cloud_lut(INPUTS, "meteosat-10_seviri_v1.inst", "liquid-water_stg.mm",
                                             "liquid-water-cloud_test.lut", tmp_path, 2, channelid=[9], srf_quad=1,
                                             version=25, work_path=tmp_path / "work",
                                             cloud_vertical_profile="cirrostratus")
   with pytest.raises(ValueError, match="cloud_vertical_profile"):
      create_orac_luts.create_orac_cloud_lut(INPUTS, "meteosat-10_seviri_v1.inst", "water-ice_sph.mm",
                                             "ice-cloud_test.lut", tmp_path, 2, channelid=[9], srf_quad=1,
                                             version=25, work_path=tmp_path / "work",
                                             cloud_vertical_profile="adiabatic")


def test_v25_production_run_files_select_the_profile_by_phase():
   runs = sorted((ROOT / "runs").glob("*_v25.run"))
   assert len(runs) == 40
   for run in runs:
      settings = create_orac_luts.read_runfile(run)
      expected = "cirrostratus" if "water-ice" in run.name else "isothermal"
      assert settings["cloud_vertical_profile"] == expected, run.name
      assert settings["version"] == 25, run.name
      twin = run.with_name(run.name.replace("_v25.run", "_v24.run"))
      v24 = create_orac_luts.read_runfile(twin)
      for key, value in v24.items():
         if key not in ("version", "cloud_vertical_profile"):
            assert settings[key] == value, (run.name, key)


# ---------------------------------------------------------------------------
# End to end on the compact ice test grid: V24 unchanged, V25 changes E_md only
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_generator_isothermal_is_unchanged_and_cirrostratus_changes_only_e_md(tmp_path):
   from oraclut.io.lut import read_lut

   products = {}
   for mode, options in (("default", {}), ("isothermal", {"cloud_vertical_profile": "isothermal"}),
                         ("cirrostratus", {"cloud_vertical_profile": "cirrostratus"})):
      out = tmp_path / mode
      out.mkdir()
      status = create_orac_luts.create_orac_cloud_lut(
         INPUTS, "meteosat-10_seviri_v1.inst", "water-ice_sph.mm", "ice-cloud_test.lut", out, 2,
         channelid=[4, 9], srf_quad=1, version=25, nstreams=60, work_path=out / "work", **options)
      assert status == 0
      products[mode] = read_lut(next(out.glob("*.nc")))
   default, iso, new = products["default"], products["isothermal"], products["cirrostratus"]
   for name in default.variable_names:
      assert np.array_equal(default.variables[name], iso.variables[name]), name
      if name != "E_md":
         assert np.array_equal(iso.variables[name], new.variables[name]), name
   assert not np.array_equal(iso.variables["E_md"], new.variables["E_md"])
   assert np.all(new.variables["E_md"] >= iso.variables["E_md"])          # a warmer lower cloud emits more
   assert np.all(np.isfinite(new.variables["E_md"]))
   assert [int(c) for c in new.variables["thermal_channel_id"]] == [4, 9]   # mixed (3.9 um) and thermal (10.8 um)
   assert "cloud_vertical_profile" not in iso.global_attributes
   assert new.global_attributes["cloud_vertical_profile"] == "cirrostratus"
   assert new.global_attributes["cloud_top_temperature_K"] == 240.0
   assert new.global_attributes["cloud_vertical_profile_md5"] == "7d7b472f3ff9681901c03afe85c3202b"
   assert new.global_attributes["cloud_emission_layers"] == ct.EMISSION_LAYERS
   assert tuple(iso.variable_attributes["E_md"]["valid_range"]) == (0.0, 1.0)
   assert new.variable_attributes["E_md"]["valid_range"][1] > 1.0
