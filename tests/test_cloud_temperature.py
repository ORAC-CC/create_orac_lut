"""V25 cloud temperature profile (src/oraclut/cloud_temperature.py) and its use by create_orac_luts.py.

The V25 formulation: T_top = 240 K (reference pressure 628 hPa); the path
length below the cloud top z = x / beta_ext_055 for cumulative 0.55-um
optical depth x, not limited and independent of any altitude; a phase-pure
saturated adiabat (liquid water or ice) throughout, with the pressure evolved
hydrostatically along it.
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
GRID_B_TAU_MAX = 256.0


@pytest.fixture(scope="module")
def liquid():
   return ct.cloud_temperature_model("liquid-water", GRID_B_TAU_MAX)


@pytest.fixture(scope="module")
def ice():
   return ct.cloud_temperature_model("water-ice", GRID_B_TAU_MAX)


# ---------------------------------------------------------------------------
# Thermodynamics
# ---------------------------------------------------------------------------

def test_saturation_vapour_pressures_match_murphy_and_koop_reference_values():
   # Murphy and Koop (2005): e_sw(273.16) = 611.657 Pa; e_si(273.16) = 611.657 Pa (triple point)
   assert ct.saturation_vapour_pressure_liquid(273.16) == pytest.approx(611.657, rel=2e-4)
   assert ct.saturation_vapour_pressure_ice(273.16) == pytest.approx(611.657, rel=2e-4)
   # supercooled water is more volatile than ice below the triple point
   assert ct.saturation_vapour_pressure_liquid(240.0) > ct.saturation_vapour_pressure_ice(240.0)


def test_latent_heats_are_in_the_textbook_range():
   assert ct.latent_heat_vaporisation(273.15) == pytest.approx(2.501e6)
   assert ct.latent_heat_sublimation(273.15) == pytest.approx(2.834e6, rel=2e-3)
   assert ct.latent_heat_sublimation(240.0) > ct.latent_heat_vaporisation(240.0)


def test_saturated_lapse_rates_lie_between_known_limits():
   dry = ct.G0 / ct.C_P_DRY
   for phase in ("liquid", "ice"):
      cold = ct.saturated_lapse_rate(240.0, 62800.0, phase)
      warm = ct.saturated_lapse_rate(290.0, 90000.0, phase)
      assert 0.85 * dry < cold < dry                      # little vapour at 240 K: close to the dry rate
      assert 0.003 < warm < 0.0045                        # about 3-4.5 K/km at 290 K
   with pytest.raises(ValueError):
      ct.saturated_lapse_rate(240.0, 62800.0, "mixed")


# ---------------------------------------------------------------------------
# The reference state and the path-length coordinate
# ---------------------------------------------------------------------------

def test_cloud_top_is_exactly_240_k_for_both_phases(liquid, ice):
   assert ct.T_TOP_K == 240.0
   for model in (liquid, ice):
      assert model.t_top_k == 240.0
      assert model.temperature_k[0] == 240.0
      assert model.temperature_at_depth_k(0.0) == 240.0
      assert model.boundary_temperatures_k(16.0, [0.0, 0.5, 1.0])[0] == 240.0
      assert model.cloud_base_temperature_k(0.0) == 240.0
      assert model.p_top_hpa == 628.0 and model.pressure_hpa[0] == pytest.approx(628.0)


def test_extinction_coefficients_are_exactly_the_model_constants(liquid, ice):
   assert ct.CLOUD_MODELS == {"liquid-water": ("liquid", 20.0), "water-ice": ("ice", 1.0)}
   assert liquid.beta_ext_055_per_km == 20.0
   assert ice.beta_ext_055_per_km == 1.0


def test_distance_is_exactly_optical_depth_over_beta(liquid, ice):
   assert liquid.depth_of_optical_depth_km(10.0) == 0.5            # liquid: x = 10 -> 0.5 km
   assert ice.depth_of_optical_depth_km(1.0) == 1.0                 # ice: x = 1 -> 1 km
   assert liquid.depth_of_optical_depth_km(1.0) == 0.05
   assert np.array_equal(liquid.depth_of_optical_depth_km([0.0, 2.0, 10.0, 256.0]), [0.0, 0.1, 0.5, 12.8])
   assert np.array_equal(ice.depth_of_optical_depth_km([0.0, 4.0, 64.0, 256.0]), [0.0, 4.0, 64.0, 256.0])
   with pytest.raises(ValueError):
      liquid.depth_of_optical_depth_km(-1.0)


def test_there_is_no_depth_limit_in_the_model(liquid):
   source = inspect.getsource(ct).lower()
   for word in ("h_max", "hmax", "geometric_depth", "maximum_geometric", "capped", "cap "):
      assert word not in source, word
   assert set(inspect.signature(ct.cloud_temperature_model).parameters) == {"substance", "max_tau_055", "t_top_k", "p_top_hpa"}
   for name in ("h_max_km", "geometric_depth_km", "cloud_top_km"):
      assert not hasattr(liquid, name)


def test_distance_keeps_increasing_with_optical_depth(liquid, ice):
   taus = np.array([0.5, 1.0, 4.0, 16.0, 64.0, 128.0, 256.0])
   for model in (liquid, ice):
      z = model.depth_of_optical_depth_km(taus)
      assert np.all(np.diff(z) > 0.0)
      assert np.allclose(np.diff(z) / np.diff(taus), 1.0 / model.beta_ext_055_per_km)   # never saturates
      t = np.array([model.cloud_base_temperature_k(tau) for tau in taus])
      assert np.all(np.diff(t) > 0.0)


def test_liquid_uses_the_liquid_adiabat_and_ice_the_ice_adiabat_throughout():
   # Each tabulated profile is the integration of its own phase's lapse rate
   # at every point, including far below the freezing level.
   for substance, own, other in (("liquid-water", "liquid", "ice"), ("water-ice", "ice", "liquid")):
      model = ct.cloud_temperature_model(substance, GRID_B_TAU_MAX)
      assert model.phase == own
      n = len(model.depth_km)
      for i in (0, n // 4, n // 2, n - 2):
         t, p = model.temperature_k[i], 100.0 * model.pressure_hpa[i]
         h = 1000.0 * (model.depth_km[i + 1] - model.depth_km[i])
         step = (model.temperature_k[i + 1] - t) / h
         assert step == pytest.approx(ct.saturated_lapse_rate(t, p, own), rel=2e-3)
      # at the freezing level the two adiabats differ by more than the integration tolerance
      i = int(np.searchsorted(model.temperature_k, 273.15))
      t, p = model.temperature_k[i], 100.0 * model.pressure_hpa[i]
      assert abs(ct.saturated_lapse_rate(t, p, own) - ct.saturated_lapse_rate(t, p, other)) > 1e-5


def test_no_freezing_ceiling_or_phase_switch(liquid, ice):
   # Both adiabats pass through the freezing level without any change of
   # behaviour: strictly increasing temperature and a smoothly varying lapse rate.
   for model in (liquid, ice):
      t = model.temperature_k
      assert t[-1] > 273.15 and np.all(np.diff(t) > 0.0)
      crossing = int(np.searchsorted(t, 273.15))
      rates = np.diff(t[crossing - 50:crossing + 50])                # 10 m steps either side of 273.15 K
      assert np.all(rates > 0.0)
      curvature = np.abs(np.diff(rates))                            # a phase switch or ceiling would be an outlier
      assert curvature.max() < 3.0 * np.median(curvature)
   source = inspect.getsource(ct)
   assert "273.15" not in source.replace("T_FREEZE = 273.15", "")       # the constant serves the latent-heat fit only


def test_absolute_altitude_does_not_enter(liquid):
   # The model is built from the substance and the grid's largest optical
   # depth only: no atmosphere, no layer profile, no cloud height.
   parameters = inspect.signature(ct.cloud_temperature_model).parameters
   assert "atmstr" not in parameters and "scatreltau" not in parameters
   assert not hasattr(liquid, "cloud_top_km")
   source = inspect.getsource(ct).lower()
   assert "load_atmstr" not in source and "atmstr" not in source and "height" not in source and "altitude" in source
   # the pressure along the adiabat is hydrostatic in the adiabat's own temperature
   z, t, p = liquid.depth_km, liquid.temperature_k, 100.0 * liquid.pressure_hpa
   i = len(z) // 2
   expected = p[i] * ct.G0 / (ct.R_DRY * t[i])
   assert (p[i + 1] - p[i - 1]) / (2000.0 * (z[i + 1] - z[i])) == pytest.approx(expected, rel=1e-3)


def test_all_channels_share_one_thermodynamic_profile(liquid):
   # The temperatures depend on the node's 0.55-um optical depth and the
   # particle layer fractions only: two channels with different layer optical
   # depths and optical properties receive the same boundary temperatures.
   pmom_a = np.ones((5, 1), np.float32)
   pmom_b = np.full((5, 1), 0.5, np.float32)
   _, _, _, ta = ct.emission_layers([3.0], [0.5], pmom_a, [1.0], 16.0, liquid)
   _, _, _, tb = ct.emission_layers([0.7], [0.99], pmom_b, [1.0], 16.0, liquid)
   assert np.array_equal(ta, tb)
   assert ta.size == ct.EMISSION_SUBLAYERS + 1
   assert ta[0] == 240.0 and ta[-1] == pytest.approx(liquid.cloud_base_temperature_k(16.0), abs=1e-3)


def test_no_microphysical_extinction_enters_the_distance():
   source = inspect.getsource(ct)
   assert "generate_scattering_properties" not in source.replace("of\ngenerate_scattering_properties", "")
   assert "import" not in source.split("def saturation_vapour_pressure_liquid")[0].split('"""', 2)[2].replace(
      "from dataclasses import dataclass", "").replace("import numpy as np", "")
   with pytest.raises(ValueError):
      ct.cloud_temperature_model("sulphuric-acid", 1.0)


def test_zero_optical_depth_is_safe(liquid, ice):
   for model in (liquid, ice):
      with np.errstate(all="raise"):
         assert model.depth_of_optical_depth_km(0.0) == 0.0
         t = model.boundary_temperatures_k(0.0, np.linspace(0.0, 1.0, 13))
      assert np.all(t == 240.0)
      assert model.cloud_base_temperature_k(0.0) == 240.0
   small = ct.cloud_temperature_model("water-ice", 0.0)                   # a grid whose largest node is zero
   assert small.cloud_base_temperature_k(0.0) == 240.0


def test_emission_sublayers_preserve_the_layer_optical_depth_and_properties(liquid):
   dtau = np.asarray([2.0, 1.0], np.float32)
   ssalb = np.asarray([0.5, 0.6], np.float32)
   pmom = np.asarray([[1.0, 1.0], [0.8, 0.7]], np.float32)
   weights = np.asarray([0.75, 0.25], np.float32)
   emtau, emssa, empmo, temper = ct.emission_layers(dtau, ssalb, pmom, weights, 64.0, liquid, nsub=4)
   assert emtau.shape == (8,) and emssa.shape == (8,) and empmo.shape == (2, 8) and temper.shape == (9,)
   assert np.allclose(emtau[:4].sum(), 2.0) and np.allclose(emtau[4:].sum(), 1.0)
   assert np.all(emssa[:4] == 0.5) and np.all(emssa[4:] == 0.6)
   assert np.all(empmo[1, :4] == 0.8) and np.all(empmo[1, 4:] == 0.7)
   # the boundary after the first layer sits at 75 % of the node's optical depth: z = 0.75 x 64 / 20 km
   assert temper[4] == pytest.approx(liquid.temperature_at_depth_k(0.75 * 64.0 / 20.0), abs=1e-3)
   assert temper[-1] == pytest.approx(liquid.temperature_at_depth_k(64.0 / 20.0), abs=1e-3)
   assert ct.EMISSION_SUBLAYERS == 12
   with pytest.raises(ValueError, match="MXCLY"):
      ct.emission_layers(dtau, ssalb, pmom, weights, 64.0, liquid, nsub=24)


# ---------------------------------------------------------------------------
# Run-file setting and production run files
# ---------------------------------------------------------------------------

def write_run(tmp_path, extra):
   text = "\n".join([
      "platform = 'meteosat-10'", "instrument = 'seviri'", "forward_model = 'cloud'",
      "in_path = 'create_orac_lut/input_files'", "instfile = 'meteosat-10_seviri_v1.inst'",
      "mmfile = 'liquid-water_stg.mm'", "lutfile = 'liquid-water-cloud_test.lut'", "atmospheres = 2",
      "channelid = [9]", "gas = 0", "no_rayleigh = 0", "srf_quad = 1", "nstreams = 60", "version = 25",
      "out_path = 'validation/tmp/v25/runfile'"] + extra) + "\n"
   run = tmp_path / "case.run"
   run.write_text(text)
   return run


def test_run_file_defaults_to_the_legacy_isothermal_cloud(tmp_path):
   assert create_orac_luts.read_runfile(write_run(tmp_path, []))["cloud_temperature_profile"] == "isothermal"
   assert create_orac_luts.read_runfile(write_run(tmp_path, ["cloud_temperature_profile = 'adiabatic'"]))[
      "cloud_temperature_profile"] == "adiabatic"


def test_run_file_rejects_unknown_profiles_and_aerosol_use(tmp_path):
   with pytest.raises(ValueError, match="cloud_temperature_profile"):
      create_orac_luts.read_runfile(write_run(tmp_path, ["cloud_temperature_profile = 'linear'"]))
   run = write_run(tmp_path, ["cloud_temperature_profile = 'adiabatic'"])
   run.write_text(run.read_text().replace("forward_model = 'cloud'", "forward_model = 'aerosol'"))
   with pytest.raises(ValueError, match="cloud forward model only"):
      create_orac_luts.read_runfile(run)


def test_generator_rejects_an_unknown_profile(tmp_path):
   with pytest.raises(ValueError, match="cloud_temperature_profile"):
      create_orac_luts.create_orac_cloud_lut(INPUTS, "meteosat-10_seviri_v1.inst", "liquid-water_stg.mm",
                                             "liquid-water-cloud_test.lut", tmp_path, 2, channelid=[9], srf_quad=1,
                                             version=25, work_path=tmp_path / "work",
                                             cloud_temperature_profile="linear")


def test_v25_production_run_files_select_the_adiabatic_profile_and_version_25():
   runs = sorted((ROOT / "runs").glob("*_v25.run"))
   assert len(runs) == 40
   for run in runs:
      settings = create_orac_luts.read_runfile(run)
      assert settings["cloud_temperature_profile"] == "adiabatic", run.name
      assert settings["version"] == 25, run.name
      twin = run.with_name(run.name.replace("_v25.run", "_v24.run"))
      v24 = create_orac_luts.read_runfile(twin)
      for key, value in v24.items():
         if key not in ("version", "cloud_temperature_profile"):
            assert settings[key] == value, (run.name, key)


def test_production_grid_extremes_are_finite_and_monotonic():
   # The deepest Grid B node, tau_055 = 256, maps to 12.8 km (liquid) and 256 km (ice).
   for lutfile, substance in (("liquid-water-cloud-grid-b.lut", "liquid-water"), ("ice-cloud-grid-b.lut", "water-ice")):
      lutstr = load_lutstr(INPUTS / "lut" / lutfile, 75.0)
      tau_max = float(np.max(lutstr.opd))
      assert tau_max == GRID_B_TAU_MAX
      model = ct.cloud_temperature_model(substance, tau_max)
      assert model.depth_km[-1] == pytest.approx(tau_max / model.beta_ext_055_per_km)
      assert np.all(np.isfinite(model.temperature_k)) and np.all(np.isfinite(model.pressure_hpa))
      assert np.all(np.diff(model.temperature_k) > 0.0) and np.all(np.diff(model.pressure_hpa) > 0.0)
      for tau in lutstr.opd:
         assert np.isfinite(model.cloud_base_temperature_k(float(tau)))


# ---------------------------------------------------------------------------
# End to end on the compact test grid: V24 unchanged, V25 changes E_md only
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_generator_isothermal_is_unchanged_and_adiabatic_changes_only_e_md(tmp_path):
   from oraclut.io.lut import read_lut

   products = {}
   for mode, options in (("default", {}), ("isothermal", {"cloud_temperature_profile": "isothermal"}),
                         ("adiabatic", {"cloud_temperature_profile": "adiabatic"})):
      out = tmp_path / mode
      out.mkdir()
      status = create_orac_luts.create_orac_cloud_lut(
         INPUTS, "meteosat-10_seviri_v1.inst", "liquid-water_stg.mm", "liquid-water-cloud_test.lut", out, 2,
         channelid=[4, 9], srf_quad=1, version=25, nstreams=60, work_path=out / "work", **options)
      assert status == 0
      products[mode] = read_lut(next(out.glob("*.nc")))
   default, iso, adi = products["default"], products["isothermal"], products["adiabatic"]
   for name in default.variable_names:
      assert np.array_equal(default.variables[name], iso.variables[name]), name
      if name != "E_md":
         assert np.array_equal(iso.variables[name], adi.variables[name]), name
   assert not np.array_equal(iso.variables["E_md"], adi.variables["E_md"])
   assert np.all(adi.variables["E_md"] >= iso.variables["E_md"])          # a warmer lower cloud emits more
   assert np.all(np.isfinite(adi.variables["E_md"]))
   assert [int(c) for c in adi.variables["thermal_channel_id"]] == [4, 9]   # mixed (3.9 um) and thermal (10.8 um)
   assert "cloud_temperature_profile" not in iso.global_attributes
   assert adi.global_attributes["cloud_temperature_profile"] == "adiabatic"
   assert adi.global_attributes["cloud_top_temperature_K"] == 240.0
   assert adi.global_attributes["cloud_top_reference_pressure_hPa"] == 628.0
   assert adi.global_attributes["cloud_extinction_coefficient_055um_per_km"] == 20.0
   assert not any("geometric" in k or "height" in k for k in adi.global_attributes)
   assert tuple(iso.variable_attributes["E_md"]["valid_range"]) == (0.0, 1.0)
   assert adi.variable_attributes["E_md"]["valid_range"][1] > 1.0
