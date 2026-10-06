"""V25 cloud temperature profile (src/oraclut/cloud_temperature.py) and its use by create_orac_luts.py."""

import inspect
from pathlib import Path

import numpy as np
import pytest

import create_orac_luts
from oraclut import cloud_temperature as ct
from oraclut.idl_mirror.interpol import interpol
from oraclut.idl_mirror.load_atmstr import load_atmstr
from oraclut.idl_mirror.load_mmdat import load_mmdat

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"


def production_profile(mmfile):
   """(mmstr, atmstr, scatreltau) exactly as create_orac_cloud_lut forms them for atmosphere 2."""

   mmstr = load_mmdat(INPUTS / "microphysics" / mmfile, INPUTS)
   atmstr = load_atmstr(INPUTS / "atm" / "mls.atm", 2)
   nlayers = atmstr.nlevels - 1
   hlayers = (atmstr.height[0:nlayers] + atmstr.height[1:nlayers + 1]) / np.float32(2)
   scatreltau = interpol(mmstr.rext, mmstr.height, hlayers)
   return mmstr, atmstr, scatreltau / np.sum(scatreltau, dtype=np.float32)


@pytest.fixture(scope="module")
def liquid():
   mmstr, atmstr, scatreltau = production_profile("liquid-water_stg.mm")
   return ct.cloud_temperature_model(mmstr.substance, atmstr, scatreltau)


@pytest.fixture(scope="module")
def ice():
   mmstr, atmstr, scatreltau = production_profile("water-ice_sph.mm")
   return ct.cloud_temperature_model(mmstr.substance, atmstr, scatreltau)


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


def test_liquid_and_ice_thermodynamics_differ(liquid, ice):
   # same cloud top, different adiabats
   assert liquid.phase == "liquid" and ice.phase == "ice"
   assert liquid.cloud_top_km == ice.cloud_top_km == 4.0
   assert ct.saturated_lapse_rate(260.0, 70000.0, "liquid") != ct.saturated_lapse_rate(260.0, 70000.0, "ice")
   assert not np.allclose(liquid.temperature_at_depth_k([1.0, 2.0, 2.5]), ice.temperature_at_depth_k([1.0, 2.0, 2.5]), atol=1e-3)


# ---------------------------------------------------------------------------
# Profile and the optical-depth to geometric-depth model
# ---------------------------------------------------------------------------

def test_cloud_top_is_exactly_240_k(liquid, ice):
   for model in (liquid, ice):
      assert model.t_top_k == 240.0
      assert model.temperature_k[0] == 240.0
      assert model.temperature_at_depth_k(0.0) == 240.0
      assert model.boundary_temperatures_k(16.0, [0.0, 0.5, 1.0])[0] == 240.0
      assert model.cloud_base_temperature_k(0.0) == 240.0


def test_temperature_increases_monotonically_downward(liquid, ice):
   for model in (liquid, ice):
      assert np.all(np.diff(model.temperature_k) > 0.0)
      for tau in (0.25, 1.0, 4.0, 16.0, 64.0, 256.0):
         t = model.boundary_temperatures_k(tau, np.linspace(0.0, 1.0, 13))
         assert np.all(np.diff(t) > 0.0), tau


def test_cloud_top_state_comes_from_the_production_atmosphere(liquid):
   assert liquid.cloud_top_km == 4.0                          # top of the 3-4 km layer holding the particles
   assert liquid.cloud_top_hpa == pytest.approx(628.0)


def test_zero_optical_depth_is_handled_without_division_or_gradient(liquid, ice):
   for model in (liquid, ice):
      with np.errstate(all="raise"):
         assert model.geometric_depth_km(0.0) == 0.0
         assert np.all(model.depth_of_optical_depth_km(0.0, [0.0, 0.0]) == 0.0)
         t = model.boundary_temperatures_k(0.0, np.linspace(0.0, 1.0, 13))
      assert np.all(t == 240.0)
   with pytest.raises(ValueError):
      liquid.geometric_depth_km(-1.0)


def test_nominal_thickness_follows_the_representative_extinction(liquid, ice):
   assert liquid.beta_ext_055_per_km == 20.0 and ice.beta_ext_055_per_km == 1.0
   assert liquid.geometric_depth_km(10.0) == pytest.approx(0.5)      # liquid: tau 10 -> 0.5 km
   assert ice.geometric_depth_km(1.0) == pytest.approx(1.0)          # ice: tau 1 -> 1 km
   assert liquid.geometric_depth_km(1.0) == pytest.approx(0.05)
   # uncapped mapping z = x / beta
   assert np.allclose(liquid.depth_of_optical_depth_km(10.0, [0.0, 2.0, 10.0]), [0.0, 0.1, 0.5])


def test_thickness_caps(liquid, ice):
   assert liquid.h_max_km == 2.5 and ice.h_max_km == 6.0
   assert liquid.geometric_depth_km(50.0) == 2.5                     # nominal 2.5 km: at the cap
   assert liquid.geometric_depth_km(256.0) == 2.5
   assert ice.geometric_depth_km(6.0) == 6.0
   assert ice.geometric_depth_km(256.0) == 6.0
   assert liquid.geometric_depth_km(49.0) < 2.5 and ice.geometric_depth_km(5.9) < 6.0


def test_capped_cloud_is_mapped_continuously_over_the_whole_optical_depth(liquid, ice):
   for model, tau in ((liquid, 256.0), (ice, 64.0)):
      fractions = np.linspace(0.0, 1.0, 25)
      z = model.depth_of_optical_depth_km(tau, fractions * tau)
      assert np.allclose(z, model.h_max_km * fractions)             # linear across 0..H_max
      t = model.boundary_temperatures_k(tau, fractions)
      assert np.all(np.diff(t) > 0.0)                               # no isothermal remainder
      assert t[-1] == pytest.approx(model.temperature_at_depth_k(model.h_max_km))
      steps = np.diff(t)
      assert np.all(steps[1:] / steps[:-1] > 0.8)                   # smooth: the lapse rate only falls slowly with depth
      assert np.all(steps[1:] / steps[:-1] < 1.05)


def test_all_channels_share_one_physical_cloud(liquid):
   # The temperatures depend on the 0.55-um optical depth and the particle
   # layer fractions only: two channels with different layer optical depths
   # and optical properties receive the same boundary temperatures.
   pmom_a = np.ones((5, 1), np.float32)
   pmom_b = np.full((5, 1), 0.5, np.float32)
   _, _, _, ta = ct.emission_layers([3.0], [0.5], pmom_a, [1.0], 16.0, liquid)
   _, _, _, tb = ct.emission_layers([0.7], [0.99], pmom_b, [1.0], 16.0, liquid)
   assert np.array_equal(ta, tb)
   assert ta.size == ct.EMISSION_SUBLAYERS + 1
   assert ta[0] == 240.0 and ta[-1] == pytest.approx(liquid.cloud_base_temperature_k(16.0), abs=1e-3)


def test_no_microphysical_extinction_enters_the_physical_depth():
   # The model is built from the substance, atmosphere and layer profile alone;
   # the per-particle / normalised extinction of generate_scattering_properties
   # is not an input.
   parameters = inspect.signature(ct.cloud_temperature_model).parameters
   assert set(parameters) == {"substance", "atmstr", "scatreltau", "t_top_k"}
   source = inspect.getsource(ct)
   assert "generate_scattering_properties" not in source.replace("of\ngenerate_scattering_properties", "")  # not imported
   assert "import" not in source.split("def saturation_vapour_pressure_liquid")[0].split('"""', 2)[2].replace(
      "from dataclasses import dataclass", "").replace("import numpy as np", "")
   with pytest.raises(ValueError):
      ct.cloud_temperature_model("sulphuric-acid", None, [1.0])


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
   # the boundary after the first layer sits at 75 % of the (capped) depth
   assert temper[4] == pytest.approx(liquid.temperature_at_depth_k(0.75 * 2.5), abs=1e-3)
   assert temper[-1] == pytest.approx(liquid.temperature_at_depth_k(2.5), abs=1e-3)
   assert np.all(np.abs(np.diff(temper)) < 10.0)                   # within DISORT's accuracy guidance


def test_dynamic_range_of_the_sublayer_temperature_steps(liquid, ice):
   for model in (liquid, ice):
      t = model.boundary_temperatures_k(256.0, np.linspace(0.0, 1.0, ct.EMISSION_SUBLAYERS + 1))
      assert np.diff(t).max() < 5.0                                 # ice: 4.6 K in the first 0.5 km; liquid: 1.9 K


# ---------------------------------------------------------------------------
# Run-file setting
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
   assert [int(c) for c in adi.variables["thermal_channel_id"]] == [4, 9]   # mixed (3.9 um) and thermal (10.8 um)
   assert np.all(np.isfinite(adi.variables["E_md"]))
   assert "cloud_temperature_profile" not in iso.global_attributes
   assert adi.global_attributes["cloud_temperature_profile"] == "adiabatic"
   assert adi.global_attributes["cloud_top_temperature_K"] == 240.0
   assert adi.global_attributes["cloud_extinction_coefficient_055um_per_km"] == 20.0
   assert adi.global_attributes["cloud_maximum_geometric_depth_km"] == 2.5
   assert tuple(iso.variable_attributes["E_md"]["valid_range"]) == (0.0, 1.0)
   assert adi.variable_attributes["E_md"]["valid_range"][1] > 1.0
