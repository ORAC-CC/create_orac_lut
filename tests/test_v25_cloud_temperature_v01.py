"""LUT release V25: the cloud temperature profile of the thermal-emission call (oraclut.idl_mirror.cloud_temperature).

The profile data, the interpolation to the LUT optical depths, the
liquid-water adiabat (against the historical V25 values), the sub-layered
emission-call inputs, the run-file setting, and the DISORT interface: with a
non-scattering cloud the emitted radiance at the top must equal the analytic
integral of the Planck source that is linear in optical depth between the
boundary temperatures, top down.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import create_orac_luts   # noqa: E402
from oraclut.idl_mirror import cloud_temperature as ct   # noqa: E402

PROFILE_SHA256 = "8e9e68bba7f65ef0231dd309d1c129a0292fed64b9a49f53ec70458d23042aa6"   # md5 7d7b472f3ff9681901c03afe85c3202b
SUPPLIED_DEPTH = {0.0625: 2.0, 0.125: 2.0774, 0.25: 2.2323, 0.5: 2.5419, 1.0: 3.1613, 2.0: 4.4, 4.0: 5.2, 8.0: 7.0,
                  16.0: 9.72, 32.0: 10.8067, 64.0: 10.91, 128.0: 11.0, 256.0: 11.0}


def test_supplied_cirrostratus_profile_is_read_and_verified():
   profstr = ct.read_cloud_profile(ct.CIRROSTRATUS_FILE)
   assert profstr.sha256 == PROFILE_SHA256
   assert profstr.cloud_type == "Cirrostratus"
   assert profstr.origin == "/tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf"
   assert profstr.extinction.size == ct.PROFILE_ROWS and abs(profstr.extinction.sum() - 1.0) < 1e-5
   assert np.array_equal(profstr.cot, list(SUPPLIED_DEPTH))
   assert np.allclose(profstr.depth, list(SUPPLIED_DEPTH.values()), atol=1e-4)


def test_corrupted_profile_file_is_refused(tmp_path):
   lines = ct.CIRROSTRATUS_FILE.read_text().splitlines()
   row = lines[7].split()
   row[4] = f"{float(row[4]) + 0.01:.6f}"                    # temperature no longer 8 K/km x depth
   lines[7] = "   ".join(row)
   bad = tmp_path / "bad.dat"
   bad.write_text("\n".join(lines) + "\n")
   with pytest.raises(ValueError, match="8.0 K/km"):
      ct.read_cloud_profile(bad)


def test_cirrostratus_temperatures_follow_the_supplied_profile():
   cloudprof = ct.setup_cloud_profile("cirrostratus", 256.0, 250.0)
   fraction = np.linspace(0.0, 1.0, 41)
   for cot, depth in SUPPLIED_DEPTH.items():
      t = ct.cloud_temperature(cloudprof, cot, fraction)
      assert t[0] == 250.0 and abs(t[-1] - (250.0 + 8.0 * depth)) < 1e-3
      assert np.all(np.diff(t) > 0.0)
   # linear in COT between supplied profiles, held below the smallest, never extrapolated
   depth = cloudprof.profstr.depth                            # as tabulated (7.000001 km at COT 8)
   assert abs(ct.profile_depth(cloudprof.profstr, 5.65685424949238) - (depth[6] + (depth[7] - depth[6]) * (5.65685424949238 - 4.0) / 4.0)) < 1e-12
   assert ct.profile_depth(cloudprof.profstr, 1e-10) == ct.profile_depth(cloudprof.profstr, 0.0625)
   with pytest.raises(ValueError, match="exceeds the largest supplied profile"):
      ct.setup_cloud_profile("cirrostratus", 300.0, 250.0)


def test_wet_adiabat_reproduces_the_historical_v25_reference_state():
   # validation/v25/results/REPORT_v25_cloud_temperature_profile.md section 7 (240 K, 628 hPa, tau 256 -> 12.8 km)
   historical = ct.setup_cloud_profile("wet_adiabat", 256.0, 240.0)
   assert abs(historical.z[-1] - 12.8) < 1e-12
   assert abs(historical.t[-1] - 317.49) < 0.005 and abs(historical.p[-1] - 2957.2) < 0.05
   for tau, t in ((16.0, 247.09), (64.0, 266.22), (128.0, 286.60)):
      assert abs(ct.cloud_temperature(historical, tau, [1.0])[0] - t) < 0.005
   # V25: the same adiabat from the V24 cloud-top temperature 250 K
   v25 = ct.setup_cloud_profile("wet_adiabat", 256.0, 250.0)
   t = ct.cloud_temperature(v25, 256.0, np.linspace(0.0, 1.0, 41))
   assert t[0] == 250.0 and np.all(np.diff(t) > 0.0) and 315.0 < t[-1] < 325.0


def test_wet_adiabatic_lapse_rate_limits():
   dry = ct.G0 / ct.C_P_DRY * 1000.0                          # K/km
   assert 0.9 * dry < 1000.0 * ct.wet_adiabatic_lapse_rate(220.0, 30000.0) < dry
   assert 3.0 < 1000.0 * ct.wet_adiabatic_lapse_rate(293.15, 100000.0) < 5.0


def test_emission_layers_conserve_the_layer_and_follow_the_profile():
   cloudprof = ct.setup_cloud_profile("cirrostratus", 256.0, 250.0)
   f32 = np.float32
   emtau, emssa = np.asarray([3.2], f32), np.asarray([0.6], f32)
   empmo = np.asfortranarray(np.linspace(1.0, 0.0, 11, dtype=f32)[:, None])
   subtau, subssa, subpmo, subtemp = ct.emission_layers(emtau, emssa, empmo, np.asarray([1.0], f32), 16.0, cloudprof, 250.0)
   n = ct.EMISSION_LAYERS
   assert subtau.size == n and subssa.size == n and subpmo.shape == (11, n) and subtemp.size == n + 1
   assert subtau.dtype == subssa.dtype == subpmo.dtype == subtemp.dtype == np.float32
   assert abs(float(subtau.sum()) - 3.2) < 1e-5 and np.all(subssa == f32(0.6)) and np.all(subpmo == empmo)
   assert np.array_equal(subtemp, ct.cloud_temperature(cloudprof, 16.0, np.arange(n + 1) / n).astype(f32))
   with pytest.raises(ValueError, match="differs from the cloud-top temperature"):
      ct.emission_layers(emtau, emssa, empmo, np.asarray([1.0], f32), 16.0, cloudprof, 240.0)
   with pytest.raises(ValueError, match="MXCLY"):
      ct.emission_layers(np.asarray([1.0, 2.0], f32), np.asarray([0.6, 0.6], f32), np.asfortranarray(np.ones((11, 2), f32)),
                         np.asarray([0.5, 0.5], f32), 16.0, cloudprof, 250.0)


def test_run_file_setting(tmp_path):
   text = (ROOT / "runs" / "template.run").read_text()
   assert create_orac_luts.read_runfile(ROOT / "runs" / "template.run")["cloud_vertical_profile"] == "isothermal"
   for value in ("cirrostratus", "wet_adiabat", "isothermal"):
      path = tmp_path / f"{value}.run"
      path.write_text(text + f"\ncloud_vertical_profile = '{value}'\n")
      assert create_orac_luts.read_runfile(path)["cloud_vertical_profile"] == value
   path = tmp_path / "bad.run"
   path.write_text(text + "\ncloud_vertical_profile = 'adiabat'\n")
   with pytest.raises(ValueError, match="cloud_vertical_profile must be one of"):
      create_orac_luts.read_runfile(path)


def test_production_run_files_select_the_profile_of_their_substance():
   for name, profile in (("liquid-water_old", "wet_adiabat"), ("liquid-water_240", "wet_adiabat"), ("water-ice_agg", "cirrostratus")):
      settings = create_orac_luts.read_runfile(ROOT / "runs" / f"terra_modis_cloud_{name}_v25.run")
      assert settings["cloud_vertical_profile"] == profile and settings["version"] == 25
      assert settings["out_path"] == "staging/lut_v25" and settings["nmom"] == 1000 and settings["nstreams"] == 60
      v24 = create_orac_luts.read_runfile(ROOT / "runs" / f"terra_modis_cloud_{name}_v24.run")
      for key in set(settings) - {"cloud_vertical_profile", "version", "out_path"}:
         assert settings[key] == v24[key], key


def _analytic_top_radiance(tau, planck, mu):
   """Radiance leaving the top of a non-scattering stack with a Planck source linear in optical depth in each layer."""

   radiance, top = 0.0, 0.0
   for k, dtau in enumerate(tau):
      e = np.exp(-dtau / mu)
      slope = (planck[k + 1] - planck[k]) / dtau
      radiance += np.exp(-top / mu) * (planck[k] * (1.0 - e) + slope * (mu * (1.0 - e) - dtau * e))
      top += dtau
   return radiance


def test_disort_emission_uses_each_boundary_temperature_top_down():
   from oraclut.idl_mirror import call_disort, setup_disort
   from oraclut.radiative_transfer.legacy_disort import plkavg

   f32 = np.float32
   cloudprof = ct.setup_cloud_profile("cirrostratus", 256.0, 250.0)
   n = ct.EMISSION_LAYERS
   subtau, subssa, subpmo, subtemp = ct.emission_layers(
      np.asarray([4.0], f32), np.asarray([0.0], f32), np.asfortranarray(np.zeros((61, 1), f32)), np.asarray([1.0], f32),
      16.0, cloudprof, 250.0)
   subpmo[0, :] = 1.0
   saz = np.asarray([0.0, 30.0, 60.0], f32)
   umu = np.concatenate([-np.cos(np.radians(saz)), np.cos(np.radians(saz))[::-1]]).astype(f32)
   umu[np.abs(umu) == 1.0] *= f32(0.99999)
   disort_vars = setup_disort(60, n, saz.size, 1, 61)
   wn = 1e4 / 11.0
   wnlo, wnhi = f32(0.995 * wn), f32(1.005 * wn)
   uu = call_disort(disort_vars, subtau, subssa, subpmo, np.asarray([0.0, subtau.sum()], f32), umu, np.asarray([0.0], f32),
                    0.0, f32(1.0), 0.0, plank=True, wnlo=wnlo, wnhi=wnhi, temp=subtemp, nlayer=n)[5]
   planck = np.asarray([plkavg(float(wnlo), float(wnhi), float(t)) for t in subtemp])
   for i, mu in enumerate(umu[saz.size:]):
      expected = _analytic_top_radiance(subtau.astype(np.float64), planck, float(mu))
      reversed_profile = _analytic_top_radiance(subtau.astype(np.float64), planck[::-1], float(mu))
      assert abs(uu[saz.size + i, 0, 0] / expected - 1.0) < 1e-4
      assert abs(reversed_profile / expected - 1.0) > 0.1     # the check distinguishes the order
