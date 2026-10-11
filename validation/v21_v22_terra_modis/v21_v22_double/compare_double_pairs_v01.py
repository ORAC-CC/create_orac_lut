"""V21 (IDL) / V22 (Python) equivalence benchmark for the double-precision Terra MODIS LUTs.

   python compare_double_pairs.py --v21-dir staging/lut_v21_double --v22-dir staging/lut_v22_double \
                                  --out validation/v21_v22_terra_modis/v21_v22_double/results [--samples 4000]

For each of the three products (pold, p240, pagg):

1. Node level: compare_luts.compare over every variable and coordinate of the
   full LUT grid (including the SZA 80 nodes used by interpolation).
2. Operational domain: the validated measurement-impact calculation of
   operational_domain.py (solar black and Lambertian surfaces, thermal and
   mixed brightness temperatures; SZA < 75 with the 72.5 and 75- slices that
   carry the 80 degree support node; production MODIS uncertainty).
3. ORAC interpolation and Jacobians: at random operational retrieval states
   (SZA in [0, 75), VZA in [0, 65], RAA in [0, 180], log10 tau in
   [-2.1, 2.408], r_e within the grid) the LUT coefficients are interpolated
   as ORAC does (linear in the three angles, src/set_gzero.F90 and
   int_lut_*.F90; bicubic in log10 tau and r_e with finite-difference corner
   derivatives, Int_LUT_Common, mirrored by validation/v24/disort_options/
   orac_fm.bicubic), giving F, dF/dlog10(tau) and dF/dr_e for
     solar channels:   black-surface reflectance R_0v
     thermal channels: brightness temperature (mid-latitude summer atmosphere)
     mixed channels:   brightness temperature of thermal + solar radiance
   Differences V22 - V21 are expressed in units of the production
   measurement uncertainty sigma: forward model |dF|/sigma, Jacobians
   |dJ| * step / sigma with steps 0.1 in log10(tau) and 1 um in r_e.

Acceptance (fixed before the comparison, not tuned to the result): every
forward-model and step-scaled Jacobian difference <= 0.1 sigma.  This keeps
the numerical V21/V22 difference an order of magnitude below the measurement
noise, so it changes a retrieval cost function by < 1 % per channel.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
TERRA = HERE.parent
ROOT = TERRA.parents[1]
sys.path.insert(0, str(TERRA))
sys.path.insert(0, str(ROOT / "validation" / "v24" / "disort_options"))
import compare_luts                                        # noqa: E402
import operational_domain as od                             # noqa: E402
import orac_fm as fm                                        # noqa: E402
from atmosphere import terms as atmosphere_terms            # noqa: E402
from radiance_impact import _values                         # noqa: E402

PRODUCTS = {
   "pold": ("terra_modis_m_liquid-water_a01_pold", "liquid"),
   "p240": ("terra_modis_m_liquid-water_a01_p240", "liquid"),
   "pagg": ("terra_modis_m_water-ice_a01_pagg", "ice"),
}
TOLERANCE_SIGMA = 0.1
STEP_LOG10_TAU = 0.1
STEP_RE = 1.0


def angle_weights(grid, x):
   i0 = int(np.clip(np.searchsorted(grid, x, side="right") - 1, 0, grid.size - 2))
   w = (x - grid[i0]) / (grid[i0 + 1] - grid[i0])
   return i0, w


def interpolate_angles(table, axes, point):
   """Linear interpolation of ``table`` along the named angle axes at ``point`` (ORAC GZero weights)."""
   out = table
   for axis, (grid, x) in sorted(axes.items(), reverse=True):
      i0, w = angle_weights(grid, x)
      out = (1.0 - w) * np.take(out, i0, axis=axis) + w * np.take(out, i0 + 1, axis=axis)
   return out


class Lut:
   def __init__(self, path):
      self.path = path
      with Dataset(path) as d:
         self.tau = _values(d, "optical_depth")
         self.re = _values(d, "effective_radius")
         self.sza = _values(d, "solar_zenith")
         self.vza = _values(d, "satellite_zenith")
         self.raa = _values(d, "relative_azimuth")
         self.channels = _values(d, "channel_id").astype(int)
         self.solar = _values(d, "solar_channel_id").astype(int)
         self.thermal = _values(d, "thermal_channel_id").astype(int)
         self.mixed = set(_values(d, "mixed_channel_id").astype(int))
         self.snr = _values(d, "snr")
         self.nedt = _values(d, "nedt")
         self.refbt = _values(d, "refbt")
         self.f0 = _values(d, "F0")
         self.planck = {int(c): tuple(float(_values(d, n)[i]) for n in ("B1", "B2", "T1", "T2"))
                        for i, c in enumerate(self.thermal)}
         # Tables with (tau, re) first, angles next, channel last.
         self.r0v = _values(d, "R_0v", ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith",
                                        "relative_azimuth", "solar_channels"))
         self.tdv = _values(d, "T_dv", ("optical_depth", "effective_radius", "satellite_zenith", "channels"))
         self.rdv = _values(d, "R_dv", ("optical_depth", "effective_radius", "satellite_zenith", "channels"))
         self.em = _values(d, "E_md", ("optical_depth", "effective_radius", "satellite_zenith", "thermal_channels"))
      self.xt = np.log10(self.tau)

   def solar_reflectance(self, sza, vza, raa, lt, re):
      g = interpolate_angles(self.r0v, {2: (self.sza, sza), 3: (self.vza, vza), 4: (self.raa, raa)}, None)
      return fm.bicubic(g, self.xt, self.re, lt, re)                 # each (n_solar,)

   def thermal_terms(self, vza, lt, re):
      tdv = fm.bicubic(interpolate_angles(self.tdv, {2: (self.vza, vza)}, None), self.xt, self.re, lt, re)
      rdv = fm.bicubic(interpolate_angles(self.rdv, {2: (self.vza, vza)}, None), self.xt, self.re, lt, re)
      em = fm.bicubic(interpolate_angles(self.em, {2: (self.vza, vza)}, None), self.xt, self.re, lt, re)
      return tdv, rdv, em


def sample_states(lut, n, seed):
   rng = np.random.default_rng(seed)
   return np.column_stack([rng.uniform(0.0, 74.999, n), rng.uniform(0.0, min(65.0, lut.vza.max()), n),
                           rng.uniform(0.0, 180.0, n), rng.uniform(-2.1, np.log10(256.0), n),
                           rng.uniform(lut.re.min(), lut.re.max(), n)])


def interpolated_comparison(v21, v22, cloud, n, seed):
   states = sample_states(v21, n, seed)
   solar_only = [c for c in v21.solar if c not in v21.mixed]
   thermal_index = {int(c): i for i, c in enumerate(v21.thermal)}
   all_index = {int(c): i for i, c in enumerate(v21.channels)}
   solar_index = {int(c): i for i, c in enumerate(v21.solar)}
   records = {"solar": [], "thermal": [], "mixed": []}
   atmospheres = {}
   for sza, vza, raa, lt, re in states:
      f21, t21, r21 = v21.solar_reflectance(sza, vza, raa, lt, re)
      f22, t22, r22 = v22.solar_reflectance(sza, vza, raa, lt, re)
      for c in solar_only:
         i = solar_index[c]
         noise = 1.0 / v21.snr[i] if v21.snr[i] > 0 else 0.0
         sigma = float(fm.sigma_solar(abs(f21[i]), noise, production=True))
         if abs(f21[i]) < 0.01 or sigma <= 0:
            continue
         records["solar"].append((c, sza, vza, raa, lt, re, f21[i], (f22[i] - f21[i]) / sigma,
                                  (t22[i] - t21[i]) * STEP_LOG10_TAU / sigma, (r22[i] - r21[i]) * STEP_RE / sigma))
      key = round(float(vza), 1)
      if key not in atmospheres:
         atmospheres[key] = atmosphere_terms("terra", "modis", list(v21.channels), v21.planck, cloud, float(vza))
      atm = atmospheres[key]
      th21 = v21.thermal_terms(vza, lt, re)
      th22 = v22.thermal_terms(vza, lt, re)
      for c in v21.thermal:
         c = int(c)
         e = atm[c]
         ai, ti = all_index[c], thermal_index[c]

         def radiance(terms, solar=None):
            (tdv, dtdv_t, dtdv_r), (rdv, drdv_t, drdv_r), (em, dem_t, dem_r) = terms
            values = []
            for k, (a, b, m) in enumerate(((tdv[ai], rdv[ai], em[ti]), (dtdv_t[ai], drdv_t[ai], dem_t[ti]),
                                            (dtdv_r[ai], drdv_r[ai], dem_r[ti]))):
               r = (e["Rbc_up"] * a + e["B_c"] * m + e["Rac_dwn"] * b) * e["Tac_lw"] + (e["Rac_up"] if k == 0 else 0.0)
               values.append(r)
            if solar is not None:
               transmission = e["Tac"] ** (1.0 / np.cos(np.deg2rad(sza))) * e["Tac"] ** (1.0 / np.cos(np.deg2rad(vza)))
               f0 = v21.f0[solar_index[c]]
               values = [values[k] + f0 * transmission * solar[k] for k in range(3)]
            return values

         is_mixed = c in v21.mixed
         s21 = (f21[solar_index[c]], t21[solar_index[c]], r21[solar_index[c]]) if is_mixed else None
         s22 = (f22[solar_index[c]], t22[solar_index[c]], r22[solar_index[c]]) if is_mixed else None
         rad21, rad22 = radiance(th21, s21), radiance(th22, s22)
         bt21, dbdr21 = fm.r2t(rad21[0], v21.planck[c])
         bt22, dbdr22 = fm.r2t(rad22[0], v21.planck[c])
         sigma = float(fm.sigma_bt(bt21, v21.nedt[ti], v21.refbt[ti], v21.planck[c], mixed=is_mixed, production=True))
         jt21, jr21 = dbdr21 * rad21[1], dbdr21 * rad21[2]
         jt22, jr22 = dbdr22 * rad22[1], dbdr22 * rad22[2]
         records["mixed" if is_mixed else "thermal"].append(
            (c, sza, vza, raa, lt, re, bt21, (bt22 - bt21) / sigma, (jt22 - jt21) * STEP_LOG10_TAU / sigma,
             (jr22 - jr21) * STEP_RE / sigma))
   summary = {}
   names = ("channel", "sza", "vza", "raa", "log10_tau", "re", "reference", "dF_sigma", "dJtau_sigma", "dJre_sigma")
   for kind, rows in records.items():
      if not rows:
         continue
      a = np.asarray(rows, dtype=float)
      kind_summary = {"cases": int(a.shape[0])}
      for col, label in ((7, "forward_model"), (8, "jacobian_log10_tau"), (9, "jacobian_re")):
         v = np.abs(a[:, col])
         worst = int(np.argmax(v))
         kind_summary[label] = {
            "max_sigma": float(v.max()), "rms_sigma": float(np.sqrt(np.mean(v ** 2))),
            "p99_sigma": float(np.percentile(v, 99)), "p999_sigma": float(np.percentile(v, 99.9)),
            "count_over_tolerance": int(np.sum(v > TOLERANCE_SIGMA)),
            "worst_case": {n: float(a[worst, i]) for i, n in enumerate(names)}}
      summary[kind] = kind_summary
   return summary


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("--v21-dir", type=Path, required=True)
   parser.add_argument("--v22-dir", type=Path, required=True)
   parser.add_argument("--out", type=Path, required=True)
   parser.add_argument("--samples", type=int, default=4000)
   parser.add_argument("--v21-suffix", default="_v21.nc")
   parser.add_argument("--v22-suffix", default="_v22.nc")
   args = parser.parse_args()
   args.out.mkdir(parents=True, exist_ok=True)
   report = {"tolerance_sigma": TOLERANCE_SIGMA, "step_log10_tau": STEP_LOG10_TAU, "step_re_um": STEP_RE, "products": {}}
   all_pass = True
   for name, (stem, cloud) in PRODUCTS.items():
      v21_path = args.v21_dir / (stem + args.v21_suffix)
      v22_path = args.v22_dir / (stem + args.v22_suffix)
      if not (v21_path.is_file() and v22_path.is_file()):
         report["products"][name] = {"status": "missing", "v21": str(v21_path), "v22": str(v22_path)}
         all_pass = False
         continue
      result, coordinate_rows, call_rows = compare_luts.compare(v21_path, v22_path)
      (args.out / f"{name}_node_comparison.json").write_text(json.dumps(result, indent=1, default=str) + "\n")
      compare_luts._write_csv(args.out / f"{name}_direct_beam_calls.csv", call_rows[:1000])
      with Dataset(v21_path) as reference, Dataset(v22_path) as candidate:
         sza = _values(reference, "solar_zenith")
         geometries, weights = od.operational_sza_weights(sza, od.MAX_SOL_ZEN)
         measurements = od.measurement_rows(reference, candidate, cloud, geometries, weights)
      od.write_csv(args.out / f"{name}_operational_measurements.csv", measurements)
      interpolated = interpolated_comparison(Lut(v21_path), Lut(v22_path), cloud, args.samples, seed=21)
      worst = max(k[m]["max_sigma"] for k in interpolated.values() for m in ("forward_model", "jacobian_log10_tau", "jacobian_re"))
      node_max = {v: result["variables"][v]["maximum_absolute_difference"] for v in result.get("variables", {})
                  if isinstance(result["variables"][v], dict) and "maximum_absolute_difference" in result["variables"][v]}
      passed = worst <= TOLERANCE_SIGMA
      all_pass &= passed
      report["products"][name] = {"v21": str(v21_path), "v22": str(v22_path), "node_maximum_absolute_difference": node_max,
                                  "operational_max_sigma": max((r.get("maximum_sigma", 0.0) for r in measurements
                                                                if r.get("uncertainty") == "production"), default=None),
                                  "interpolated": interpolated, "worst_sigma": worst, "passed": passed}
      print(f"{name}: worst interpolated difference {worst:.3g} sigma -> {'PASS' if passed else 'FAIL'}", flush=True)
   report["all_passed"] = all_pass
   (args.out / "equivalence_summary.json").write_text(json.dumps(report, indent=1, default=float) + "\n")
   print(("ALL THREE PAIRS PASS" if all_pass else "EQUIVALENCE NOT ESTABLISHED") + f"; summary {args.out / 'equivalence_summary.json'}")
   return 0 if all_pass else 1


if __name__ == "__main__":
   raise SystemExit(main())
