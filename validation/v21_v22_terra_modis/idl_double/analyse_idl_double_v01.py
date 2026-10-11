"""Analyse a double-precision DISORT validation run made through IDL.

Reads results/<run>/ written by run_idl_double_validation.sh and reports, per
case and per DISORT call:

* interface test results;
* that the compact single-precision (validation) LUT reproduces the
  production V21 LUT, and that LUT coefficients recomputed here from the
  captured DISORT outputs reproduce the LUTs IDL wrote;
* precision effect: historical single precision versus DISORT_DP at NSTR 60;
* stream-count convergence of DISORT_DP at NSTR 48, 60, 80, 100;
* arithmetic-precision error: DISORT_DP versus a gfortran REAL*16 build;
* cross-build agreement: DISORT_DP in IDL versus the ifort kernel called
  directly and versus a gfortran double build;
* repeatability, physical boundary conditions and energy balance;
* reflectance errors relative to the production MODIS uncertainty, and the
  effect on finite-difference Jacobians in log10(tau) and r_e;
* timing of single- and double-precision calls.

Usage: python analyse_idl_double.py [RUN_DIRECTORY] [--quad-all]
"""

from __future__ import annotations

import argparse
import ctypes
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset
from scipy.io import readsav


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import disort_reference as ref                           # noqa: E402
import radiance_impact as impact                         # noqa: E402
from radiance_impact import fm                           # noqa: E402

PRODUCTION_V21 = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/"
                      "terra_modis_m_liquid-water_a01_p240_v21.nc")
CASES = {
   "A": {"channel": 6, "centre": (256.0, 29.0, 80.0), "geometry": (0.0, 180.0)},
   "B": {"channel": 20, "centre": (1.0, 11.0, 30.0), "geometry": (36.111, 126.0)},
   "C": {"channel": 6, "centre": (0.125, 21.0, 30.0), "geometry": (43.333, 0.0)},
   "D": {"channel": 1, "centre": (8.0, 11.0, 40.0), "geometry": (28.889, 90.0)},
}
OUT = ref.OUTPUTS


# ---------------------------------------------------------------------------
# Captures
# ---------------------------------------------------------------------------

def load_capture(path):
   c = readsav(str(path))
   n = int(c["nlyr1"])
   temper = np.zeros(n + 1)
   given = np.asarray(c["temper_in"], np.float64).ravel()
   temper[:min(n + 1, given.size)] = given[:n + 1]
   args = {
      "nlyr": n, "nmom": int(c["nmom"]),
      "dtauc": np.asarray(c["dtauc1_in"], np.float64)[:n], "ssalb": np.asarray(c["ssalb1_in"], np.float64)[:n],
      # readsav gives (layer, moment), or (moment,) for a single-layer call.
      "pmom": np.atleast_2d(np.asarray(c["pmom1_in"], np.float64))[:n].T, "temper": temper,
      "utau": np.asarray(c["utau_in"], np.float64), "umu": np.asarray(c["umu_in"], np.float64),
      "phi": np.asarray(c["phi_in"], np.float64), "plank": int(c["plank"]),
      **{name: float(c[name]) for name in ref.SCALARS},
   }
   numu, ntau, nphi = len(args["umu"]), len(args["utau"]), len(args["phi"])
   nstr_list = [int(x) for x in np.atleast_1d(c["nstr_list"])]
   # readsav reverses IDL dimension order; restore Fortran order (numu, ntau, nphi[, k]).
   sp = {k: np.asarray(c[k], np.float64) for k in OUT}
   sp["uu"] = np.asarray(c["uu"], np.float64).T
   dp = {}
   for k_index, nstr in enumerate(nstr_list):
      dp[nstr] = {k: np.asarray(c["dp_" + k], np.float64).T[..., k_index] if k == "uu"
                  else np.asarray(c["dp_" + k], np.float64)[k_index] for k in OUT}
   repeat = {k: np.asarray(c["r_" + k], np.float64).T if k == "uu" else np.asarray(c["r_" + k], np.float64)
             for k in OUT}
   for k in OUT:
      if k == "uu":
         for d in [sp, repeat] + list(dp.values()):
            assert d["uu"].shape == (numu, ntau, nphi), d["uu"].shape
   kind = "emission" if args["plank"] else ("direct" if args["fbeam"] > 0 else "diffuse")
   return {"args": args, "sp": sp, "dp": dp, "repeat": repeat, "kind": kind, "nstr_v21": int(c["nstr_v21"]),
           "time_sp": float(c["time_sp"]), "time_dp": dict(zip(nstr_list, map(float, c["dp_time"]))),
           "type_dp_uu": int(c["type_dp_uu"]), "wvnmlo": float(c["wvnmlo"]), "wvnmhi": float(c["wvnmhi"])}


def map_calls(captures, lut_path, channel):
   """Assign (tau, r_e, SZA) to each capture from the generator loop order."""
   with Dataset(lut_path) as d:
      taus = np.asarray(d["optical_depth"][:], np.float64)
      res = np.asarray(d["effective_radius"][:], np.float64)
      szas = np.asarray(d["solar_zenith"][:], np.float64)
   thermal = channel >= 20
   i = 0
   for a, tau in enumerate(taus):
      for r, re in enumerate(res):
         order = ["diffuse"] + (["emission"] if thermal else []) + ["direct"] * len(szas)
         s = 0
         for kind in order:
            cap = captures[i]
            if cap["kind"] != kind:
               raise AssertionError(f"capture {i + 1}: expected {kind}, found {cap['kind']}")
            cap.update({"tau": float(tau), "re": float(re), "tau_index": a, "re_index": r})
            if kind == "direct":
               expected = np.float32(np.cos(np.float32(szas[s]) * np.float32(np.pi / 180.0)))
               if abs(cap["args"]["umu0"] - float(expected)) > 2e-7:
                  raise AssertionError(f"capture {i + 1}: UMU0 {cap['args']['umu0']} not cos({szas[s]})")
               cap.update({"sza": float(szas[s]), "sza_index": s})
               s += 1
            i += 1
   if i != len(captures):
      raise AssertionError(f"{len(captures)} captures, {i} mapped")
   return taus, res, szas


# ---------------------------------------------------------------------------
# LUT coefficients from DISORT outputs (create_orac_cloud_lut / write_v2_lut)
# ---------------------------------------------------------------------------

_PLKAVG = {}


def plkavg(variant, wnlo, wnhi, temperature=250.0):
   if variant not in _PLKAVG:
      library = ctypes.CDLL(str(ref.LIBRARIES[variant]), mode=os.RTLD_LAZY)
      real = ctypes.c_double if variant == "ifort_dp" else ctypes.c_float
      library.plkavg_.restype = real
      _PLKAVG[variant] = (library, real)
   library, real = _PLKAVG[variant]
   return float(library.plkavg_(ctypes.byref(real(wnlo)), ctypes.byref(real(wnhi)), ctypes.byref(real(temperature))))


def coefficients(cap, out):
   """LUT coefficients for one call (already divided by 100 as written)."""
   args = cap["args"]
   umu = args["umu"]
   nsaz = len(umu) // 2
   up = 2 * nsaz - np.arange(nsaz) - 1
   uu = out["uu"]
   if cap["kind"] == "direct":
      r0v = (uu[up, 0, :] * math.pi)[:, ::-1] / 100.0              # (vza, raa)
      return {"R_0v": r0v, "R_0d": out["flup"][0] / out["rfldir"][0], "T_0d": out["rfldn"][1] / out["rfldir"][0],
              "T_00": out["rfldir"][1] / out["rfldir"][0]}
   if cap["kind"] == "diffuse":
      fisot = args["fisot"]
      tdv = (uu[np.arange(nsaz), 1, 0] - 100.0 * np.exp(args["utau"][1] / umu[:nsaz])) / 100.0
      return {"R_dd": out["flup"][0] / (fisot * math.pi), "T_dd": out["rfldn"][1] / (fisot * math.pi),
              "R_dv": uu[up, 0, 0] / 100.0, "T_dv": tdv}
   bbe = plkavg("sp_dlm", cap["wvnmlo"], cap["wvnmhi"])           # IDL uses the single-precision PLKAVG
   return {"E_md": uu[up, 0, 0] / bbe}


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def physical_checks(cap, out):
   """Boundary conditions, sign and energy balance (black surface)."""
   args = cap["args"]
   umu, uu = args["umu"], out["uu"]
   down = umu < 0
   r = {}
   if cap["kind"] == "direct":
      incident = args["fbeam"] * args["umu0"]
      r["toa_downward_diffuse_flux"] = float(out["rfldn"][0])
      r["toa_downward_intensity_max_abs"] = float(np.max(np.abs(uu[down, 0, :])))
      r["direct_flux_relative_error"] = float(abs(out["rfldir"][0] - incident) / incident)
      absorbed = incident - out["flup"][0] - out["rfldir"][1] - out["rfldn"][1]
      r["absorbed_fraction"] = float(absorbed / incident)
      r["reflected_fraction"] = float(out["flup"][0] / incident)
   elif cap["kind"] == "diffuse":
      incident = math.pi * args["fisot"]
      r["toa_downward_diffuse_flux_relative_error"] = float(abs(out["rfldn"][0] - incident) / incident)
      r["toa_downward_intensity_max_abs_error"] = float(np.max(np.abs(uu[down, 0, :] - args["fisot"])))
      absorbed = incident - out["flup"][0] - out["rfldn"][1]
      r["absorbed_fraction"] = float(absorbed / incident)
   else:
      r["toa_downward_intensity_max_abs"] = float(np.max(np.abs(uu[down, 0, :])))
      r["toa_downward_diffuse_flux"] = float(out["rfldn"][0])
   r["minimum_intensity"] = float(np.min(uu))
   r["minimum_flux"] = float(min(np.min(out["flup"]), np.min(out["rfldn"]), np.min(out["uavg"])))
   return r


def max_abs(a, b):
   return float(np.max(np.abs(np.asarray(a, np.float64) - np.asarray(b, np.float64))))


def coefficient_differences(ca, cb):
   return {k: max_abs(ca[k], cb[k]) for k in ca}


def compare_luts(path_a, path_b, names):
   result = {}
   with Dataset(path_a) as a, Dataset(path_b) as b:
      for name in names:
         if name not in a.variables or name not in b.variables:
            continue
         x = np.asarray(a[name][:], np.float64)
         y = np.asarray(b[name][:], np.float64)
         result[name] = {"bitwise_equal": bool(np.array_equal(x, y)), "max_abs": float(np.max(np.abs(x - y)))}
   return result


def compact_vs_production(compact_path, production_path, names):
   """Compare a compact LUT with the production LUT at the same coordinate values."""
   result = {}
   with Dataset(compact_path) as c, Dataset(production_path) as p:
      for name in names:
         if name not in c.variables:
            continue
         dims = c[name].dimensions
         index = []
         for dim in dims:
            coordinate = {"solar_channels": "solar_channel_id", "channels": "channel_id",
                          "thermal_channels": "thermal_channel_id"}.get(dim, dim)
            cv = np.asarray(c[coordinate][:], np.float64)
            pv = np.asarray(p[coordinate][:], np.float64)
            idx = [int(np.flatnonzero(np.isclose(pv, v, rtol=0, atol=1e-6))[0]) for v in cv]
            index.append(idx)
         pvals = np.asarray(p[name][:], np.float64)[np.ix_(*index)]
         cvals = np.asarray(c[name][:], np.float64)
         result[name] = {"bitwise_equal": bool(np.array_equal(cvals, pvals)), "max_abs": float(np.max(np.abs(cvals - pvals))),
                         "elements": int(cvals.size)}
   return result


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyse_case(name, case_dir, quad_all):
   spec = CASES[name]
   vdir = case_dir / "validation"
   gdir = case_dir / "dp_generator"
   captures = [load_capture(p) for p in sorted((vdir / "captures").glob("call_*.sav"))]
   lut_sp = next((vdir / "luts" / f"case_{name}").glob("*_v21.nc"))
   lut_dp = next((gdir / "luts" / f"case_{name}").glob("*_v90.nc"))
   taus, res, szas = map_calls(captures, lut_sp, spec["channel"])
   names = ("R_0v", "R_0d", "T_0d", "T_00", "R_dv", "T_dv", "R_dd", "T_dd", "E_md")
   report = {"channel": spec["channel"], "taus": taus.tolist(), "res": res.tolist(), "szas": szas.tolist(),
             "calls": len(captures), "compact_v21_vs_production_v21": compact_vs_production(lut_sp, PRODUCTION_V21, names),
             "v90_minus_v21_compact": compare_luts(lut_sp, lut_dp, names)}
   nstrs = sorted(captures[0]["dp"])
   rows = []
   for i, cap in enumerate(captures):
      sp_c = coefficients(cap, cap["sp"])
      dp_c = {n: coefficients(cap, cap["dp"][n]) for n in nstrs}
      row = {"call": i + 1, "kind": cap["kind"], "tau": cap["tau"], "re": cap["re"], "sza": cap.get("sza"),
             "output_type_double": cap["type_dp_uu"] == 5,
             "repeat_bitwise": all(np.array_equal(cap["repeat"][k], cap["dp"][cap["nstr_v21"]][k]) for k in OUT),
             "sp_minus_dp60": coefficient_differences(sp_c, dp_c[60]),
             "stream_convergence_vs_100": {str(n): coefficient_differences(dp_c[n], dp_c[100]) for n in nstrs if n != 100},
             "stream_successive": {f"{a}-{b}": coefficient_differences(dp_c[a], dp_c[b]) for a, b in zip(nstrs[:-1], nstrs[1:])},
             "physical_sp": physical_checks(cap, cap["sp"]),
             "physical_dp": {str(n): physical_checks(cap, cap["dp"][n]) for n in nstrs},
             "time_sp": cap["time_sp"], "time_dp": {str(k): v for k, v in cap["time_dp"].items()}}
      # Cross-build: same kernel called directly must be bitwise identical to the IDL DLM.
      status, direct = ref.run("ifort_dp", cap["args"], 60)
      row["idl_dp_vs_ifort_kernel_bitwise"] = bool(status == 0 and all(np.array_equal(direct[k], cap["dp"][60][k]) for k in OUT))
      status, gdp = ref.run("gfort_dp", cap["args"], 60)
      row["idl_dp_vs_gfortran_dp"] = {k: max_abs(gdp[k], cap["dp"][60][k]) for k in OUT}
      row["idl_dp_vs_gfortran_dp_coefficients"] = coefficient_differences(coefficients(cap, gdp), dp_c[60])
      centre = (cap["tau"], cap["re"]) == spec["centre"][:2] and cap.get("sza", spec["centre"][2]) == spec["centre"][2]
      if quad_all or centre:
         row["quad"] = {}
         for n in (nstrs if centre else [60]):
            status, q = ref.run("gfort_qp", cap["args"], n)
            row["quad"][str(n)] = {"status": status, "dp_minus_quad": coefficient_differences(dp_c[n], coefficients(cap, q)),
                                   "dp_minus_quad_outputs": {k: max_abs(cap["dp"][n][k], q[k]) for k in OUT}}
      rows.append(row)
   report["calls_detail"] = rows

   # Central node: reflectance error relative to production uncertainty, and Jacobians.
   direct = {(c["tau_index"], c["re_index"], c["sza_index"]): c for c in captures if c["kind"] == "direct"}
   s_index = int(np.flatnonzero(szas == spec["centre"][2])[0])
   a0 = int(np.flatnonzero(np.isclose(taus, spec["centre"][0]))[0])
   r0 = int(np.flatnonzero(np.isclose(res, spec["centre"][1]))[0])
   with Dataset(PRODUCTION_V21) as p:
      solar = list(np.asarray(p["solar_channel_id"][:]).astype(int))
      snr = float(np.asarray(p["snr"][:])[solar.index(spec["channel"])])
      vza_grid = np.asarray(p["satellite_zenith"][:], np.float64)
      raa_grid = np.asarray(p["relative_azimuth"][:], np.float64)
   noise = 1.0 / snr if snr > 0 else 0.0

   def r0v(cap, which):
      out = cap["sp"] if which == "sp" else cap["dp"][which]
      return coefficients(cap, out)["R_0v"]
   centre_cap = direct[(a0, r0, s_index)]
   values = {w: r0v(centre_cap, w) for w in ["sp"] + nstrs}
   reference_values = values[100]
   sigma = fm.sigma_solar(np.abs(reference_values), noise, production=True)
   valid = np.abs(reference_values) >= 0.01
   jv = int(np.argmin(np.abs(vza_grid - spec["geometry"][0])))
   jr = int(np.argmin(np.abs(raa_grid - spec["geometry"][1])))
   central = {"snr": snr, "noise_used": noise, "geometry_vza_raa": [float(vza_grid[jv]), float(raa_grid[jr])],
              "R_0v_at_geometry": {str(w): float(values[w][jv, jr]) for w in values},
              "max_sigma_sp_minus_dp100": float(np.max(np.abs(values["sp"] - reference_values)[valid] / sigma[valid])) if valid.any() else None,
              "max_sigma_dp60_minus_dp100": float(np.max(np.abs(values[60] - reference_values)[valid] / sigma[valid])) if valid.any() else None,
              "max_sigma_dp48_minus_dp100": float(np.max(np.abs(values[48] - reference_values)[valid] / sigma[valid])) if valid.any() else None,
              "max_abs_sp_minus_dp100": max_abs(values["sp"], reference_values),
              "max_abs_dp60_minus_dp100": max_abs(values[60], reference_values),
              "minimum_R_0v": {str(w): float(np.min(values[w])) for w in values}}
   # Finite-difference Jacobians at the centre (central where neighbours exist), as ORAC's bicubic scheme.
   jac = {}
   lt = np.log10(taus)
   for w in ["sp"] + nstrs:
      def g(a, r):
         return r0v(direct[(a, r, s_index)], w)
      am, ap = max(a0 - 1, 0), min(a0 + 1, len(taus) - 1)
      rm, rp = max(r0 - 1, 0), min(r0 + 1, len(res) - 1)
      jac[str(w)] = {"dR_dlog10tau": (g(ap, r0) - g(am, r0)) / (lt[ap] - lt[am]),
                     "dR_dre": (g(a0, rp) - g(a0, rm)) / (res[rp] - res[rm])}
   central["jacobian_max_abs_sp_minus_dp100"] = {k: max_abs(jac["sp"][k], jac["100"][k]) for k in jac["sp"]}
   central["jacobian_max_abs_dp60_minus_dp100"] = {k: max_abs(jac["60"][k], jac["100"][k]) for k in jac["sp"]}
   central["jacobian_max_abs_dp100"] = {k: float(np.max(np.abs(jac["100"][k]))) for k in jac["sp"]}
   central["jacobian_at_geometry"] = {w: {k: float(v[jv, jr]) for k, v in jac[w].items()} for w in jac}
   report["central"] = central

   # Generator LUT (version 90) against coefficients recomputed from captured DISORT_DP NSTR 60.
   with Dataset(lut_dp) as d:
      lut_r0v = np.asarray(d["R_0v"][:], np.float64)   # (raa, vza, sza, tau, re, chan) in file order
      dims = d["R_0v"].dimensions
   worst = 0.0
   for (a, r, s), cap in direct.items():
      mine = coefficients(cap, cap["dp"][60])["R_0v"]                        # (vza, raa)
      order = [dims.index(x) for x in ("solar_channels", "optical_depth", "effective_radius", "solar_zenith",
                                       "satellite_zenith", "relative_azimuth")]
      written = np.transpose(lut_r0v, order)[0, a, r, s]
      worst = max(worst, float(np.max(np.abs(written - np.float32(mine)))))
   report["v90_lut_vs_captured_dp60_R_0v_max_abs"] = worst
   return report


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("run", nargs="?", type=Path, default=HERE / "results" / "latest")
   parser.add_argument("--quad-all", action="store_true", help="REAL*16 reference for every call at NSTR 60")
   args = parser.parse_args()
   run = args.run.resolve()
   summary = {"run": str(run)}
   interface = readsav(str(run / "interface" / "dp_interface_test.sav"))
   summary["interface"] = {"info": interface["info"].decode() if isinstance(interface["info"], bytes) else str(interface["info"]),
                           "ok": int(interface["ok"]), "sp_bitwise": int(interface["sp_bitwise"]),
                           "sp_maxdiff": float(interface["sp_maxdiff"])}
   summary["cases"] = {}
   for name in CASES:
      summary["cases"][name] = analyse_case(name, run / f"case_{name}", args.quad_all)
   out = run / "analysis"
   out.mkdir(exist_ok=True)
   (out / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
   print(f"wrote {out / 'summary.json'}")
   return 0


if __name__ == "__main__":
   raise SystemExit(main())
