"""Operator-level comparison of the DISORT-option experiments (validation only).

For every case and experiment (run_variant.py products) this compares the LUT
with the production baseline and with the converged reference
(double-precision kernel, NSTR = 236) at every compact-grid node.

Operator level only.  R_0v is the LUT's bidirectional reflectance of the
atmosphere-plus-cloud column over a black surface; it is ONE input of the
ORAC forward model, not the TOA measurement (that is propagate.py, through
orac_fm.py).  E_md is the cloud emissivity operator; a change dE changes
the TOA brightness temperature by at most dE lambda T^2 / c2 (Planck
linearisation at the instrument file's reference brightness temperature).
These operator-level numbers are retained as diagnostics.

Measurement uncertainty (the instrument files that the LUTs carry to ORAC):
   solar, MODIS        sigma_L / L = 1 / snr
   solar, SLSTR        sigma_L / L = rua (rub = ruc = 0), and, separately, the
                       legacy SAD noise-equivalent reflectance oldnefr (absolute)
   thermal             NEdT (K)
Normalised difference r = |dL| / sigma_L.  Classification (stated criterion):
   r < 0.1 negligible; 0.1 <= r < 0.5 potentially material but below the
   uncertainty; 0.5 <= r < 2 comparable; r >= 2 clearly larger.

   python analyse.py            -> results/*.csv and printed summary
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
from run_variant import CASES, EXPERIMENTS, INPUTS, TMP   # noqa: E402
from oraclut.io.lut import read_lut                       # noqa: E402
from oraclut.idl_mirror import load_inststr               # noqa: E402

C2 = 1.4387769e-2           # m K
RESULTS = ROOT / "validation" / "v24" / "results" / "disort_options"
REFERENCE = "double_nstr236"
OPERATORS = ("R_0d", "R_dd", "R_dv", "T_dd", "T_dv", "T_0d", "T_00")


def classify(r):
   return "negligible" if r < 0.1 else "below" if r < 0.5 else "comparable" if r < 2.0 else "larger"


def load(case, experiment):
   path = TMP / "runs" / case / experiment
   files = sorted(path.glob("*.nc"))
   if not files or not (path / "record.json").exists():
      return None, None
   return read_lut(files[0]), json.loads((path / "record.json").read_text())


def channel_info(case, lut):
   instfile = CASES[case][0]
   ids = [int(c) for c in np.asarray(lut.variables["channel_id"])]
   # nadir channels (the oblique-view channels of a dual-view LUT are copies of them)
   inst = load_inststr(INPUTS / "inst" / instfile, requestedchannelid=sorted(set(ids) & set(CASES[case][3])))
   info = {}
   for i, ch in enumerate(inst.channelid):
      ch = int(ch)
      info[ch] = {"wavelength_um": float(np.asarray(lut.variables["channel_wl_abs" if "channel_wl_abs" in lut.variables else "wavelength"])[ids.index(ch)])
                  if ("channel_wl_abs" in lut.variables or "wavelength" in lut.variables) else None,
                  "sigma_rel": (float(inst.rua[i]) if int(getattr(inst, "view", 0)) > 0 and float(inst.rua[i]) > 0
                                else (1.0 / float(inst.snr[i]) if float(inst.snr[i]) > 0 else None)),
                  "sigma_rel_source": "rua" if int(getattr(inst, "view", 0)) > 0 and float(inst.rua[i]) > 0 else "1/snr",
                  "nefr": float(inst.oldnefr[i]) if float(inst.oldnefr[i]) > 0 else None,
                  "nedt": float(inst.nedt[i]) if float(inst.nedt[i]) > 0 else None,
                  "refbt": float(inst.refbt[i]) if float(inst.refbt[i]) > 0 else None}
   return info


def geometry_theta(lut, dims, shape):
   sza, vza, raa = (np.deg2rad(np.asarray(lut.variables[d], dtype=np.float64))
                    for d in ("solar_zenith", "satellite_zenith", "relative_azimuth"))
   cos_theta = (-np.cos(sza)[:, None, None] * np.cos(vza)[None, :, None]
                - np.sin(sza)[:, None, None] * np.sin(vza)[None, :, None] * np.cos(raa)[None, None, :])
   theta = np.rad2deg(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
   order = [dims.index(d) for d in ("solar_zenith", "satellite_zenith", "relative_azimuth")]
   return np.broadcast_to(np.moveaxis(theta.reshape(theta.shape + (1,) * (len(shape) - 3)), [0, 1, 2], order), shape)


def compare(case, experiment, against, lut, base, info):
   """Rows (one per channel) of TOA-radiance differences of lut from base."""

   rows = []
   dims = base.variable_dimensions["R_0v"]
   axis = dims.index("solar_channels")
   solar = [int(c) for c in np.asarray(base.variables["solar_channel_id"])]
   a = np.asarray(base.variables["R_0v"], dtype=np.float64)
   b = np.asarray(lut.variables["R_0v"], dtype=np.float64)
   theta = geometry_theta(base, dims, a.shape)
   coords = {d: np.asarray(base.variables[d], dtype=np.float64) for d in dims if d != "solar_channels"}
   for c, ch in enumerate(solar):
      if ch not in info:
         continue
      x, y, t = (np.take(v, c, axis=axis) for v in (a, b, theta))
      d = np.abs(y - x)
      rel = d / np.maximum(np.abs(x), 1e-6)
      d = np.where(np.isfinite(d), d, 0.0)                 # reference breakdown point (see the report)
      rel = np.where(np.isfinite(rel), rel, 0.0)
      s = info[ch]["sigma_rel"]
      r = rel / s if s else np.full_like(rel, np.nan)
      where = np.unravel_index(int(np.nanargmax(r)) if s else int(np.argmax(d)), d.shape)
      sub = [dd for dd in dims if dd != "solar_channels"]
      row = {"case": case, "experiment": experiment, "against": against, "quantity": "R_0v operator", "channel": ch,
             "sigma_definition": f"{info[ch]['sigma_rel_source']} = {s:.4g}" if s else "none",
             "max_abs": float(d.max()), "max_rel": float(rel.max()), "rms_rel": float(np.sqrt(np.mean(rel**2))),
             "max_dBT_K": None, "max_r": float(np.nanmax(r)) if s else None, "p99_r": float(np.nanpercentile(r, 99)) if s else None,
             "max_r_theta_lt170": float(np.nanmax(r[t < 170.0])) if s else None,
             "max_r_nefr": float(d.max() / info[ch]["nefr"]) if info[ch]["nefr"] else None,
             "class": classify(float(np.nanmax(r))) if s else None,
             "worst_at": {dd: float(coords[dd][where[k]]) for k, dd in enumerate(sub)} | {"theta": float(t[where])}}
      row["worst_at"] = json.dumps(row["worst_at"])
      rows.append(row)
   if "E_md" in base.variables:
      dims = base.variable_dimensions["E_md"]
      axis = dims.index("thermal_channels")
      thermal = [int(c) for c in np.asarray(base.variables["thermal_channel_id"])]
      a = np.asarray(base.variables["E_md"], dtype=np.float64)
      b = np.asarray(lut.variables["E_md"], dtype=np.float64)
      for c, ch in enumerate(thermal):
         if ch not in info:
            continue
         d = np.abs(np.take(b, c, axis=axis) - np.take(a, c, axis=axis))
         wl = float(np.asarray(base.variables["thermal_channel_wl"] if "thermal_channel_wl" in base.variables else [np.nan])[0])
         rows.append({"case": case, "experiment": experiment, "against": against, "quantity": "E_md operator", "channel": ch,
                      "sigma_definition": f"NEdT = {info[ch]['nedt']} K", "max_abs": float(d.max()),
                      "max_rel": float((d / np.maximum(np.abs(np.take(a, c, axis=axis)), 1e-6)).max()),
                      "rms_rel": None, "max_dBT_K": None, "max_r": None, "p99_r": None, "max_r_theta_lt170": None,
                      "max_r_nefr": None, "class": None, "worst_at": ""})
   return rows


def thermal_bt(rows, info, wavelengths):
   for row in rows:
      if row["quantity"].startswith("E_md"):
         ch = row["channel"]
         k = wavelengths[ch] * 1e-6 * info[ch]["refbt"] ** 2 / C2      # K per unit relative radiance
         row["max_dBT_K"] = row["max_abs"] * k
         row["max_r"] = row["max_dBT_K"] / info[ch]["nedt"]
         row["class"] = classify(row["max_r"])
   return rows


def operators(case, experiment, against, lut, base):
   out = []
   for name in OPERATORS:
      if name in base.variables:
         a = np.asarray(base.variables[name], dtype=np.float64)
         b = np.asarray(lut.variables[name], dtype=np.float64)
         out.append({"case": case, "experiment": experiment, "against": against, "operator": name,
                     "max_abs": float(np.nanmax(np.abs(b - a))),
                     "max_rel": float(np.nanmax(np.abs(b - a) / np.maximum(np.abs(a), 1e-3)))})
   return out


def wavelengths_of(case):
   from oraclut.idl_mirror import load_srfstrarr
   instfile, _, _, channels, _ = CASES[case]
   inst = load_inststr(INPUTS / "inst" / instfile, requestedchannelid=channels)
   srf, _ = load_srfstrarr(inst, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
   return {int(ch): float(s.wvl_centre) for ch, s in zip(inst.channelid, srf)}


def main():
   RESULTS.mkdir(parents=True, exist_ok=True)
   radiance, operator_rows, runtime = [], [], []
   for case in [c for c in CASES if not (c.endswith("_bs") or c.startswith("defect"))]:
      base, base_record = load(case, "baseline")
      if base is None:
         continue
      info = channel_info(case, base)
      wl = wavelengths_of(case)
      for ch in info:
         info[ch]["wavelength_um"] = wl[ch]
      ref, _ = load(case, REFERENCE)
      for experiment in EXPERIMENTS:
         lut, record = load(case, experiment)
         if lut is None:
            continue
         runtime.append({"case": case, "experiment": experiment, "nstreams": record["nstreams"], "kernel": record["kernel"],
                         "accur": record["accur"], "deltam": record["deltam"], "corint": record["corint"],
                         "moments": record["moments"], "disort_cpu_s": record["disort_cpu_total_s"],
                         "ratio_to_baseline": record["disort_cpu_total_s"] / base_record["disort_cpu_total_s"],
                         "direct_calls": record["disort_calls"]["direct"], "diffuse_calls": record["disort_calls"]["diffuse"],
                         "emission_calls": record["disort_calls"]["emission"]})
         if experiment != "baseline":
            radiance += thermal_bt(compare(case, experiment, "baseline", lut, base, info), info, wl)
            operator_rows += operators(case, experiment, "baseline", lut, base)
         if ref is not None and experiment != REFERENCE:
            radiance += thermal_bt(compare(case, experiment, REFERENCE, lut, ref, info), info, wl)
            operator_rows += operators(case, experiment, REFERENCE, lut, ref)
      with open(RESULTS / f"channels_{case}.json", "w") as handle:
         json.dump(info, handle, indent=2)
   for name, rows in (("operator_r0v_emd.csv", radiance), ("lut_operators.csv", operator_rows), ("runtime_single.csv", runtime)):
      if rows:
         with open(RESULTS / name, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
   # printed summary: worst normalised difference over channels per case/experiment
   for case in [c for c in CASES if not (c.endswith("_bs") or c.startswith("defect"))]:
      for against in ("baseline", REFERENCE):
         print(f"\n{case}: worst |dL|/sigma over channels, against {against}")
         for experiment in EXPERIMENTS:
            rows = [r for r in radiance if r["case"] == case and r["experiment"] == experiment and r["against"] == against]
            if not rows:
               continue
            solar = [r for r in rows if r["quantity"].startswith("R_0v") and r["max_r"] is not None]
            thermal = [r for r in rows if r["quantity"].startswith("E_md")]
            ws = max(solar, key=lambda r: r["max_r"]) if solar else None
            wt = max(thermal, key=lambda r: r["max_r"]) if thermal else None
            rt = next((r for r in runtime if r["case"] == case and r["experiment"] == experiment), None)
            print(f"  {experiment:24s} solar {ws['max_r']:9.3f} (ch {ws['channel']}, rel {ws['max_rel']:.1e}, "
                  f"off-glory {ws['max_r_theta_lt170']:.3f}) {ws['class']:10s} | thermal "
                  + (f"{wt['max_dBT_K']:.4f} K = {wt['max_r']:.3f} NEdT (ch {wt['channel']})" if wt else "-")
                  + (f" | DISORT CPU x{rt['ratio_to_baseline']:.2f}" if rt else ""))


if __name__ == "__main__":
   main()
