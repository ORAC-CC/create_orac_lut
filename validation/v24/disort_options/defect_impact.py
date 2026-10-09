"""Effect of the reproduced corrupted direct-beam DISORT calls on the ORAC measurement (validation only).

For each reproduction case (run_variant.py defect_*), compares the production
kernel (baseline, NSTR 60 single precision) with double precision at NSTR 60
in every direct-beam output (R_0v, R_0d, T_0d, T_00) of every
(tau, re, sza) cell.  For a black surface, ORAC equation form 3 reduces to
Ref = Tac_0v R_0v (check_fm.py), and the production solar uncertainty is
relative (sigma = s Ref, s^2 = noise^2 + 0.01^2 + 0.02^2, get_measurements.F90),
so r = |dR_0v| / (s R_0v,double) exactly, independent of the clear-sky
transmittance.  Over a reflecting surface the R_0d/T_0d/T_00 errors of the
same call add the surface terms.

   python defect_impact.py  ->  results/disort_options/defect_impact.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import propagate as pg                                 # noqa: E402

CASES = {"defect_slstr240": "slstr", "defect_modis_agg": "modis", "defect_modis_ghm5": "modis", "defect_modis_ghm17": "modis"}
DIRECT = ("R_0v", "R_0d", "T_0d", "T_00")


def arrays(l):
   out = {}
   for name in DIRECT:
      dims = l.variable_dimensions[name]
      order = [dims.index(d) for d in ("solar_channels", "optical_depth", "effective_radius", "solar_zenith") if d in dims]
      rest = [i for i in range(len(dims)) if i not in order]
      out[name] = np.moveaxis(np.asarray(l.variables[name], dtype=np.float64), order + rest, range(len(dims)))
   return out


def main():
   rows = []
   for case, instrument in CASES.items():
      base, dbl = pg.lut(case, "baseline"), pg.lut(case, "double")
      ch_ids = np.asarray(base.variables["solar_channel_id"]).astype(int)
      if instrument == "slstr":
         noise = np.asarray(base.variables["rua"], dtype=np.float64)
      else:
         snr = np.asarray(base.variables["snr"], dtype=np.float64)
         noise = np.where(snr > 0, 1.0 / np.where(snr > 0, snr, 1.0), np.nan)   # snr = 0: no ORAC solar sigma (r = nan)
      s = np.sqrt(noise ** 2 + 0.01 ** 2 + 0.02 ** 2)
      a, b = arrays(base), arrays(dbl)
      ax = {n: np.asarray(base.variables[n]) for n in ("optical_depth", "effective_radius", "solar_zenith")}
      for idx in np.ndindex(a["R_0v"].shape[:4]):
         d = {n: float(np.abs(a[n][idx] - b[n][idx]).max()) for n in DIRECT}
         if d["R_0v"] < 1e-3:
            continue
         ref = b["R_0v"][idx]
         r = np.abs(a["R_0v"][idx] - ref) / (s[idx[0]] * np.maximum(ref, 1e-6))
         rows.append({"case": case, "channel": int(ch_ids[idx[0]]), "tau": float(ax["optical_depth"][idx[1]]),
                      "re": float(ax["effective_radius"][idx[2]]), "sza": float(ax["solar_zenith"][idx[3]]),
                      "relative_sigma": float(s[idx[0]]), "max_abs_dR_0v": d["R_0v"], "max_abs_dR_0d": d["R_0d"],
                      "max_abs_dT_0d": d["T_0d"], "max_abs_dT_00": d["T_00"],
                      "share_views_r_gt_2": float((r > 2).mean()), "median_r_black": float(np.median(r)),
                      "max_r_black": float(r.max()), "baseline_range": f"{a['R_0v'][idx].min():.3f} {a['R_0v'][idx].max():.3f}",
                      "double_range": f"{ref.min():.3f} {ref.max():.3f}"})
         print(rows[-1], flush=True)
   with open(pg.RESULTS / "defect_impact.csv", "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)


if __name__ == "__main__":
   main()
