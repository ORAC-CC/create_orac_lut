"""Operator cancellation or amplification in the ORAC forward model (validation only).

For a candidate configuration C (default: the production baseline) and the
numerical reference R, at every node state, geometry, surface and channel:

   dF_total  = F(all operators of C) - F(all operators of R)
   dF_k      = F(operator k of C, all others of R) - F(R)      for each operator k
   ratio     = |dF_total| / sum_k |dF_k|

ratio ~ 1: the operator errors add (or one dominates); ratio << 1: they
cancel in the measurement.  "dominant" is the operator with the largest
|dF_k|; "R_0v-only" compares dF_total with dF for R_0v alone (the quantity a
black-surface, operator-level study looks at).

Operators substituted (as the forward model uses them): R_0v, T_00 (both its
uses, T_00 and T_vv), T_0d, T_dv, R_dd, R_dv, E_md.

   python cancellation.py [CANDIDATE]   -> results/disort_options/fm_cancellation_<CANDIDATE>.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import propagate as pg                 # noqa: E402
from atmosphere import terms           # noqa: E402

OPERATORS = ("R_0v", "T_00", "T_0d", "T_dv", "R_dd", "R_dv", "E_md")


def measurement(t, const, atm_by_vza, surf):
   per_vza = [pg.evaluate(t, const, atm, surf) for atm in atm_by_vza]
   out = {}
   for ch in per_vza[0]:
      out[ch] = np.stack([per_vza[iv][ch]["y"][:, :, iv, :] for iv in range(len(per_vza))], axis=2)
   return out


def main():
   candidate = sys.argv[1] if len(sys.argv) > 1 else "baseline"
   rows = []
   for case in [c for c in pg.CASES if not (c.endswith("_bs") or c.startswith("defect"))]:
      c_lut, r_lut = pg.lut(case, candidate), pg.lut(case, pg.REFERENCE)
      if c_lut is None or r_lut is None:
         continue
      platform, instrument, cloud = pg.PLATFORM[case]
      pg.CURRENT_CASE[0] = case
      const = pg.channel_constants(r_lut, instrument)
      planck = {ch: const[ch]["planck"] for ch in const if "planck" in const[ch]}
      tc, tr = pg.tables(c_lut), pg.tables(r_lut)
      axes = tr["axes"]
      nodes = np.array([k == "node" for _, _, k in pg.states(axes)])
      atm_by_vza = [terms(platform, instrument, pg.CASES[case][3], planck, cloud, float(v)) for v in axes["satellite_zenith"]]
      for surf in ("black", "dark", "vegetation", "snow"):
         f_ref = measurement(tr, const, atm_by_vza, surf)
         f_all = measurement(tc, const, atm_by_vza, surf)
         single = {}
         for op in OPERATORS:
            hybrid = dict(tr)
            hybrid[op] = tc[op]
            single[op] = measurement(hybrid, const, atm_by_vza, surf)
         for ch in f_ref:
            total = (f_all[ch] - f_ref[ch])[nodes]
            parts = np.stack([(single[op][ch] - f_ref[ch])[nodes] for op in OPERATORS])
            finite = np.isfinite(total) & np.all(np.isfinite(parts), axis=0)     # reference breakdown excluded
            total = np.where(finite, total, 0.0)
            parts = np.where(finite, parts, 0.0)
            mag = np.abs(parts).sum(axis=0)
            if np.abs(total).max() == 0.0:
               rows.append({"case": case, "candidate": candidate, "surface": surf, "channel": ch,
                            "max_abs_dF_total": 0.0, "max_abs_dF_R0v_only": float(np.abs(parts[0]).max()), "points": 0})
               continue
            significant = np.abs(total) >= 0.1 * np.abs(total).max()      # where the effect matters
            ratio = np.abs(total[significant]) / np.maximum(mag[significant], 1e-30)
            dominant = np.argmax(np.abs(parts), axis=0)[significant]
            r0v = np.abs(parts[0][significant])
            rows.append({"case": case, "candidate": candidate, "surface": surf, "channel": ch,
                         "max_abs_dF_total": float(np.abs(total).max()),
                         "max_abs_dF_R0v_only": float(np.abs(parts[0]).max()),
                         "points": int(significant.sum()),
                         "median_ratio_total_to_sum": float(np.median(ratio)),
                         "min_ratio_total_to_sum": float(ratio.min()),
                         "median_total_over_R0v_only": float(np.median(np.abs(total[significant]) / np.maximum(r0v, 1e-30))),
                         **{f"share_dominant_{op}": float((dominant == k).mean()) for k, op in enumerate(OPERATORS)},
                         **{f"max_abs_dF_{op}": float(np.abs(parts[k]).max()) for k, op in enumerate(OPERATORS)}})
         print(case, surf, "done", flush=True)
   out = pg.RESULTS / f"fm_cancellation_{candidate}.csv"
   with open(out, "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=max((list(r) for r in rows), key=len), restval="")
      writer.writeheader()
      writer.writerows(rows)
   print("written", out)


if __name__ == "__main__":
   main()
