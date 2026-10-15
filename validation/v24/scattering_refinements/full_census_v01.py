"""Isolated single-precision DISORT call failures in the full V23 and V24 Terra MODIS products (read-only).

   python validation/v24/scattering_refinements/full_census.py   -> results/full_census.json

The established detector of validation/v24/disort_options/corrupted_calls.py
(isolated_score, local_max, SCORE_MIN, ABS_MIN, imported unchanged) is applied
to every direct-beam call (solar channel, optical depth, effective radius,
solar zenith) of R_0v (all solar zeniths, as in corrupted_calls.py) and of
R_0d and T_0d (interior solar zeniths only: at the 75 deg edge the linear
extrapolation flags smooth curvature), for the V23 (staging/lut_v23) and V24
(staging/lut_v24) liquid-water products, which share the Grid B nodes.  Calls
flagged in one product only identify which product holds the failure at a
node where the two differ.  Partial failures (some views only) can escape this
median-over-views detector, so the census is a lower bound; the largest
remaining differences are attributed call by call in full_inspection.py.
"""

import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "validation" / "v24" / "disort_options"))
from corrupted_calls import ABS_MIN, SCORE_MIN, isolated_score, local_max   # noqa: E402

CALL_DIMS = ("solar_channels", "optical_depth", "effective_radius", "solar_zenith")


def census(path):
   cells = []
   with Dataset(path) as d:
      d.set_auto_maskandscale(False)
      for name in ("R_0v", "R_0d", "T_0d"):
         dims = list(d[name].dimensions)
         order = [dims.index(k) for k in CALL_DIMS] + [i for i, k in enumerate(dims) if k not in CALL_DIMS]
         x = np.transpose(np.asarray(d[name][:], np.float64), order)
         while x.ndim < 6:
            x = x[..., None]
         score, anomaly = isolated_score(x)
         hit = (score > SCORE_MIN) & (anomaly > ABS_MIN) & local_max(score)
         if name != "R_0v":
            hit[..., [0, -1]] = False
         for i in np.argwhere(hit):
            cells.append({"variable": name, "channel": int(d["solar_channel_id"][i[0]]), "tau": float(d["optical_depth"][i[1]]),
                          "re": float(d["effective_radius"][i[2]]), "sza": float(d["solar_zenith"][i[3]]),
                          "score": round(float(score[tuple(i)]), 3), "median_abs_anomaly": float(anomaly[tuple(i)])})
   return cells


def calls(cells):
   return {(c["channel"], c["tau"], c["re"], c["sza"]) for c in cells}


out = {}
for stem in ("liquid-water_a01_pold", "liquid-water_a01_p240"):
   v23 = census(ROOT / "staging" / "lut_v23" / f"terra_modis_m_{stem}_v23.nc")
   v24 = census(ROOT / "staging" / "lut_v24" / f"terra_modis_m_{stem}_v24.nc")
   a, b = calls(v23), calls(v24)
   out[stem] = {"v23_cells": v23, "v24_cells": v24, "v23_calls": len(a), "v24_calls": len(b),
                "v23_only": sorted(a - b), "v24_only": sorted(b - a), "both": sorted(a & b)}
   print(f"== {stem}: flagged calls V23 {len(a)}, V24 {len(b)}; V23 only {len(a - b)}, V24 only {len(b - a)}, both {len(a & b)}")
   for label, cells in (("V23", v23), ("V24", v24)):
      for c in cells:
         print(f"   {label} {c['variable']:4s} ch {c['channel']:2d} tau {c['tau']:9.4g} r_e {c['re']:4g} SZA {c['sza']:6.2f}  "
               f"score {c['score']:.3g}  median anomaly {c['median_abs_anomaly']:.3g}")
(HERE / "results" / "full_census.json").write_text(json.dumps(out, indent=1, default=list) + "\n")
