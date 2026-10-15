"""Attribute the largest V24 - V23 differences of the full liquid-water products (read-only).

   python validation/v24/scattering_refinements/full_inspection.py   (after full_census.py) -> results/full_inspection.json

For pold and p240 (staging/lut_v23, staging/lut_v24):
* the largest |V24 - V23| of R_0v, R_0d and T_0d once the direct-beam calls
  flagged in either product by full_census.py are excluded;
* the number of R_0v values below -0.01 in each product and the calls they
  belong to;
* energy conservation of every flagged call, R_0d + T_0d + T_00 at the call
  and at its solar-zenith neighbours (a non-absorbing call must stay near 1);
* the TOP_CALLS unflagged direct-beam calls with the largest |V24 - V23| in
  R_0v, each with its roughness in either product: the largest deviation over
  the views from the mean of the two solar-zenith neighbours.  A call much
  rougher in one product than in the other (and than its neighbours) is a
  partial DISORT failure of that product, which the median-over-views
  detector of full_census.py can miss; similar roughness means a change of
  the optics.
"""

from collections import Counter
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

TOP_CALLS = 8
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
census = json.loads((HERE / "results" / "full_census.json").read_text())
out = {}
for stem in ("liquid-water_a01_pold", "liquid-water_a01_p240"):
   a = Dataset(ROOT / "staging" / "lut_v23" / f"terra_modis_m_{stem}_v23.nc")
   b = Dataset(ROOT / "staging" / "lut_v24" / f"terra_modis_m_{stem}_v24.nc")
   for d in (a, b):
      d.set_auto_maskandscale(False)
   ch = [int(v) for v in a["solar_channel_id"][:]]
   tau = np.asarray(a["optical_depth"][:], np.float64)
   re = np.asarray(a["effective_radius"][:], np.float64)
   sza = np.asarray(a["solar_zenith"][:], np.float64)
   flagged = sorted({tuple(c) for key in ("v23_only", "v24_only", "both") for c in census[stem][key]})
   index = lambda c, t, r, s: (ch.index(c), int(np.argmin(abs(tau - t))), int(np.argmin(abs(re - r))), int(np.argmin(abs(sza - s))))
   record = {"flagged_calls": flagged, "excluding_flagged": {}, "energy_of_flagged_calls": []}
   for name in ("R_0v", "R_0d", "T_0d"):
      dims = list(a[name].dimensions)
      x, y = np.asarray(a[name][:], np.float64), np.asarray(b[name][:], np.float64)
      mask = np.zeros(x.shape, bool)
      for call in flagged:
         c, t, r, s = index(*call)
         idx = [slice(None)] * x.ndim
         for dim, i in (("solar_channels", c), ("optical_depth", t), ("effective_radius", r), ("solar_zenith", s)):
            idx[dims.index(dim)] = i
         mask[tuple(idx)] = True
      diff = np.where(mask, 0.0, np.abs(y - x))
      i = np.unravel_index(np.argmax(diff), diff.shape)
      at = {dim: (ch[j] if dim == "solar_channels" else float(a[dim][j])) for dim, j in zip(dims, i)}
      record["excluding_flagged"][name] = {"max": float(diff.max()), "at": at, "v23": float(x[i]), "v24": float(y[i])}
   # the largest unflagged R_0v changes, call by call, with each product's roughness
   dims = list(a["R_0v"].dimensions)
   order = [dims.index(k) for k in ("solar_channels", "optical_depth", "effective_radius", "solar_zenith")] + \
           [i for i, k in enumerate(dims) if k not in ("solar_channels", "optical_depth", "effective_radius", "solar_zenith")]
   xa = np.transpose(np.asarray(a["R_0v"][:], np.float64), order)
   xb = np.transpose(np.asarray(b["R_0v"][:], np.float64), order)
   xa, xb = xa.reshape(xa.shape[:4] + (-1,)), xb.reshape(xb.shape[:4] + (-1,))
   change = np.max(np.abs(xb - xa), axis=-1)
   for call in flagged:
      change[index(*call)] = 0.0
   rough = lambda x, j: float(np.max(np.abs(x[..., j, :] - 0.5 * (x[..., j - 1, :] + x[..., j + 1, :])))) if 0 < j < sza.size - 1 else float("nan")
   record["top_unflagged_R_0v_calls"] = []
   for flat in np.argsort(change, axis=None)[::-1][:TOP_CALLS]:
      c, t, r, j = np.unravel_index(flat, change.shape)
      record["top_unflagged_R_0v_calls"].append({
         "channel": ch[c], "tau": float(tau[t]), "re": float(re[r]), "sza": float(sza[j]), "max_change": float(change[c, t, r, j]),
         "roughness_v23": rough(xa[c, t, r], j), "roughness_v24": rough(xb[c, t, r], j),
         "roughness_neighbour_calls_v23": [rough(xa[c, t, r], k) for k in (j - 1, j + 1)],
         "roughness_neighbour_calls_v24": [rough(xb[c, t, r], k) for k in (j - 1, j + 1)]})
   for label, d in (("v23", a), ("v24", b)):
      r0v = np.asarray(d["R_0v"][:])
      neg = np.argwhere(r0v < -0.01)
      record[f"{label}_R_0v_below_-0.01"] = int(len(neg))
      record[f"{label}_negative_calls"] = [[*key, n] for key, n in Counter(
         (ch[k[5]], float(tau[k[3]]), float(re[k[4]]), float(sza[k[2]])) for k in neg).most_common()]
   for call in flagged:
      c, t, r, s = index(*call)
      rows = []
      for j in range(max(s - 1, 0), min(s + 2, sza.size)):
         rows.append({"sza": float(sza[j]), **{label: float(sum(d[n][j, t, r, c] for n in ("R_0d", "T_0d", "T_00")))
                                              for label, d in (("v23", a), ("v24", b))}})
      record["energy_of_flagged_calls"].append({"call": list(call), "R_0d+T_0d+T_00": rows})
   out[stem] = record
   print(f"== {stem}: {len(flagged)} flagged calls excluded")
   for name, e in record["excluding_flagged"].items():
      print(f"   {name}: max {e['max']:.3g} at {e['at']} (V23 {e['v23']:.4f}, V24 {e['v24']:.4f})")
   print(f"   R_0v < -0.01: V23 {record['v23_R_0v_below_-0.01']} {record['v23_negative_calls']}, "
         f"V24 {record['v24_R_0v_below_-0.01']} {record['v24_negative_calls']}")
   for e in record["top_unflagged_R_0v_calls"]:
      print(f"   top ch {e['channel']:2d} tau {e['tau']:8.4g} r_e {e['re']:4g} SZA {e['sza']:6.2f}: change {e['max_change']:.3g}; "
            f"roughness V23 {e['roughness_v23']:.3g} (neighbours {e['roughness_neighbour_calls_v23'][0]:.2g}, "
            f"{e['roughness_neighbour_calls_v23'][1]:.2g}), V24 {e['roughness_v24']:.3g} "
            f"(neighbours {e['roughness_neighbour_calls_v24'][0]:.2g}, {e['roughness_neighbour_calls_v24'][1]:.2g})")
   for e in record["energy_of_flagged_calls"]:
      print(f"   energy {e['call']}: " + ", ".join(f"SZA {row['sza']:g}: {row['v23']:.4f}/{row['v24']:.4f}" for row in e["R_0d+T_0d+T_00"]))
(HERE / "results" / "full_inspection.json").write_text(json.dumps(out, indent=1) + "\n")
