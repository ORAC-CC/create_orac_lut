"""Census of corrupted direct-beam DISORT calls in ORAC cloud LUT products (validation only).

R_0v (= pi L / F0, including cos(sza)) is computed by one direct-beam DISORT
call per (channel, optical depth, radius, solar zenith) cell; the call fills
every view zenith and relative azimuth.  Two kinds of fault are counted:

  isolated   the cell is inconsistent with its two solar-zenith neighbours
             (same channel, tau, re): score = median over (vza, raa) of
             |R_0v - neighbour prediction| / (median |prediction| + 1e-3) above
             SCORE_MIN, with the median absolute anomaly above ABS_MIN, and the
             score a local maximum in solar zenith (the neighbours of a spike
             inherit about half of it).  The prediction is the neighbours' mean
             (linear extrapolation at the solar-zenith edges).  SCORE_MIN = 0.2:
             the reproduced faults score 0.25-2.8, smooth large-radius ringing
             at most ~0.13.  Reproductions
             (run_variant.py cases defect_*) show these are single-precision
             DISORT breakdowns: double precision, NSTR 56, 68 and 100 agree with
             each other and differ from the production value only in that cell.
  range      any R_0v < CATASTROPHIC_MIN or > CATASTROPHIC_MAX ("catastrophic"),
             or minimum in [CATASTROPHIC_MIN, -0.01) ("mild negative"); the
             mild class is dominated by smooth negative ringing in strongly
             absorbing bands (e.g. Baum GHM, 4-4.5 um, large radius), which is
             not an isolated-call fault (run_variant.py case ringing_modis_ghm24).

SLSTR oblique-view channels whose R_0v equals a nadir channel's exactly (the
same DISORT computation) are left out of the counts.  Read-only.

   python corrupted_calls.py [PRODUCT.nc ...]
      -> results/disort_options/corrupted_calls.csv          (one row per product)
         results/disort_options/corrupted_calls_isolated.csv (one row per candidate cell, score > LISTED_MIN)
   (default: every *_v23.nc and *_v24.nc cloud product in the permanent LUT directory)
"""

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.io.lut import read_lut   # noqa: E402

LUT_DIR = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
RESULTS = ROOT / "validation" / "v24" / "results" / "disort_options"
CALL_DIMS = ("solar_channels", "optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")
SCORE_MIN = 0.2
LISTED_MIN = 0.15                       # cells listed in the detail table (sensitivity of SCORE_MIN)
ABS_MIN = 0.002
CATASTROPHIC_MIN = -0.1
CATASTROPHIC_MAX = 1.5


def isolated_score(x):
   """x: (..., sza, vza, raa) -> (relative score, median absolute anomaly) per (..., sza) cell."""
   p = np.empty_like(x)
   p[..., 1:-1, :, :] = 0.5 * (x[..., :-2, :, :] + x[..., 2:, :, :])
   p[..., 0, :, :] = 2.0 * x[..., 1, :, :] - x[..., 2, :, :]
   p[..., -1, :, :] = 2.0 * x[..., -2, :, :] - x[..., -3, :, :]
   anomaly = np.median(np.abs(x - p).reshape(x.shape[:-2] + (-1,)), axis=-1)
   scale = np.median(np.abs(p).reshape(x.shape[:-2] + (-1,)), axis=-1)
   return anomaly / (scale + 1e-3), anomaly


def local_max(score):
   """True where the score is not exceeded by either solar-zenith neighbour (last axis)."""
   left = np.concatenate([np.full(score.shape[:-1] + (1,), -np.inf), score[..., :-1]], axis=-1)
   right = np.concatenate([score[..., 1:], np.full(score.shape[:-1] + (1,), -np.inf)], axis=-1)
   return (score >= left) & (score >= right)


def main():
   paths = [Path(p) for p in sys.argv[1:]] or sorted(p for p in LUT_DIR.glob("*_v2[34].nc")
                                                     if "liquid-water" in p.name or "water-ice" in p.name)
   rows, cells = [], []
   for path in paths:
      lut = read_lut(path)
      dims = lut.variable_dimensions["R_0v"]
      x = np.moveaxis(np.asarray(lut.variables["R_0v"], dtype=np.float64), [dims.index(d) for d in CALL_DIMS], range(6))
      ch = np.asarray(lut.variables["solar_channel_id"]).astype(int)
      tau = np.asarray(lut.variables["optical_depth"])
      re = np.asarray(lut.variables["effective_radius"])
      sza = np.asarray(lut.variables["solar_zenith"])
      keep = [j for j in range(x.shape[0]) if not any(np.array_equal(x[i], x[j]) for i in range(j))]
      x = x[keep]
      ch = ch[keep]
      score, anomaly = isolated_score(x)
      candidate = (score > LISTED_MIN) & (anomaly > ABS_MIN) & local_max(score)
      isolated = candidate & (score > SCORE_MIN)
      lo, hi = x.min(axis=(4, 5)), x.max(axis=(4, 5))
      catastrophic = (lo < CATASTROPHIC_MIN) | (hi > CATASTROPHIC_MAX)
      mild = ~catastrophic & (lo < -0.01)
      for i in np.argwhere(candidate):
         i = tuple(i)
         cells.append({"product": path.name, "channel": int(ch[i[0]]), "tau": float(tau[i[1]]), "re": float(re[i[2]]),
                       "sza": float(sza[i[3]]), "score": float(score[i]), "median_abs_anomaly": float(anomaly[i]),
                       "min_r0v": float(lo[i]), "max_r0v": float(hi[i]), "catastrophic_range": bool(catastrophic[i]),
                       "isolated": bool(isolated[i])})
      total = int(np.prod(score.shape))
      rows.append({"product": path.name, "channels_counted": int(len(ch)), "direct_cells": total,
                   "isolated": int(isolated.sum()), "isolated_share": float(isolated.sum()) / total,
                   "isolated_channels": " ".join(map(str, sorted(set(ch[np.nonzero(isolated)[0]].tolist())))),
                   "isolated_distinct_re_sza": int(len({(c[0], c[2], c[3]) for c in map(tuple, np.argwhere(isolated))})),
                   "catastrophic_range": int(catastrophic.sum()), "catastrophic_not_isolated": int((catastrophic & ~isolated).sum()),
                   "mild_negative": int(mild.sum()), "min_r0v": float(lo.min()), "max_r0v": float(hi.max())})
      print(rows[-1], flush=True)
   RESULTS.mkdir(parents=True, exist_ok=True)
   for name, table in (("corrupted_calls.csv", rows), ("corrupted_calls_isolated.csv", cells)):
      if table:
         with open(RESULTS / name, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)


if __name__ == "__main__":
   main()
