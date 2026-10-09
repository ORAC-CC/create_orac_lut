"""Baseline reproduction check (validation only).

The harness baseline (run_variant.py CASE baseline) uses the production
generator on a compact grid whose every coordinate is a Grid B node.  This
compares each of its LUT variables with the finished V24 production product
of the same platform/model at the matching nodes, channels and views.

   python check_baseline.py CASE V24_PRODUCT.nc
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.io.lut import read_lut   # noqa: E402

COORDS = ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")


def main():
   case, product = sys.argv[1], Path(sys.argv[2])
   base = ROOT / "validation" / "tmp" / "disort_options" / "runs" / case / "baseline"
   small = read_lut(next(base.glob("*.nc")))
   full = read_lut(product)
   index = {}
   for name in COORDS:
      s = np.asarray(small.variables[name], dtype=np.float64)
      f = np.asarray(full.variables[name], dtype=np.float64)
      index[name] = np.array([int(np.argmin(np.abs(f - v))) for v in s])
      assert np.allclose(f[index[name]], s, rtol=1e-6, atol=1e-6), name
   cs = [int(c) for c in np.asarray(small.variables["channel_id"])]
   cf = [int(c) for c in np.asarray(full.variables["channel_id"])]
   index["channels"] = np.array([cf.index(c) for c in cs])
   for dim, var in (("solar_channels", "solar_channel_id"), ("thermal_channels", "thermal_channel_id")):
      if var in small.variables:
         a = [int(c) for c in np.asarray(small.variables[var])]
         b = [int(c) for c in np.asarray(full.variables[var])]
         index[dim] = np.array([b.index(c) for c in a])
   worst = 0.0
   compared = 0
   for name in small.variable_names:
      x = np.asarray(small.variables[name])
      if x.dtype.kind != "f" or name in COORDS:
         continue
      dims = small.variable_dimensions[name]
      if not all(d in index for d in dims):
         continue
      y = np.asarray(full.variables[name])[np.ix_(*[index[d] for d in dims])]
      d = float(np.nanmax(np.abs(x.astype(np.float64) - y.astype(np.float64)))) if x.size else 0.0
      worst = max(worst, d)
      compared += 1
      print(f"{name:32s} {str(dims):95s} max |diff| {d:.3e}{'  (bitwise)' if np.array_equal(x, y, equal_nan=True) else ''}")
   print(f"{case}: {compared} variables compared with {product.name}; largest difference {worst:.3e}")


if __name__ == "__main__":
   main()
