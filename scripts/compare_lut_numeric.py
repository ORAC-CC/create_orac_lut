"""Compare the numeric variables of two ORAC LUT files, with a float32 ulp tolerance.

   python scripts/compare_lut_numeric.py A.nc B.nc [--max-ulp N] [--record OUT.json]

Every numeric variable present in both files is compared on its raw stored
values (fill values included, so a fill in one file and a value in the other
is a difference).  Reports per variable whether the values are bitwise
identical and the largest difference in float32 units in the last place.
Exit status 1 if any variable differs by more than --max-ulp (default 0,
i.e. bitwise) or if the shapes or variable sets differ.  Read-only.
"""

import argparse
import json
import sys

import numpy as np
from netCDF4 import Dataset


def ulp_difference(a, b):
   a = np.asarray(a, np.float32).ravel()
   b = np.asarray(b, np.float32).ravel()
   ia = a.view(np.int32).astype(np.int64)
   ib = b.view(np.int32).astype(np.int64)
   ia = np.where(ia < 0, np.int64(-2 ** 31) - ia, ia)          # monotone integer ordering of floats
   ib = np.where(ib < 0, np.int64(-2 ** 31) - ib, ib)
   return int(np.max(np.abs(ia - ib))) if a.size else 0


def compare(path_a, path_b, max_ulp):
   result, failed = {}, False
   with Dataset(path_a) as a, Dataset(path_b) as b:
      a.set_auto_maskandscale(False)
      b.set_auto_maskandscale(False)
      names_a = {n for n, v in a.variables.items() if v.dtype.kind in "iuf"}
      names_b = {n for n, v in b.variables.items() if v.dtype.kind in "iuf"}
      if names_a != names_b:
         failed = True
         result["_variable_sets"] = {"only_a": sorted(names_a - names_b), "only_b": sorted(names_b - names_a)}
      for name in sorted(names_a & names_b):
         x, y = np.asarray(a[name][:]), np.asarray(b[name][:])
         if x.shape != y.shape:
            result[name] = {"shape_a": list(x.shape), "shape_b": list(y.shape)}
            failed = True
            continue
         bitwise = bool(np.array_equal(x, y))
         ulp = 0 if bitwise else (ulp_difference(x, y) if x.dtype.kind == "f" else None)
         result[name] = {"bitwise": bitwise, "max_float32_ulp": ulp,
                         "max_abs": float(np.max(np.abs(x.astype(np.float64) - y.astype(np.float64)))) if x.size else 0.0}
         if not bitwise and (ulp is None or ulp > max_ulp):
            failed = True
   return result, not failed


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("a")
   parser.add_argument("b")
   parser.add_argument("--max-ulp", type=int, default=0)
   parser.add_argument("--record")
   args = parser.parse_args()
   result, ok = compare(args.a, args.b, args.max_ulp)
   if args.record:
      with open(args.record, "w") as handle:
         json.dump({"a": args.a, "b": args.b, "max_ulp": args.max_ulp, "passed": ok, "variables": result}, handle, indent=1)
   differing = {k: v for k, v in result.items() if not v.get("bitwise", False)}
   print(f"{args.a}\n  vs {args.b}: {len(result)} variables, {len(differing)} not bitwise; "
         + ("PASS" if ok else "FAIL") + f" (tolerance {args.max_ulp} float32 ulp)")
   for name, value in differing.items():
      print(f"  {name}: {value}")
   return 0 if ok else 1


if __name__ == "__main__":
   raise SystemExit(main())
