"""Validate the full V24 Terra MODIS products against the V23 products on the same Grid B (read-only).

   python validation/v24/scattering_refinements/compare_full_products.py   -> results/full_<product>.json

For each of pold, p240 and pagg in staging/lut_v24/:
* structure (validation/v23_terra_modis/compare_with_reference.structure): the
  V23 variable set, Grid B coordinates bitwise (liquid-water-cloud_v02.lut /
  ice-cloud_v02.lut), 36 channels, no NaN/Inf, fill pattern as V23;
* provenance: LUT release V24 made by releases/lut_v24.release, scientific code
  lut-code-v24.1;
* against staging/lut_v23/ (identical grid nodes): pagg must be identical in
  every scientific variable (the Baum class is calculated as in V23); for the
  liquid-water products the median, 99th percentile and maximum |V24 - V23| of
  every RT operator and of the bulk optics, with the location of the maximum.
"""

import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "validation" / "v23_terra_modis"))
from compare_with_reference import FILL, PRODUCTS, RT, structure   # noqa: E402

OPTICS = ("extinction_coefficient", "extinction_coefficient_ratio", "single_scatter_albedo", "asymmetry_parameter")
PROVENANCE_ATTRIBUTES = ("lut_release", "lut_release_version")


def main():
   for product, (stem, grid) in PRODUCTS.items():
      v24, v23 = ROOT / "staging" / "lut_v24" / f"{stem}_v24.nc", ROOT / "staging" / "lut_v23" / f"{stem}_v23.nc"
      if not v24.is_file():
         print(product, "not staged")
         continue
      record = {"product": str(v24), "v23": str(v23), "structure": structure(v24, grid), "operators": {}}
      with Dataset(v23) as a, Dataset(v24) as b:
         for d in (a, b):
            d.set_auto_maskandscale(False)
            d.set_auto_chartostring(False)
         a_fill = {n: int(np.count_nonzero(np.asarray(a[n][:]) == FILL)) for n in a.variables if np.issubdtype(a[n].dtype, np.floating)}
         record["fill_as_v23"] = record["structure"]["fill"] == {n: c for n, c in a_fill.items() if c}
         record["identical_variables"] = [n for n in a.variables if np.asarray(a[n][:]).tobytes() == np.asarray(b[n][:]).tobytes()]
         record["all_scientific_variables_identical"] = len(record["identical_variables"]) == len(a.variables)
         for name in RT + OPTICS:
            x, y = np.asarray(a[name][:], np.float64), np.asarray(b[name][:], np.float64)
            valid = (x != FILL) & (y != FILL)
            d = np.where(valid, np.abs(y - x), 0.0)
            imax = np.unravel_index(np.argmax(d), d.shape)
            record["operators"][name] = {"median": float(np.median(d[valid])), "p99": float(np.percentile(d[valid], 99)),
                                         "max": float(d.max()), "max_index": [int(i) for i in imax],
                                         "dimensions": list(a[name].dimensions), "values_at_max": [float(x[imax]), float(y[imax])]}
      out = HERE / "results" / f"full_{product}.json"
      out.write_text(json.dumps(record, indent=1, default=str) + "\n")
      s = record["structure"]
      print(f"{product}: coordinates bitwise {all(v['bitwise'] for v in s.values() if isinstance(v, dict) and 'bitwise' in v)}, "
            f"NaN/Inf {s['nan_inf'] or 'none'}, fill as V23 {record['fill_as_v23']}, release {s['provenance']['lut_release']} "
            f"{s['provenance']['scientific_code']}; identical variables {len(record['identical_variables'])}")
      for name, e in record["operators"].items():
         print(f"   {name:30s} median {e['median']:.1e} p99 {e['p99']:.1e} max {e['max']:.2e} at {e['max_index']}")


if __name__ == "__main__":
   main()
