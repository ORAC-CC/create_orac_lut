"""Stage 9: compare the corrected single-precision V22 Terra MODIS LUTs with the production V21 IDL LUTs.

Reference: the shared V21 products (releases/idl_lut_v21_single.release).
Candidate: the definitive V22 products (releases/lut_v22_single-r2.release), promoted
to the same shared directory; set ORAC_LUT_V22_CANDIDATES to compare another directory
(e.g. a regeneration into staging/lut_v22_single-r2/).

For every product:
* dimensions, variable set and order, attributes (compare_luts.py);
* every variable: raw-byte equality, elements differing, maximum absolute and
  relative difference, RMS, p95/p99 of |difference|, coordinates of the
  largest difference (compare_luts.py, decoded values; plus a raw-byte check
  of the stored arrays);
* inside ORAC's operational angular domain, including the SZA 80 support
  nodes through interpolation (operational_domain.py): coefficient
  differences, direct-beam calls, and forward-model measurement differences
  in units of the measurement uncertainty.

If every RT operator is raw-byte identical, the reconstructed radiances and
the retrieval Jacobians (which are functions of the operators only) are
identical as well; this is stated in the result rather than recomputed.

Usage: python compare_full_luts.py   -> results/full_luts/<product>.json, summary.json
"""

import json
import os
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import compare_luts                                          # noqa: E402
import operational_domain as od                              # noqa: E402

V21 = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
STAGE = Path(os.environ.get("ORAC_LUT_V22_CANDIDATES", "/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS"))
OUT = HERE / "results" / "full_luts"


def raw_equality(reference_path, candidate_path):
   rows = {}
   with Dataset(reference_path) as r, Dataset(candidate_path) as c:
      for d in (r, c):
         d.set_auto_maskandscale(False)
         d.set_auto_chartostring(False)
      for name in r.variables:
         if name not in c.variables:
            rows[name] = None
            continue
         a, b = np.asarray(r[name][:]), np.asarray(c[name][:])
         rows[name] = bool(a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes())
   return rows


def main():
   OUT.mkdir(parents=True, exist_ok=True)
   summary = {}
   for case, (stem, cloud) in od.CASES.items():
      reference_path = V21 / f"{stem}_v21.nc"
      candidate_path = STAGE / f"{stem}_v22.nc"
      if not candidate_path.is_file():
         summary[case] = {"status": "not staged", "candidate": str(candidate_path)}
         print(case, "not staged yet")
         continue
      result, coordinate_rows, calls = compare_luts.compare(reference_path, candidate_path)
      raw = raw_equality(reference_path, candidate_path)
      variables = {}
      for c in result["comparisons"]:
         keys = ("status", "elements", "differing_elements", "maximum_absolute_difference", "maximum_relative_difference",
                 "rms_difference", "absolute_difference_percentiles", "maximum_absolute_index", "reference_minimum", "candidate_minimum")
         variables[c["variable"]] = {k: c.get(k) for k in keys} | {"raw_bitwise": raw.get(c["variable"])}
      rt_bitwise = all(raw.get(name) for name in od.RT_VARIABLES if name in raw)
      record = {"reference": str(reference_path), "candidate": str(candidate_path),
                "dimensions_equal": result["dimensions"]["same_names_and_lengths"],
                "variable_order_equal": result["variable_order"]["same"],
                "attribute_differences": result["attribute_differences"], "variables": variables,
                "raw_bitwise_variables": sum(1 for v in raw.values() if v), "variables_total": len(raw),
                "rt_operators_raw_bitwise": rt_bitwise}
      with Dataset(reference_path) as reference, Dataset(candidate_path) as candidate:
         sza = od.impact._values(reference, "solar_zenith")
         geometries, weights = od.operational_sza_weights(sza, od.MAX_SOL_ZEN)
         sza_classes = od.classify_nodes(sza, od.MAX_SOL_ZEN, inclusive=False)
         record["negative_values"] = {name: {"v21": int(np.count_nonzero(np.asarray(reference[name][:]) < 0)),
                                             "v22": int(np.count_nonzero(np.asarray(candidate[name][:]) < 0))}
                                      for name in od.RT_VARIABLES if name in reference.variables}
         if rt_bitwise:
            record["operational_domain"] = ("every RT operator is raw-byte identical: reconstructed radiances and "
                                            "retrieval Jacobians are identical in the whole domain")
         else:
            coefficient_rows, _ = od.coefficient_analysis(reference, candidate, sza, geometries, weights, sza_classes)
            call_rows = od.classify_calls(od.call_analysis(reference, candidate, sza_classes))
            record["operational_domain"] = {
               "coefficients": coefficient_rows,
               "call_classes": {cls: {s: sum(1 for r in call_rows if r["call_class"] == cls and r["sza_class"] == s)
                                      for s in ("in_domain", "support", "unused")}
                                for cls in ("isolated_significant", "significant_not_isolated", "below_1_sigma")},
               "measurements": od.measurement_rows(reference, candidate, cloud, geometries, weights)}
      (OUT / f"{case}.json").write_text(json.dumps(record, indent=1, default=str) + "\n")
      differing = {k: v["differing_elements"] for k, v in variables.items() if v.get("differing_elements")}
      not_raw = [k for k, v in raw.items() if not v]
      summary[case] = {"raw_bitwise": f"{record['raw_bitwise_variables']}/{record['variables_total']}", "not_raw_bitwise": not_raw,
                       "decoded_differing": differing, "rt_operators_raw_bitwise": rt_bitwise}
      print(case, json.dumps(summary[case]))
   (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
   main()
