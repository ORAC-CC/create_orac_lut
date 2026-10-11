"""Stage 6: compare compact V22 Python LUTs with the V21 IDL compact LUTs, every NetCDF variable.

Reference: idl_double/results/run_20261010T214143/case_<X>/validation/luts/case_<X>/
terra_modis_m_liquid-water_a01_p240_v21.nc (single-precision V21 IDL, equal to the
production V21 products in every RT coefficient).  Candidates: work_<tag>/case_<X>/luts/
(and the untagged 17de4c0 baseline work/ for comparison).  Uses
validation/v21_v22_terra_modis/compare_luts.py.

Usage: python compare_compact_luts.py <tag>   -> results/<tag>/compare_compact_luts.json
"""

import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import compare_luts                                        # noqa: E402

RUN = HERE.parent / "idl_double" / "results" / "run_20261010T214143"
NAME = "terra_modis_m_liquid-water_a01_p240_{}.nc"


def negatives(path):
   counts = {}
   with Dataset(path) as d:
      d.set_auto_mask(False)
      for name in compare_luts.RT_OPERATORS:
         if name in d.variables:
            counts[name] = int(np.count_nonzero(np.asarray(d[name][:]) < 0))
   return counts


def summary(reference, candidate):
   result, rows, calls = compare_luts.compare(reference, candidate)
   table = {}
   for c in result["comparisons"]:
      table[c["variable"]] = {k: c.get(k) for k in ("status", "elements", "differing_elements", "exactly_equal",
                                                     "maximum_absolute_difference", "maximum_relative_difference",
                                                     "rms_difference", "maximum_absolute_index")}
   return {"dimensions_equal": result["dimensions"]["same_names_and_lengths"], "variable_order_equal": result["variable_order"]["same"],
           "attribute_differences": result["attribute_differences"], "variables": table,
           "negative_values": {"reference": negatives(reference), "candidate": negatives(candidate)}}


def main():
   tag = sys.argv[1]
   out = {}
   for case in "ABCD":
      reference = RUN / f"case_{case}" / "validation" / "luts" / f"case_{case}" / NAME.format("v21")
      out[case] = {"corrected": summary(reference, HERE / f"work_{tag}" / f"case_{case}" / "luts" / NAME.format("v22")),
                   "baseline_17de4c0": summary(reference, HERE / "work" / f"case_{case}" / "luts" / NAME.format("v22"))}
      for label, s in out[case].items():
         differing = {k: v["differing_elements"] for k, v in s["variables"].items() if v.get("differing_elements")}
         print(case, label, "dims equal", s["dimensions_equal"], "| variables differing:", differing or "none",
               "| attribute differences:", len(s["attribute_differences"]), "| negatives ref/cand:",
               sum(s["negative_values"]["reference"].values()), sum(s["negative_values"]["candidate"].values()))
   (HERE / "results" / tag / "compare_compact_luts.json").write_text(json.dumps(out, indent=1, default=str) + "\n")


if __name__ == "__main__":
   main()
