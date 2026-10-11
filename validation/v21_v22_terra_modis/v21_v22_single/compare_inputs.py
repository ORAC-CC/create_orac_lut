"""Stage 5: compare every DISORT argument and output of a Python capture with the V21 IDL captures.

Reference: reference/v21_idl_single_calls.npz (90 V21 IDL calls).
Python: results/<tag>/python_single_case_<X>.npz from capture_python_single.py <tag>.

For each argument: calls differing, elements differing, max |difference|,
max ulp, first difference.  Integer, logical and dimension arguments are
compared exactly; outputs bitwise (the Python capture holds the outputs of the
kernel actually used).

Usage: python compare_inputs.py <tag>   -> results/<tag>/compare_inputs.json
"""

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference" / "v21_idl_single_calls.npz"
ARRAYS = ("dtauc", "ssalb", "pmom", "utau", "umu", "phi", "temper")
SCALARS = ("umu0", "fbeam", "fisot", "wvnmlo", "wvnmhi")
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")


def ulps(a, b):
   a = np.asarray(a, np.float32).ravel().view(np.int32).astype(np.int64)
   b = np.asarray(b, np.float32).ravel().view(np.int32).astype(np.int64)
   a = np.where(a < 0, -(a & 0x7FFFFFFF), a)
   b = np.where(b < 0, -(b & 0x7FFFFFFF), b)
   return np.abs(a - b)


def calls(tag):
   ref = np.load(REFERENCE)
   out = []
   index = 0
   for case in "ABCD":
      z = np.load(HERE / "results" / tag / f"python_single_case_{case}.npz")
      for j in range(int(z["ncalls"])):
         k, p = f"c{index:03d}_", f"c{j:03d}_"
         assert str(ref[k + "case"]) == case
         idl = {name: ref[k + name] for name in ARRAYS + SCALARS}
         idl.update({"out_" + o: ref[k + "out_" + o] for o in OUTPUTS})
         idl.update(nlyr=int(ref[k + "nlyr"]), nmom=int(ref[k + "nmom"]), nstr=int(ref[k + "nstr_v21"]), plank=int(ref[k + "plank"]))
         py = {name: z[p + name] for name in ("dtauc", "ssalb", "pmom", "utau", "umu", "phi")}
         n = py["dtauc"].size
         temper = np.zeros(n + 1, np.float32)
         temper[:z[p + "temper"].size] = z[p + "temper"]
         py["temper"] = temper
         py.update({name: z[p + name] for name in SCALARS})
         py.update({"out_" + o: z[p + "out_" + o] for o in OUTPUTS})
         py.update(nlyr=n, nmom=py["pmom"].shape[0] - 1, nstr=int(z[p + "nstr"]), plank=int(z[p + "plank"]))
         py["dtypes"] = str(z[p + "dtypes"]) if p + "dtypes" in z.files else ""
         kind = "emission" if idl["plank"] else ("direct" if idl["fbeam"] > 0 else "diffuse")
         out.append((case, index, kind, idl, py))
         index += 1
   assert index == int(ref["ncalls"])
   return out


def main():
   tag = sys.argv[1]
   rows = {}
   config = {"nlyr": 0, "nmom": 0, "nstr": 0, "plank": 0, "shape": 0}
   bitwise_calls = []
   dtypes = set()
   for case, index, kind, idl, py in calls(tag):
      dtypes.add(py["dtypes"])
      for name in ("nlyr", "nmom", "nstr", "plank"):
         config[name] += idl[name] != py[name]
      for name in ARRAYS + SCALARS + tuple("out_" + o for o in OUTPUTS):
         a, b = np.asarray(idl[name], np.float32), np.asarray(py[name], np.float32)
         row = rows.setdefault(name, dict(calls_differing=0, elements_differing=0, max_abs=0.0, max_ulp=0, first=None, shape_mismatch=0))
         if a.shape != b.shape:
            row["shape_mismatch"] += 1
            continue
         d = a.ravel().view(np.uint32) != b.ravel().view(np.uint32)
         if d.any():
            row["calls_differing"] += 1
            row["elements_differing"] += int(d.sum())
            row["max_abs"] = max(row["max_abs"], float(np.max(np.abs(a.astype(np.float64) - b.astype(np.float64)))))
            row["max_ulp"] = max(row["max_ulp"], int(ulps(a, b).max()))
            if row["first"] is None:
               e = int(np.argmax(d))
               row["first"] = dict(case=case, call=index, kind=kind, element=int(e), idl=float(a.ravel()[e]), python=float(b.ravel()[e]))
      if all(np.array_equal(np.asarray(idl["out_" + o], np.float32).view(np.uint32), np.asarray(py["out_" + o], np.float32).view(np.uint32)) for o in OUTPUTS):
         bitwise_calls.append(index)
   summary = {"tag": tag, "arguments": rows, "configuration_mismatches": config, "bitwise_output_calls": len(bitwise_calls),
              "calls": 90, "non_bitwise_calls": sorted(set(range(90)) - set(bitwise_calls)), "backend_input_dtypes": sorted(dtypes)}
   (HERE / "results" / tag / "compare_inputs.json").write_text(json.dumps(summary, indent=1) + "\n")
   for name, row in rows.items():
      print(f"{name:12s} calls differing {row['calls_differing']:3d}  elements {row['elements_differing']:6d}  max abs {row['max_abs']:.3g}  max ulp {row['max_ulp']}  shape mismatch {row['shape_mismatch']}")
   print("configuration mismatches", config)
   print("bitwise-identical solver calls", len(bitwise_calls), "/ 90")
   print("backend input dtypes", sorted(dtypes))


if __name__ == "__main__":
   main()
