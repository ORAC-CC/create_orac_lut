"""Single-precision V21 (IDL) / V22 (Python) DISORT reproduction experiments.

Inputs:
  V21: IDL captures of the historical single-precision calls (pre-call copies),
       idl_double/results/run_20261010T214143/case_*/validation/captures
  V22: results/python_single_case_*.npz (capture_python_single.py; active
       TOTAL-corrected Python, ORAC_LUT_DISORT_PRECISION=single)

1. Every DISORT argument compared call by call (IDL vs Python).
2. Four-way experiment per call: A = V21 inputs / Intel DLM kernel (IDL memory
   layout), B = V22 inputs / Intel DLM, C = V21 inputs / gfortran kernel,
   D = V22 inputs / gfortran kernel; plus A vs the IDL-captured outputs and D
   vs the Python-captured outputs.
3. Input-substitution experiments on the failing p240 call (case A: channel 6,
   tau 256, r_e 29, SZA 80), forward (V22 base, V21 groups substituted) and
   reverse (V21 base, V22 groups), all with the Intel DLM kernel.

Usage: python four_way.py      (intel-compilers/2022 runtime on LD_LIBRARY_PATH)
Writes results/four_way.json.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "idl_double"))
import analyse_idl_double as idl                       # noqa: E402
import sp_kernels as K                                  # noqa: E402

IDL_RUN = HERE.parent / "idl_double" / "results" / "run_20261010T214143"
RESULTS = HERE / "results"
ARRAYS = ("dtauc", "ssalb", "pmom", "utau", "umu", "phi", "temper")
SCALARS = ("fbeam", "umu0", "fisot", "wvnmlo", "wvnmhi")
# Fixed arguments of the Python C adapter (legacy_disort.c) versus the IDL DISORT_VARS values.
PY_FIXED = {"usrtau": 1, "usrang": 1, "ibcnd": 0, "lamber": 1, "onlyfl": 0, "phi0": 0.0, "albedo": 0.0,
            "btemp": 0.0, "ttemp": 0.0, "temis": 0.0, "accur": float(np.float32(1e-8)), "nstr": 60}


def ordered(x):
   """Monotone integer image of float32 values (ulp distance = difference)."""
   i = np.asarray(x, np.float32).view(np.int32).astype(np.int64)
   return np.where(i < 0, np.int64(-2 ** 31) - i, i)


def compare_array(a, b):
   a = np.asarray(a, np.float32)
   b = np.asarray(b, np.float32)
   if a.shape != b.shape:
      return {"shape_idl": list(a.shape), "shape_python": list(b.shape)}
   diff = a != b
   d = np.abs(a.astype(np.float64) - b)
   rel = np.where(a != 0, d / np.abs(np.where(a != 0, a, 1)), np.where(d > 0, np.inf, 0))
   record = {"dtype": "float32", "shape": list(a.shape), "differing": int(diff.sum()), "size": int(a.size),
             "max_abs": float(d.max()) if a.size else 0.0, "max_rel": float(rel.max()) if a.size else 0.0,
             "max_ulp": int(np.max(np.abs(ordered(a) - ordered(b)))) if a.size else 0}
   if diff.any():
      first = tuple(int(v) for v in np.argwhere(diff)[0])
      record["first_differing_index"] = list(first)
      record["first_values_idl_python"] = [float(a[first]), float(b[first])]
   return record


def load_idl():
   calls = {}
   for case in "ABCD":
      calls[case] = []
      for path in sorted((IDL_RUN / f"case_{case}" / "validation" / "captures").glob("call_*.sav")):
         c = idl.load_capture(path)
         a = dict(c["args"])
         a["temper"] = np.asarray(a["temper"], np.float32)[:a["nlyr"] + 1] if a["plank"] else np.zeros(a["nlyr"], np.float32)
         c["args"] = a
         calls[case].append(c)
   return calls


def load_python():
   calls = {}
   for case in "ABCD":
      z = np.load(RESULTS / f"python_single_case_{case}.npz")
      calls[case] = []
      for i in range(int(z["ncalls"])):
         g = lambda k: z[f"c{i:03d}_{k}"]
         pm = np.asarray(g("pmom"), np.float32)
         args = {"nlyr": g("dtauc").size, "nmom": pm.shape[0] - 1, "dtauc": g("dtauc"), "ssalb": g("ssalb"), "pmom": pm,
                 "utau": g("utau"), "umu": g("umu"), "phi": g("phi"), "fbeam": float(g("fbeam")), "umu0": float(g("umu0")),
                 "fisot": float(g("fisot")), "plank": int(g("plank")), "wvnmlo": float(g("wvnmlo")), "wvnmhi": float(g("wvnmhi")),
                 "temper": np.asarray(g("temper"), np.float32), "accur": float(np.float32(1e-8)), "phi0": 0.0, "nstr": int(g("nstr"))}
         out = {k: np.asarray(g("out_" + k)) for k in K.OUTPUTS}
         calls[case].append({"args": args, "out": out, "kind": "emission" if args["plank"] else ("direct" if args["fbeam"] > 0 else "diffuse")})
   return calls


def output_difference(x, y):
   return {k: {"differing": int((np.float32(x[k]) != np.float32(y[k])).sum()),
               "max_abs": float(np.max(np.abs(np.float64(x[k]) - np.float64(y[k]))))} for k in K.OUTPUTS}


def bitwise(x, y):
   return all(np.array_equal(np.float32(x[k]), np.float32(y[k])) for k in K.OUTPUTS)


def main():
   RESULTS.mkdir(exist_ok=True)
   v21, v22 = load_idl(), load_python()
   report = {"inputs": {}, "four_way": {}, "fixed_arguments": {}}

   # 1. Arguments.
   aggregate = {name: {"calls_differing": 0, "max_ulp": 0, "max_abs": 0.0, "first": None} for name in ARRAYS + SCALARS}
   integer_mismatch = []
   for case in "ABCD":
      if len(v21[case]) != len(v22[case]):
         raise AssertionError(f"case {case}: {len(v21[case])} IDL calls, {len(v22[case])} Python calls")
      for index, (c, p) in enumerate(zip(v21[case], v22[case])):
         a, b = c["args"], p["args"]
         if c["kind"] != p["kind"] or a["nlyr"] != b["nlyr"] or a["nmom"] != b["nmom"]:
            integer_mismatch.append((case, index, c["kind"], p["kind"], a["nlyr"], b["nlyr"]))
         per_call = {}
         for name in ARRAYS:
            r = compare_array(a[name][:a["nlyr"] + (1 if name == "temper" and a["plank"] else 0)] if name in ("dtauc", "ssalb", "temper")
                              else (np.asarray(a[name])[:, :a["nlyr"]] if name == "pmom" else a[name]),
                              b[name] if name != "temper" else b[name][:b["nlyr"] + (1 if b["plank"] else 0)])
            per_call[name] = r
         for name in SCALARS:
            per_call[name] = compare_array(np.array([a.get(name, 0.0)]), np.array([b.get(name, 0.0)]))
         for name, r in per_call.items():
            g = aggregate[name]
            if r.get("differing", 0):
               g["calls_differing"] += 1
               g["max_ulp"] = max(g["max_ulp"], r["max_ulp"])
               g["max_abs"] = max(g["max_abs"], r["max_abs"])
               if g["first"] is None:
                  g["first"] = {"case": case, "call": index + 1, "kind": c["kind"], **{k: r[k] for k in ("first_differing_index", "first_values_idl_python", "max_rel")}}
   report["inputs"] = {"per_argument": aggregate, "integer_or_kind_mismatches": integer_mismatch}
   c0 = v21["A"][0]
   idl_fixed = {"usrtau": 1, "usrang": 1, "ibcnd": 0, "lamber": 1, "onlyfl": 0, "phi0": c0["args"]["phi0"], "albedo": c0["args"]["albedo"],
                "btemp": c0["args"]["btemp"], "ttemp": c0["args"]["ttemp"], "temis": c0["args"]["temis"], "accur": c0["args"]["accur"],
                "nstr": c0["nstr_v21"]}
   report["fixed_arguments"] = {k: {"idl": idl_fixed[k], "python": PY_FIXED[k], "equal": idl_fixed[k] == PY_FIXED[k]} for k in PY_FIXED}

   # 2. Four-way.
   tally = {"A_eq_idl": 0, "D_eq_python": 0, "A_eq_C": 0, "B_eq_D": 0, "A_eq_B": 0, "B_eq_idl": 0, "calls": 0}
   worst = {"A_vs_C": 0.0, "B_vs_D": 0.0, "A_vs_B": 0.0, "D_vs_A": 0.0}
   by_kind = {}
   for case in "ABCD":
      for c, p in zip(v21[case], v22[case]):
         A = K.call("intel_dlm", c["args"])[1]
         B = K.call("intel_dlm", p["args"])[1]
         C = K.call("gfortran", c["args"])[1]
         D = K.call("gfortran", p["args"])[1]
         tally["calls"] += 1
         tally["A_eq_idl"] += bitwise(A, c["sp"])
         tally["D_eq_python"] += bitwise(D, p["out"])
         tally["A_eq_C"] += bitwise(A, C)
         tally["B_eq_D"] += bitwise(B, D)
         tally["A_eq_B"] += bitwise(A, B)
         tally["B_eq_idl"] += bitwise(B, c["sp"])
         for label, (x, y) in {"A_vs_C": (A, C), "B_vs_D": (B, D), "A_vs_B": (A, B), "D_vs_A": (D, A)}.items():
            m = max(v["max_abs"] for v in output_difference(x, y).values())
            worst[label] = max(worst[label], m)
            k = by_kind.setdefault(c["kind"], {})
            k[label] = max(k.get(label, 0.0), m)
   report["four_way"] = {"bitwise_counts": tally, "max_abs_output_difference": worst, "by_call_kind": by_kind}

   # 3. Substitution on the failing call.
   # Case A loop order: tau (64, 128, 256) outer, r_e (27, 29, 31), then diffuse, SZA 70, SZA 80 calls;
   # tau 256 (a = 2), r_e 29 (r = 1), SZA 80 is capture (2*3 + 1)*3 + 2 + 1 = 24.
   idx = (2 * 3 + 1) * 3 + 2
   assert v21["A"][idx]["kind"] == "direct" and abs(v21["A"][idx]["args"]["umu0"] - np.cos(np.deg2rad(80))) < 1e-6
   base21, base22 = v21["A"][idx]["args"], v22["A"][idx]["args"]
   groups = {"DTAUC+UTAU": ("dtauc", "utau"), "SSALB": ("ssalb",), "PMOM": ("pmom",), "geometry (UMU, PHI, UMU0)": ("umu", "phi", "umu0"),
             "scalars (FBEAM, FISOT, WVNMLO/HI)": ("fbeam", "fisot", "wvnmlo", "wvnmhi")}
   ref21 = K.call("intel_dlm", base21)[1]
   ref22 = K.call("intel_dlm", base22)[1]
   toa_up = lambda o: np.float64(o["uu"][10:, 0, :])
   def run(base, donor, keys):
      args = dict(base)
      for key in keys:
         args[key] = donor[key]
      return K.call("intel_dlm", args)[1]
   substitution = {"call": f"case A capture {idx + 1}", "v21_vs_v22_max_abs_toa_up": float(np.max(np.abs(toa_up(ref21) - toa_up(ref22)))),
                   "forward_from_v22": {}, "reverse_from_v21": {}, "cumulative_from_v22": {}}
   for name, keys in groups.items():
      f = run(base22, base21, keys)
      substitution["forward_from_v22"][name] = {"max_abs_toa_up_vs_v21": float(np.max(np.abs(toa_up(f) - toa_up(ref21)))),
                                                "max_abs_toa_up_vs_v22": float(np.max(np.abs(toa_up(f) - toa_up(ref22)))),
                                                "bitwise_v21": bitwise(f, ref21)}
      r = run(base21, base22, keys)
      substitution["reverse_from_v21"][name] = {"max_abs_toa_up_vs_v21": float(np.max(np.abs(toa_up(r) - toa_up(ref21)))),
                                                "max_abs_toa_up_vs_v22": float(np.max(np.abs(toa_up(r) - toa_up(ref22)))),
                                                "bitwise_v22": bitwise(r, ref22)}
   keys = []
   for name, group in groups.items():
      keys += list(group)
      f = run(base22, base21, keys)
      substitution["cumulative_from_v22"][name] = {"max_abs_toa_up_vs_v21": float(np.max(np.abs(toa_up(f) - toa_up(ref21)))),
                                                   "bitwise_v21": bitwise(f, ref21)}
   report["substitution"] = substitution
   (RESULTS / "four_way.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
   print(json.dumps({k: report[k] for k in ("four_way", "fixed_arguments")}, indent=1, default=str))
   print(json.dumps(report["inputs"], indent=1, default=str))
   print(json.dumps(report["substitution"], indent=1, default=str))


if __name__ == "__main__":
   main()
