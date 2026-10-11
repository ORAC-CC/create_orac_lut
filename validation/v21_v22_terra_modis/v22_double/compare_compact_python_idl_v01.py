"""Compact V22 (Python, double-precision DISORT) versus IDL double-precision results.

Runs the active Python generator (create_orac_luts.create_orac_cloud_lut) on
the four compact Terra MODIS p240 cases of the IDL double-precision
validation (validation/v21_v22_terra_modis/idl_double/cases), capturing every
DISORT call, and compares with the IDL run
idl_double/results/run_20261010T214143:

* DISORT inputs, call by call, against the IDL captures (pre-call copies);
* DISORT outputs against IDL DISORT_DP at NSTR 60;
* the compact Python LUT against the compact IDL LUT made with DISORT_DP at
  NSTR 60 (version 90 label; identical numerics to the integrated V21 path).

Usage: python compare_compact_python_idl.py [OUTPUT_JSON]
Precision is taken from ORAC_LUT_DISORT_PRECISION (default double).
"""

from __future__ import annotations

import importlib
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
IDL_RUN = HERE.parent / "idl_double" / "results" / "run_20261010T214143"
CASES_DIR = HERE.parent / "idl_double" / "cases"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE.parent / "idl_double"))
import analyse_idl_double as idl_analysis     # noqa: E402

CASES = {"A": 6, "B": 20, "C": 6, "D": 1}
RT = ("R_0v", "R_0d", "T_0d", "T_00", "R_dv", "T_dv", "R_dd", "T_dd", "E_md")
OPTICS = ("extinction_coefficient", "extinction_coefficient_ratio", "single_scatter_albedo", "asymmetry_parameter")
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")


def stage(case):
   in_path = WORK / f"case_{case}" / "input_files"
   if in_path.exists():
      shutil.rmtree(in_path.parent)
   (in_path / "lut").mkdir(parents=True)
   for entry in (ROOT / "create_orac_lut" / "input_files").iterdir():
      if entry.name not in ("lut", "driver"):
         (in_path / entry.name).symlink_to(entry)
   shutil.copy(CASES_DIR / f"case_{case}.lut", in_path / "lut" / f"case_{case}.lut")
   return in_path


def run_python(case, channel):
   production = importlib.import_module("create_orac_luts")
   call_module = importlib.import_module("oraclut.idl_mirror.call_disort")
   original = call_module.disort
   calls = []

   def capture(dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kw):
      out = original(dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kw)
      calls.append({"dtauc": np.array(dtauc, np.float32), "ssalb": np.array(ssalb, np.float32),
                    "pmom": np.array(pmom, np.float32), "utau": np.array(utau, np.float32),
                    "umu": np.array(umu, np.float32), "phi": np.array(phi, np.float32),
                    "fbeam": float(np.float32(fbeam)), "umu0": float(np.float32(umu0)),
                    "fisot": float(np.float32(fisot)), "plank": bool(kw.get("plank")),
                    "out": {k: np.array(v, np.float64) for k, v in out.items()}})
      return out

   in_path = stage(case)
   out_dir = WORK / f"case_{case}" / "luts"
   out_dir.mkdir(parents=True, exist_ok=True)
   call_module.disort = capture
   try:
      production.create_orac_cloud_lut(in_path, "terra_modis_v1.inst", "liquid-water_240.mm", f"case_{case}.lut",
                                       out_dir, "2", channelid=[channel], gas=0, no_rayleigh=0, srf_quad=1,
                                       reuse_scat=0, scat_only=0, tmatrix_path=None, version=22, driver=None,
                                       nstreams=60, nmom=1000, work_path=WORK / f"case_{case}" / "work")
   finally:
      call_module.disort = original
   lut = next(out_dir.rglob("*_v22.nc"))
   return calls, lut


def ulp32(a, b):
   a = np.asarray(a, np.float32).ravel()
   b = np.asarray(b, np.float32).ravel()
   ia = a.view(np.int32).astype(np.int64)
   ib = b.view(np.int32).astype(np.int64)
   return int(np.max(np.abs(ia - ib))) if a.size else 0


def compare_case(case, channel):
   calls, py_lut = run_python(case, channel)
   captures = [idl_analysis.load_capture(p) for p in sorted((IDL_RUN / f"case_{case}" / "validation" / "captures").glob("call_*.sav"))]
   if len(calls) != len(captures):
      raise AssertionError(f"case {case}: {len(calls)} Python calls, {len(captures)} IDL calls")
   inputs = {k: 0 for k in ("dtauc", "ssalb", "pmom", "utau", "umu", "phi")}
   scalars_equal = True
   outputs = {k: 0.0 for k in OUTPUTS}
   for py, cap in zip(calls, captures):
      a = cap["args"]
      n = a["nlyr"]
      if py["plank"] != bool(a["plank"]):
         raise AssertionError(f"case {case}: call kinds differ")
      inputs["dtauc"] = max(inputs["dtauc"], ulp32(py["dtauc"], a["dtauc"][:n]))
      inputs["ssalb"] = max(inputs["ssalb"], ulp32(py["ssalb"], a["ssalb"][:n]))
      inputs["pmom"] = max(inputs["pmom"], ulp32(py["pmom"], a["pmom"][:, :n]))
      for k in ("utau", "umu", "phi"):
         inputs[k] = max(inputs[k], ulp32(py[k], a[k]))
      scalars_equal &= all(py[k] == a[k] for k in ("fbeam", "umu0", "fisot"))
      for k in OUTPUTS:
         outputs[k] = max(outputs[k], float(np.max(np.abs(py["out"][k] - cap["dp"][60][k]))))
   idl_lut = next((IDL_RUN / f"case_{case}" / "dp_generator" / "luts" / f"case_{case}").glob("*_v90.nc"))
   lut = {}
   with Dataset(py_lut) as p, Dataset(idl_lut) as q:
      for name in RT + OPTICS + ("F0",):
         if name not in p.variables or name not in q.variables:
            continue
         x = np.asarray(p[name][:], np.float64)
         y = np.asarray(q[name][:], np.float64)
         if x.shape != y.shape:
            raise AssertionError(f"case {case}: {name} shapes {x.shape} {y.shape}")
         d = np.abs(x - y)
         rel = np.where(np.abs(y) >= 1e-3, d / np.maximum(np.abs(y), 1e-30), 0.0)
         lut[name] = {"bitwise_fraction": float(np.mean(x == y)), "max_abs": float(d.max()), "max_rel_where_ge_1e-3": float(rel.max()),
                      "max_ulp": ulp32(x, y)}
   return {"calls": len(calls), "python_lut": str(py_lut), "idl_lut": str(idl_lut),
           "disort_input_max_ulp": inputs, "disort_scalars_identical": bool(scalars_equal),
           "disort_output_max_abs": outputs, "lut": lut}


def main():
   output = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "compact_python_vs_idl.json"
   precision = os.environ.get("ORAC_LUT_DISORT_PRECISION", "double") or "double"
   report = {"disort_precision": precision, "idl_reference_run": str(IDL_RUN), "cases": {}}
   for case, channel in CASES.items():
      report["cases"][case] = compare_case(case, channel)
      print(case, json.dumps(report["cases"][case]["lut"]["R_0v"]), flush=True)
   output.write_text(json.dumps(report, indent=1) + "\n")
   print(f"wrote {output}")


if __name__ == "__main__":
   main()
