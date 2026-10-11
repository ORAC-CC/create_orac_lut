"""Capture the DISORT calls of the V22 Python generator in single precision (pre-call copies).

Runs the active Python generator (create_orac_luts.py, TOTAL-corrected
revisions; legacy_disort.py with ORAC_LUT_DISORT_PRECISION=single, i.e. the
historical gfortran single-precision kernel) on the four compact Terra MODIS
p240 cases of validation/v21_v22_terra_modis/idl_double/cases and saves, for
every DISORT call, copies of all arguments taken *before* the call (the
historical kernel dithers SSALB = 1 in place) and the outputs.

Usage: ORAC_LUT_DISORT_PRECISION=single python capture_python_single.py
Writes results/python_single_case_<X>.npz and results/python_single_sources.json.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASES_DIR = HERE.parent / "idl_double" / "cases"
WORK = HERE / "work"
RESULTS = HERE / "results"
CASES = {"A": 6, "B": 20, "C": 6, "D": 1}
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))


def sha256(path):
   return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stage(case):
   in_path = WORK / f"case_{case}" / "input_files"
   if in_path.parent.exists():
      shutil.rmtree(in_path.parent)
   (in_path / "lut").mkdir(parents=True)
   for entry in (ROOT / "create_orac_lut" / "input_files").iterdir():
      if entry.name != "lut":
         (in_path / entry.name).symlink_to(entry)
   shutil.copy(CASES_DIR / f"case_{case}.lut", in_path / "lut" / f"case_{case}.lut")
   return in_path


def main():
   if os.environ.get("ORAC_LUT_DISORT_PRECISION") != "single":
      sys.exit("set ORAC_LUT_DISORT_PRECISION=single explicitly (the production default is double)")
   RESULTS.mkdir(parents=True, exist_ok=True)
   production = importlib.import_module("create_orac_luts")
   call_module = importlib.import_module("oraclut.idl_mirror.call_disort")
   legacy = importlib.import_module("oraclut.radiative_transfer.legacy_disort")
   assert legacy.disort_precision() == "single"
   original = call_module.disort
   sources = {name: {"path": str(Path(m.__file__).relative_to(ROOT)), "sha256": sha256(m.__file__)}
              for name, m in sorted(sys.modules.items()) if name.startswith("oraclut") and getattr(m, "__file__", None)}
   sources["create_orac_luts"] = {"path": "create_orac_luts.py", "sha256": sha256(ROOT / "create_orac_luts.py")}
   sources["single_precision_library"] = {"path": str(legacy._library_path().relative_to(ROOT)), "sha256": sha256(legacy._library_path())}
   (RESULTS / "python_single_sources.json").write_text(json.dumps(sources, indent=1) + "\n")
   for case, channel in CASES.items():
      calls = []

      def capture(dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kw):
         record = {"dtauc": np.array(dtauc, np.float32), "ssalb": np.array(ssalb, np.float32),
                   "pmom": np.array(pmom, np.float32), "utau": np.array(utau, np.float32),
                   "umu": np.array(umu, np.float32), "phi": np.array(phi, np.float32),
                   "fbeam": np.float32(fbeam), "umu0": np.float32(umu0), "fisot": np.float32(fisot),
                   "plank": int(bool(kw.get("plank"))), "nstr": int(kw["nstreams"]),
                   "wvnmlo": np.float32(kw.get("wavenumber_low") or 0.0),
                   "wvnmhi": np.float32(kw.get("wavenumber_high") or 0.0),
                   "temper": np.array(kw["temperature"] if kw.get("temperature") is not None else np.zeros(np.size(dtauc)), np.float32)}
         out = original(dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kw)
         record.update({"out_" + k: np.array(v) for k, v in out.items()})
         calls.append(record)
         return out

      in_path = stage(case)
      call_module.disort = capture
      try:
         production.create_orac_cloud_lut(in_path, "terra_modis_v1.inst", "liquid-water_240.mm", f"case_{case}.lut",
                                          WORK / f"case_{case}" / "luts", "2", channelid=[channel], gas=0, no_rayleigh=0,
                                          srf_quad=1, reuse_scat=0, scat_only=0, tmatrix_path=None, version=22,
                                          driver=None, nstreams=60, nmom=1000, work_path=WORK / f"case_{case}" / "work")
      finally:
         call_module.disort = original
      flat = {f"c{i:03d}_{k}": v for i, record in enumerate(calls) for k, v in record.items()}
      np.savez(RESULTS / f"python_single_case_{case}.npz", ncalls=len(calls), **flat)
      print(f"case {case}: {len(calls)} calls captured", flush=True)


if __name__ == "__main__":
   main()
