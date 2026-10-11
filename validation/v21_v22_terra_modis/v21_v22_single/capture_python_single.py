"""Capture the DISORT calls of the V22 Python generator in single precision (pre-call copies).

Runs the active Python generator (create_orac_luts.py, TOTAL-corrected
revisions; legacy_disort.py with ORAC_LUT_DISORT_PRECISION=single, i.e. the
historical gfortran single-precision kernel) on the four compact Terra MODIS
p240 cases of validation/v21_v22_terra_modis/idl_double/cases and saves, for
every DISORT call, copies of all arguments taken *before* the call (the
historical kernel dithers SSALB = 1 in place) and the outputs.

Usage: ORAC_LUT_DISORT_PRECISION=single python capture_python_single.py [tag]
Writes results/python_single_case_<X>.npz and results/python_single_sources.json,
or, with a tag, results/<tag>/... and work_<tag>/case_<X>/ (the untagged files of
the 17de4c0 baseline are kept).  Arguments are recorded at the entry of the
selected single-precision kernel backend (legacy_disort._disort_single_intel
or _disort_single), after the float32 conversion the backend applies; the
dtypes the backend received are recorded too.
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
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
WORK = HERE / (f"work_{TAG}" if TAG else "work")
RESULTS = HERE / "results" / TAG if TAG else HERE / "results"
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
   legacy = importlib.import_module("oraclut.radiative_transfer.legacy_disort")
   assert legacy.disort_precision() == "single"
   backend_name = "_disort_single_intel" if legacy.single_kernel() == "intel_dlm" else "_disort_single"
   original = getattr(legacy, backend_name)
   sources = {name: {"path": str(Path(m.__file__).relative_to(ROOT)), "sha256": sha256(m.__file__)}
              for name, m in sorted(sys.modules.items()) if name.startswith("oraclut") and getattr(m, "__file__", None)}
   sources["create_orac_luts"] = {"path": "create_orac_luts.py", "sha256": sha256(ROOT / "create_orac_luts.py")}
   sources["single_precision_library"] = {"path": str(legacy.single_library().relative_to(ROOT)), "sha256": sha256(legacy.single_library())}
   mie = importlib.import_module("oraclut.optics.legacy_mie")
   sources["mie_library"] = {"path": str(mie._library_path().relative_to(ROOT)), "sha256": sha256(mie._library_path())}
   sources["kernels"] = {"single": legacy.single_kernel(), "backend": backend_name}
   (RESULTS / "python_single_sources.json").write_text(json.dumps(sources, indent=1) + "\n")
   for case, channel in CASES.items():
      calls = []

      def capture(dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kw):
         record = {"dtypes": np.array(",".join(str(np.asarray(v).dtype) for v in (dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot))),"dtauc": np.array(dtauc, np.float32), "ssalb": np.array(ssalb, np.float32),
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
      setattr(legacy, backend_name, capture)
      try:
         production.create_orac_cloud_lut(in_path, "terra_modis_v1.inst", "liquid-water_240.mm", f"case_{case}.lut",
                                          WORK / f"case_{case}" / "luts", "2", channelid=[channel], gas=0, no_rayleigh=0,
                                          srf_quad=1, reuse_scat=0, scat_only=0, tmatrix_path=None, version=22,
                                          driver=None, nstreams=60, nmom=1000, work_path=WORK / f"case_{case}" / "work")
      finally:
         setattr(legacy, backend_name, original)
      flat = {f"c{i:03d}_{k}": v for i, record in enumerate(calls) for k, v in record.items()}
      np.savez(RESULTS / f"python_single_case_{case}.npz", ncalls=len(calls), **flat)
      print(f"case {case}: {len(calls)} calls captured", flush=True)


if __name__ == "__main__":
   main()
