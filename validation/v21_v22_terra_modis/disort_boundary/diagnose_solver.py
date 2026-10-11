"""Separate DISORT input differences from solver-build differences.

The captured V21 IDL and corrected V22 Python DISORT arguments for the selected
p240 direct-beam call are passed directly to the Fortran ``disort_`` entry point
of one or more compiled libraries, with every one of the 44 arguments under
explicit control (including MAXCLY padding as used by IDL).  Each result is
compared with the captured V21 and V22 outputs.

Usage: python diagnose_solver.py OUTPUT_JSON LIBRARY [LIBRARY ...]

DISORT warnings are written by the Fortran to stdout and are left there.
"""

from __future__ import annotations

import ctypes
import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.io import readsav


CASE = Path(__file__).resolve().parent
RESULTS = CASE / "results"
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")


def load_v21():
   capture = RESULTS / "v21_capture_latest"
   i = readsav(str(capture / "v21_disort_boundary_input.sav"))
   o = readsav(str(capture / "v21_disort_boundary_output.sav"))
   n = int(i["nlyr1"])
   maxcly = int(i["maxcly1"])
   inputs = {
      "nlyr": n, "maxcly_padded": maxcly,
      "dtauc": np.asarray(i["dtauc1"], dtype=np.float32),
      "ssalb": np.asarray(i["ssalb1"], dtype=np.float32),
      # readsav reverses IDL dimension order: (layer, moment) -> (moment, layer)
      "pmom": np.asfortranarray(np.asarray(i["pmom1"], dtype=np.float32).T),
      "temper": np.asarray(i["temper"], dtype=np.float32),
      "utau": np.asarray(i["utau"], dtype=np.float32),
      "umu": np.asarray(i["umu"], dtype=np.float32),
      "phi": np.asarray(i["phi"], dtype=np.float32),
      "nmom": int(i["nmom"]), "nstr": int(i["nstr"]),
      "fbeam": np.float32(i["fbeam"]), "umu0": np.float32(i["umu0"]),
      "fisot": np.float32(i["fisot"]), "accur": np.float32(i["accur"]),
   }
   # IDL UU is (MAXUMU, MAXULV, MAXPHI); readsav returns it reversed.
   outputs = {k: np.asarray(o[k], dtype=np.float32) for k in OUTPUTS}
   outputs["uu"] = np.asfortranarray(outputs["uu"].T)
   return inputs, outputs


def load_v22():
   p = np.load(RESULTS / "corrected_v22_disort_boundary.npz")
   n = int(p["nlyr"])
   inputs = {
      "nlyr": n, "maxcly_padded": n,
      "dtauc": p["dtauc"], "ssalb": p["ssalb"], "pmom": np.asfortranarray(p["pmom"]),
      "temper": p["temper"], "utau": p["utau"], "umu": p["umu"], "phi": p["phi"],
      "nmom": int(p["nmom"]), "nstr": int(p["nstr"]),
      "fbeam": p["fbeam"][()], "umu0": p["umu0"][()], "fisot": p["fisot"][()],
      "accur": p["accur"][()],
   }
   # The capture recorded SSALB after DISORT had dithered SSALB = 1 in place
   # (DISORTfunctions.f: SSALB(LC) = 1.0 - DITHER); restore the pre-call value.
   dither = np.float32(1.0) - np.float32(10.0) * np.finfo(np.float32).eps
   inputs["ssalb"] = np.where(inputs["ssalb"] == dither, np.float32(1.0), inputs["ssalb"]).astype(np.float32)
   outputs = {k: np.asarray(p[k], dtype=np.float32) for k in OUTPUTS}
   return inputs, outputs


def call(library, inputs, padded):
   """Call disort_ with IDL-equivalent scalar settings and chosen MAXCLY."""
   n = inputs["nlyr"]
   maxcly = inputs["maxcly_padded"] if padded else n
   nmom = inputs["nmom"]
   ntau = inputs["utau"].size
   numu = inputs["umu"].size
   nphi = inputs["phi"].size

   dtauc = np.zeros(maxcly, dtype=np.float32)
   ssalb = np.zeros(maxcly, dtype=np.float32)
   pmom = np.zeros((nmom + 1, maxcly), dtype=np.float32, order="F")
   dtauc[:n] = inputs["dtauc"][:n]
   ssalb[:n] = inputs["ssalb"][:n]
   pmom[:, :n] = inputs["pmom"][:, :n]
   # TEMPER is only read when PLANK is set; size it generously.
   temper = np.zeros(maxcly + 1, dtype=np.float32)

   i32 = lambda v: ctypes.byref(ctypes.c_int32(int(v)))
   f32 = lambda v: ctypes.byref(ctypes.c_float(float(v)))
   ptr = lambda a: a.ctypes.data_as(ctypes.c_void_p)

   rfldir = np.zeros(ntau, dtype=np.float32)
   rfldn = np.zeros(ntau, dtype=np.float32)
   flup = np.zeros(ntau, dtype=np.float32)
   dfdt = np.zeros(ntau, dtype=np.float32)
   uavg = np.zeros(ntau, dtype=np.float32)
   uu = np.zeros((numu, ntau, nphi), dtype=np.float32, order="F")
   albmed = np.zeros(numu, dtype=np.float32)
   trnmed = np.zeros(numu, dtype=np.float32)
   utau = np.ascontiguousarray(inputs["utau"], dtype=np.float32)
   umu = np.ascontiguousarray(inputs["umu"], dtype=np.float32)
   phi = np.ascontiguousarray(inputs["phi"], dtype=np.float32)
   prnt = np.zeros(5, dtype=np.int32)

   status = library.disort_(
      i32(n), ptr(dtauc), ptr(ssalb), i32(nmom), ptr(pmom), ptr(temper),
      f32(0.0), f32(0.0), i32(1), i32(ntau), ptr(utau), i32(inputs["nstr"]),
      i32(1), i32(numu), ptr(umu), i32(nphi), ptr(phi), i32(0),
      f32(inputs["fbeam"]), f32(inputs["umu0"]), f32(0.0), f32(inputs["fisot"]),
      i32(1), f32(0.0), f32(0.0), f32(0.0), f32(0.0), i32(0), i32(0),
      f32(inputs["accur"]), ptr(prnt), i32(maxcly), i32(ntau), i32(numu),
      i32(nphi), i32(nmom), ptr(rfldir), ptr(rfldn), ptr(flup), ptr(dfdt),
      ptr(uavg), ptr(uu), ptr(albmed), ptr(trnmed),
   )
   return status, {"rfldir": rfldir, "rfldn": rfldn, "flup": flup, "dfdt": dfdt,
                   "uavg": uavg, "uu": uu, "albmed": albmed, "trnmed": trnmed}


def call_double(library_path, inputs):
   """Double-precision reference through the validation adapter (production options)."""
   library = ctypes.CDLL(str(library_path), mode=os.RTLD_LAZY)
   library.oraclut_set_options.argtypes = [ctypes.c_double, ctypes.c_int32, ctypes.c_int32]
   library.oraclut_set_options(float(inputs["accur"]), 1, 1)
   n = inputs["nlyr"]
   nmom = inputs["nmom"]
   f64 = lambda a: np.ascontiguousarray(a, dtype=np.float64)
   dtauc, ssalb = f64(inputs["dtauc"][:n]), f64(inputs["ssalb"][:n])
   pmom = np.asfortranarray(inputs["pmom"][:, :n], dtype=np.float64)
   utau, umu, phi = f64(inputs["utau"]), f64(inputs["umu"]), f64(inputs["phi"])
   ntau, numu, nphi = utau.size, umu.size, phi.size
   temper = np.zeros(n + 1)
   out = {k: np.zeros(ntau) for k in ("rfldir", "rfldn", "flup", "dfdt", "uavg")}
   out["uu"] = np.zeros((numu, ntau, nphi), order="F")
   out["albmed"], out["trnmed"] = np.zeros(numu), np.zeros(numu)
   ptr = lambda a: a.ctypes.data_as(ctypes.c_void_p)
   d = ctypes.c_double
   i = ctypes.c_int32
   status = library.oraclut_disort(
      i(n), ptr(dtauc), ptr(ssalb), i(nmom), ptr(pmom), ptr(temper), d(0.0), d(0.0),
      i(0), i(ntau), ptr(utau), i(inputs["nstr"]), i(numu), ptr(umu), i(nphi), ptr(phi),
      d(float(inputs["fbeam"])), d(float(inputs["umu0"])), d(float(inputs["fisot"])),
      *(ptr(out[k]) for k in OUTPUTS),
   )
   return status, out


def compare(a, b):
   """Bitwise count and maximum absolute difference over all outputs."""
   summary = {}
   for k in OUTPUTS:
      x = np.asarray(a[k], dtype=np.float32)
      y = np.asarray(b[k], dtype=np.float32)
      same = x.view(np.int32) == y.view(np.int32)
      summary[k] = {
         "bitwise_equal": int(same.sum()), "size": int(x.size),
         "max_abs_diff": float(np.max(np.abs(x.astype(np.float64) - y))),
      }
   return summary


def toa_upwelling(out):
   """UU at UTAU = 0 for upward directions (UMU > 0): the quantity behind RBD."""
   return np.asarray(out["uu"], dtype=np.float64)[10:, 0, :]


def attribution(library, v21_in, v22_in, v21_out, reference):
   """Swap one V21 input group at a time for its V22 value."""
   groups = {"dtauc+utau": ("dtauc", "utau"), "pmom": ("pmom",), "umu": ("umu",),
             "all": ("dtauc", "utau", "pmom", "umu", "ssalb")}
   rows = []
   for name, keys in groups.items():
      inputs = dict(v21_in)
      for k in keys:
         inputs[k] = v22_in[k]
      _, out = call(library, inputs, False)
      rows.append({
         "swapped_to_v22": name,
         "max_abs_toa_up_vs_v21": float(np.max(np.abs(toa_upwelling(out) - toa_upwelling(v21_out)))),
         "max_abs_toa_up_vs_double": float(np.max(np.abs(toa_upwelling(out) - toa_upwelling(reference)))),
         "flup_top": float(out["flup"][0]),
      })
   return rows


def main(output, libraries):
   v21_in, v21_out = load_v21()
   v22_in, v22_out = load_v22()
   report = {"v21_vs_v22_captured": compare(v21_out, v22_out), "runs": [], "double": []}
   double = os.environ.get("ORACLUT_DISORT_DOUBLE")
   if double:
      refs = {}
      for name, inputs in (("v21_inputs", v21_in), ("v22_inputs", v22_in)):
         status, out = call_double(double, inputs)
         refs[name] = out
         report["double"].append({
            "library": double, "inputs": name, "status": int(status),
            "top": {"rfldn": float(out["rfldn"][0]), "flup": float(out["flup"][0]),
                    "uavg": float(out["uavg"][0])},
            "vs_v21_captured": compare(out, v21_out),
            "vs_v22_captured": compare(out, v22_out),
         })
      report["double_v21_vs_v22_inputs"] = {
         k: float(np.max(np.abs(refs["v21_inputs"][k] - refs["v22_inputs"][k]))) for k in OUTPUTS}
      np.savez(RESULTS / "double_reference_disort_boundary.npz",
               **{f"{n}_{k}": v for n, r in refs.items() for k, v in r.items()})
   for lib_path in libraries:
      library = ctypes.CDLL(str(lib_path), mode=os.RTLD_LAZY)
      library.disort_.restype = ctypes.c_int32
      for name, inputs in (("v21_inputs", v21_in), ("v22_inputs", v22_in)):
         for padded in (False, True):
            if padded and inputs["maxcly_padded"] == inputs["nlyr"]:
               continue
            status, out = call(library, inputs, padded)
            run = {
               "library": str(lib_path), "inputs": name,
               "maxcly": inputs["maxcly_padded"] if padded else inputs["nlyr"],
               "status": int(status),
               "top": {"rfldn": float(out["rfldn"][0]), "flup": float(out["flup"][0]),
                       "uavg": float(out["uavg"][0])},
               "vs_v21_captured": compare(out, v21_out),
               "vs_v22_captured": compare(out, v22_out),
            }
            if double:
               run["max_abs_toa_up_vs_double"] = float(np.max(np.abs(
                  toa_upwelling(out) - toa_upwelling(refs["v21_inputs"]))))
            report["runs"].append(run)
      if double:
         report.setdefault("attribution", {})[str(lib_path)] = attribution(
            library, v21_in, v22_in, v21_out, refs["v21_inputs"])
   Path(output).write_text(json.dumps(report, indent=1))
   return 0


if __name__ == "__main__":
   raise SystemExit(main(sys.argv[1], sys.argv[2:]))
