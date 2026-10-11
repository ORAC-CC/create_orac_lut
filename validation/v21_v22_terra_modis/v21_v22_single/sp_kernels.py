"""Call the historical single-precision DISORT kernels from Python with controlled memory layout.

Kernels (neither is rebuilt or modified):
  intel_dlm  create_orac_lut/disort2/src/DISORTtoIDL.so  -- the original V21 IDL
             DLM (ifort -O3, unknown ifort version); its Fortran entry disort_
             is called directly (the IDL wrapper only forwards pointers)
  gfortran   build/oraclut/liboraclut_disort_gfortran.so -- the historical V22
             Python kernel (gfortran 12.2 -O2), entry disort_

Every array argument is placed in a fresh buffer at a chosen byte offset
modulo 64.  IDL array data starts 8 bytes past a 16-byte boundary; the
Intel -O3 kernel's results depend on the alignment of UMU (vectorised loop
peeling), so IDL behaviour is reproduced with offset 8.
Requires the intel-compilers/2022 runtime on LD_LIBRARY_PATH for intel_dlm.
"""

import ctypes
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
KERNELS = {"intel_dlm": ROOT / "create_orac_lut/disort2/src/DISORTtoIDL.so",
           "gfortran": ROOT / "build/oraclut/liboraclut_disort_gfortran.so"}
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")
_LOADED = {}


def _library(name):
   if name not in _LOADED:
      library = ctypes.CDLL(str(KERNELS[name]), mode=os.RTLD_LAZY)
      library.disort_.restype = ctypes.c_int32
      _LOADED[name] = library
   return _LOADED[name]


def _placed(values, offset, keep, order="C"):
   values = np.asarray(values, np.float32)
   buffer = np.zeros(values.size * 4 + 128, np.uint8)
   start = (64 - buffer.ctypes.data % 64) % 64 + offset
   array = np.ndarray(values.shape, np.float32, buffer, start, order=order)
   array[...] = values
   keep.append(buffer)
   return array


def call(kernel, args, nstr=60, offset=8, offsets=None):
   """Single-precision DISORT with the fixed ORAC settings; args as in the capture loaders."""
   keep = []
   n, nmom = int(args["nlyr"]), int(args["nmom"])
   offsets = offsets or {}
   p = lambda name, values, order="C": _placed(values, offsets.get(name, offset), keep, order)
   dtauc = p("dtauc", np.ravel(args["dtauc"])[:n])
   ssalb = p("ssalb", np.ravel(args["ssalb"])[:n])
   pmom = p("pmom", np.asarray(args["pmom"])[:nmom + 1, :n], "F")
   temper = p("temper", np.asarray(args.get("temper", np.zeros(n + 1)), np.float32).ravel()[:n + 1])
   utau, umu, phi = p("utau", args["utau"]), p("umu", args["umu"]), p("phi", args["phi"])
   ntau, numu, nphi = utau.size, umu.size, phi.size
   out = {k: p(k, np.zeros(ntau)) for k in ("rfldir", "rfldn", "flup", "dfdt", "uavg")}
   out["uu"] = p("uu", np.zeros((numu, ntau, nphi)), "F")
   out["albmed"], out["trnmed"] = p("albmed", np.zeros(numu)), p("trnmed", np.zeros(numu))
   i = lambda v: ctypes.byref(ctypes.c_int32(int(v)))
   r = lambda v: ctypes.byref(ctypes.c_float(float(np.float32(v))))
   q = lambda a: a.ctypes.data_as(ctypes.c_void_p)
   prnt = np.zeros(5, np.int32)
   status = _library(kernel).disort_(
      i(n), q(dtauc), q(ssalb), i(nmom), q(pmom), q(temper), r(args.get("wvnmlo", 0)), r(args.get("wvnmhi", 0)),
      i(1), i(ntau), q(utau), i(nstr), i(1), i(numu), q(umu), i(nphi), q(phi), i(0), r(args["fbeam"]), r(args["umu0"]),
      r(args.get("phi0", 0)), r(args["fisot"]), i(1), r(0), r(0), r(0), r(0), i(int(args.get("plank", 0))), i(0),
      r(args.get("accur", 1e-8)), q(prnt), i(n), i(ntau), i(numu), i(nphi), i(nmom), *(q(out[k]) for k in OUTPUTS))
   result = {k: np.array(v) for k, v in out.items()}
   return int(status), result
