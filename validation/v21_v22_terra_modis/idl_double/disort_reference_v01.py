"""Run one captured set of DISORT arguments through several solver builds.

Builds (all from the unmodified production DISORTfunctions.f):
  ifort_dp    build/libdisort_dp_ifort.so  - the kernel inside DISORTtoIDL_DP.so
  gfort_dp    build/ref_double/disort_ref_driver
  gfort_qp    build/ref_quad/disort_ref_driver   (REAL*16)
  sp_dlm      create_orac_lut/disort2/src/DISORTtoIDL.so - historical V21 kernel

``args`` is a dict with nlyr, nmom, dtauc, ssalb, pmom (nmom+1, nlyr), temper
(nlyr+1), utau, umu, phi, fbeam, umu0, phi0, fisot, accur, wvnmlo, wvnmhi,
btemp, ttemp, temis, albedo and plank.  Only the stream count is varied.
"""

from __future__ import annotations

import ctypes
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BUILD = HERE / "build"
LIBRARIES = {
   "ifort_dp": BUILD / "libdisort_dp_ifort.so",
   "sp_dlm": ROOT / "create_orac_lut" / "disort2" / "src" / "DISORTtoIDL.so",
}
DRIVERS = {"gfort_dp": BUILD / "ref_double" / "disort_ref_driver",
           "gfort_qp": BUILD / "ref_quad" / "disort_ref_driver"}
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")
SCALARS = ("fbeam", "umu0", "phi0", "fisot", "accur", "wvnmlo", "wvnmhi", "btemp", "ttemp", "temis", "albedo")


def _mantissa_exponent(value):
   value = float(value)
   if value == 0.0:
      return 0, 0
   fraction, exponent = math.frexp(value)
   mantissa = int(fraction * 2 ** 53)
   assert math.ldexp(mantissa, exponent - 53) == value
   return mantissa, exponent - 53


def _write_input(path, args, nstr):
   nlyr, nmom = int(args["nlyr"]), int(args["nmom"])
   lines = [f"{nlyr} {nmom} {len(args['utau'])} {len(args['umu'])} {len(args['phi'])} {nstr} {int(bool(args['plank']))}"]
   values = [args[name] for name in SCALARS]
   values += list(np.ravel(args["dtauc"])[:nlyr]) + list(np.ravel(args["ssalb"])[:nlyr])
   values += list(np.asarray(args["pmom"])[:nmom + 1, :nlyr].ravel(order="F"))
   values += list(np.ravel(args["temper"])[:nlyr + 1]) + list(args["utau"]) + list(args["umu"]) + list(args["phi"])
   lines += ["%d %d" % _mantissa_exponent(v) for v in values]
   Path(path).write_text("\n".join(lines) + "\n")


def _read_output(path, args):
   lines = Path(path).read_text().split("\n")
   header = lines[0].split()
   status = int(header[2])
   out, i = {}, 1
   while i < len(lines) and lines[i].strip():
      name, n = lines[i].split()
      n = int(n)
      out[name] = np.array([float(x) for x in lines[i + 1:i + 1 + n]])
      i += 1 + n
   numu, ntau, nphi = len(args["umu"]), len(args["utau"]), len(args["phi"])
   out["uu"] = out["uu"].reshape((numu, ntau, nphi), order="F")
   return status, out, int(header[4])


def run_driver(variant, args, nstr):
   """Run a gfortran reference driver; returns (status, outputs, real_bytes)."""
   with tempfile.TemporaryDirectory() as work:
      infile, outfile = Path(work) / "in.txt", Path(work) / "out.txt"
      _write_input(infile, args, nstr)
      completed = subprocess.run([str(DRIVERS[variant]), str(infile), str(outfile)], capture_output=True, text=True,
                                 env={"PATH": "/usr/bin:/bin"}, timeout=3600)
      if completed.returncode != 0:
         raise RuntimeError(f"{variant} failed: {completed.stdout} {completed.stderr}")
      status, out, real_bytes = _read_output(outfile, args)
      out["_stdout"] = completed.stdout
      return status, out, real_bytes


_LOADED = {}


def run_library(variant, args, nstr):
   """Call disort_ of a shared library directly; double or float per variant."""
   if variant not in _LOADED:
      library = ctypes.CDLL(str(LIBRARIES[variant]), mode=os.RTLD_LAZY)
      library.disort_.restype = ctypes.c_int32
      _LOADED[variant] = library
   library = _LOADED[variant]
   real = np.float64 if variant == "ifort_dp" else np.float32
   creal = ctypes.c_double if variant == "ifort_dp" else ctypes.c_float
   nlyr, nmom = int(args["nlyr"]), int(args["nmom"])
   ntau, numu, nphi = len(args["utau"]), len(args["umu"]), len(args["phi"])
   arr = lambda a: np.ascontiguousarray(a, dtype=real)
   dtauc, ssalb = arr(np.ravel(args["dtauc"])[:nlyr]), arr(np.ravel(args["ssalb"])[:nlyr])
   pmom = np.asfortranarray(np.asarray(args["pmom"])[:nmom + 1, :nlyr], dtype=real)
   temper = arr(np.ravel(args["temper"])[:nlyr + 1])
   utau, umu, phi = arr(args["utau"]), arr(args["umu"]), arr(args["phi"])
   out = {k: np.zeros(ntau, dtype=real) for k in ("rfldir", "rfldn", "flup", "dfdt", "uavg")}
   out["uu"] = np.zeros((numu, ntau, nphi), dtype=real, order="F")
   out["albmed"], out["trnmed"] = np.zeros(numu, dtype=real), np.zeros(numu, dtype=real)
   i32 = lambda v: ctypes.byref(ctypes.c_int32(int(v)))
   rv = lambda v: ctypes.byref(creal(float(v)))
   p = lambda a: a.ctypes.data_as(ctypes.c_void_p)
   prnt = np.zeros(5, dtype=np.int32)
   status = library.disort_(
      i32(nlyr), p(dtauc), p(ssalb), i32(nmom), p(pmom), p(temper), rv(args["wvnmlo"]), rv(args["wvnmhi"]),
      i32(1), i32(ntau), p(utau), i32(nstr), i32(1), i32(numu), p(umu), i32(nphi), p(phi), i32(0),
      rv(args["fbeam"]), rv(args["umu0"]), rv(args["phi0"]), rv(args["fisot"]), i32(1), rv(args["albedo"]),
      rv(args["btemp"]), rv(args["ttemp"]), rv(args["temis"]), i32(int(bool(args["plank"]))), i32(0),
      rv(args["accur"]), p(prnt), i32(nlyr), i32(ntau), i32(numu), i32(nphi), i32(nmom),
      *(p(out[k]) for k in OUTPUTS))
   return int(status), {k: np.asarray(v, dtype=np.float64) for k, v in out.items()}


def run(variant, args, nstr):
   if variant in DRIVERS:
      status, out, _ = run_driver(variant, args, nstr)
      return status, out
   return run_library(variant, args, nstr)


def args_from_v21_capture(sav_path):
   """Arguments of the existing V21 single-precision capture (disort_boundary)."""
   from scipy.io import readsav
   i = readsav(str(sav_path))
   n = int(i["nlyr1"])
   temper = np.zeros(n + 1)
   given = np.asarray(i["temper"], dtype=np.float64).ravel()
   temper[:min(n + 1, given.size)] = given[:n + 1]
   return {
      "nlyr": n, "nmom": int(i["nmom"]),
      "dtauc": np.asarray(i["dtauc1"], np.float64)[:n], "ssalb": np.asarray(i["ssalb1"], np.float64)[:n],
      "pmom": np.asarray(i["pmom1"], np.float64)[:n].T, "temper": temper,
      "utau": np.asarray(i["utau"], np.float64), "umu": np.asarray(i["umu"], np.float64),
      "phi": np.asarray(i["phi"], np.float64), "plank": int(i["plank"]),
      **{name: float(i[name]) for name in SCALARS},
   }
