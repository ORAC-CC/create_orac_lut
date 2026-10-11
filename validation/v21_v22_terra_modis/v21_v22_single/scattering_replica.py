"""Diagnostic replica of the V21 IDL scattering-property arithmetic, with switchable alternatives.

Reproduces, operation by operation, the IDL chain

   load_mmdat (CM = complex(n, -k))  ->  INTERPOL (refractive index)
   generate_scattering_properties.pro  ->  create_bwgp.pro  ->  mie_size_dist_new.pro
   (quadrature.pro 'T' / 'G', shift_quadrature.pro, Mie_dlm_single -> mieint)
   ->  mixing  ->  legpexp.pro

for a single modified-gamma Mie component (the liquid-water p240 / pold
cases), so that each candidate source of an IDL/Python difference can be
switched on or off and the result compared with the IDL scatfile.sav.  Every
option has an "idl" value (the IDL arithmetic as read from the .pro sources and
the IDL 8.9 library) and a "python" value (what the production Python
oraclut.idl_mirror code does at the time of writing).

This is validation code: production never imports it.

Usage: python scattering_replica.py            (runs the option survey on cases A-D)
"""

import ctypes
import itertools
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from scipy.io import readsav

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from oraclut.idl_mirror import load_inststr, load_lutstr, load_srfstrarr                  # noqa: E402
import importlib                                                                          # noqa: E402
# The "python" options are the uncorrected V22 arithmetic: the SourceVault copies
# preserved before the Stage 4 corrections (identical to the 17de4c0 sources).
production_bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp_v02")          # noqa: E402
_interpol_complex = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties_v02")._interpol_complex   # noqa: E402
load_mmdat = importlib.import_module("oraclut.idl_mirror.load_mmdat_v01").load_mmdat      # noqa: E402
from oraclut.idl_total import idl_total                                                  # noqa: E402
from oraclut.optics.legacy_mie import mie_single_batch                                   # noqa: E402
from oraclut.radiative_transfer import legacy_disort                                     # noqa: E402

f32 = np.float32
IDL_RUN = HERE.parent / "idl_double" / "results" / "run_20261010T214143"
CASES = {"A": 6, "B": 20, "C": 6, "D": 1}
INPUTS = ROOT / "create_orac_lut" / "input_files"
MIE_DLM = ROOT / "create_orac_lut" / "mie" / "dlm-code" / "mie_dlm_single.so"

IDL = dict(ri="idl", wn="f32", literals="f32", quad_t="idl", quad_g="idl", w1p="idl", math="libm", gamma="math",
           mie="intel", store="fltarr", mie_dqv_offset=0, mie_f_offset=8)
PYTHON = dict(ri="python", wn="f64", literals="f64", quad_t="python", quad_g="python", w1p="python", math="numpy",
              gamma="math", mie="gfortran", store="f64", mie_dqv_offset=0, mie_f_offset=0)


# ---------------------------------------------------------------------------
# Kernels and library functions
# ---------------------------------------------------------------------------

_LIBM = ctypes.CDLL("libm.so.6")
_LIBM.tgamma.restype = ctypes.c_double
_LIBM.tgamma.argtypes = [ctypes.c_double]
_MIE = None


def _intel_mie():
   global _MIE
   if _MIE is None:
      legacy_disort._load_intel_library()                 # loads the Intel 2022 run time (RTLD_GLOBAL)
      runtime = legacy_disort.intel_runtime_directory()
      ctypes.CDLL(str(runtime / "libiomp5.so"), mode=os.RTLD_LAZY | os.RTLD_GLOBAL)
      _MIE = ctypes.CDLL(str(MIE_DLM), mode=os.RTLD_LAZY | os.RTLD_LOCAL)
   return _MIE


def mie_gfortran(dx, cm, dqv):
   """The gfortran build of mieint used by the production generator (oraclut.optics.legacy_mie)."""
   return mie_single_batch(dx, complex(cm), dqv)


def _placed(shape, offset, dtype=np.float64):
   """Zeroed array whose data start ``offset`` bytes past a 64-byte boundary."""
   nbytes = int(np.prod(shape)) * np.dtype(dtype).itemsize
   buffer = np.zeros(nbytes + 128, np.uint8)
   start = (64 - buffer.ctypes.data % 64) % 64 + offset
   return np.ndarray(shape, dtype, buffer, start)


def mie_intel(dx, cm, dqv, dqv_offset, f_offset):
   """mieint_ of the original Intel IDL DLM, called as mie_dlm_single.c calls it (one thread)."""
   lib = _intel_mie()
   inp = dqv.size
   dqvdata = _placed((inp,), dqv_offset)                 # the DLM's static Dqvdata[10000] (64-aligned .bss)
   dqvdata[:] = dqv
   f11 = _placed((dx.size, inp), f_offset)                # IDL temporary arrays F11ARR..F34ARR [Inp, Npts]
   f33, f12, f34 = (_placed((dx.size, inp), f_offset) for _ in range(3))
   s1 = _placed((inp,), 0, np.complex128)                 # malloc'd in the parallel region (16-aligned)
   s2 = _placed((inp,), 0, np.complex128)
   out = {k: np.zeros(dx.size) for k in ("qext", "qsca", "qbsc", "g")}
   cm_v = (ctypes.c_double * 2)(float(np.real(cm)), float(np.imag(cm)))
   inp_v = ctypes.c_int32(inp)
   q = lambda a: a.ctypes.data_as(ctypes.c_void_p)
   for i in range(dx.size):
      dx_v, qext, qsca, qbsc, g, err = (ctypes.c_double(float(dx[i])), ctypes.c_double(), ctypes.c_double(),
                                        ctypes.c_double(), ctypes.c_double(), ctypes.c_int32(0))
      row = lambda a: ctypes.c_void_p(a.ctypes.data + i * inp * 8)
      lib.mieint_(ctypes.byref(dx_v), cm_v, ctypes.byref(inp_v), q(dqvdata), ctypes.byref(qext), ctypes.byref(qsca),
                  ctypes.byref(qbsc), ctypes.byref(g), q(s1), q(s2), row(f11), row(f33), row(f12), row(f34),
                  ctypes.byref(err))
      if err.value:
         raise RuntimeError("Size parameter overflow!")
      out["qext"][i], out["qsca"][i], out["qbsc"][i], out["g"][i] = qext.value, qsca.value, qbsc.value, g.value
   out["f11"] = np.array(f11)
   return out


# ---------------------------------------------------------------------------
# IDL library and project routines
# ---------------------------------------------------------------------------

def quadrature_g_idl(npts):
   """quadrature.pro, QuadType 'G' (Newton iteration from FirstGuess), op for op."""
   n = float(npts)
   abscissa = np.zeros(npts)
   weight = np.zeros(npts)

   def legendre(nn, x):
      pl, pm, pn = 1.0, 1.0, x
      for k in range(2, int(nn) + 1):
         inv = 1.0 / float(k)
         pl = pm
         pm = pn
         pn = (2.0 - inv) * x * pm - (1.0 - inv) * pl
      return pl, pm, pn

   tolerance = 10.0 ** (-14)
   for term, zero in enumerate(range(1, npts + 1)):
      lx = 19262.0
      x = math.cos(math.pi * (zero - 0.25) / (npts + 0.5))
      while abs(x - lx) > tolerance:
         lx = x
         pl, pm, pn = legendre(npts, x)
         x = x - (1.0 - x * x) * pn / (float(npts) * (pm - x * pn))
      abscissa[term] = x
      pl, pm, pn = legendre(n, x)
      weight[term] = 2.0 * (1.0 - x * x) / (n * pm) ** 2
   return abscissa[::-1].copy(), weight[::-1].copy()


def quadrature_t_idl(npts):
   """quadrature.pro 'T': Abscissa = -1D0+2*Dindgen(NPts)/(N-1)."""
   n = float(npts)
   abscissa = -1.0 + 2.0 * np.arange(npts, dtype=np.float64) / (n - 1.0)
   weight = np.full(npts, 2.0 / (n - 1.0))
   weight[0] = 1.0 / (n - 1.0)
   weight[npts - 1] = 1.0 / (n - 1.0)
   return abscissa, weight


def interpol_ri_idl(cm, wl, xout):
   """IDL INTERPOL(CM, WL, XOUT) for COMPLEX CM, DOUBLE WL, FLOAT XOUT (linear, irregular).

   p = (xout-x[s])*diff/(x[s+1]-x[s]) + v[s]: DOUBLE * COMPLEX -> DCOMPLEX, so
   the arithmetic is double complex (multiplication and division by a real
   value act on each part), and the result is stored in a COMPLEX array.
   """
   cm = np.asarray(cm, np.complex64)
   wl = np.asarray(wl, np.float64)
   xout = np.atleast_1d(np.asarray(xout, np.float32)).astype(np.float64)
   m = wl.size
   increasing = wl[-1] > wl[0]
   result = np.zeros(xout.shape, np.complex64)
   for i, x in enumerate(xout.ravel()):
      # VALUE_LOCATE: largest s with x >= wl[s] (increasing) or x <= wl[s] (decreasing), -1 below
      if increasing:
         s = int(np.searchsorted(wl, x, side="right")) - 1
      else:
         s = int(np.searchsorted(-wl, -x, side="right")) - 1
      s = min(max(s, 0), m - 2)
      d = complex(f32(cm[s + 1].real) - f32(cm[s].real), f32(cm[s + 1].imag) - f32(cm[s].imag))   # COMPLEX diff
      t = x - wl[s]
      h = wl[s + 1] - wl[s]
      re = (t * d.real) / h + float(cm[s].real)
      im = (t * d.imag) / h + float(cm[s].imag)
      result.ravel()[i] = np.complex64(complex(f32(re), f32(im)))
   return result.reshape(xout.shape)


# ---------------------------------------------------------------------------
# create_bwgp.pro / mie_size_dist_new.pro
# ---------------------------------------------------------------------------

def gamma(x, option):
   if option == "math":
      return math.gamma(x)
   if option == "tgamma":
      return _LIBM.tgamma(x)
   if option.startswith("ulp"):                         # math.gamma shifted by k ulp (search)
      return float(np.float64(math.gamma(x)) + int(option[3:]) * np.spacing(math.gamma(x)))
   raise ValueError(option)


def mie_size_dist(rm, s, wn, cm, qv, opt):
   """mie_size_dist_new, 'modified_gamma', Nd = 1, params = [Rm, S, 0.001, 100.0], xres = 0.4, /DLM."""
   if opt["literals"] == "f32":
      rl, ru, xres = float(f32(0.001)), float(f32(100.0)), float(f32(0.4))
   else:
      rl, ru, xres = 0.001, 100.0, 0.4
   if 2.0 * math.pi * ru * wn >= 10000.0:
      ru = (10000.0 - 1.0) / (2.0 * math.pi * wn)
   npts = max(int(2.0 * math.pi * (ru - rl) * wn / xres), 200)
   if opt["quad_t"] == "idl":
      absc, wght = quadrature_t_idl(npts)
   else:
      absc, wght = production_bwgp.quadrature("T", npts)
   r = ((rl + ru) + (ru - rl) * absc) / 2.0
   w1 = (ru - rl) * wght / 2.0
   alpha = (1.0 - 3 * s) / s
   b = 1.0 / (rm * s)
   if opt["math"] == "libm":
      n = math.pow(b, -alpha - 1.0) * gamma(alpha + 1.0, opt["gamma"])
      ralpha = np.array([math.pow(x, alpha) for x in r])
      expo = np.array([math.exp(-b * x) for x in r])
   else:
      n = b ** (-alpha - 1.0) * gamma(alpha + 1.0, opt["gamma"])
      ralpha = r ** alpha
      expo = np.exp(-b * r)
   if opt["w1p"] == "idl":
      w1p = w1 * (1.0 / n) * ralpha * expo                # IDL: W1 * (Nd/n) * R^alpha * exp(-b*R)
   else:
      w1p = w1 / n * ralpha * expo
   dx = 2.0 * math.pi * r * wn
   if opt["mie"] == "intel":
      mie = mie_intel(dx, cm, qv, opt["mie_dqv_offset"], opt["mie_f_offset"])
   else:
      mie = mie_gfortran(dx, cm, qv)
   w1pa = w1p * math.pi * (r * r)                        # IDL: W1P * !dpi * R^2
   bext = idl_total(w1pa * mie["qext"])
   bsca = idl_total(w1pa * mie["qsca"])
   g = idl_total(w1pa * mie["g"] * mie["qsca"]) / bsca
   w = bsca / bext
   spm = idl_total(w1pa[:, None] * mie["f11"] * mie["qsca"][:, None], dimension=1) / bsca
   return dict(bext=bext, bsca=bsca, w=w, g=g, spm=spm, npts=npts, r=r, w1p=w1p, dx=dx, mie=mie)


def legpexp(nmom, qv, qw, phase):
   """legpexp.pro (as oraclut.idl_mirror.create_bwgp.legpexp, whose arithmetic order is the IDL's)."""
   return production_bwgp.legpexp(nmom, qv, qw, np.asarray(phase, np.float64))[1]


# ---------------------------------------------------------------------------
# generate_scattering_properties.pro for one channel, one modified-gamma component
# ---------------------------------------------------------------------------

def inputs(case):
   inststr = load_inststr(INPUTS / "inst" / "terra_modis_v1.inst", requestedchannelid=[CASES[case]])
   lutstr = load_lutstr(HERE.parent / "idl_double" / "cases" / f"case_{case}.lut", inststr.max_sat_zenith)
   srfstrarr, nwvl_max = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
   mmstr = load_mmdat(INPUTS / "microphysics" / "liquid-water_240.mm", INPUTS)
   assert nwvl_max == 1 and mmstr.ncomp == 1 and mmstr.distname[0] == "modified_gamma"
   return inststr, lutstr, srfstrarr, mmstr


def scattering(case, opt, nmom=1000):
   inststr, lutstr, srfstrarr, mmstr = inputs(case)
   comp = mmstr.comp[0]
   wvl = f32(srfstrarr[0].wvl[0])
   if opt["ri"] == "idl":
      cm_idl = np.asarray(comp.cm, np.complex128)
      cm_idl = (cm_idl.real.astype(f32) + 1j * cm_idl.imag.astype(f32)).astype(np.complex64)   # complex(n, -k)
      aerm550 = interpol_ri_idl(cm_idl, comp.wl, f32(0.55))[0]
      aerm = interpol_ri_idl(cm_idl, comp.wl, wvl)[0]
   else:
      aerm550 = _interpol_complex(comp.cm, comp.wl, [0.55])[0]
      aerm = _interpol_complex(comp.cm, comp.wl, [float(wvl)])[0]
   if opt["quad_g"] == "idl":
      absc, weights = quadrature_g_idl(nmom)
   else:
      absc, weights = production_bwgp.quadrature("g", nmom)
   qv = ((-1.0 - 1.0) * absc + (1.0 + -1.0)) / 2.0
   s = float(f32(mmstr.s[0])) if opt["literals"] == "f32" else float(mmstr.s[0])
   mrat = float(mmstr.mrat[0])
   store = (lambda v: np.asarray(v, f32)) if opt["store"] == "fltarr" else (lambda v: np.asarray(v, np.float64))
   if opt["wn"] == "f32":
      wn550, wn = float(f32(1.0) / f32(0.55)), float(f32(1.0) / wvl)
   else:
      wn550, wn = 1.0 / 0.55, 1.0 / float(wvl)
   two_n_plus_one = (2.0 * np.arange(nmom, dtype=f32) + f32(1.0))
   out = {k: [] for k in ("bext550", "w550", "g550", "phs550", "amom550", "bextrat", "bext", "w", "g", "phs", "amom")}
   detail = []
   for r in range(lutstr.efr_n):
      rm = float(lutstr.efr[r])
      a = mie_size_dist(rm, s, wn550, aerm550, qv, opt)
      c = mie_size_dist(rm, s, wn, aerm, qv, opt)
      detail.append((a, c))
      # create_bwgp: Bext, w, g, Phi are FLTARR; then the *_c arrays (FLTARR)
      b550c, w550c, g550c, p550c = store(a["bext"]), store(a["w"]), store(a["g"]), store(a["spm"])
      bc, wc, gc, pc = store(c["bext"]), store(c["w"]), store(c["g"]), store(c["spm"])
      # mixing: lut_MRat is DBLARR, so the products and totals are double
      for name, bb, ww, gg, pp in (("550", b550c, w550c, g550c, p550c), ("", bc, wc, gc, pc)):
         mb = mrat * float(bb)
         mbw = mb * float(ww)
         bext_v = store(mb / mrat)
         w_v = store(mbw / mb)
         g_v = store((mbw * float(gg)) / mbw)
         phs_v = store((mbw * np.asarray(pp, np.float64)) / mbw)
         amom_v = store(legpexp(nmom, qv, weights, phs_v) / two_n_plus_one)
         out["bext" + name].append(bext_v)
         out["w" + name].append(w_v)
         out["g" + name].append(g_v)
         out["phs" + name].append(phs_v)
         out["amom" + name].append(amom_v)
      if opt["store"] == "fltarr":
         out["bextrat"].append(f32(out["bext"][-1]) / f32(out["bext550"][-1]))     # FLOAT / FLOAT
      else:
         out["bextrat"].append(f32(float(c["bext"]) / float(a["bext"])))
   return {k: np.asarray(v, f32) for k, v in out.items()}, detail, dict(aerm550=aerm550, aerm=aerm, qv=qv, weights=weights)


def idl_scatfile(case):
   c = readsav(str(IDL_RUN / f"case_{case}" / "validation" / "luts" / f"case_{case}" / "scatfile.sav"))
   return {"bext550": c["bext550"], "w550": c["w550"], "g550": c["g550"], "phs550": c["phs550"], "amom550": c["amom550"],
           "bextrat": c["bextrat"].reshape(-1), "bext": c["bext"].reshape(-1), "w": c["w"].reshape(-1), "g": c["g"].reshape(-1),
           "phs": c["phs"].reshape(c["phs"].shape[0], -1), "amom": c["amom"].reshape(c["amom"].shape[0], -1)}


def python_scatfile(case):
   z = np.load(HERE / "work" / f"case_{case}" / "work" / "scatfile.npz")
   e = z["phs"].shape[-1]
   return {"bext550": z["bext550"], "w550": z["w550"], "g550": z["g550"], "phs550": z["phs550"].T, "amom550": z["amom550"].T,
           "bextrat": z["bextrat"].reshape(-1), "bext": z["bext"].reshape(-1), "w": z["w"].reshape(-1), "g": z["g"].reshape(-1),
           "phs": z["phs"].reshape(-1, e).T, "amom": z["amom"].reshape(-1, e).T}


def ulp(a, b):
   a = np.asarray(a, f32).ravel().view(np.int32).astype(np.int64)
   b = np.asarray(b, f32).ravel().view(np.int32).astype(np.int64)
   a = np.where(a < 0, -(a & 0x7FFFFFFF), a)
   b = np.where(b < 0, -(b & 0x7FFFFFFF), b)
   return np.abs(a - b)


def compare(mine, ref):
   rows = {}
   for k in ref:
      x, y = np.asarray(mine[k], f32).ravel(), np.asarray(ref[k], f32).ravel()
      d = x.view(np.uint32) != y.view(np.uint32)
      rows[k] = dict(n=int(x.size), differing=int(d.sum()), max_abs=float(np.max(np.abs(x.astype(np.float64) - y))) if d.any() else 0.0,
                     max_ulp=int(ulp(x, y).max()), first=int(np.argmax(d)) if d.any() else None)
   return rows


def main():
   results = {}
   for case in CASES:
      ref = idl_scatfile(case)
      results[case] = {"python_scatfile": compare(python_scatfile(case), ref)}
      results[case]["replica_python"] = compare(scattering(case, PYTHON)[0], ref)
      results[case]["replica_idl"] = compare(scattering(case, IDL)[0], ref)
      print(case, json.dumps({k: {q: v["differing"] for q, v in r.items()} for k, r in results[case].items()}), flush=True)
   (HERE / "results").mkdir(exist_ok=True)
   (HERE / "results" / "scattering_replica.json").write_text(json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
   main()
