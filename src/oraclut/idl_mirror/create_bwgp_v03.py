"""Python counterparts of create_bwgp.pro and the Mie routines it calls.

IDL:  pro create_bwgp, distname, Rm, S, RI, wl, Dqv, Bext, w, g, Phi, scode=, tmatrix_path=, eps=, neps=, Vavg=
      pro mie_size_dist_new, distname, Nd, params, wavenumber, Cm, Dqv=, xres=, /DLM, Bext, Bsca, w, g, SPM, Vavg=
      Pro Quadrature, Quadtype, NPts, Abscissa, Weight
      pro shift_quadrature, abscissa, weights, A, B, new_abscissa, new_weights
      pro legpexp, Inp, qv, qw, phase, Inlc, lc

This is the scattering path made explicit:

   Mie: size distribution -> single-particle Mie -> bulk properties;
   T-matrix: dubovik_lognormal_multiple_eps -> bulk properties;
   both paths return Bext, Bsca, w, g and phase function -> Legendre expansion.

Arithmetic reproduces the IDL (validation/v21_v22_terra_modis/
REPORT_single_precision_equivalence.md): DOUBLE where the IDL is DOUBLE, FLOAT
where the IDL operands are FLOAT (the literals 0.001, 100.0, 0.4 and 4./3., the
wavenumber 1.0/wl), the IDL evaluation order, glibc exp and pow as called by
IDL, the IDL quadrature.pro abscissae and weights, and FLTARR storage of the
create_bwgp results.  The Mie kernel is selected in oraclut.optics.legacy_mie.
"""

import math
from statistics import NormalDist

import numpy as np

from ..idl_math import idl_exp, idl_pow
from ..optics.legacy_mie import mie_single_batch
from .dubovik import dubovik_lognormal_multiple_eps
from .total import idl_total


def _legendre(nn, x):
   """IDL quadrature.pro: Pro Legendre, nn, x, Pl, Pm, Pn (nn >= 2)."""

   pl = 1.0                                            # value doesn't matter
   pm = 1.0
   pn = x
   for n in range(2, int(nn) + 1):
      inv = 1.0 / float(n)                             # IDL: in = 1d0/Double(n)
      pl = pm
      pm = pn
      pn = (2.0 - inv) * x * pm - (1.0 - inv) * pl
   return pl, pm, pn


def quadrature(quadtype, npts):
   """Abscissae and weights on [-1, 1]: 'T' trapezium or 'G' Gauss-Legendre (IDL quadrature.pro).

   IDL 'T': Abscissa = -1D0+2*Dindgen(NPts)/(N-1); Weight = 2D0/(N-1), ends 1D0/(N-1).
   IDL 'G': zeros of P_N by Newton iteration from FirstGuess (to 1d-14) and
   weights 2(1-x^2)/(N P_(N-1))^2, reversed into ascending order.  The IDL
   arithmetic is reproduced operation by operation; numpy's leggauss gives
   nodes and weights that differ in the last bits.
   """

   npts = int(npts)
   n = float(npts)                                     # IDL: N = Double(NPts)
   if quadtype.upper() == "T":
      abscissa = -1.0 + 2.0 * np.arange(npts, dtype=np.float64) / (n - 1.0)
      weight = np.full(npts, 2.0 / (n - 1.0), dtype=np.float64)
      weight[0] = 1.0 / (n - 1.0)
      weight[npts - 1] = 1.0 / (n - 1.0)
   elif quadtype.upper() == "G":
      if npts > 2000:
         raise ValueError("Error in quadrature: Too many quadrature points")
      abscissa = np.zeros(npts, dtype=np.float64)
      weight = np.zeros(npts, dtype=np.float64)
      tolerance = 10.0 ** (-14)                        # IDL: 1d1^(-SigFig), SigFig = 14
      for term, zero in enumerate(range(1, npts + 1)):
         lx = 19262.0
         x = math.cos(math.pi * (zero - 0.25) / (npts + 0.5))          # IDL: FirstGuess('G', NPts, Zero)
         while abs(x - lx) > tolerance:
            lx = x
            pl, pm, pn = _legendre(npts, x)
            x = x - (1.0 - x * x) * pn / (float(npts) * (pm - x * pn))   # IDL: x - NewtonG('G', NPts, x)
         abscissa[term] = x
         pl, pm, pn = _legendre(n, x)
         weight[term] = 2.0 * (1.0 - x * x) / (n * pm) ** 2
      abscissa = abscissa[::-1].copy()                 # IDL: Reverse(Abscissa)
      weight = weight[::-1].copy()
   else:
      raise NotImplementedError("quadrature: only trapezium ('T') and Gaussian ('G') are ported")
   return abscissa, weight


def shift_quadrature(abscissa, weights, a, b):
   """Shift abscissae and weights from [-1, 1] to [a, b]."""

   new_abscissa = ((a + b) + (b - a) * abscissa) / 2.0
   new_weights = (b - a) * weights / 2.0
   return new_abscissa, new_weights


def gauss_cvf(p):
   """IDL GAUSS_CVF(P): cutoff V such that the probability that X > V is P."""

   return -NormalDist().inv_cdf(p)


def mie_size_dist_new(distname, nd, params, wavenumber, cm, dqv, xres=0.1):
   """Scattering parameters of a size distribution of spheres.

   params: 'log_normal'     [mode radius Rm, spread S, -, -]
           'modified_gamma' [effective radius, effective variance, Rl, Ru]
   Returns (bext, bsca, w, g, spm, vavg) where spm[0, :] is F11 at cos(theta) = dqv
   and vavg is the average volume per particle.  Nd is the number density (1).
   """

   ru_max = 10000.0
   xres = float(xres)
   if np.imag(cm) > 0:
      raise ValueError("mie_size_dist_new: imaginary part of the refractive index must be negative")

   # Create vectors for size integration
   tq = gauss_cvf(0.999)
   if distname == "modified_gamma":
      rl = float(params[2])
      ru = float(params[3])
   elif distname == "log_normal":
      rl = math.exp(math.log(params[0]) + tq * math.log(params[1]))
      ru = math.exp(math.log(params[0]) - tq * math.log(params[1]) + math.log(4.0))
   else:
      raise ValueError("Invalid size distribution name: " + str(distname))

   if 2.0 * np.pi * rl * wavenumber >= ru_max:
      raise ValueError("Lower bound of integral is larger than maximum permitted size parameter.")
   if 2.0 * np.pi * ru * wavenumber >= ru_max:
      ru = (ru_max - 1.0) / (2.0 * np.pi * wavenumber)
      print("Warning: Radius upper bound truncated to avoid size parameter overflow.")

   # Accurate calculation requires 0.1 step size in x but this can take an age
   # so limit to 200 (IDL: Npts = (long(2D0*!dpi*(ru-rl)*wavenumber/xres)) > 200)
   npts = max(int(2.0 * np.pi * (ru - rl) * wavenumber / xres), 200)

   # quadrature on the radii
   absc, wght = quadrature("T", npts)
   r, w1 = shift_quadrature(absc, wght, rl, ru)

   if distname == "modified_gamma":
      # params[0] = effective radius, params[1] = effective variance
      alpha = (1.0 - 3 * params[1]) / params[1]
      b = 1.0 / (params[0] * params[1])
      # IDL: N = (b^(-alpha-1))*Factorial(alpha); FACTORIAL is GAMMA(alpha+1d)
      # for non-integer alpha (math.gamma reproduces the IDL value for the
      # validated cases) and the exact table value for integer alpha <= 20.
      n = math.pow(b, -alpha - 1.0) * math.gamma(alpha + 1.0)
      # IDL: W1P = W1 * (Nd/n) * R^alpha * exp(-b*R)   (DOUBLE; glibc pow and exp)
      w1p = w1 * (float(nd) / n) * idl_pow(r, alpha) * idl_exp(-b * r)
   else:
      # params[0] = mode radius, params[1] = S
      # IDL: W1P = W1 * Nd / (sqrt(2D0)*sqrt(!dpi)*R*alog(S)) * exp(-0.5D0*(alog(R/Rm)/alog(S))^2)
      log_s = math.log(params[1])
      w1p = (w1 * np.exp(-0.5 * (np.log(r / params[0]) / log_s) ** 2)
             / (math.sqrt(2.0) * math.sqrt(math.pi) * r * log_s))
      if nd != 1.0:
         w1p = w1p * nd

   dx = 2.0 * np.pi * r * wavenumber

   # Single-particle Mie efficiencies and F11 at every radius: the preserved
   # Fortran kernel (IDL: Mie_dlm_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, ..., F11, ...)
   dqv = np.asarray(dqv, dtype=np.float64)
   mie = mie_single_batch(dx, complex(cm), dqv)
   dqxt = mie["qext"]
   dqsc = mie["qsca"]
   dg = mie["g"]
   f11 = mie["f11"]                                    # (npts, inp)

   w1pa = w1p * np.pi * r ** 2                          # area weighted distribution
   w1pv = w1pa * 4.0 / 3.0 * r                          # IDL: W1PA * 4./3. * R, i.e. ((W1PA*4.)/3.)*R (DOUBLE)

   bext = float(idl_total(w1pa * dqxt))
   bsca = float(idl_total(w1pa * dqsc))
   g = float(idl_total(w1pa * dg * dqsc) / bsca)
   w = bsca / bext

   spm = np.zeros((4, dqv.size), dtype=np.float64)
   # IDL: SPM[0,i] = total(W1PA * DSPM[0,i,*] * Dqsc) / Bsca
   spm[0, :] = idl_total(w1pa[:, None] * f11 * dqsc[:, None], dimension=1) / bsca
   # F33, F12, F34 (spm[1:4]) are not needed for the scalar RT and are left zero.

   vavg = float(idl_total(w1pv))
   return bext, bsca, w, g, spm, vavg


def legpexp(inp, qv, qw, phase):
   """Expand ``phase`` (sampled at Legendre points qv, weights qw) in Legendre polynomials.

   Returns (inlc, lc): lc[n] = (2n+1)/2 * sum(phase * P_n * qw), the expansion
   stopping when |lc[n]| < 1e-9 (later coefficients remain zero), as in the IDL.
   """

   qv = np.asarray(qv, dtype=np.float64)
   qw = np.asarray(qw, dtype=np.float64)
   phase = np.asarray(phase, dtype=np.float64)
   if not (qv.shape == qw.shape == phase.shape):
      raise ValueError("legpexp: qv, qw and phase must have equal shapes")
   lc = np.zeros(inp, dtype=np.float64)
   lc[0] = idl_total(phase * qw) / 2.0
   if inp == 1:
      return 0, lc
   lc[1] = 3.0 * idl_total(phase * qv * qw) / 2.0
   lpnm2 = np.ones_like(qv)
   lpnm1 = qv.copy()
   n = 2
   while n < inp:
      # calculate the nth Legendre polynomial
      lpn = ((2.0 * n - 1.0) / n) * qv * lpnm1 - ((n - 1.0) / n) * lpnm2
      # integrate up the Legendre coefficient
      lc[n] = (2.0 * n + 1.0) * idl_total(phase * lpn * qw) / 2.0
      if abs(lc[n]) < 1.0e-9:
         break
      lpnm2 = lpnm1
      lpnm1 = lpn
      n = n + 1
   inlc = n - 1
   return inlc, lc


def create_bwgp(distname, rm, s, ri, wl, dqv, scode="mie", tmatrix_path=None, eps=None, neps=None):
   """Bulk extinction, single-scattering albedo, asymmetry and phase function.

   For each wavelength wl[i] (with refractive index ri[i]) the Mie size
   distribution calculation is run and Bext, w, g and Phi[:, i] (F11 at
   cos(theta) = dqv) are returned, plus Vavg (average volume per particle,
   from the last wavelength; called with the single 0.55 micron reference
   wavelength when Vavg is wanted).
   """

   wl = np.atleast_1d(np.asarray(wl, dtype=np.float32))  # IDL: FLOAT wavelengths (0.55, srfstrarr.wvl)
   ri = np.atleast_1d(np.asarray(ri))
   if ri.size != wl.size:
      raise ValueError("create_bwgp: Array size mismatch!")

   inp = np.asarray(dqv).size
   inw = wl.size
   # load_srfstrarr pads shorter QM=2 SRFs with zero wavelength entries so
   # every channel shares one rectangular array.  IDL's loop skips those
   # entries; use the same zero result without emitting a divide-by-zero
   # warning before the per-point zero guard below.
   wn = np.divide(np.float32(1.0), wl, out=np.zeros_like(wl), where=wl != 0.0)  # IDL: wn = 1.0/wl (FLOAT)
   # IDL: Bext = fltarr(Inw); w = Bext; g = Bext; Phi = fltarr(Inp,Inw)
   bext = np.zeros(inw, dtype=np.float32)
   w = np.zeros(inw, dtype=np.float32)
   g = np.zeros(inw, dtype=np.float32)
   phi = np.zeros((inp, inw), dtype=np.float32)
   vavg = 0.0

   dotmatrix = scode is not None and scode.lower() == "tmatrix"
   if dotmatrix and np.any(wl < 6.0) and (eps is None or neps is None):
      raise ValueError('Asymmetry parameters "eps" and "neps" are required for T-matrix calculations')
   if dotmatrix and np.any(wl < 6.0) and tmatrix_path is None:
      raise FileNotFoundError("T-matrix scattering requires tmatrix_path")

   for i in range(inw):                                # LOOP OVER WAVELENGTH OR CHANNEL
      if wl[i] != 0:
         if dotmatrix and wl[i] < 6.0:
            ri_i = complex(ri[i])
            if ri_i.imag > -5.0e-4:
               ri_i = complex(ri_i.real, -5.0e-4)
            elif ri_i.imag < -0.5:
               raise ValueError(f"Imaginary RI is too negative for the Dubovik LUT range: k = {ri_i.imag}")
            bexttmp, bscatmp, wtmp, gtmp, phase = dubovik_lognormal_multiple_eps(
               tmatrix_path, 1.0, rm, s, float(wn[i]), ri_i, eps, neps, dqv=dqv,
               renorm_ph=True, silent=True)
            phi[:, i] = phase
         else:
            # The legacy wrapper uses Mie for zero/long-wave points and for
            # thermal channels (wl >= 6 micron), even for a T-matrix component.
            # IDL: [Rm, S, 0.001, 100.0] is DOUBLE (Rm is DOUBLE): the FLOAT
            # S and literals are promoted from single precision; xres=0.4 is FLOAT.
            params = [float(rm), float(np.float32(s)), float(np.float32(0.001)), float(np.float32(100.0))]
            bexttmp, bscatmp, wtmp, gtmp, spm, vavgtmp = mie_size_dist_new(
               distname, 1.0, params, float(wn[i]), ri[i], dqv, xres=np.float32(0.4))
            phi[:, i] = spm[0, :]
            vavg = vavgtmp
         bext[i] = bexttmp
         w[i] = wtmp
         g[i] = gtmp
   return bext, w, g, phi, vavg
