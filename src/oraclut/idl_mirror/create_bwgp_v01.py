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

Arithmetic follows the validated Python (double precision throughout, single
precision only where the IDL result is stored in FLTARR by the caller).  Where
the IDL evaluates an expression in a different order the IDL form is quoted in
a comment; the validated order is kept so that products stay bitwise identical.
"""

import math
from statistics import NormalDist

import numpy as np

from ..optics.legacy_mie import mie_single_batch
from .dubovik import dubovik_lognormal_multiple_eps


def quadrature(quadtype, npts):
   """Abscissae and weights on [-1, 1]: 'T' trapezium or 'G' Gauss-Legendre.

   IDL 'T': Abscissa = -1D0+2*Dindgen(NPts)/(N-1); Weight = 2D0/(N-1), ends 1D0/(N-1).
   IDL 'G': Gauss-Legendre zeros by Newton iteration, returned in ascending
   order; numpy's leggauss gives the same nodes in the same order.
   """

   npts = int(npts)
   if quadtype.upper() == "T":
      abscissa = np.linspace(-1.0, 1.0, npts, dtype=np.float64)      # IDL: -1D0+2*Dindgen(NPts)/(N-1)
      weight = np.full(npts, 2.0 / (npts - 1), dtype=np.float64)
      weight[0] = 1.0 / (npts - 1)
      weight[npts - 1] = 1.0 / (npts - 1)
   elif quadtype.upper() == "G":
      if npts > 2000:
         raise ValueError("Error in quadrature: Too many quadrature points")
      abscissa, weight = np.polynomial.legendre.leggauss(npts)
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
      alpha = (1.0 - 3.0 * params[1]) / params[1]
      b = 1.0 / (params[0] * params[1])
      n = b ** (-alpha - 1.0) * math.gamma(alpha + 1.0)          # IDL: (b^(-alpha-1))*Factorial(alpha)
      w1p = w1 / n * r ** alpha * np.exp(-b * r)                 # IDL: W1 * (Nd/n) * R^alpha * exp(-b*R), Nd = 1
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
   w1pv = w1pa * (4.0 / 3.0) * r                        # IDL: W1PA * 4./3. * R  (volume weighted)

   bext = float(np.sum(w1pa * dqxt))
   bsca = float(np.sum(w1pa * dqsc))
   g = float(np.sum(w1pa * dg * dqsc) / bsca)
   w = bsca / bext

   spm = np.zeros((4, dqv.size), dtype=np.float64)
   # IDL: SPM[0,i] = total(W1PA * DSPM[0,i,*] * Dqsc) / Bsca
   spm[0, :] = np.sum(w1pa[:, None] * f11 * dqsc[:, None], axis=0) / bsca
   # F33, F12, F34 (spm[1:4]) are not needed for the scalar RT and are left zero.

   vavg = float(np.sum(w1pv))
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
   lc[0] = np.sum(phase * qw) / 2.0
   if inp == 1:
      return 0, lc
   lc[1] = 3.0 * np.sum(phase * qv * qw) / 2.0
   lpnm2 = np.ones_like(qv)
   lpnm1 = qv.copy()
   n = 2
   while n < inp:
      # calculate the nth Legendre polynomial
      lpn = ((2.0 * n - 1.0) / n) * qv * lpnm1 - ((n - 1.0) / n) * lpnm2
      # integrate up the Legendre coefficient
      lc[n] = (2.0 * n + 1.0) * np.sum(phase * lpn * qw) / 2.0
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

   wl = np.atleast_1d(np.asarray(wl, dtype=np.float64))
   ri = np.atleast_1d(np.asarray(ri))
   if ri.size != wl.size:
      raise ValueError("create_bwgp: Array size mismatch!")

   inp = np.asarray(dqv).size
   inw = wl.size
   # load_srfstrarr pads shorter QM=2 SRFs with zero wavelength entries so
   # every channel shares one rectangular array.  IDL's loop skips those
   # entries; use the same zero result without emitting a divide-by-zero
   # warning before the per-point zero guard below.
   wn = np.divide(1.0, wl, out=np.zeros_like(wl), where=wl != 0.0)  # IDL: wn = 1.0/wl
   bext = np.zeros(inw, dtype=np.float64)
   w = np.zeros(inw, dtype=np.float64)
   g = np.zeros(inw, dtype=np.float64)
   phi = np.zeros((inp, inw), dtype=np.float64)
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
               tmatrix_path, 1.0, rm, s, wn[i], ri_i, eps, neps, dqv=dqv,
               renorm_ph=True, silent=True)
            phi[:, i] = phase
         else:
            # The legacy wrapper uses Mie for zero/long-wave points and for
            # thermal channels (wl >= 6 micron), even for a T-matrix component.
            bexttmp, bscatmp, wtmp, gtmp, spm, vavgtmp = mie_size_dist_new(
               distname, 1.0, [rm, s, 0.001, 100.0], wn[i], ri[i], dqv, xres=0.4)
            phi[:, i] = spm[0, :]
            vavg = vavgtmp
         bext[i] = bexttmp
         w[i] = wtmp
         g[i] = gtmp
   return bext, w, g, phi, vavg
