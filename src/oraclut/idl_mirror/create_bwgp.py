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

Arithmetic reproduces the IDL (V21_V22_REFERENCE.md; detailed evidence in Git,
commit f9b0b3b): DOUBLE where the IDL is DOUBLE, FLOAT
where the IDL operands are FLOAT (the literals 0.001, 100.0, 0.4 and 4./3., the
wavenumber 1.0/wl), the IDL evaluation order, glibc exp and pow as called by
IDL, the IDL quadrature.pro abscissae and weights, and FLTARR storage of the
create_bwgp results.  The Mie kernel is selected in oraclut.optics.legacy_mie.

Deliberate differences from the IDL (LUT release V24), used only when the
caller asks for them (generate_scattering_properties, Mie classes), so that the
default behaviour of every routine here is that of the IDL:

1. Integration limit (historical commit c6ad545).  A modified-gamma
   integration can end at the first node of the legacy radius lattice at or
   beyond a multiple of the effective radius (beyond 100 um when necessary)
   instead of at 100 um: lattice_upper_radius.  Keeping the legacy nodes makes
   the change the removed or added tail alone.
2. Radius refinement (historical commit 96bebbe).  Every interval of that
   legacy trapezoid can be divided into 2^k equal parts (nested: the legacy
   nodes and the limits are kept), k being the smallest level giving a given
   size-parameter step and at least MIE_NODES_PER_SPREAD nodes per
   area-weighted standard deviation: radius_refinement_level.  The Mie kernel
   is then called on blocks of radius nodes and the phase-function sum is
   accumulated in the order of the single IDL TOTAL, so the result does not
   depend on the blocking.
3. The phase function can be returned in double precision (phase_dtype) for
   the adaptive Legendre expansion of legendre_expansion.py, whose 1e-9
   termination threshold is below the FLOAT rounding of a forward peak.
"""

import math
from statistics import NormalDist

import numpy as np

from ..idl_math import idl_exp, idl_pow
from ..optics.legacy_mie import mie_single_batch
from .dubovik import dubovik_lognormal_multiple_eps
from .total import idl_total

# Radius limits (microns) and size-parameter step that the IDL create_bwgp.pro
# passes to mie_size_dist_new: [Rm, S, 0.001, 100.0], xres=0.4, all FLOAT.
MIE_RADIUS_LOWER = np.float32(0.001)
MIE_RADIUS_UPPER = np.float32(100.0)
MIE_XRES = np.float32(0.4)
# Refined radius grids (V24): smallest number of nodes per area-weighted
# standard deviation re sqrt(ve) of a modified-gamma distribution, and the
# number of radius nodes x angles passed to one Mie kernel call (32 MB per
# float64 phase array), which bounds the memory of refined grids.
MIE_NODES_PER_SPREAD = 3.0
MIE_BLOCK_ELEMENTS = 4_000_000


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


def size_integration_limits(distname, params, wavenumber):
   """(rl, ru, truncated): the radius limits of the mie_size_dist_new size integration.

   Computed inline in the IDL mie_size_dist_new.pro; a separate routine here so
   that the adaptive Legendre expansion uses the largest size parameter,
   2 pi ru wavenumber, that the integration actually reaches.
   """

   ru_max = 10000.0
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
   truncated = 2.0 * np.pi * ru * wavenumber >= ru_max
   if truncated:
      ru = (ru_max - 1.0) / (2.0 * np.pi * wavenumber)
   return rl, ru, truncated


def legacy_npts(rl, ru, wavenumber, xres):
   """The IDL trapezoid node count: Npts = (long(2D0*!dpi*(ru-rl)*wavenumber/xres)) > 200."""

   return max(int(2.0 * np.pi * (ru - rl) * wavenumber / float(xres)), 200)


def mie_size_dist_new(distname, nd, params, wavenumber, cm, dqv, xres=0.1, npts=None):
   """Scattering parameters of a size distribution of spheres.

   params: 'log_normal'     [mode radius Rm, spread S, -, -]
           'modified_gamma' [effective radius, effective variance, Rl, Ru]
   Returns (bext, bsca, w, g, spm, vavg) where spm[0, :] is F11 at cos(theta) = dqv
   and vavg is the average volume per particle.  Nd is the number density (1).

   npts: number of trapezoid radius nodes on [Rl, Ru].  None is the IDL
   (legacy_npts); mie_integration_limits gives it for the V24 lattice limit
   and refined grids.
   """

   xres = float(xres)
   if np.imag(cm) > 0:
      raise ValueError("mie_size_dist_new: imaginary part of the refractive index must be negative")

   # Create vectors for size integration
   rl, ru, truncated = size_integration_limits(distname, params, wavenumber)
   if truncated:
      if npts is not None:
         raise ValueError("mie_size_dist_new: an explicit node count cannot be combined with a truncated upper radius")
      print("Warning: Radius upper bound truncated to avoid size parameter overflow.")

   # Accurate calculation requires 0.1 step size in x but this can take an age
   # so limit to 200 (IDL: Npts = (long(2D0*!dpi*(ru-rl)*wavenumber/xres)) > 200)
   if npts is None:
      npts = legacy_npts(rl, ru, wavenumber, xres)

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
   # Fortran kernel (IDL: Mie_dlm_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, ..., F11, ...).
   # The kernel is called on blocks of at most MIE_BLOCK_ELEMENTS radius nodes x
   # angles (no IDL counterpart; the IDL holds every node at once).  Each node's
   # results do not depend on the block, and the sum of W1PA * F11 * Dqsc over
   # the nodes is carried from block to block in the sequential order of the
   # single IDL TOTAL, so the result is that of one call.
   dqv = np.asarray(dqv, dtype=np.float64)
   w1pa = w1p * np.pi * r ** 2                          # area weighted distribution
   w1pv = w1pa * 4.0 / 3.0 * r                          # IDL: W1PA * 4./3. * R, i.e. ((W1PA*4.)/3.)*R (DOUBLE)
   dqxt = np.empty(npts, dtype=np.float64)
   dqsc = np.empty(npts, dtype=np.float64)
   dg = np.empty(npts, dtype=np.float64)
   phase_sum = None
   block = max(1, MIE_BLOCK_ELEMENTS // dqv.size)
   for start in range(0, npts, block):
      nodes = slice(start, min(start + block, npts))
      mie = mie_single_batch(dx[nodes], complex(cm), dqv)
      dqxt[nodes] = mie["qext"]
      dqsc[nodes] = mie["qsca"]
      dg[nodes] = mie["g"]
      terms = w1pa[nodes, None] * mie["f11"] * dqsc[nodes, None]          # f11: (nodes, inp)
      if phase_sum is not None:
         terms = np.concatenate((phase_sum[None, :], terms))
      phase_sum = idl_total(terms, dimension=1)

   bext = float(idl_total(w1pa * dqxt))
   bsca = float(idl_total(w1pa * dqsc))
   g = float(idl_total(w1pa * dg * dqsc) / bsca)
   w = bsca / bext

   spm = np.zeros((4, dqv.size), dtype=np.float64)
   # IDL: SPM[0,i] = total(W1PA * DSPM[0,i,*] * Dqsc) / Bsca
   spm[0, :] = phase_sum / bsca
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


def lattice_upper_radius(rm, factor, wavenumber):
   """(upper radius, node count) of a modified-gamma integration ending at factor * rm (V24).

   The legacy grid of the IDL (create_bwgp.pro, mie_size_dist_new.pro) is the
   trapezoid on [0.001, 100] um with legacy_npts nodes.  Moving its endpoint
   would move every node, and on this coarse grid that alone changes the
   optics by up to ~0.5% in extinction where the removed tail is negligible.
   The legacy spacing is therefore kept and the grid stops at the first
   lattice node at or beyond factor * rm (beyond 100 um when necessary).  The
   limits and xres are the FLOAT values of the IDL, as in create_bwgp.
   """

   rl, ru = float(MIE_RADIUS_LOWER), float(MIE_RADIUS_UPPER)
   spacing = (ru - rl) / (legacy_npts(rl, ru, wavenumber, MIE_XRES) - 1)
   intervals = int(math.ceil((factor * rm - rl) / spacing))
   return rl + intervals * spacing, intervals + 1


def radius_refinement_level(spacing, wavenumber, rm, s, refined_xres):
   """Smallest k >= 0 such that spacing / 2^k meets both resolution bounds (V24).

   The bounds are a size-parameter step 2 pi h wavenumber <= refined_xres and
   h <= re sqrt(ve) / MIE_NODES_PER_SPREAD (rm = re, s = ve).  The legacy step
   exceeds the nominal xres = 0.4 by up to 1% because the node count is an
   integer, so both bounds carry a 2% margin: the nominal 0.4 / 2^k counts as
   0.4 / 2^k.
   """

   limit = min(refined_xres / (2.0 * np.pi * wavenumber), rm * math.sqrt(s) / MIE_NODES_PER_SPREAD)
   level = 0
   while spacing / 2**level > 1.02 * limit:
      level += 1
   return level


def mie_integration_limits(distname, rm, s, wavenumber, radius_upper_factor=None, refined_xres=None):
   """(params, npts) that create_bwgp passes to mie_size_dist_new at one wavenumber.

   params = [Rm, S, 0.001, 100.0] as in the IDL (DOUBLE, from the FLOAT S and
   literals) and npts = None (legacy_npts), unless, for a modified-gamma
   distribution (V24):
      radius_upper_factor  ends the legacy lattice at the first node >= factor * rm
                           (lattice_upper_radius), which sets params[3] and npts;
      refined_xres         divides every interval of that legacy grid into 2^k
                           equal parts (radius_refinement_level), so
                           npts = (n - 1) 2^k + 1 for the legacy n; the limits and
                           the legacy nodes are unchanged.
   """

   # IDL: [Rm, S, 0.001, 100.0] is DOUBLE (Rm is DOUBLE): the FLOAT S and
   # literals are promoted from single precision; xres=0.4 is FLOAT.
   params = [float(rm), float(np.float32(s)), float(MIE_RADIUS_LOWER), float(MIE_RADIUS_UPPER)]
   npts = None
   if distname == "modified_gamma" and radius_upper_factor is not None:
      params[3], npts = lattice_upper_radius(rm, radius_upper_factor, wavenumber)
   if distname == "modified_gamma" and refined_xres is not None:
      if npts is None:
         npts = legacy_npts(params[2], params[3], wavenumber, MIE_XRES)
      level = radius_refinement_level((params[3] - params[2]) / (npts - 1), wavenumber, rm, s, refined_xres)
      npts = (npts - 1) * 2**level + 1
   return params, npts


def create_bwgp(distname, rm, s, ri, wl, dqv, scode="mie", tmatrix_path=None, eps=None, neps=None,
                radius_upper_factor=None, refined_xres=None, phase_dtype=np.float32):
   """Bulk extinction, single-scattering albedo, asymmetry and phase function.

   For each wavelength wl[i] (with refractive index ri[i]) the Mie size
   distribution calculation is run and Bext, w, g and Phi[:, i] (F11 at
   cos(theta) = dqv) are returned, plus Vavg (average volume per particle,
   from the last wavelength; called with the single 0.55 micron reference
   wavelength when Vavg is wanted).

   radius_upper_factor, refined_xres: None (the IDL) or the V24 modified-gamma
   integration limit and radius refinement (mie_integration_limits).
   phase_dtype: the IDL Phi is FLTARR (float32); the V24 adaptive Legendre
   expansion asks for float64.
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
   phi = np.zeros((inp, inw), dtype=phase_dtype)
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
            # IDL: mie_size_dist_new, distname, 1.0, [Rm, S, 0.001, 100.0], ..., xres=0.4
            params, npts = mie_integration_limits(distname, rm, s, float(wn[i]), radius_upper_factor, refined_xres)
            bexttmp, bscatmp, wtmp, gtmp, spm, vavgtmp = mie_size_dist_new(
               distname, 1.0, params, float(wn[i]), ri[i], dqv, xres=MIE_XRES, npts=npts)
            phi[:, i] = spm[0, :]
            vavg = vavgtmp
         bext[i] = bexttmp
         w[i] = wtmp
         g[i] = gtmp
   return bext, w, g, phi, vavg
