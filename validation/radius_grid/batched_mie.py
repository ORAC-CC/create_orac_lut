"""Memory-bounded mie_size_dist_new for refined radius grids (validation copy).

Validation-only code.  Identical in mathematics to
oraclut.idl_mirror.create_bwgp.mie_size_dist_new (same limits, trapezoid
nodes and weights, distribution, Mie kernel and sums) except that the Mie
kernel is called on blocks of radius nodes and the sums are accumulated
block by block, so that the (nodes x angles) phase arrays of a refined grid
never exist at once.  Only the summation order differs (rounding level).

   install(module)   replace module.mie_size_dist_new by the batched version
"""

import math

import numpy as np

from oraclut.optics.legacy_mie import mie_single_batch

ELEMENTS_PER_BLOCK = 4_000_000          # radius nodes x angles per Mie call (32 MB per float64 array)


def install(bwgp):
   production = bwgp.mie_size_dist_new

   def mie_size_dist_new(distname, nd, params, wavenumber, cm, dqv, xres=0.1, npts=None):
      if np.imag(cm) > 0:
         raise ValueError("mie_size_dist_new: imaginary part of the refractive index must be negative")
      rl, ru, truncated = bwgp.size_integration_limits(distname, params, wavenumber)
      if truncated:
         return production(distname, nd, params, wavenumber, cm, dqv, xres=xres, npts=npts)
      if npts is None:
         npts = max(int(2.0 * np.pi * (ru - rl) * wavenumber / xres), 200)
      absc, wght = bwgp.quadrature("T", npts)
      r, w1 = bwgp.shift_quadrature(absc, wght, rl, ru)
      if distname == "modified_gamma":
         alpha = (1.0 - 3.0 * params[1]) / params[1]
         b = 1.0 / (params[0] * params[1])
         n = b ** (-alpha - 1.0) * math.gamma(alpha + 1.0)
         w1p = w1 / n * r ** alpha * np.exp(-b * r)
      else:
         log_s = math.log(params[1])
         w1p = (w1 * np.exp(-0.5 * (np.log(r / params[0]) / log_s) ** 2)
                / (math.sqrt(2.0) * math.sqrt(math.pi) * r * log_s))
      if nd != 1.0:
         w1p = w1p * nd
      dx = 2.0 * np.pi * r * wavenumber
      dqv = np.asarray(dqv, dtype=np.float64)
      w1pa = w1p * np.pi * r ** 2
      w1pv = w1pa * (4.0 / 3.0) * r
      block = max(1, ELEMENTS_PER_BLOCK // dqv.size)
      bext = bsca = gsca = 0.0
      phase = np.zeros(dqv.size, dtype=np.float64)
      for start in range(0, r.size, block):
         sl = slice(start, start + block)
         mie = mie_single_batch(dx[sl], complex(cm), dqv)
         bext += float(np.sum(w1pa[sl] * mie["qext"]))
         bsca += float(np.sum(w1pa[sl] * mie["qsca"]))
         gsca += float(np.sum(w1pa[sl] * mie["g"] * mie["qsca"]))
         phase += np.sum(w1pa[sl, None] * mie["f11"] * mie["qsca"][:, None], axis=0)
      spm = np.zeros((4, dqv.size), dtype=np.float64)
      spm[0, :] = phase / bsca
      return bext, bsca, bsca / bext, gsca / bsca, spm, float(np.sum(w1pv))

   bwgp.mie_size_dist_new = mie_size_dist_new
   return production
