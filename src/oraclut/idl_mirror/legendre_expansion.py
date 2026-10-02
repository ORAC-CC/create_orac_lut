"""Adaptive Legendre expansion of Mie size-distribution phase functions.

There is no IDL counterpart.  The IDL, and the Python generator up to source
revision c6ad545 (including the V23 reference LUTs), used one fixed number,
NMom = 1000, as the Gauss-Legendre quadrature order, the number of Legendre
coefficients and the number of moments passed to DISORT
(generate_scattering_properties.pro:103, legpexp.pro).  Since the 2026-10
Legendre change the three quantities are separated for Mie size
distributions:

   Nq   the Gauss-Legendre quadrature order on which the size-distribution-
        averaged phase function is calculated and integrated;
   D    an upper bound on the polynomial degree of that averaged phase
        function, known before any scattering calculation;
   L    the number of Legendre coefficients retained and passed to DISORT,
        found from the averaged phase function itself.

The expansion follows Grainger (1990), DPhil thesis, sections 4.4-4.5
(documents/Grainger 1990.pdf, thesis pp. 53-57):

   p(mu) = sum_l omega_l P_l(mu),  omega_l = (2l+1)/2 int p P_l dmu   [4.4.1-2]

terminated by King's criterion, |omega_l| < 1e-9 for all l >= L, and checked
by reconstructing the phase function to six significant figures.

Why D bounds the averaged phase function.  The preserved Mie kernel
(mie/dlm-code/mieint.f) sums the amplitude series to NStop(x) terms, so the
F11 of one sphere is a polynomial in mu of degree 2*NStop(x).  The size
distribution average is a weighted sum of such polynomials, so its degree is
at most D = 2*NStop(x_max), where x_max is the largest size parameter of the
radius integration actually performed (create_bwgp.mie_integration_limits:
3.5 x effective radius on the legacy lattice for liquid water, 100 um or the
log-normal limits otherwise).  D is a safeguard on the angular resolution only: the
averaged phase function is much smoother than that of its largest particle,
and L is usually well below D.

Why Nq > D.  An Nq-point Gauss-Legendre rule integrates polynomials of degree
below 2*Nq exactly, so omega_l is exact when D + l < 2*Nq.  With Nq > D every
calculated coefficient is free of aliasing, and the coefficients with
D < l < Nq are exactly zero.  Their computed values therefore measure the
finite-precision noise of this particular calculation.

The expansion length L.  L is determined from the coefficients of the
size-distribution-averaged phase function alone (no radiative transfer is
involved): King's criterion is applied literally to the whole tail,
|omega_l| < 1e-9 for every l >= L, where the coefficients beyond D are known
to be zero.  The series is then accepted only if it reproduces the directly
calculated averaged phase function, in double precision, at independent
scattering angles to six significant figures (relative error <= 5e-7).

Finite precision.  For the most forward-peaked phase functions (p(0) ~ 1e6)
the rounding noise of double-precision quadrature, measured in the band just
above D where the exact coefficients vanish, is ~1e-8 to 1e-7.  It then
exceeds 1e-9, and coefficients between the noise level and 1e-9 cannot be
told apart from noise.  They are kept rather than cut at the noise level:
cutting there discards real band-edge structure and makes L depend on Nq,
whereas keeping them can extend L at most to D + 1 because the degree bound
caps the expansion.  The rule is therefore conservative and stable when Nq is
increased.  A noise level above 1e-6, or a failed reconstruction, means that
the adequacy of the quadrature has not been demonstrated; the caller then
increases Nq and finally fails.
"""

import math

import numpy as np


KING_THRESHOLD = 1.0e-9                 # Grainger (1990) section 4.5, after King (1983)
# "accurate to 6 significant figures" (Grainger 1990 section 4.5), read
# strictly: a relative error of 5e-7 leaves six significant figures whatever
# the leading digit.  The Mie phase functions here never fall below ~1e-8 of
# their forward value, so the relative error is not pathological.
RECONSTRUCTION_TOLERANCE = 5.0e-7
QUADRATURE_NOISE_LIMIT = 1.0e-6         # largest admissible |omega_l| where the exact value is zero (l > D)
MIE_MAX_ANGLES = 10000                  # mieint.f: Parameter (Imaxnp = 10000)
MAX_QUADRATURE_ATTEMPTS = 3             # Nq is increased by half when adequacy is not demonstrated

# Scattering angles (degrees) at which the reconstructed phase function is
# compared with the directly calculated one: the forward peak every 0.01
# degree, then every 0.5 degree, including the 30-180 degree range of the
# ORAC viewing geometries.  None coincides with a Gauss-Legendre node.
CHECK_THETA = np.concatenate((np.arange(0.0, 1.0, 0.01), np.arange(1.0, 180.0001, 0.5)))

_GAUSS_LEGENDRE = {}


def gauss_legendre(npts):
   """Gauss-Legendre abscissae (ascending) and weights on [-1, 1].

   numpy.polynomial.legendre.leggauss loses accuracy in the weights at high
   order: against a 50-digit reference its weight nearest x = 1 is wrong by
   2e-7 at 3000 points and 2e-6 at 6000 (scipy.special.roots_legendre
   behaves the same).  Multiplied by a forward peak of ~1e6 that error alone
   creates a spurious Legendre tail of ~1e-4.  The leggauss abscissae are
   accurate, so they are refined by Newton iteration on the three-term
   recurrence and the weights recomputed as 2 / ((1 - x^2) P_N'(x)^2), with
   1 - x^2 formed as (1 - x)(1 + x).  Results are cached by order.
   """

   npts = int(npts)
   if npts in _GAUSS_LEGENDRE:
      return _GAUSS_LEGENDRE[npts]
   if npts < 2:
      raise ValueError("gauss_legendre: at least two points are required")
   x = np.polynomial.legendre.leggauss(npts)[0]
   for iteration in range(3):
      pnm1 = np.ones_like(x)
      pn = x.copy()
      for n in range(2, npts + 1):
         pnm1, pn = pn, ((2.0 * n - 1.0) * x * pn - (n - 1.0) * pnm1) / n
      dpn = npts * (pnm1 - x * pn) / ((1.0 - x) * (1.0 + x))
      x = x - pn / dpn
   pnm1 = np.ones_like(x)
   pn = x.copy()
   for n in range(2, npts + 1):
      pnm1, pn = pn, ((2.0 * n - 1.0) * x * pn - (n - 1.0) * pnm1) / n
   dpn = npts * (pnm1 - x * pn) / ((1.0 - x) * (1.0 + x))
   weights = 2.0 / ((1.0 - x) * (1.0 + x) * dpn**2)
   _GAUSS_LEGENDRE[npts] = (x, weights)
   return x, weights


def mie_nstop(x):
   """Number of Mie amplitude terms summed by mieint.f for size parameter x.

   mieint.f uses Dx + 4.00 Dx^(1/3) + 2 (Dx <= 8 or Dx >= 4200) or
   Dx + 4.05 Dx^(1/3) + 2 otherwise, truncated to an integer.  The larger
   4.05 form plus one is returned so that the result is never below the
   kernel's value whatever the rounding of its single-precision exponent.
   """

   x = float(x)
   if x < 0.02:
      return 2
   return int(x + 4.05 * x**(1.0 / 3.0) + 2.0) + 1


def mie_phase_degree(x_max):
   """Upper bound D on the polynomial degree in mu of a Mie F11 averaged over sizes up to x_max."""

   return 2 * mie_nstop(x_max)


def noise_band(degree):
   """Number of coefficients beyond the degree bound used to measure the rounding noise."""

   return 32 + int(degree) // 32


QUADRATURE_ORDER_STEP = 64              # orders are rounded up to a multiple of this, so few distinct rules are computed


def initial_quadrature_order(degree):
   """Gauss-Legendre order Nq: above the degree bound plus the noise band, rounded up to QUADRATURE_ORDER_STEP.

   Rounding up only adds nodes (the coefficients stay exact) and keeps the
   number of distinct Gauss-Legendre rules, which gauss_legendre caches, small.
   """

   minimum = int(degree) + 1 + noise_band(degree)
   return QUADRATURE_ORDER_STEP * int(math.ceil(minimum / QUADRATURE_ORDER_STEP))


def legendre_coefficients(qv, qw, phase):
   """All Legendre coefficients omega_l, l = 0 ... npts-1, of phase sampled at the nodes qv.

   omega_l = (2l+1)/2 sum(phase * P_l(qv) * qw), the legpexp convention
   (create_bwgp.py) without its single-coefficient stop.  phase may hold one
   function per column, shape (npts, nfun); the result has the same shape.
   """

   qv = np.asarray(qv, dtype=np.float64)
   qw = np.asarray(qw, dtype=np.float64)
   phase = np.asarray(phase, dtype=np.float64)
   single = phase.ndim == 1
   if single:
      phase = phase[:, None]
   if phase.shape[0] != qv.size or qw.size != qv.size:
      raise ValueError("legendre_coefficients: qv, qw and phase must have the same number of points")
   npts = qv.size
   weighted = phase * qw[:, None]
   omega = np.zeros((npts, phase.shape[1]), dtype=np.float64)
   lpnm1 = np.ones_like(qv)
   lpn = qv.copy()
   omega[0, :] = lpnm1 @ weighted / 2.0
   if npts > 1:
      omega[1, :] = 3.0 * (lpn @ weighted) / 2.0
   for n in range(2, npts):
      # calculate the nth Legendre polynomial and integrate up the coefficient
      lpnm1, lpn = lpn, ((2.0 * n - 1.0) / n) * qv * lpn - ((n - 1.0) / n) * lpnm1
      omega[n, :] = (2.0 * n + 1.0) * (lpn @ weighted) / 2.0
   return omega[:, 0] if single else omega


def king_expansion_length(omega, degree, threshold):
   """Number of coefficients L such that |omega_l| < threshold for every l >= L.

   This is King's criterion applied to the whole tail, not to a single
   coefficient: L is one more than the highest order l <= degree with
   |omega_l| >= threshold.  omega must come from a quadrature order
   Nq > degree, so that the coefficients beyond the degree bound are exactly
   zero apart from rounding.
   """

   omega = np.asarray(omega, dtype=np.float64)
   significant = np.flatnonzero(np.abs(omega[:int(degree) + 1]) >= threshold)
   if significant.size == 0:
      raise ValueError("king_expansion_length: no coefficient exceeds the termination threshold")
   return int(significant[-1]) + 1


def rounding_noise(omega, degree):
   """Largest |omega_l| in the band degree < l <= degree + noise_band(degree).

   The exact coefficients vanish there, so the computed values measure the
   rounding noise of the quadrature near the end of the expansion.  The band
   has a fixed width so that the measure does not grow with Nq (the noise
   grows roughly as 2l + 1).
   """

   omega = np.asarray(omega, dtype=np.float64)
   degree = int(degree)
   if omega.ndim != 1 or omega.size <= degree + 1:
      raise ValueError("rounding_noise: the coefficients must extend beyond the degree bound")
   return float(np.max(np.abs(omega[degree + 1:degree + 1 + noise_band(degree)])))


def reconstruct_phase_function(omega, length, mu):
   """Phase function sum_{l < length} omega_l P_l(mu) (Grainger 1990 eq. 4.4.1)."""

   return np.polynomial.legendre.legval(np.asarray(mu, dtype=np.float64), np.asarray(omega[:length], dtype=np.float64))


def expansion_length(omega, degree, check_mu, check_phase):
   """Return (L, noise, reconstruction_error) for one size-distribution-averaged phase function.

   omega are the Legendre coefficients from a quadrature order Nq > degree,
   and check_phase is the directly calculated averaged phase function at the
   independent angles check_mu.

   L follows King's criterion literally, |omega_l| < 1e-9 for all l >= L,
   with the coefficients beyond the degree bound known to be zero.  When the
   rounding noise exceeds 1e-9, coefficients down to the noise level are
   kept, so L <= degree + 1 (conservative; see the module docstring).  The
   reconstruction error is the largest relative difference, in double
   precision, between the series of L terms and check_phase.  The caller
   accepts the expansion only if the error is within RECONSTRUCTION_TOLERANCE
   and the noise within QUADRATURE_NOISE_LIMIT.
   """

   noise = rounding_noise(omega, degree)
   length = king_expansion_length(omega, degree, KING_THRESHOLD)
   reconstructed = reconstruct_phase_function(omega, length, check_mu)
   error = float(np.max(np.abs(reconstructed / check_phase - 1.0)))
   return length, noise, error
