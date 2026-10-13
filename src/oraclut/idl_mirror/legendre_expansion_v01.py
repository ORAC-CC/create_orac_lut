"""Adaptive Legendre expansion of size-distribution-averaged Mie phase functions (V24).

There is no IDL counterpart.  The IDL generate_scattering_properties.pro, and
the Python generator up to LUT release V23, used one fixed number NMom = 1000
for three different things: the Gauss-Legendre quadrature order on which the
phase function is sampled, the number of Legendre coefficients computed by
legpexp, and the number of moments passed to DISORT.  At visible wavelengths
and moderate to large effective radii a 1000-point rule aliases the
coefficients of the strongly forward-peaked Mie phase function and the
expansion is cut short.  V24 separates the three quantities for Mie size
distributions (historical method: commit dbcc42c and
validation/REPORT_lut_numerics_development.md, sections 7-9):

   D    an upper bound on the polynomial degree in mu = cos(theta) of the
        averaged phase function, known before the scattering calculation;
   Nq   a Gauss-Legendre order above D, so that no coefficient is aliased;
   L    the number of Legendre coefficients kept, from King's criterion
        applied to the coefficients of the averaged phase function itself.

The expansion is that of Grainger (1990), DPhil thesis, sections 4.4-4.5
(documents/Grainger 1990.pdf):

   p(mu) = sum_l omega_l P_l(mu),   omega_l = (2l+1)/2 int p(mu) P_l(mu) dmu

terminated by King's criterion |omega_l| < 1e-9, here applied to the whole
tail (every l >= L), and accepted only if the series reproduces the directly
calculated phase function to six significant figures.  omega_l includes the
factor 2l+1, as legpexp's coefficients do; the DISORT moments are
chi_l = omega_l / (2l+1).

Degree bound.  The preserved Mie kernel (mie/dlm-code/mieint.f) sums the
amplitude series for n = 1 ... NStop(x), with
NStop = x + 4.05 x^(1/3) + 2 for 8 < x < 4200 (4.00 otherwise, 2 for
x < 0.02), so S1 and S2 are polynomials in mu of degree NStop and the F11 of
one sphere has degree 2 NStop.  A size-distribution average is a weighted sum
of such polynomials, so its degree is at most D = 2 NStop(x_max), x_max being
the largest size parameter of the radius integration actually performed.

Quadrature.  An Nq-point Gauss-Legendre rule integrates polynomials of degree
below 2 Nq exactly, so omega_l is exact for D + l < 2 Nq.  With Nq > D every
computed coefficient is free of aliasing and those with D < l < Nq are exactly
zero; their computed values measure the rounding noise of the calculation.
The IDL quadrature.pro refuses more than 2000 points and numpy's
leggauss weights lose accuracy near mu = 1 at these orders (about 2e-6 at
6000 points, which a forward peak of 1e6 turns into a spurious tail of 1e-4),
so the nodes are Newton-refined here and the weights recomputed.

Expansion length.  L is one more than the highest order l <= D with
|omega_l| >= 1e-9.  For the most strongly forward-peaked phase functions the
double-precision rounding noise just above D can exceed 1e-9; coefficients
down to that noise level are then kept (L <= D + 1) rather than cut at the
noise, which keeps L stable when Nq is increased.  The calculation is
rejected if the noise exceeds 1e-6 or the reconstruction error exceeds 5e-7;
the caller then raises Nq and finally fails.
"""

import math

import numpy as np


KING_THRESHOLD = 1.0e-9                 # Grainger (1990) section 4.5, after King (1983)
RECONSTRUCTION_TOLERANCE = 5.0e-7       # "six significant figures" (Grainger 1990 section 4.5), relative
QUADRATURE_NOISE_LIMIT = 1.0e-6         # largest admissible |omega_l| where the exact value is zero (l > D)
MIE_MAX_ANGLES = 10000                  # mieint.f: Parameter (Imaxnp = 10000)
MAX_QUADRATURE_ATTEMPTS = 3             # Nq is raised by half when adequacy is not demonstrated
QUADRATURE_ORDER_STEP = 64              # Nq is rounded up to a multiple of this (few distinct rules)

# Scattering angles (degrees) at which the series is compared with the
# directly calculated phase function: every 0.01 degree across the forward
# peak, then every 0.5 degree to 180 degrees.
CHECK_THETA = np.concatenate((np.arange(0.0, 1.0, 0.01), np.arange(1.0, 180.0001, 0.5)))

_RULES = {}


def gauss_legendre(npts):
   """Gauss-Legendre nodes (ascending) and weights on [-1, 1], Newton-refined; cached by order."""

   npts = int(npts)
   if npts < 2:
      raise ValueError("gauss_legendre: at least two points are required")
   if npts in _RULES:
      return _RULES[npts]
   x = np.polynomial.legendre.leggauss(npts)[0]
   for iteration in range(4):
      # P_N and P_(N-1) by the three-term recurrence; the last pass only
      # evaluates them for the weights
      p0 = np.ones_like(x)
      p1 = x.copy()
      for n in range(2, npts + 1):
         p0, p1 = p1, ((2.0 * n - 1.0) * x * p1 - (n - 1.0) * p0) / n
      derivative = npts * (p0 - x * p1) / ((1.0 - x) * (1.0 + x))
      if iteration == 3:
         break
      x = x - p1 / derivative
   weights = 2.0 / ((1.0 - x) * (1.0 + x) * derivative**2)
   _RULES[npts] = (x, weights)
   return x, weights


def mie_nstop(x):
   """An upper bound on the number of amplitude terms mieint.f sums for size parameter x.

   The 4.05 form is used throughout (it is never below the 4.00 form) and one
   is added so that the bound holds whatever the rounding of the kernel's
   single-precision exponent 1./3.
   """

   x = float(x)
   if x < 0.02:
      return 2
   return int(x + 4.05 * x**(1.0 / 3.0) + 2.0) + 1


def degree_bound(x_max):
   """D: an upper bound on the degree in mu of a Mie F11 averaged over size parameters up to x_max."""

   return 2 * mie_nstop(x_max)


def noise_band(degree):
   """Number of coefficients beyond the degree bound in which the rounding noise is measured."""

   return 32 + int(degree) // 32


def quadrature_order(degree):
   """Initial Nq: above D plus the noise band, rounded up to QUADRATURE_ORDER_STEP."""

   minimum = int(degree) + 1 + noise_band(degree)
   return QUADRATURE_ORDER_STEP * int(math.ceil(minimum / QUADRATURE_ORDER_STEP))


def legendre_coefficients(qv, qw, phase):
   """omega_l = (2l+1)/2 sum(phase * P_l(qv) * qw) for l = 0 ... len(qv)-1 (double precision)."""

   qv = np.asarray(qv, dtype=np.float64)
   weighted = np.asarray(phase, dtype=np.float64) * np.asarray(qw, dtype=np.float64)
   if weighted.shape != qv.shape:
      raise ValueError("legendre_coefficients: qv, qw and phase must have equal shapes")
   npts = qv.size
   omega = np.zeros(npts, dtype=np.float64)
   pnm1 = np.ones_like(qv)
   pn = qv.copy()
   omega[0] = np.dot(pnm1, weighted) / 2.0
   omega[1] = 3.0 * np.dot(pn, weighted) / 2.0
   for n in range(2, npts):
      # the nth Legendre polynomial at the nodes, and its coefficient
      pnm1, pn = pn, ((2.0 * n - 1.0) / n) * qv * pn - ((n - 1.0) / n) * pnm1
      omega[n] = (2.0 * n + 1.0) * np.dot(pn, weighted) / 2.0
   return omega


def rounding_noise(omega, degree):
   """Largest |omega_l| for degree < l <= degree + noise_band(degree), where the exact values vanish."""

   degree = int(degree)
   band = omega[degree + 1:degree + 1 + noise_band(degree)]
   if band.size < noise_band(degree):
      raise ValueError("rounding_noise: the coefficients do not extend through the noise band")
   return float(np.max(np.abs(band)))


def king_length(omega, degree):
   """L such that |omega_l| < KING_THRESHOLD for every l >= L (coefficients beyond D are zero)."""

   significant = np.flatnonzero(np.abs(omega[:int(degree) + 1]) >= KING_THRESHOLD)
   if significant.size == 0:
      raise ValueError("king_length: no Legendre coefficient reaches the termination threshold")
   return int(significant[-1]) + 1


def series(omega, length, mu):
   """sum_{l < length} omega_l P_l(mu) (Grainger 1990 eq. 4.4.1), double precision."""

   return np.polynomial.legendre.legval(np.asarray(mu, dtype=np.float64), omega[:length])


def expand(omega, degree, check_mu, check_phase):
   """(L, noise, reconstruction error) of one averaged phase function.

   omega from a Gauss-Legendre order Nq > degree + noise_band(degree);
   check_phase is the directly calculated phase function at check_mu.  The
   expansion is acceptable when noise <= QUADRATURE_NOISE_LIMIT and the error
   <= RECONSTRUCTION_TOLERANCE.
   """

   noise = rounding_noise(omega, degree)
   length = king_length(omega, degree)
   error = float(np.max(np.abs(series(omega, length, check_mu) / np.asarray(check_phase, dtype=np.float64) - 1.0)))
   return length, noise, error
