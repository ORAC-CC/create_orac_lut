"""Tests of the adaptive Legendre expansion of Mie phase functions.

See src/oraclut/idl_mirror/legendre_expansion.py and
validation/REPORT_lut_numerics_development.md.  The slow test is the known
difficult case: ice spheres, effective radius 93 microns, MODIS channel 8
(0.412 microns), where the fixed NMom = 1000 expansion (the IDL, and the
Python generator up to c6ad545) aliases even the normalisation and asymmetry
moments.  The liquid-water cases use the production size integration
(3.5 x effective radius on the legacy lattice).
"""

import importlib
from pathlib import Path

import numpy as np
import pytest

from oraclut.idl_mirror import legendre_expansion as le
from oraclut.idl_mirror.create_bwgp import create_bwgp, legpexp, quadrature
from oraclut.idl_mirror.load_mmdat import load_mmdat


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
MODIS_CHANNEL_8 = float(np.float32(0.4118525))     # terra_modis_v1.inst centre wavelength, microns
# the module (oraclut.idl_mirror exports a function of the same name)
gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")


def _legendre_squares(x, w, orders):
   """sum(w * P_k(x)^2) for each k in orders (exactly 2 / (2k + 1) for a Gauss rule of order > k)."""

   values = {}
   pnm1 = np.ones_like(x)
   pn = x.copy()
   for k in range(max(orders) + 1):
      if k == 0:
         current = pnm1
      elif k == 1:
         current = pn
      else:
         pnm1, pn = pn, ((2.0 * k - 1.0) * x * pn - (k - 1.0) * pnm1) / k
         current = pn
      if k in orders:
         values[k] = float(np.sum(w * current**2))
   return values


def test_gauss_legendre_is_accurate_at_high_order():
   npts = 3300
   x, w = le.gauss_legendre(npts)
   assert np.all(np.diff(x) > 0.0) and np.allclose(x, -x[::-1], rtol=0.0, atol=1e-15)
   assert abs(np.sum(w) - 2.0) < 1e-13
   for k, value in _legendre_squares(x, w, [0, 1, 1000, npts - 1]).items():
      assert abs(value * (2 * k + 1) / 2.0 - 1.0) < 1e-10, k
   # at low order it reproduces numpy's rule (whose weights are good to ~1e-12 at 64 points)
   xn, wn = np.polynomial.legendre.leggauss(64)
   x64, w64 = le.gauss_legendre(64)
   assert np.allclose(x64, xn, rtol=0.0, atol=1e-15) and np.allclose(w64, wn, rtol=1e-11, atol=0.0)


def test_accurate_weights_remove_the_spurious_tail_of_a_forward_peaked_polynomial():
   # A strongly forward-peaked band-limited phase function (truncated
   # Henyey-Greenstein, g = 0.995, degree 3000): its coefficients beyond the
   # degree are exactly zero.
   degree = 3000
   omega_true = (2.0 * np.arange(degree + 1) + 1.0) * 0.995**np.arange(degree + 1)
   npts = le.initial_quadrature_order(degree)
   errors = {}
   tails = {}
   for name, (x, w) in {"numpy": np.polynomial.legendre.leggauss(npts), "refined": le.gauss_legendre(npts)}.items():
      qv = -x
      phase = np.polynomial.legendre.legval(qv, omega_true)
      omega = le.legendre_coefficients(qv, w, phase)
      errors[name] = float(np.max(np.abs(omega[:degree + 1] - omega_true)))
      tails[name] = le.rounding_noise(omega, degree)
   # double-precision noise for a peak p(1) ~ 8e4; numpy's weights are ~200-500 times worse
   assert errors["refined"] < 3e-8 and tails["refined"] < 3e-8
   assert errors["numpy"] > 100.0 * errors["refined"] and tails["numpy"] > 100.0 * tails["refined"]


def test_king_criterion_is_applied_to_the_whole_tail():
   # Coefficients decaying as exp(-l/10) with one exact zero at l = 100.  A
   # single-coefficient test (legpexp) stops at the zero; King's criterion,
   # |omega_l| < 1e-9 for every l >= L, does not.
   omega = np.exp(-np.arange(400) / 10.0)
   omega[100] = 0.0
   first_small = int(np.flatnonzero(np.abs(omega[2:]) < le.KING_THRESHOLD)[0]) + 2
   assert first_small == 100
   length = le.king_expansion_length(omega, 399, le.KING_THRESHOLD)
   assert length == int(np.flatnonzero(omega >= le.KING_THRESHOLD)[-1]) + 1 == 208


def test_degree_bound_contains_the_averaged_mie_phase_function():
   # Liquid water at 10.8 microns: every coefficient beyond the degree bound
   # of the size-distribution-averaged phase function vanishes.
   mmstr = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS)
   wl = 10.8
   ri = gsp._interpol_complex(mmstr.comp[0].cm, mmstr.comp[0].wl, [wl])
   degree = le.mie_phase_degree(2.0 * np.pi * 100.0 / wl)
   npts = le.initial_quadrature_order(degree)
   x, w = le.gauss_legendre(npts)
   _, _, _, phase, _ = create_bwgp("modified_gamma", 10.0, mmstr.s[0], ri, [wl], -x)
   omega = le.legendre_coefficients(-x, w, phase[:, 0])
   assert abs(omega[0] - 1.0) < 1e-12
   assert le.rounding_noise(omega, degree) < 1e-11
   assert le.king_expansion_length(omega, degree, le.KING_THRESHOLD) < degree


@pytest.mark.slow
def test_difficult_case_needs_and_gets_far_more_than_1000_moments(monkeypatch):
   mmstr = load_mmdat(INPUTS / "microphysics" / "water-ice_sph.mm", INPUTS)
   ri = gsp._interpol_complex(mmstr.comp[0].cm, mmstr.comp[0].wl, [MODIS_CHANNEL_8])
   lut_mrat = np.ones((1, 1))
   lut_rm = np.full((1, 1), 93.0)
   check_mu = np.cos(np.deg2rad(le.CHECK_THETA))

   # the fixed expansion: 1000 Gauss points, 1000 moments
   abscissas, weights = quadrature("g", 1000)
   qv = -abscissas
   _, _, g, phase, _ = create_bwgp("modified_gamma", 93.0, mmstr.s[0], ri, [MODIS_CHANNEL_8], np.concatenate((qv, check_mu)))
   inlc, omega_fixed = legpexp(1000, qv, weights, phase[:1000, 0])
   assert omega_fixed[0] < 0.8                                   # normalisation aliased (thesis omega_0 test fails)
   assert omega_fixed[1] / 3.0 < g[0] - 0.2                      # asymmetry aliased
   fixed_error = np.max(np.abs(le.reconstruct_phase_function(omega_fixed, 1000, check_mu) / phase[1000:, 0] - 1.0))
   assert fixed_error > 0.1

   # adaptive quadrature order and expansion length
   bext_c, w_c, g_c, vavg_c, phs, omega, lmom = gsp._adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, ri, MODIS_CHANNEL_8, "test")
   length = int(lmom[0])
   assert length > 3000
   assert abs(omega[0, 0] - 1.0) < 1e-9 and abs(omega[1, 0] / 3.0 - g_c[0, 0]) < 1e-9
   _, _, _, check_phase, _ = create_bwgp("modified_gamma", 93.0, mmstr.s[0], ri, [MODIS_CHANNEL_8], check_mu)
   error = np.max(np.abs(le.reconstruct_phase_function(omega[:, 0], length, check_mu) / check_phase[:, 0] - 1.0))
   assert error <= le.RECONSTRUCTION_TOLERANCE

   # the result is stable when the quadrature order is increased by a quarter
   original = le.initial_quadrature_order
   monkeypatch.setattr(gsp, "initial_quadrature_order", lambda degree: original(degree) + original(degree) // 4)
   _, _, _, _, phs2, omega2, lmom2 = gsp._adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, ri, MODIS_CHANNEL_8, "test")
   assert phs2.shape[0] > phs.shape[0]
   assert abs(int(lmom2[0]) - length) <= 0.01 * length
   common = min(length, int(lmom2[0]))
   chi = omega[:common, 0] / (2.0 * np.arange(common) + 1.0)
   chi2 = omega2[:common, 0] / (2.0 * np.arange(common) + 1.0)
   assert np.max(np.abs(chi - chi2)) < 1e-7


def _liquid_expansion(effective_radius, wavelength):
   mmstr = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS)
   ri = gsp._interpol_complex(mmstr.comp[0].cm, mmstr.comp[0].wl, [wavelength])
   return gsp._adaptive_mie_wavelength(mmstr, np.ones((1, 1)), np.full((1, 1), effective_radius), ri, wavelength, "test")


@pytest.mark.parametrize("effective_radius, wavelength, low, high", [
   (10.0, 2.1142, 150, 400),             # substantially fewer than 1000 moments (L = 220)
   (15.0, 0.6462791, 900, 1300),         # close to 1000 (L = 1093)
   (40.0, 0.6462791, 2500, 3000),        # more than 1000 (L = 2819)
])
def test_liquid_expansion_lengths_follow_the_averaged_phase_function(effective_radius, wavelength, low, high):
   bext_c, w_c, g_c, vavg_c, phs, omega, lmom = _liquid_expansion(effective_radius, wavelength)
   length = int(lmom[0])
   assert low <= length <= high
   assert abs(omega[0, 0] - 1.0) < 1e-9 and abs(omega[1, 0] / 3.0 - g_c[0, 0]) < 1e-9
   check_mu = np.cos(np.deg2rad(le.CHECK_THETA))
   _, _, _, direct, _ = create_bwgp("modified_gamma", effective_radius, 0.1111111,
                                    gsp._interpol_complex(*_liquid_tables(), [wavelength]), [wavelength], check_mu,
                                    radius_upper_factor=gsp.LIQUID_UPPER_RADIUS_FACTOR)
   error = np.max(np.abs(le.reconstruct_phase_function(omega[:, 0], length, check_mu) / direct[:, 0] - 1.0))
   assert error <= le.RECONSTRUCTION_TOLERANCE


def _liquid_tables():
   component = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS).comp[0]
   return component.cm, component.wl


def test_adaptive_stopping_is_reproducible():
   first = _liquid_expansion(15.0, 0.6462791)
   second = _liquid_expansion(15.0, 0.6462791)
   assert np.array_equal(first[6], second[6])
   assert np.array_equal(first[5], second[5]) and np.array_equal(first[4], second[4])
