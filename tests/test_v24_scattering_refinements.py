"""LUT release V24: the three scattering-property refinements of Mie classes.

1. liquid-water size integration to 3.5 effective radii on the legacy radius lattice;
2. adaptive Legendre expansion of every averaged phase function (legendre_expansion.py);
3. nested refinement of the liquid-water and ice-sphere radius integration;
and the driver rule that DISORT receives max(nmom, L) moments without truncation.
"""

import importlib
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

import create_orac_luts
from oraclut.idl_mirror import load_mmdat

cb = importlib.import_module("oraclut.idl_mirror.create_bwgp")
gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
le = importlib.import_module("oraclut.idl_mirror.legendre_expansion")

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"


def _pinned_v23_create_bwgp():
   """create_bwgp_v03.py, the file of LUT release V23 (releases/lut_v23.release)."""
   path = ROOT / "src" / "oraclut" / "idl_mirror" / "create_bwgp_v03.py"
   spec = importlib.util.spec_from_file_location("oraclut.idl_mirror._create_bwgp_v23", path)
   module = importlib.util.module_from_spec(spec)
   module.__package__ = "oraclut.idl_mirror"
   spec.loader.exec_module(module)
   return module


def _identical(a, b):
   return all(np.asarray(x).dtype == np.asarray(y).dtype and np.atleast_1d(x).tobytes() == np.atleast_1d(y).tobytes()
              for x, y in zip(a, b))


# ---------------------------------------------------------------------------
# legendre_expansion.py
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("npts", [64, 1000, 2048])
def test_gauss_legendre_integrates_polynomials_below_degree_2n_exactly(npts):
   x, w = le.gauss_legendre(npts)
   assert np.all(np.diff(x) > 0) and abs(np.sum(w) - 2.0) < 1e-13
   # orthogonality of P_k and P_j (degree k + j < 2N) through the rule
   for k, j in ((0, 2 * npts - 2), (npts - 1, npts - 1), (npts // 2, npts + 3)):
      if k + j < 2 * npts:
         pk = np.polynomial.legendre.legval(x, np.eye(k + 1)[k])
         pj = np.polynomial.legendre.legval(x, np.eye(j + 1)[j])
         exact = 2.0 / (2 * k + 1) if k == j else 0.0
         assert abs(np.sum(w * pk * pj) - exact) < 1e-12


def test_degree_bound_is_twice_the_mieint_series_length():
   # mieint.f: NStop = Dx + 4.00 Dx^(1/3) + 2 (Dx <= 8), + 4.05 (8 < Dx < 4200), 2 for Dx < 0.02
   for x in (0.01, 3.0, 8.0, 100.0, 2082.0, 4199.0):
      nstop = 2 if x < 0.02 else int(x + (4.00 if x <= 8.0 or x >= 4200.0 else 4.05) * x ** (1.0 / 3.0) + 2.0)
      assert le.mie_nstop(x) >= nstop and le.degree_bound(x) >= 2 * nstop
   degree = le.degree_bound(500.0)
   order = le.quadrature_order(degree)
   assert order % le.QUADRATURE_ORDER_STEP == 0 and order > degree + le.noise_band(degree)


def test_henyey_greenstein_coefficients_and_king_length():
   g = 0.85
   degree = 400                                       # g^l < 1e-9 / (2l + 1) well before l = 400
   order = le.quadrature_order(degree)
   x, w = le.gauss_legendre(order)
   mu = -x
   phase = (1.0 - g * g) / (1.0 + g * g - 2.0 * g * mu) ** 1.5
   omega = le.legendre_coefficients(mu, w, phase)
   l = np.arange(150)
   assert np.max(np.abs(omega[:150] - (2 * l + 1) * g ** l)) < 1e-10
   length = le.king_length(omega, degree)
   exact = next(n for n in range(1, 400) if all((2 * k + 1) * g ** k < le.KING_THRESHOLD for k in range(n, 400)))
   assert length == exact


# ---------------------------------------------------------------------------
# create_bwgp.py: the IDL path is unchanged; lattice limit and nested refinement
# ---------------------------------------------------------------------------

def test_default_create_bwgp_is_bitwise_the_v23_release_including_blocked_mie_sums(monkeypatch):
   old = _pinned_v23_create_bwgp()
   qv = -old.quadrature("G", 400)[0]
   cases = (("modified_gamma", 11.0, 0.1111111, 0.646, 1.332 - 1e-8j), ("modified_gamma", 1.0, 0.1111111, 11.03, 1.15 - 0.09j),
            ("log_normal", 0.12, 1.8, 0.55, 1.45 - 0.001j))
   for dist, rm, s, wl, ri in cases:
      a = old.create_bwgp(dist, rm, s, np.asarray([ri], np.complex64), [wl], qv)
      b = cb.create_bwgp(dist, rm, s, np.asarray([ri], np.complex64), [wl], qv)
      assert _identical(a, b), (dist, wl)
   monkeypatch.setattr(cb, "MIE_BLOCK_ELEMENTS", 400 * 37)     # 37 radius nodes per Mie call
   a = old.create_bwgp("modified_gamma", 21.0, 0.1111111, np.asarray([1.332 - 1e-8j], np.complex64), [0.646], qv)
   b = cb.create_bwgp("modified_gamma", 21.0, 0.1111111, np.asarray([1.332 - 1e-8j], np.complex64), [0.646], qv)
   assert _identical(a, b)


def test_lattice_limit_keeps_the_legacy_nodes_and_ends_at_the_first_node_beyond_the_factor():
   wavenumber = float(np.float32(1.0) / np.float32(0.646))
   rl, ru = float(cb.MIE_RADIUS_LOWER), float(cb.MIE_RADIUS_UPPER)
   n_legacy = cb.legacy_npts(rl, ru, wavenumber, cb.MIE_XRES)
   spacing = (ru - rl) / (n_legacy - 1)
   for re in (5.0, 21.0, 39.0):
      upper, npts = cb.lattice_upper_radius(re, 3.5, wavenumber)
      assert upper >= 3.5 * re > upper - spacing
      assert abs((upper - rl) / (npts - 1) - spacing) < 1e-12 * spacing
   assert cb.lattice_upper_radius(39.0, 3.5, wavenumber)[0] > 100.0       # beyond the IDL's 100 um when needed


def test_nested_refinement_keeps_the_legacy_nodes_and_meets_both_bounds():
   for wl, re, s, xres in ((0.646, 10.0, 0.1111111, 0.025), (11.03, 1.0, 0.1111111, 0.025), (0.646, 41.0, 0.1111111, 0.05)):
      wavenumber = float(np.float32(1.0) / np.float32(wl))
      factor = 3.5 if xres == 0.025 else None
      params, npts = cb.mie_integration_limits("modified_gamma", re, s, wavenumber, factor, xres)
      legacy_params, legacy_npts = cb.mie_integration_limits("modified_gamma", re, s, wavenumber, factor, None)
      legacy_npts = legacy_npts or cb.legacy_npts(legacy_params[2], legacy_params[3], wavenumber, cb.MIE_XRES)
      assert params == legacy_params                                     # unchanged limits
      level = int(round(math.log2((npts - 1) / (legacy_npts - 1))))
      assert npts == (legacy_npts - 1) * 2 ** level + 1                  # every legacy node is kept
      h = (params[3] - params[2]) / (npts - 1)
      assert 2.0 * math.pi * h * wavenumber <= 1.02 * xres and h <= 1.02 * re * math.sqrt(s) / cb.MIE_NODES_PER_SPREAD
      assert level == 0 or 2.0 * h * 2.0 * math.pi * wavenumber > 1.02 * xres or 2.0 * h > 1.02 * re * math.sqrt(s) / 3.0


def test_the_v24_rules_select_their_classes():
   liquid = load_mmdat(INPUTS / "microphysics" / "liquid-water_240.mm", INPUTS)
   spheres = load_mmdat(INPUTS / "microphysics" / "water-ice_sph.mm", INPUTS)
   aerosol = load_mmdat(INPUTS / "microphysics" / "aerosol_a79.mm", INPUTS)
   assert gsp.radius_upper_factor(liquid, 0) == 3.5 and gsp.refined_xres(liquid, 0) == 0.025
   assert gsp.radius_upper_factor(spheres, 0) is None and gsp.refined_xres(spheres, 0) == 0.05
   assert all(gsp.radius_upper_factor(aerosol, c) is None and gsp.refined_xres(aerosol, c) is None for c in range(aerosol.ncomp))
   assert gsp.uses_adaptive_legendre(liquid) and gsp.uses_adaptive_legendre(spheres)
   baum = load_mmdat(INPUTS / "microphysics" / "water-ice_agg.mm", INPUTS)
   assert not gsp.uses_adaptive_legendre(baum)
   liquid.s[0] = 0.2
   with pytest.raises(ValueError, match="effective variance"):
      gsp.radius_upper_factor(liquid, 0)


# ---------------------------------------------------------------------------
# The adaptive expansion of a real averaged phase function (cheap case)
# ---------------------------------------------------------------------------

def test_adaptive_expansion_of_a_thermal_liquid_phase_function():
   mm = load_mmdat(INPUTS / "microphysics" / "liquid-water_240.mm", INPUTS)
   ri = gsp._interpol_complex(mm.comp[0].cm, mm.comp[0].wl, [11.03])
   bext, w, g, vavg, phs, omega, lmom, diagnostics = gsp._adaptive_mie_wavelength(
      mm, np.ones((1, 2)), np.asarray([[2.0, 5.0]]), ri, 11.03, "test")
   for r in range(2):
      assert diagnostics["error"][r] <= le.RECONSTRUCTION_TOLERANCE and diagnostics["noise"][r] <= le.QUADRATURE_NOISE_LIMIT
      assert 2 < lmom[r] <= diagnostics["degree"][r] + 1 < diagnostics["order"][r]
      assert abs(omega[0, r] - 1.0) < 1e-6 and abs(omega[1, r] / 3.0 - g[0, r]) < 1e-6
      assert np.all(np.abs(omega[lmom[r]:int(diagnostics["degree"][r]) + 1, r]) < le.KING_THRESHOLD)
   assert lmom[1] > lmom[0]                                              # larger particles need more moments


# ---------------------------------------------------------------------------
# Driver: DISORT receives max(nmom, L) moments, never a truncation
# ---------------------------------------------------------------------------

def test_disort_moments_are_max_of_nmom_and_l_and_never_truncated():
   moments = np.zeros(3000, np.float32)
   moments[:2500] = np.linspace(1.0, 1e-6, 2500, dtype=np.float32)
   n, pm = create_orac_luts._disort_moments(1000, 2500, moments)
   assert n == 2500 and np.array_equal(pm, moments[:2500])
   short = np.zeros(1000, np.float32)
   short[:300] = moments[:300]
   n, pm = create_orac_luts._disort_moments(1000, 300, short)
   assert n == 1000 and np.array_equal(pm, short)                        # the IDL's 1000-moment array
   n, pm = create_orac_luts._disort_moments(1000, 300, moments[:400] * 0 + np.pad(moments[:300], (0, 100)))
   assert n == 1000 and np.array_equal(pm[:300], moments[:300]) and not pm[300:].any()
   bad = moments.copy()
   with pytest.raises(ValueError, match="beyond"):
      create_orac_luts._disort_moments(1000, 1000, bad)                 # a non-zero moment beyond those passed
