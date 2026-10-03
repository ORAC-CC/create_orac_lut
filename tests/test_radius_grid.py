"""Refined radius integration of liquid-water and ice-sphere Mie size distributions (V24).

Validation evidence: validation/radius_grid/ and
validation/REPORT_lut_numerics_development.md.
"""

import importlib
from pathlib import Path

import numpy as np
import pytest

from oraclut.idl_mirror.load_mmdat import load_mmdat

bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")

ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
VARIANCE = 0.1111111
LIQUID = gsp.LIQUID_REFINED_XRES
ICE = gsp.ICE_SPHERE_REFINED_XRES


def _nodes(params, npts):
   absc, wght = bwgp.quadrature("T", npts)
   return bwgp.shift_quadrature(absc, wght, params[2], params[3])[0]


def test_blocked_mie_sum_is_bitwise_the_single_call_sum(monkeypatch):
   mu = np.cos(np.deg2rad(np.linspace(0.0, 180.0, 181)))
   for distname, params, wavelength, ri in (("modified_gamma", [10.0, VARIANCE, 0.001, 35.0], 0.6462791, complex(1.332, -1e-8)),
                                            ("log_normal", [0.5, 1.8, 0.001, 100.0], 0.55, complex(1.5, -0.01))):
      single = bwgp.mie_size_dist_new(distname, 1.0, params, 1.0 / wavelength, ri, mu, xres=0.4)
      monkeypatch.setattr(bwgp, "MIE_BLOCK_ELEMENTS", 97 * mu.size)
      blocked = bwgp.mie_size_dist_new(distname, 1.0, params, 1.0 / wavelength, ri, mu, xres=0.4)
      monkeypatch.undo()
      for a, b in zip(single, blocked):
         assert np.array_equal(a, b)


@pytest.mark.parametrize("wavelength, effective_radius, xres, level", [
   (0.4118525, 10.0, LIQUID, 4),       # legacy step 0.4 in x -> 0.025
   (0.6462791, 1.0, LIQUID, 4),
   (2.1142, 40.0, LIQUID, 4),
   (0.6462791, 30.0, ICE, 3),          # -> 0.05
   (11.0262, 10.0, LIQUID, 4),         # 200-node floor, h = 0.5 um (x step 0.285) -> 0.018
   (11.0262, 10.0, ICE, 3),            # -> 0.036
   (11.0262, 0.5, ICE, 4),             # distribution spread: re sqrt(ve) / 3 = 0.056 um
])
def test_refinement_level(wavelength, effective_radius, xres, level):
   wavenumber = 1.0 / wavelength
   spacing = bwgp.legacy_radius_spacing(wavenumber)
   assert bwgp.radius_refinement_level(spacing, wavenumber, effective_radius, VARIANCE, xres) == level
   h = spacing / 2**level
   spread = effective_radius * np.sqrt(VARIANCE) / bwgp.MIE_NODES_PER_SPREAD
   assert 2.0 * np.pi * h * wavenumber <= 1.02 * xres and h <= 1.02 * spread
   coarser = 2.0 * h
   assert 2.0 * np.pi * coarser * wavenumber > 1.02 * xres or coarser > 1.02 * spread


@pytest.mark.parametrize("factor", [gsp.LIQUID_UPPER_RADIUS_FACTOR, None])      # liquid lattice, ice 0.001-100 um
@pytest.mark.parametrize("wavelength", [0.6462791, 11.0262])
def test_refined_grid_nests_the_legacy_grid(factor, wavelength):
   wavenumber = 1.0 / wavelength
   legacy_params, legacy_npts = bwgp.mie_integration_limits("modified_gamma", 20.0, VARIANCE, wavenumber, factor)
   params, npts = bwgp.mie_integration_limits("modified_gamma", 20.0, VARIANCE, wavenumber, factor, refined_xres=LIQUID)
   assert params == legacy_params                      # limits unchanged
   if legacy_npts is None:
      legacy_npts = max(int(2.0 * np.pi * (100.0 - 0.001) * wavenumber / 0.4), 200)
   level = bwgp.radius_refinement_level(bwgp.legacy_radius_spacing(wavenumber), wavenumber, 20.0, VARIANCE, LIQUID)
   assert level > 0 and npts - 1 == (legacy_npts - 1) * 2**level
   legacy = _nodes(params, legacy_npts)
   refined = _nodes(params, npts)
   assert np.allclose(refined[::2**level], legacy, rtol=0.0, atol=1e-10)


def test_log_normal_ignores_the_refinement():
   mu = -np.polynomial.legendre.leggauss(32)[0]
   ri = [complex(1.45, -1e-3)]
   plain = bwgp.create_bwgp("log_normal", 0.07, 1.7, ri, [0.865], mu)
   refined = bwgp.create_bwgp("log_normal", 0.07, 1.7, ri, [0.865], mu, refined_xres=LIQUID)
   for a, b in zip(plain[:4], refined[:4]):
      assert np.array_equal(a, b)


def test_refinement_converges_towards_a_finer_grid(monkeypatch):
   # liquid water near 2.1 um, r_e = 10 um: the legacy grid is in error by
   # ~2e-4 in extinction and ~1e-1 in back-scattering (validation/radius_grid/)
   ri = [complex(np.float32(1.2856), -np.float32(3.3e-4))]
   theta = np.array([0.0, 30.0, 90.0, 140.0, 170.0, 180.0])
   mu = np.cos(np.deg2rad(theta))
   factor = gsp.LIQUID_UPPER_RADIUS_FACTOR
   legacy = bwgp.create_bwgp("modified_gamma", 10.0, VARIANCE, ri, [2.1142], mu, radius_upper_factor=factor)
   refined = bwgp.create_bwgp("modified_gamma", 10.0, VARIANCE, ri, [2.1142], mu, radius_upper_factor=factor,
                              refined_xres=LIQUID)
   monkeypatch.setattr(bwgp, "radius_refinement_level", lambda *args: 6)
   finest = bwgp.create_bwgp("modified_gamma", 10.0, VARIANCE, ri, [2.1142], mu, radius_upper_factor=factor,
                             refined_xres=LIQUID)
   error = lambda a: (abs(a[0][0] / finest[0][0] - 1.0), abs(a[2][0] - finest[2][0]),
                      np.max(np.abs(a[3][:, 0] / finest[3][:, 0] - 1.0)))
   legacy_error, refined_error = error(legacy), error(refined)
   assert refined_error[0] < 1e-5 and refined_error[1] < 1e-5 and refined_error[2] < 1e-2
   assert legacy_error[2] > 10.0 * refined_error[2]


def test_refinement_applies_to_liquid_water_and_ice_spheres_only():
   for name, expected in (("liquid-water_stg.mm", [LIQUID]), ("liquid-water_old.mm", [LIQUID]), ("water-ice_sph.mm", [ICE]),
                          ("water-ice_agg.mm", [None]), ("aerosol_a79.mm", [None, None]), ("aerosol_a75.mm", [None, None, None])):
      mmstr = load_mmdat(INPUTS / "microphysics" / name, INPUTS)
      assert [gsp.refined_xres(mmstr, c) for c in range(mmstr.ncomp)] == expected
