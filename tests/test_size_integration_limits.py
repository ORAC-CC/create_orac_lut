"""Liquid-water modified-gamma size-integration limit (3.5 x effective radius on the legacy lattice).

Validation evidence: validation/size_distribution_limits/ and
validation/REPORT_lut_numerics_development.md.
"""

import copy
import importlib
from pathlib import Path

import numpy as np
import pytest

from oraclut.idl_mirror.load_mmdat import load_mmdat

bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")

ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"


def _legacy_nodes(wavelength):
   npts = max(int(2.0 * np.pi * (100.0 - 0.001) / wavelength / 0.4), 200)
   absc, wght = bwgp.quadrature("T", npts)
   return bwgp.shift_quadrature(absc, wght, 0.001, 100.0)[0]


@pytest.mark.parametrize("wavelength", [0.6462791, 2.1142, 11.0262])
@pytest.mark.parametrize("effective_radius", [5.0, 10.0, 30.0, 40.0])
def test_lattice_keeps_the_legacy_nodes(wavelength, effective_radius):
   factor = gsp.LIQUID_UPPER_RADIUS_FACTOR
   upper, npts = bwgp.lattice_upper_radius(effective_radius, factor, 1.0 / wavelength)
   legacy = _legacy_nodes(wavelength)
   spacing = legacy[1] - legacy[0]
   assert upper >= factor * effective_radius and upper - spacing < factor * effective_radius
   absc, wght = bwgp.quadrature("T", npts)
   nodes = bwgp.shift_quadrature(absc, wght, 0.001, upper)[0]
   assert np.allclose(np.diff(nodes), spacing, rtol=1e-9, atol=0.0)
   common = min(nodes.size, legacy.size)
   assert np.allclose(nodes[:common], legacy[:common], rtol=0.0, atol=1e-10)
   if effective_radius >= 30.0:
      assert upper > 100.0 and nodes.size > legacy.size        # r_e = 30 and 40 um integrate beyond 100 um


def test_liquid_limit_changes_only_the_tail():
   ri = complex(np.float32(1.2856), -np.float32(3.3e-4))       # liquid water near 2.1 um
   mu = -np.polynomial.legendre.leggauss(64)[0]
   legacy = bwgp.create_bwgp("modified_gamma", 10.0, 0.1111111, [ri], [2.1142], mu)
   adopted = bwgp.create_bwgp("modified_gamma", 10.0, 0.1111111, [ri], [2.1142], mu,
                              radius_upper_factor=gsp.LIQUID_UPPER_RADIUS_FACTOR)
   rel_ext = abs(adopted[0][0] / legacy[0][0] - 1.0)
   assert 0.0 < rel_ext < 1e-5
   assert abs(adopted[1][0] - legacy[1][0]) < 2e-6 and abs(adopted[2][0] - legacy[2][0]) < 2e-6
   assert np.max(np.abs(adopted[3][:, 0] / legacy[3][:, 0] - 1.0)) < 1e-4


def test_other_distributions_keep_the_legacy_integration():
   mu = -np.polynomial.legendre.leggauss(32)[0]
   ri = [complex(1.45, -1e-3)]
   for distname, rm, s in (("log_normal", 0.07, 1.7), ("modified_gamma", 20.0, 0.1111111)):
      plain = bwgp.create_bwgp(distname, rm, s, ri, [0.865], mu)
      if distname == "log_normal":     # the factor is ignored for log-normal limits
         factored = bwgp.create_bwgp(distname, rm, s, ri, [0.865], mu, radius_upper_factor=3.0)
         for a, b in zip(plain[:4], factored[:4]):
            assert np.array_equal(a, b)
      direct = bwgp.mie_size_dist_new(distname, 1.0, [rm, s, 0.001, 100.0], 1.0 / 0.865, ri[0], mu, xres=0.4)
      assert plain[0][0] == direct[0] and np.array_equal(plain[3][:, 0], direct[4][0, :])


def test_rule_applies_to_liquid_water_modified_gamma_only():
   liquid = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS)
   narrow = load_mmdat(INPUTS / "microphysics" / "liquid-water_ap02.mm", INPUTS)
   ice = load_mmdat(INPUTS / "microphysics" / "water-ice_sph.mm", INPUTS)
   aerosol = load_mmdat(INPUTS / "microphysics" / "aerosol_a79.mm", INPUTS)
   assert gsp.radius_upper_factor(liquid, 0) == 3.5
   assert gsp.radius_upper_factor(narrow, 0) == 3.5
   assert gsp.radius_upper_factor(ice, 0) is None
   assert all(gsp.radius_upper_factor(aerosol, c) is None for c in range(aerosol.ncomp))
   wide = copy.deepcopy(liquid)
   wide.s = np.array([0.2])
   with pytest.raises(ValueError, match="validated for effective variance"):
      gsp.radius_upper_factor(wide, 0)


def test_explicit_node_count_refuses_a_truncated_upper_radius():
   with pytest.raises(ValueError, match="explicit node count"):
      bwgp.mie_size_dist_new("modified_gamma", 1.0, [40.0, 0.1111111, 0.001, 2000.0], 1.0 / 0.4, complex(1.33, -1e-8),
                             np.array([1.0, 0.0, -1.0]), xres=0.4, npts=300)
