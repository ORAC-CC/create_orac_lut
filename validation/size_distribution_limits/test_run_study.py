"""Focused checks for the standalone integration-limit study."""

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_study import CASES, REFERENCE_BOUNDS, distribution_density, radius_quadrature


def test_modified_gamma_parameter_is_effective_radius_on_full_domain():
    case = CASES[1]
    radius, weights = radius_quadrature(*REFERENCE_BOUNDS, wavelength=0.55)
    density = distribution_density(case, radius)
    calculated = np.sum(weights * density * radius**3) / np.sum(weights * density * radius**2)
    np.testing.assert_allclose(calculated, case.effective_radius, rtol=2e-7)


def test_lognormal_legacy_bounds_scale_with_mode_radius():
    small = next(case for case in CASES if case.name == "aerosol_small")
    lower, upper = small.legacy_bounds
    assert 0.0 < lower < small.radius_parameter < upper < 5.0


def test_radius_quadrature_matches_production_point_rule_and_weights():
    radius, weights = radius_quadrature(0.001, 100.0, wavelength=0.55)
    assert radius.size == max(int(2.0 * np.pi * (100.0 - 0.001) / 0.55 / 0.4), 200)
    np.testing.assert_allclose(np.sum(weights), 100.0 - 0.001)
    assert weights[0] == weights[-1] == 0.5 * weights[1]
