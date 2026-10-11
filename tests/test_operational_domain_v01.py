"""Tests for the ORAC operational-angular-domain logic used in the V21/V22 reanalysis."""

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "validation" / "v21_v22_terra_modis" / "operational_domain.py"
SPEC = importlib.util.spec_from_file_location("operational_domain", MODULE_PATH)
DOMAIN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DOMAIN)

SZA = np.arange(0.0, 91.0, 10.0)


def test_bracket_matches_orac_locate_plus_one():
   # ORAC: lower node from locate(bounded), upper node = lower + 1.
   assert DOMAIN.orac_bracket(SZA, 74.999) == (7, np.float64(74.999 - 70.0) / 10.0)
   low, w = DOMAIN.orac_bracket(SZA, 70.0)
   assert (low, w) == (7, 0.0)
   low, w = DOMAIN.orac_bracket(SZA, 90.0)          # bounded: last interval
   assert (low, w) == (8, 1.0)
   low, w = DOMAIN.orac_bracket(SZA, 0.0)
   assert (low, w) == (0, 0.0)


def test_strict_solar_limit_classes():
   classes = {c["node"]: c for c in DOMAIN.classify_nodes(SZA, 75.0, inclusive=False)}
   assert all(classes[n]["class"] == "in_domain" for n in SZA[SZA <= 70.0])
   assert classes[80.0]["class"] == "support"
   assert 0.49 < classes[80.0]["maximum_weight"] < 0.5
   assert classes[90.0]["class"] == "unused"
   assert classes[90.0]["maximum_weight"] == 0.0


def test_limit_on_a_node():
   grid = np.array([0.0, 25.0, 50.0, 75.0, 90.0])
   strict = {c["node"]: c["class"] for c in DOMAIN.classify_nodes(grid, 75.0, inclusive=False)}
   inclusive = {c["node"]: c["class"] for c in DOMAIN.classify_nodes(grid, 75.0, inclusive=True)}
   assert strict[75.0] == "support" and strict[90.0] == "unused"
   assert inclusive[75.0] == "in_domain" and inclusive[90.0] == "unused"


def test_operational_weights_and_interpolation():
   geometries, weights = DOMAIN.operational_sza_weights(SZA, 75.0)
   assert geometries.tolist() == [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 72.5, 75.0]
   np.testing.assert_allclose(weights.sum(axis=1), 1.0)
   assert weights[:, 9].max() == 0.0                       # 90 deg never used
   np.testing.assert_allclose(weights[-1, 7:9], [0.5, 0.5])
   values = np.random.default_rng(1).normal(size=(3, SZA.size, 4))
   interpolated = DOMAIN.apply_sza(values, 1, weights)
   np.testing.assert_allclose(interpolated[:, -1], 0.5 * (values[:, 7] + values[:, 8]))
   np.testing.assert_allclose(interpolated[:, :8], values[:, :8])
