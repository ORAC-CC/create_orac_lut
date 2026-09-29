"""Inexpensive tests for the quadrature-convergence analysis tools (no forward model)."""

import json

import numpy as np
import pytest

from oraclut.validation.quadrature import (
    ErrorMetrics, assert_subset, classify, error_metrics, interpolate, linear_nodes, log10_nodes,
    log_nodes, measurement_uncertainty, normalised_error, subset, to_jsonable, transform, write_json,
)


def test_grid_constructors_include_endpoints_and_production_nodes():
    assert linear_nodes(1.0, 39.0, 20).tolist() == list(range(1, 40, 2))
    octaves = log_nodes(0.0078125, 256.0, 1)
    assert octaves[0] == 0.0078125 and octaves[-1] == 256.0 and octaves.size == 16
    fine = log_nodes(0.0078125, 256.0, 8)
    assert fine.size == 121 and np.all(np.isin(np.round(octaves, 12), np.round(fine, 12)))
    aerosol = log10_nodes(0.01, 10.0, 20)
    assert aerosol[0] == pytest.approx(0.01) and aerosol[-1] == pytest.approx(10.0)


def test_subset_keeps_first_and_last_nodes():
    truth = np.arange(1.0, 39.5, 0.5)
    production = subset(truth, 4)
    assert production.tolist() == list(range(1, 40, 2))
    coarse = subset(truth, 8)
    assert coarse[0] == 1.0 and coarse[-1] == 39.0 and coarse.size == 11
    assert_subset(production, truth)
    with pytest.raises(ValueError):
        assert_subset([1.0, 2.25], truth)


def test_transforms():
    x = np.array([1.0, 10.0, 100.0])
    assert transform(x, "linear").tolist() == x.tolist()
    assert transform(x, "log10").tolist() == [0.0, 1.0, 2.0]
    np.testing.assert_allclose(transform(x, "log1p"), np.log1p(x))
    np.testing.assert_allclose(transform(np.array([0.0, 60.0]), "cos"), [1.0, 0.5])
    with pytest.raises(ValueError):
        transform(np.array([0.0]), "log10")
    with pytest.raises(ValueError):
        transform(x, "spline")


def test_interpolation_is_exact_on_nodes_and_linear_in_the_chosen_coordinate():
    truth_x = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    truth_y = np.log10(truth_x) * 3.0 + 1.0            # linear in log10(x)
    candidate_x = np.array([1.0, 4.0, 16.0])
    candidate_y = np.interp(candidate_x, truth_x, truth_y)
    in_log = interpolate(candidate_x, candidate_y, truth_x, "log10")
    np.testing.assert_allclose(in_log, truth_y, atol=1e-12)     # exact: truth is linear in log10
    in_lin = interpolate(candidate_x, candidate_y, truth_x, "linear")
    assert not np.allclose(in_lin, truth_y)                      # not exact in the wrong coordinate
    # trailing axes are carried through
    stacked = np.stack([candidate_y, 2 * candidate_y], axis=1)
    out = interpolate(candidate_x, stacked, truth_x, "log10")
    assert out.shape == (5, 2) and np.allclose(out[:, 1], 2 * truth_y)


def test_classification_of_off_node_points():
    truth_x = np.arange(0.0, 8.5, 0.5)
    candidate_x = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
    truth_y = truth_x ** 3                                        # curvature grows with x
    classes = classify(truth_x, candidate_x, truth_y, "linear")
    assert classes.on_node.sum() == 5 and set(truth_x[classes.on_node]) == {0, 2, 4, 6, 8}
    assert set(truth_x[classes.midpoint]) == {1.0, 3.0, 5.0, 7.0}
    assert set(truth_x[classes.edge]) == {0.5, 1.0, 1.5, 6.5, 7.0, 7.5}
    assert not np.any(classes.high_curvature & classes.on_node)
    assert truth_x[classes.high_curvature].min() >= 4.0          # highest curvature at large x
    assert classes.off_node.sum() == truth_x.size - 5


def test_error_metrics_and_relative_floor():
    truth = np.array([0.0, 0.5, 1.0, 2.0])
    interp = np.array([0.1, 0.5, 1.1, 2.0])
    m = error_metrics(truth, interp, truth_x=np.array([10.0, 20.0, 30.0, 40.0]))
    assert (m.points, m.max_abs) == (4, pytest.approx(0.1))
    assert m.mean_abs == pytest.approx(0.05) and m.rms_abs == pytest.approx(np.sqrt(0.02 / 4))
    assert m.max_rel == pytest.approx(0.1)            # the zero reference is excluded by the floor
    assert m.max_abs_location in (10.0, 30.0)
    masked = error_metrics(truth, interp, mask=np.array([False, True, False, True]))
    assert masked.max_abs == 0.0 and masked.points == 2
    assert error_metrics(truth, interp, mask=np.zeros(4, bool)).max_abs is None


def test_measurement_uncertainty_and_normalisation():
    ch1 = measurement_uncertainty(1, wavelength_um=0.638, solar=True, thermal=False,
                                  oldnefr=0.000274, snr=48.8, nedt=None, refbt=None)
    assert ch1.solar_reflectance_sigma == 0.000274 and ch1.solar_relative_sigma == pytest.approx(1 / 48.8)
    assert ch1.thermal_sigma_operator() is None
    ch9 = measurement_uncertainty(9, wavelength_um=10.7963, solar=False, thermal=True,
                                  oldnefr=None, snr=None, nedt=0.06, refbt=300.0)
    # dT/dE = lambda T^2 / c2 = 10.7963e-6 * 9e4 / 1.4388e-2 ~ 67.5 K per unit emissivity
    assert ch9.thermal_kelvin_per_unit_operator == pytest.approx(67.54, rel=1e-3)
    assert ch9.thermal_sigma_operator() == pytest.approx(0.06 / 67.54, rel=1e-3)
    metrics = ErrorMetrics(10, 1e-3, 5e-4, 4e-4, 0.01, 5.0)
    thermal = normalised_error(metrics, ch9, "thermal")
    assert thermal["max_over_sigma_thermal"] == pytest.approx(1e-3 * 67.54 / 0.06, rel=1e-3)
    assert thermal["max_equivalent_kelvin"] == pytest.approx(0.0675, rel=1e-3)
    solar = normalised_error(metrics, ch1, "solar")
    assert solar["max_over_sigma_reflectance"] == pytest.approx(1e-3 / 0.000274)
    assert solar["max_rel_over_inverse_snr"] == pytest.approx(0.01 * 48.8)
    assert normalised_error(ErrorMetrics(0, None, None, None, None, None), ch1, "solar") == {"max_over_sigma": None}


def test_serialisation_round_trip(tmp_path):
    payload = {"metrics": ErrorMetrics(3, 1.0, 0.5, 0.25, None, 2.0), "nodes": np.array([1.0, 2.0]),
               "flag": np.bool_(True), "nested": {"x": np.float32(1.5)}}
    path = write_json(tmp_path / "out" / "summary.json", payload)
    data = json.loads(path.read_text())
    assert data["metrics"]["max_abs"] == 1.0 and data["nodes"] == [1.0, 2.0]
    assert data["flag"] is True and data["nested"]["x"] == 1.5
    assert to_jsonable([np.int64(3)]) == [3]


# ---------------------------------------------------------------------------
# Refinement-experiment tools
# ---------------------------------------------------------------------------

from oraclut.validation.quadrature import (  # noqa: E402
    Acceptance, acceptance, cache_key, cell_fractions, cell_midpoints, equidistributed_nodes,
    interpolate_2d, merge_nodes, sectioned_nodes,
)


def test_sectioned_truth_grid_has_each_boundary_once_and_required_spacings():
    r = sectioned_nodes([(1.0, 2.0, 0.05), (2.0, 5.0, 0.1), (5.0, 10.0, 0.25), (10.0, 20.0, 0.5), (20.0, 39.0, 1.0)])
    assert r.size == 110 and r[0] == 1.0 and r[-1] == 39.0
    assert np.count_nonzero(r == 2.0) == 1 and np.count_nonzero(r == 5.0) == 1 and np.count_nonzero(r == 10.0) == 1
    assert np.allclose(np.diff(r[r <= 2.0]), 0.05) and np.allclose(np.diff(r[(r >= 20.0)]), 1.0)
    assert np.allclose(np.diff(r[(r >= 5.0) & (r <= 10.0)]), 0.25)
    with pytest.raises(ValueError):
        sectioned_nodes([(1.0, 2.0, 0.3)])          # not an integer number of steps


def test_low_tau_nodes_merge_with_the_main_log_grid_without_duplicates():
    low = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
    tau = merge_nodes(low, log_nodes(0.0078125, 256.0, 12))
    assert tau.size == 13 + 15 * 12 + 1
    assert tau[0] == 1e-10 and tau[-1] == 256.0 and np.all(np.diff(tau) > 0)
    assert np.all(np.isin(np.round([1e-10, 1e-3, 5e-3, 0.0078125, 1.0, 256.0], 12), np.round(tau, 12)))
    # production nodes 2^k are exact subsets of the dense truth
    assert_subset([1e-10] + [2.0 ** k for k in range(-7, 9)], tau)


def test_small_radius_candidates_are_truth_subsets():
    truth = sectioned_nodes([(1.0, 2.0, 0.05), (2.0, 5.0, 0.1), (5.0, 10.0, 0.25), (10.0, 20.0, 0.5), (20.0, 39.0, 1.0)])
    production = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39]
    refined = [1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 35, 39]
    assert_subset(production, truth) is not None
    assert_subset(refined, truth) is not None
    assert_subset([1.0, 2.3], truth)                 # 2.3 is a valid 0.1 node
    with pytest.raises(ValueError):
        assert_subset([1.0, 1.33], truth)            # 1.33 is not


def test_equidistributed_design_concentrates_nodes_where_curvature_is_high():
    x = np.linspace(1.0, 39.0, 381)
    y = 1.0 / x                                          # curvature falls steeply with x
    nodes = equidistributed_nodes(x, y, 12)
    assert nodes[0] == 1.0 and nodes[-1] == 39.0 and nodes.size == 12
    assert np.all(np.isin(nodes, x))
    assert np.sum(nodes < 5.0) >= 5                     # most nodes at small x
    assert np.diff(nodes)[-1] > np.diff(nodes)[0] * 5   # cells widen where the function is flat
    with pytest.raises(ValueError):
        equidistributed_nodes(x, y, 1)


def test_midpoints_and_quarter_points_are_never_nodes_and_are_disjoint():
    nodes = np.array([1.0, 2.0, 4.0, 8.0])
    mids_lin = cell_midpoints(nodes, "linear"); mids_log = cell_midpoints(nodes, "log10")
    quarters = cell_fractions(nodes, 0.25, "log10")
    assert mids_lin.tolist() == [1.5, 3.0, 6.0]
    np.testing.assert_allclose(mids_log, [np.sqrt(2), np.sqrt(8), np.sqrt(32)])
    assert not np.any(np.isin(np.round(mids_log, 9), np.round(nodes, 9)))
    assert not np.any(np.isclose(quarters[:, None], mids_log[None, :]))   # independent set disjoint from midpoints
    with pytest.raises(ValueError):
        cell_fractions(nodes, 1.0)


def test_2d_interpolation_is_bilinear_in_the_chosen_coordinates():
    rx = np.array([1.0, 3.0, 5.0]); ty = np.array([0.5, 2.0, 8.0])
    # f linear in r and in log10(tau): bilinear interpolation must be exact
    f = lambda r, t: 2.0 * r[:, None] + 3.0 * np.log10(t)[None, :]
    node_values = f(rx, ty)
    qx = np.array([1.5, 2.0, 4.5]); qy = np.array([1.0, 4.0])
    out = interpolate_2d(rx, ty, node_values, qx, qy, "linear", "log10")
    np.testing.assert_allclose(out, f(qx, qy), atol=1e-12)
    assert out.shape == (3, 2)
    wrong = interpolate_2d(rx, ty, node_values, qx, qy, "linear", "linear")
    assert not np.allclose(wrong, f(qx, qy))
    # trailing axes carried through
    stacked = np.stack([node_values, -node_values], axis=-1)
    assert interpolate_2d(rx, ty, stacked, qx, qy, "linear", "log10").shape == (3, 2, 2)


def test_acceptance_requires_every_channel_within_half_sigma():
    ok = acceptance({1: 0.3, 4: 0.45, 9: 0.5})
    assert isinstance(ok, Acceptance) and ok.accepted and ok.worst_channel == 9 and ok.worst_ratio == 0.5
    bad = acceptance({1: 0.3, 4: 0.51, 9: 0.2})
    assert not bad.accepted and bad.worst_channel == 4
    unjudged = acceptance({1: 0.3, 4: None, 9: 0.2})
    assert unjudged.accepted and unjudged.ratios["4"] is None
    assert not acceptance({4: None}).accepted
    assert not acceptance({1: 0.6}, limit_fraction=0.5).accepted and acceptance({1: 0.6}, limit_fraction=1.0).accepted


def test_cache_key_is_deterministic_and_separates_configurations_and_states():
    base = dict(instrument="i", microphysics="m", atmosphere=2, streams=60, phase_order=1000, rayleigh=True, gas=False)
    kw = dict(channel_set=[1, 4, 9], solar_zenith=[30.0, 80.0], satellite_zenith=[40.0], relative_azimuth=[90.0])
    k1 = cache_key(base, effective_radius=2.0, optical_depth=0.5, **kw)
    assert k1 == cache_key(dict(reversed(list(base.items()))), effective_radius=2.0, optical_depth=0.5, **kw)
    assert k1 != cache_key(base, effective_radius=2.05, optical_depth=0.5, **kw)
    assert k1 != cache_key(base, effective_radius=2.0, optical_depth=0.50001, **kw)
    assert k1 != cache_key({**base, "streams": 32}, effective_radius=2.0, optical_depth=0.5, **kw)
    assert k1 != cache_key(base, effective_radius=2.0, optical_depth=0.5, **{**kw, "solar_zenith": [30.0, 70.0]})
    assert len(k1) == 64
