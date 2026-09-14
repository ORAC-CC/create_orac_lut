"""The particle vertical profile must be interpolated exactly as IDL INTERPOL does.

Expected values were established by running IDL 8.9 INTERPOL directly (see
docs/validation_matrix.md): linear inside the tabulated range, linear
extrapolation from the two nearest end nodes on either side, identical for
ascending and descending abscissae.  Every extrapolation case here fails under
the previous ``np.interp`` implementation, which clamps to the end values.
"""

from pathlib import Path

import numpy as np
import pytest

from oraclut.config import read_atmosphere, read_microphysics
from oraclut.pipeline import _atmosphere_path, _interpolate_profile


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"

# aerosol_a79.mm profile as loaded (descending heights), and the IDL results for
# these query heights: INTERPOL(v, h, x) -> 1085.23 932.015 778.8 625.585 472.37 0 0 0 0 0
A79_HEIGHT = np.array([5.5, 4.5, 3.5, 2.5, 1.5], dtype=np.float32)
A79_AMOUNT = np.array([0.0, 0.0, 0.0, 472.37, 778.80], dtype=np.float32)
QUERY = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 5.5, 6.5, 26.25, 97.5], dtype=np.float32)
IDL_RESULT = np.array([1085.23, 932.015, 778.8, 625.585, 472.37, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32)


def test_matches_idl_interpol_for_the_a79_profile_descending_input():
    np.testing.assert_allclose(_interpolate_profile(A79_HEIGHT, A79_AMOUNT, QUERY), IDL_RESULT, rtol=1e-6, atol=1e-3)


def test_ascending_input_gives_the_same_result_as_idl():
    result = _interpolate_profile(A79_HEIGHT[::-1], A79_AMOUNT[::-1], QUERY)
    np.testing.assert_allclose(result, IDL_RESULT, rtol=1e-6, atol=1e-3)


def test_interior_and_exact_end_nodes():
    x = np.array([1.0, 2.0, 3.0, 4.0]); y = np.array([1.0, 2.0, 4.0, 8.0])
    result = _interpolate_profile(x, y, np.array([1.0, 1.5, 2.5, 4.0]))
    np.testing.assert_allclose(result, [1.0, 1.5, 3.0, 8.0], rtol=1e-6)


def test_lower_and_upper_extrapolation_use_the_nearest_end_segments():
    # IDL: INTERPOL([1,2,4,8], [1,2,3,4], [0, 5]) -> 0.0, 12.0  (asymmetric slopes)
    x = np.array([1.0, 2.0, 3.0, 4.0]); y = np.array([1.0, 2.0, 4.0, 8.0])
    np.testing.assert_allclose(_interpolate_profile(x, y, np.array([0.0, 5.0])), [0.0, 12.0], rtol=1e-6)


def test_two_point_profile_extrapolates_linearly_either_way():
    # IDL: INTERPOL([10,20],[1,2],[0,1,1.5,2,3]) -> 0, 10, 15, 20, 30 for ascending and descending x
    expected = [0.0, 10.0, 15.0, 20.0, 30.0]
    query = np.array([0.0, 1.0, 1.5, 2.0, 3.0])
    np.testing.assert_allclose(_interpolate_profile(np.array([1.0, 2.0]), np.array([10.0, 20.0]), query), expected, rtol=1e-6)
    np.testing.assert_allclose(_interpolate_profile(np.array([2.0, 1.0]), np.array([20.0, 10.0]), query), expected, rtol=1e-6)


def test_result_is_single_precision():
    assert _interpolate_profile(A79_HEIGHT, A79_AMOUNT, QUERY).dtype == np.float32


def test_malformed_profiles_are_rejected():
    with pytest.raises(ValueError):
        _interpolate_profile(np.array([1.0]), np.array([1.0]), np.array([0.5]))
    with pytest.raises(ValueError):
        _interpolate_profile(np.array([1.0, 1.0, 2.0]), np.array([1.0, 2.0, 3.0]), np.array([0.5]))
    with pytest.raises(ValueError):
        _interpolate_profile(np.array([1.0, 2.0]), np.array([1.0, 2.0, 3.0]), np.array([0.5]))


def test_a79_lowest_layer_weights_reproduce_legacy():
    """Legacy layer fractions for the 0.5, 1.5, 2.5 km layers: 46.449 %, 33.333 %, 20.218 %."""

    model = read_microphysics(INPUTS / "microphysics" / "aerosol_a79.mm")
    atmosphere = read_atmosphere(_atmosphere_path(INPUTS, 2), 2)
    layer_height = (atmosphere.height_km[:-1] + atmosphere.height_km[1:]) / np.float32(2.0)
    raw = _interpolate_profile(model.profile_height_km, model.profile_relative_amount, layer_height)
    weights = raw / np.sum(raw, dtype=np.float32)
    assert raw[-1] == pytest.approx(1085.23, abs=1e-2)          # 0.5 km, extrapolated
    np.testing.assert_allclose(weights[-3:], [0.20218, 0.33333, 0.46449], atol=5e-5)
    assert np.all(weights[:-3] == 0.0)


def test_cloud_stg_profile_is_unaffected_by_extrapolation():
    """Both end-node pairs of the STG profile are zero, so extrapolation gives exactly zero."""

    model = read_microphysics(INPUTS / "microphysics" / "liquid-water_stg.mm")
    atmosphere = read_atmosphere(_atmosphere_path(INPUTS, 2), 2)
    layer_height = (atmosphere.height_km[:-1] + atmosphere.height_km[1:]) / np.float32(2.0)
    raw = _interpolate_profile(model.profile_height_km, model.profile_relative_amount, layer_height)
    clamped = np.interp(layer_height.astype(np.float64), model.profile_height_km[::-1].astype(np.float64),
                        model.profile_relative_amount[::-1].astype(np.float64)).astype(np.float32)
    np.testing.assert_array_equal(raw, clamped)
