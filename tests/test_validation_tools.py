from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from oraclut.validation.compare import compare_arrays, compare_lut_files
from oraclut.validation.test_case import enforce_small_case, load_existing_test_case


ROOT = Path(__file__).parents[1]


def test_existing_liquid_water_test_grid_is_small_and_exact():
    case = load_existing_test_case(
        ROOT,
        ROOT / "create_orac_lut/input_files/lut/liquid-water-cloud_test.lut",
        (1,),
    )
    enforce_small_case(case)
    assert case.dimensions.optical_depth == 2
    assert case.dimensions.effective_radius == 3
    assert case.dimensions.selected_channels == 1
    assert case.optical_depth == (0.5, 10.0)
    assert case.effective_radius == (1.0, 5.5, 10.0)


def test_array_comparison_reports_scientific_residuals():
    result = compare_arrays("quantity", np.array([1.0, 2.0]), np.array([1.0, 2.5]), units="u")
    assert result.status == "ok"
    assert result.maximum_absolute_difference == 0.5
    assert result.maximum_absolute_index == (1,)
    assert result.units == "u"


def test_array_comparison_reports_shape_mismatch():
    result = compare_arrays("quantity", np.zeros((2,)), np.zeros((1,)))
    assert result.status == "shape mismatch"


def test_frozen_tiny_v2_reference_has_exact_structure_and_small_operator_residual():
    reference = ROOT / "validation/generated/meteosat10_seviri_liquid_water_stg_cloud_test/meteosat10_seviri_liquid_water_stg_cloud_test_legacy_reference_v21.nc"
    candidate = ROOT / "validation/generated/meteosat10_seviri_liquid_water_stg_cloud_test/python_legacy_equivalent_ifort_v21.nc"
    if not candidate.exists():
        pytest.skip("generated Python tiny candidate is unavailable")
    comparison = compare_lut_files(reference, candidate)
    assert comparison["structural"]["status"] == "exact"
    assert comparison["serialization"]["byte_equal"] is False
    residuals = {
        item["quantity"]: item["maximum_absolute_difference"]
        for item in comparison["comparisons"]
        if item["quantity"] in {"R_0d", "R_0v", "R_dd", "R_dv", "T_00", "T_0d", "T_dd", "T_dv"}
    }
    assert max(residuals.values()) < 2.0e-4
    assert sum(
        item.get("nan_inf_difference", {}).get("fill_mask_difference", 0)
        for item in comparison["comparisons"]
    ) == 12


def test_full_r0v_forensic_summary_preserves_sensitive_cluster_separately():
    summary_path = ROOT / "validation/diagnostics/full_r0v_residual_summary.json"
    if not summary_path.exists():
        pytest.skip("forensic full R_0v summary is unavailable")
    import json

    summary = json.loads(summary_path.read_text())
    assert summary["top_50"][0]["index"] == [0, 8, 5, 9, 11, 0]
    assert summary["top_50"][0]["absolute_difference"] > 9.7e-2
    assert summary["suspected_sensitive_cluster"]["count"] == 110
    assert summary["outside_suspected_sensitive_cluster"]["max_abs"] < 1.4e-2
