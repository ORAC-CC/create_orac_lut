"""Inexpensive tests for the legacy-versus-Python comparison infrastructure.

These exercise the comparison machinery only. They run no scientific
calculation, no IDL, no Mie and no DISORT.
"""

from pathlib import Path

import numpy as np
import pytest

from oraclut.io.lut import read_lut
from oraclut.validation.compare import ArrayComparison, compare_arrays
from oraclut.validation.product_comparison import (
    NUMERICAL_BANDS,
    aerosol_checks,
    classify_numerical,
    collect_warnings,
    compare_products,
    render_markdown,
    structural_comparison,
    write_report,
)


ROOT = Path(__file__).parents[1]
PYTHON_AEROSOL = (
    ROOT / "create_orac_lut" / "luts" / "meteosat-10_seviri_aerosol"
    / "meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc"
)
LEGACY_CLOUD = (
    ROOT / "validation" / "generated" / "meteosat10_seviri_liquid_water_stg_cloud_test"
    / "meteosat10_seviri_liquid_water_stg_cloud_test_legacy_reference_v21.nc"
)
PYTHON_CLOUD = (
    ROOT / "create_orac_lut" / "luts" / "meteosat-10_seviri_cloud"
    / "meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"
)
LEGACY_RUNNER = ROOT / "scripts" / "generate_aerosol_legacy_reference.sh"
PROVENANCE_DIAGNOSTIC = ROOT / "validation" / "diagnostics" / "aerosol_legacy_provenance.pro"


def _operator(quantity: str, maximum: float | None, *, compared: int = 96,
              differing: int = 96, status: str = "ok", nan_inf=None) -> ArrayComparison:
    return ArrayComparison(
        quantity, None, (96,), (96,), "float32", "float32", 0.0, 0.0, 1.0, 1.0,
        maximum, maximum, maximum, (0,), (0,), nan_inf or {}, status,
        mean_absolute_difference=maximum, compared_elements=compared,
        differing_elements=differing,
    )


def test_compare_arrays_reports_mean_and_differing_element_counts():
    result = compare_arrays(
        "quantity", np.array([1.0, 2.0, 3.0, 4.0]), np.array([1.0, 2.0, 3.5, 4.5])
    )
    assert result.status == "ok"
    assert result.maximum_absolute_difference == pytest.approx(0.5)
    assert result.mean_absolute_difference == pytest.approx(0.25)
    assert (result.compared_elements, result.differing_elements) == (4, 2)


def test_classification_is_green_for_floating_point_scale_residuals():
    verdict = classify_numerical([_operator("R_0v", 1.2e-4), _operator("T_dv", 7.0e-5)])
    assert verdict["result"] == "GREEN"
    assert verdict["worst_operator_absolute_difference"] == pytest.approx(1.2e-4)


def test_classification_is_amber_for_modest_discrepancies():
    verdict = classify_numerical([_operator("R_0v", 3.0e-3, differing=4)])
    assert verdict["result"] == "AMBER"


def test_classification_is_red_for_large_or_systematic_discrepancies():
    assert classify_numerical([_operator("R_0v", 0.4, differing=4)])["result"] == "RED"
    systematic = classify_numerical([_operator("R_dd", 5.0e-3, compared=96, differing=96)])
    assert systematic["result"] == "RED"
    assert any("of elements differ" in reason for reason in systematic["reasons"])


def test_classification_is_red_for_shape_mismatch_or_nonfinite_candidate():
    mismatch = classify_numerical([_operator("R_0v", None, status="shape mismatch")])
    assert mismatch["result"] == "RED"
    nonfinite = classify_numerical(
        [_operator("R_0v", 1.0e-6, nan_inf={"candidate_nan": 3})]
    )
    assert nonfinite["result"] == "RED"
    assert any("nonfinite" in reason for reason in nonfinite["reasons"])


def test_classification_is_red_without_operator_variables():
    assert classify_numerical([_operator("central_wavelength", 0.0)])["result"] == "RED"


def test_bands_are_documented_as_case_specific():
    assert "not a universal scientific tolerance" in NUMERICAL_BANDS["rationale"]


def test_warning_collection_records_the_disort_near_singular_message(tmp_path):
    log = tmp_path / "stdout.log"
    log.write_text(
        "starting\n"
        " UPBEAM--SGECO says matrix near singular\n"
        " UPBEAM--SGECO says matrix near singular\n"
        "done\n"
    )
    warnings = collect_warnings(log, None)
    assert warnings["upbeam_sgeco_near_singular_occurrences"] == 2
    assert len(warnings["legacy"]["lines"]) == 2
    assert warnings["python"]["log"] is None
    assert "not treated as failure" in warnings["interpretation"]


def test_missing_log_is_reported_rather_than_raising(tmp_path):
    warnings = collect_warnings(tmp_path / "absent.log", None)
    assert warnings["legacy"]["note"] == "log not found"


@pytest.mark.skipif(not PYTHON_AEROSOL.is_file(), reason="Python aerosol product not generated")
def test_aerosol_checks_confirm_pressure_grid_and_state_count():
    product = read_lut(PYTHON_AEROSOL, validate_v21=False)
    checks = aerosol_checks(product, product, expected_channels=[1])
    assert checks["result"] == "PASS"
    details = {check["check"]: check for check in checks["checks"]}
    assert details["reference: surface_pressure values"]["ok"] is True
    assert "[950.0, 1013.0, 1050.0]" in details["reference: surface_pressure values"]["detail"]
    assert details["reference: total aerosol grid states"]["ok"] is True
    assert "96" in details["reference: total aerosol grid states"]["detail"]
    position = details["reference: pressure dimension position in R_0v"]
    assert position["ok"] is True and "surface_pressure at index 5" in position["detail"]
    assert details["reference: channel selection"]["ok"] is True


@pytest.mark.skipif(not PYTHON_AEROSOL.is_file(), reason="Python aerosol product not generated")
def test_aerosol_checks_fail_when_the_pressure_grid_is_wrong():
    product = read_lut(PYTHON_AEROSOL, validate_v21=False)
    checks = aerosol_checks(product, product, expected_surface_pressure=(900.0, 1000.0, 1100.0))
    assert checks["result"] == "FAIL"
    assert checks["failures"]


@pytest.mark.skipif(not PYTHON_AEROSOL.is_file(), reason="Python aerosol product not generated")
def test_self_comparison_is_structurally_identical_and_green(tmp_path):
    report = compare_products(PYTHON_AEROSOL, PYTHON_AEROSOL, formulation="aerosol")
    assert report["structural"]["result"] == "PASS"
    assert report["structural"]["mismatches"] == []
    assert report["numerical"]["result"] == "GREEN"
    assert report["reference"]["sha256"] == report["candidate"]["sha256"]
    assert set(report["variable_groups"]["radiative_transfer_operators"]) == {
        "R_0v", "R_0d", "R_dv", "R_dd", "T_00", "T_0d", "T_dv", "T_dd",
    }
    json_path, markdown_path = write_report(report, tmp_path / "out")
    assert json_path.is_file() and markdown_path.is_file()
    text = markdown_path.read_text()
    for heading in ("## Structure", "## Numerical residuals", "## Aerosol-specific checks",
                    "## Warnings", "## Classification"):
        assert heading in text
    assert "surface_pressure" in text


@pytest.mark.skipif(
    not (LEGACY_CLOUD.is_file() and PYTHON_CLOUD.is_file()),
    reason="captured legacy cloud reference or Python cloud product missing",
)
def test_validated_cloud_pair_reproduces_its_published_result():
    """The cloud case is the calibration context for the aerosol comparison."""

    report = compare_products(LEGACY_CLOUD, PYTHON_CLOUD, formulation="cloud")
    assert report["structural"]["result"] == "PASS"
    assert report["numerical"]["result"] == "GREEN"
    assert report["numerical"]["worst_operator_absolute_difference"] < 2.0e-4
    # The legacy writer emits fill values for singleton-channel optics; the
    # Python product writes finite values. Reported, never silently equated.
    fill = report["structural"]["fill_and_finite"]
    assert fill["single_scatter_albedo"]["fill_pattern_differs"] is True
    assert fill["single_scatter_albedo"]["reference"]["fill"] == 3
    assert fill["single_scatter_albedo"]["candidate"]["finite"] == 3


@pytest.mark.skipif(not PYTHON_AEROSOL.is_file(), reason="Python aerosol product not generated")
def test_structural_comparison_reports_a_dimension_mismatch():
    product = read_lut(PYTHON_AEROSOL, validate_v21=False)

    class Altered:
        path = product.path
        file_format = product.file_format
        dimensions = dict(product.dimensions) | {"surface_pressure": 4}
        coordinates = product.coordinates
        variables = product.variables
        variable_dimensions = product.variable_dimensions
        variable_dtypes = product.variable_dtypes
        variable_attributes = product.variable_attributes
        global_attributes = product.global_attributes
        variable_names = product.variable_names

    result = structural_comparison(product, Altered())
    assert result["result"] == "FAIL"
    assert any("surface_pressure" in item for item in result["mismatches"])


def test_legacy_runner_and_diagnostic_are_present_and_executable():
    assert LEGACY_RUNNER.is_file()
    assert LEGACY_RUNNER.stat().st_mode & 0o111, "the legacy runner must be executable"
    text = (ROOT / "scripts" / "generate_legacy_reference.sh").read_text()
    # The runner must never silently accept a mixed-generation source tree.
    assert "Routine provenance verified" in text
    assert "not from the selected legacy source tree" in text
    assert "resolved from the active V2.1 tree" in text
    assert PROVENANCE_DIAGNOSTIC.is_file()
    assert "PROVENANCE " in PROVENANCE_DIAGNOSTIC.read_text()


def test_legacy_runner_points_at_the_extracted_working_tree():
    text = (ROOT / "scripts" / "generate_legacy_reference.sh").read_text()
    assert (
        "LEGACY_SRC_DEFAULT=/network/scratch/grainger/oraclut-reconciliation/"
        "legacy-reference-source/current_local_create_orac_lut" in text
    )
    # The stale preservation build copy is incomplete and must not be used.
    assert "preservation/git_preservation" not in text
    # The source root is defined once; everything else derives from it.
    assert text.count("/network/scratch") == 1


def test_legacy_runner_preflights_every_routine_the_provenance_gate_needs():
    text = (ROOT / "scripts" / "generate_legacy_reference.sh").read_text()
    diagnostic = PROVENANCE_DIAGNOSTIC.read_text()
    required = text.split("COHERENT_REQUIRED=(", 1)[1].split(")", 1)[0].split()
    listed = text.split("LEGACY_SOURCE_FILES=(", 1)[1].split(")", 1)[0].split()
    assert "makerunfile_v2.pro" in listed
    assert "create_orac_aerosol_lut.pro" in listed
    assert "create_orac_aerosol_lut_wrapper.pro" in listed
    # Every gated routine must be one the diagnostic actually reports, and must
    # be backed either by its own source file or by a file that defines it.
    embedded = {"lut_quadrature": "load_lutstr.pro",
                "read_srfstr": "input_files/read_srfstr.pro",
                "bbconstants": "input_files/bbconstants.pro"}
    for routine in required:
        assert f"'{routine}'" in diagnostic, f"{routine} is not reported by the diagnostic"
        expected = embedded.get(routine, f"{routine}.pro")
        assert expected in listed, f"{routine} has no preflighted source file"


def test_render_markdown_survives_a_report_without_optional_sections():
    report = {
        "generated_utc": "2026-09-13T00:00:00+00:00",
        "reference": {"path": "a.nc", "sha256": "0" * 64, "size_bytes": 1, "role": "legacy"},
        "candidate": {"path": "b.nc", "sha256": "1" * 64, "size_bytes": 2, "role": "python"},
        "formulation": "aerosol",
        "structural": {"result": "PASS", "mismatches": [], "dimensions": {"reference": {}, "candidate": {}},
                       "coordinates": {}, "fill_and_finite": {}},
        "numerical": {"result": "GREEN", "reasons": [], "bands": NUMERICAL_BANDS},
        "comparisons": [],
        "variable_groups": {},
        "warnings": {"upbeam_sgeco_near_singular_occurrences": 0, "interpretation": "n/a"},
    }
    text = render_markdown(report)
    assert "**Structural: PASS**" in text and "**Numerical: GREEN**" in text


# ---------------------------------------------------------------------------
# Visible + infrared validation matrix
# ---------------------------------------------------------------------------

from oraclut.config import read_instrument  # noqa: E402
from oraclut.generate import main, read_configuration_file  # noqa: E402
from oraclut.validation.product_comparison import (  # noqa: E402
    channel_diagnostics, per_channel_comparison, systematic_relative_fraction,
)

GENERAL_RUNNER = ROOT / "scripts" / "generate_legacy_reference.sh"
IR_CHANNEL = 9
TWO_CHANNEL_CLOUD = (
    ROOT / "validation" / "generated" / "cloud_visible_ir_comparison" / "python"
    / "meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"
)
TWO_CHANNEL_AEROSOL = (
    ROOT / "validation" / "generated" / "aerosol_visible_ir_comparison" / "python"
    / "meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc"
)


def test_seviri_channel_1_is_solar_and_channel_9_is_thermal_only():
    instrument = read_instrument(ROOT / "create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst")
    assert 1 in instrument.solar_channels and 1 not in instrument.thermal_channels
    assert IR_CHANNEL in instrument.thermal_channels and IR_CHANNEL not in instrument.solar_channels
    # Channel 4 (3.9 um) is the mixed solar+thermal channel; it is deliberately
    # not the IR test because it does not isolate the thermal route.
    assert 4 in instrument.solar_channels and 4 in instrument.thermal_channels
    assert instrument.srf_files[IR_CHANNEL] == "rtcoef_msg_3_seviri_srf_ch09.txt"


def test_ir_channel_inputs_exist_for_gas_and_srf():
    inputs = ROOT / "create_orac_lut" / "input_files"
    assert (inputs / "srf" / "rtcoef_msg_3_seviri_srf_ch09.txt").is_file()
    for atmosphere in range(7):
        assert (inputs / "gas" / f"ModtranGasOpd_A{atmosphere}_meteosat-10_seviri_ch09.gas").is_file()
    text = (inputs / "bbconstants.pro").read_text()
    assert "'rtcoef_msg_3_seviri_srf_ch09.txt'" in text


@pytest.mark.parametrize(
    ("config", "channels", "formulation", "gas"),
    [
        ("validation_meteosat-10_seviri_cloud_liquid-water_stg_test_ch09.driver", (9,), "cloud", False),
        ("validation_meteosat-10_seviri_cloud_liquid-water_stg_test_ch01_ch09.driver", (1, 9), "cloud", False),
        ("validation_meteosat-10_seviri_aerosol_a79_test_ch09.driver", (9,), "aerosol", True),
        ("validation_meteosat-10_seviri_aerosol_a79_test_ch01_ch09.driver", (1, 9), "aerosol", True),
    ],
)
def test_validation_configurations_select_the_intended_channels(config, channels, formulation, gas):
    settings = read_configuration_file(ROOT / "configs" / config)
    assert tuple(settings["channels"]) == channels
    assert settings["forward_model"] == formulation
    assert settings.get("gas", False) is gas
    assert settings["srf_quad"] == 1 and settings["atmosphere"] == 2
    # Validation products must never land on the default legacy path.
    assert str(settings["output"]).startswith("validation/generated/")


def test_validation_configurations_dry_run_and_report_solar_thermal_split(capsys):
    assert main(["--config", "configs/validation_meteosat-10_seviri_cloud_liquid-water_stg_test_ch01_ch09.driver", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "channels:            2 of 11 available: 1, 9" in out
    assert "solar [1]; thermal [9]" in out
    assert main(["--config", "configs/validation_meteosat-10_seviri_aerosol_a79_test_ch09.driver", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "solar none; thermal [9]" in out and "gas absorption: on" in out


def test_gas_is_rejected_for_the_cloud_ir_case(capsys):
    config = ROOT / "validation" / "generated" / "test_cloud_ir_gas.driver"
    config.write_text(
        "input_files\nmeteosat-10_seviri_v1.inst\nliquid-water_stg.mm\nliquid-water-cloud_test.lut\n2\n"
        "channelid=[9]\ngas=1\n"
    )
    try:
        with pytest.raises(SystemExit):
            main(["--config", str(config), "--dry-run"])
        assert "aerosol formulation only" in capsys.readouterr().err
    finally:
        config.unlink(missing_ok=True)


def test_systematic_relative_fraction_ignores_near_zero_references():
    reference = np.array([0.5, 0.5, 0.5, 0.5, 1e-4])
    candidate = np.array([0.505, 0.505, 0.5, 0.5, 2e-4])   # 1% bias on two of four usable elements
    assert systematic_relative_fraction(reference, candidate) == pytest.approx(0.5)
    assert systematic_relative_fraction(np.array([1e-3, 1e-3]), np.array([2e-3, 3e-3])) is None


def test_relative_bias_holds_classification_at_amber():
    verdict = classify_numerical([_operator("R_dd", 4.0e-4)], systematic_relative={"R_dd": 0.75})
    assert verdict["result"] == "AMBER"
    assert any("systematic bias" in reason for reason in verdict["reasons"])
    verdict = classify_numerical([_operator("R_dd", 4.0e-4)], systematic_relative={"R_dd": 0.25})
    assert verdict["result"] == "GREEN"


def test_underflow_reports_are_counted_separately_and_not_failures(tmp_path):
    err = tmp_path / "stderr.log"
    err.write_text("% Program caused arithmetic error: Floating underflow\n")
    warnings = collect_warnings(None, None, legacy_stderr=err)
    assert warnings["floating_underflow_occurrences"] == 1
    assert warnings["upbeam_sgeco_near_singular_occurrences"] == 0
    assert "not treated as failure" in warnings["interpretation"]


@pytest.mark.skipif(not TWO_CHANNEL_CLOUD.is_file(), reason="two-channel cloud product not generated")
def test_two_channel_cloud_product_is_split_per_channel_with_thermal_diagnostics():
    product = read_lut(TWO_CHANNEL_CLOUD, validate_v21=False)
    rows = per_channel_comparison(product, product)
    by_quantity = {}
    for row in rows:
        by_quantity.setdefault(row["quantity"], []).append(row["channel_id"])
    assert by_quantity["T_dv"] == [1, 9]          # channels dimension
    assert by_quantity["R_0v"] == [1]             # solar_channels dimension
    assert by_quantity["E_md"] == [9]             # thermal_channels dimension
    assert all(row["maximum_absolute_difference"] == 0.0 for row in rows if row["status"] == "ok")
    diagnostics = channel_diagnostics(product, product)
    kinds = {c["channel_id"]: c["kind"] for c in diagnostics["reference"]["channels"]}
    assert kinds == {1: "solar", 9: "thermal"}
    assert diagnostics["reference"]["has_E_md"] is True
    assert diagnostics["metadata_mismatches"] == []
    assert diagnostics["metadata"]["B1"]["reference"] == pytest.approx([9577.68])
    assert diagnostics["reference"]["product_name_codes"]["atmospheric_model_code"] == "01"
    assert diagnostics["reference"]["product_name_codes"]["gas_absorption_in_product_name"] is False


@pytest.mark.skipif(not TWO_CHANNEL_AEROSOL.is_file(), reason="two-channel aerosol product not generated")
def test_two_channel_aerosol_product_keeps_pressure_and_gas_code():
    product = read_lut(TWO_CHANNEL_AEROSOL, validate_v21=False)
    report = compare_products(TWO_CHANNEL_AEROSOL, TWO_CHANNEL_AEROSOL, formulation="aerosol",
                              expected_channels=[1, 9])
    assert report["structural"]["result"] == "PASS" and report["numerical"]["result"] == "GREEN"
    assert report["aerosol_checks"]["result"] == "PASS"
    assert "surface_pressure" in product.variable_dimensions["E_md"]
    codes = report["channel_diagnostics"]["reference"]["product_name_codes"]
    assert codes["atmospheric_model_code"] == "12" and codes["gas_absorption_in_product_name"] is True
    text = render_markdown(report)
    assert "## Per-channel residuals" in text and "## Channel and thermal-path diagnostics" in text


def test_general_legacy_runner_covers_both_formulations_and_is_delegated_to():
    text = GENERAL_RUNNER.read_text()
    assert GENERAL_RUNNER.stat().st_mode & 0o111
    assert "create_orac_cloud_lut_wrapper" in text and "create_orac_aerosol_lut_wrapper" in text
    assert "LEGACY_SRC_DEFAULT=/network/scratch/grainger/oraclut-reconciliation/legacy-reference-source/current_local_create_orac_lut" in text
    assert text.count("/network/scratch") == 1
    assert "Routine provenance verified" in text and "resolved from the active V2.1 tree" in text
    wrapper = LEGACY_RUNNER.read_text()
    assert "generate_legacy_reference.sh" in wrapper and "--formulation aerosol" in wrapper
    diagnostic = PROVENANCE_DIAGNOSTIC.read_text()
    assert "'create_orac_cloud_lut'" in diagnostic and "'create_orac_cloud_lut_wrapper'" in diagnostic


# ---------------------------------------------------------------------------
# Regression: aerosol thermal emissivity convention
# ---------------------------------------------------------------------------

LEGACY_AEROSOL_IR = (
    ROOT / "validation" / "generated" / "aerosol_ir_comparison" / "legacy"
    / "meteosat-10_seviri_m_aerosol_a12_pa79_v21_legacy.nc"
)


@pytest.mark.slow
@pytest.mark.skipif(not LEGACY_AEROSOL_IR.is_file(), reason="captured legacy aerosol channel-9 reference missing")
def test_aerosol_thermal_emissivity_uses_the_legacy_fraction_convention(tmp_path):
    """E_md must be a fraction in 0-1 like the legacy V2 product, not a percentage.

    The visible+IR validation matrix found the aerosol emission block scaling
    E_md by 100 (the cloud block did not); this test compares the aerosol
    channel-9 emissivity directly with the captured legacy reference and fails
    for any factor-100 convention slip in either direction.
    """

    from oraclut.pipeline import generate_aerosol, write_generation

    inputs = ROOT / "create_orac_lut" / "input_files"
    result = generate_aerosol(
        ROOT,
        lut_path=inputs / "lut" / "aerosol_test.lut",
        microphysics_path=inputs / "microphysics" / "aerosol_a79.mm",
        instrument_path=inputs / "inst" / "meteosat-10_seviri_v1.inst",
        channels=(9,), atmosphere_code=2, srf_quad=1, rayleigh=True, gas=True,
    )
    output = tmp_path / "aerosol_ch09.nc"
    write_generation(output, result)
    python = read_lut(output, validate_v21=False)
    legacy = read_lut(LEGACY_AEROSOL_IR, validate_v21=False)
    em_python = np.asarray(python.variables["E_md"], dtype=np.float64)
    em_legacy = np.asarray(legacy.variables["E_md"], dtype=np.float64)
    assert em_python.shape == em_legacy.shape
    assert python.variable_dimensions["E_md"] == legacy.variable_dimensions["E_md"]
    # Same physical convention: a fraction, with legacy-level agreement, not 100x.
    assert 0.0 < em_legacy.max() <= 1.0 and 0.0 < em_python.max() <= 1.0
    assert np.max(np.abs(em_python - em_legacy)) < 1.0e-3
    assert np.median(em_python / em_legacy) == pytest.approx(1.0, abs=1.0e-3)


def test_aerosol_and_cloud_emission_blocks_share_one_normalisation():
    """Static guard: neither emission block may carry a stray percentage factor."""

    import inspect
    from oraclut import pipeline

    source = inspect.getsource(pipeline)
    emission_lines = [line for line in source.splitlines() if 'arrays["em"]' in line and "=" in line]
    assert len(emission_lines) == 2, emission_lines
    for line in emission_lines:
        assert line.strip().endswith('/ np.float32(bbe)'), line
        assert "100" not in line, line
