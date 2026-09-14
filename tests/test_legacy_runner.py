from pathlib import Path

from oraclut.config import read_driver


ROOT = Path(__file__).parents[1]
RUNNER = ROOT / "scripts/generate_liquid_water_cloud_test_reference.sh"
AEROSOL_RUNNER = ROOT / "scripts/generate_aerosol_test_reference.sh"
ENV_SCRIPT = ROOT / "scripts/legacy_idl_env.sh"
DRIVER = ROOT / "validation/reference_cases/meteosat10_seviri_liquid_water_stg_cloud_test/meteosat10_seviri_liquid_water_stg_cloud_test.driver"


def test_runner_resolves_legacy_cwd_and_transactional_output():
    text = RUNNER.read_text()
    assert 'LEGACY_ROOT="${REPO_ROOT}/create_orac_lut"' in text
    assert 'cd "${LEGACY_ROOT}"' in text
    assert 'STAGING_ROOT="${OUTPUT_BASE}/.${CASE_ID}.tmp"' in text
    assert 'mv -- "${STAGING_ROOT}" "${FINAL_ROOT}"' in text
    assert 'Refusing to overwrite completed captured reference' in text
    assert 'rm -rf -- "${STAGING_ROOT}"' in text
    assert 'source_before.json' in text
    assert 'python3 -m oraclut.validation.capture capture' in text


def test_validation_driver_keeps_legacy_relative_input_root_and_atmosphere():
    driver = read_driver(DRIVER)
    assert driver.input_root == "input_files"
    assert driver.instrument_file == "meteosat-10_seviri_v1.inst"
    assert driver.atmosphere == 2
    assert driver.options["channelid"] == "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]"


def test_runner_preflights_all_current_cloud_path_routines():
    text = RUNNER.read_text()
    for routine in (
        "create_orac_cloud_lut_wrapper", "create_orac_cloud_lut",
        "generate_scattering_properties", "create_bwgp", "load_inststr",
        "load_lutstr", "load_srfstrarr", "load_mmdat", "load_atmstr",
        "setup_disort", "call_disort", "write_v2_lut",
    ):
        assert f"resolve_routine, '{routine}'" in text


def test_runner_preflights_relative_srf_and_production_dlms():
    text = RUNNER.read_text()
    assert 'SRF="${LEGACY_ROOT}/input_files/srf/rtcoef_msg_3_seviri_srf_ch01.txt"' in text
    assert 'MIE_DLM="${REPO_ROOT}/mie/dlm-code"' in text
    assert 'DISORT_DLM="${LEGACY_ROOT}/disort2/src"' in text
    assert 'CANONICAL_OUTPUT_DIR="${LEGACY_ROOT}/luts/meteosat-10_seviri_cloud"' in text
    assert 'SOURCE_NC="${CANONICAL_OUTPUT_DIR}/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"' in text


def test_legacy_environment_keeps_recursive_source_and_nonrecursive_dlms():
    text = ENV_SCRIPT.read_text()
    assert 'export IDL_PATH="<IDL_DEFAULT_PATH>:+${oraclut_root}/"' in text
    assert 'export IDL_DLM_PATH="<IDL_DEFAULT_DLM>:${oraclut_root}/mie/dlm-code:${oraclut_root}/create_orac_lut/disort2/src"' in text
    assert "unset IDL_STARTUP" in text


def test_aerosol_reference_runner_is_repository_local_and_uses_production_disort2():
    text = AEROSOL_RUNNER.read_text()
    assert 'LEGACY_ROOT="${REPO_ROOT}/create_orac_lut"' in text
    assert 'cd "${LEGACY_ROOT}"' in text
    assert "CREATE_ORAC_AEROSOL_LUT_DRIVER" in text
    assert "aerosol_a79.mm" in text
    assert "aerosol_test.lut" in text
    assert "create_orac_lut/disort2/src" in text
    assert "DISORT4" not in text
    assert "create_orac_aerosol_lut_wrapper" in text
