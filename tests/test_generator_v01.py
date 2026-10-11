from pathlib import Path

import pytest

from oraclut.config import read_lut_grid, read_microphysics
from oraclut.generate import legacy_lut_filename, main
from oraclut.io.lut import read_lut
from oraclut.pipeline import generate_aerosol, generate_cloud, write_generation


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
INSTRUMENT = INPUTS / "inst" / "meteosat-10_seviri_v1.inst"
CONFIGS = ROOT / "configs"
CLOUD_TEST_CONFIG = "configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver"
AEROSOL_TEST_CONFIG = "configs/meteosat-10_seviri_aerosol_a79_test.driver"


@pytest.fixture
def temporary_configuration():
    """Write a configuration inside the repository; CLI paths stay repo-local."""

    written = []

    def write(text: str, name: str = "test_configuration.driver") -> str:
        path = ROOT / "validation" / "generated" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        written.append(path)
        return str(path)

    yield write
    for path in written:
        path.unlink(missing_ok=True)


def test_template_configuration_parses_and_selects_a_development_grid(capsys):
    assert main(["--config", "configs/template.driver", "--dry-run"]) == 0
    captured = capsys.readouterr().out
    assert "configuration file:  configs/template.driver" in captured
    assert "forward model:       cloud (discrete particle layer)" in captured
    assert "platform/instrument: meteosat-10 / seviri" in captured
    assert "liquid-water-cloud_test.lut  [compact development/test grid]" in captured
    assert "dry-run: no Mie or DISORT calculation was executed" in captured


@pytest.mark.parametrize("config", sorted(path.name for path in CONFIGS.glob("*.driver")))
def test_every_shipped_configuration_resolves(config, capsys):
    assert main(["--config", f"configs/{config}", "--dry-run"]) == 0
    captured = capsys.readouterr().out
    # A configuration always names its grid, and the summary says which kind.
    expected = "[compact development/test grid]" if "_test" in config else "[full grid]"
    if config != "template.driver":
        assert expected in captured
    if config.startswith("validation_"):
        # Validation cases set output= explicitly so they can never collide with
        # the legacy-named default product.
        assert "output:              validation/generated/" in captured
    else:
        assert "output:              create_orac_lut/luts/" in captured


def test_cloud_configuration_reports_dimensions_workload_and_output(capsys):
    assert main(["--config", CLOUD_TEST_CONFIG, "--dry-run"]) == 0
    captured = capsys.readouterr().out
    assert "forward model:       cloud (discrete particle layer)" in captured
    assert "microphysics:        create_orac_lut/input_files/microphysics/liquid-water_stg.mm" in captured
    assert "liquid-water / stg, 1 component(s), scattering code mie" in captured
    assert "LUT grid definition: create_orac_lut/input_files/lut/liquid-water-cloud_test.lut" in captured
    assert "channels:            1 of 11 available: 1" in captured
    assert "atmosphere:          2 (midlatitude summer)" in captured
    assert "Rayleigh scattering: on; gas absorption: off" in captured
    assert "SRF treatment:       srf_quad=1 (monochromatic at the effective channel centre)" in captured
    assert "DISORT streams:      60; phase-function order: 1000" in captured
    for line in (
        "channels             1", "optical_depth        2", "effective_radius     3",
        "solar_zenith         2", "satellite_zenith     2", "relative_azimuth     2",
        "grid points/channel  48", "grid points total    48",
    ):
        assert line in captured
    assert "surface_pressure     - (cloud formulation: no pressure dimension)" in captured
    assert "estimated workload:  18 DISORT states" in captured
    assert "meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc" in captured


def test_aerosol_configuration_reports_the_pressure_dimension(capsys):
    assert main(["--config", AEROSOL_TEST_CONFIG, "--dry-run"]) == 0
    captured = capsys.readouterr().out
    assert "forward model:       aerosol (particles distributed through the atmospheric column)" in captured
    assert "aerosol / a79, 2 component(s), scattering code mie" in captured
    assert "LUT grid definition: create_orac_lut/input_files/lut/aerosol_test.lut" in captured
    assert "Rayleigh scattering: on; gas absorption: on" in captured
    assert "surface_pressure     3" in captured
    assert "grid points/channel  96" in captured
    assert "estimated workload:  36 DISORT states" in captured
    assert "meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc" in captured


def test_command_line_overrides_the_configuration_file(capsys):
    assert main([
        "--config", AEROSOL_TEST_CONFIG,
        "--channels", "2,3", "--no-gas", "--version", "22", "--dry-run",
    ]) == 0
    captured = capsys.readouterr().out
    assert "channels:            2 of 11 available: 2, 3" in captured
    assert "gas absorption: off" in captured
    assert "_a01_pa79_v22.nc" in captured


def test_dry_run_warns_when_the_output_already_exists_but_still_succeeds(capsys):
    output = ROOT / "validation" / "generated" / "test_existing_dry_run.nc"
    output.write_bytes(b"existing product")
    try:
        assert main([
            "--config", CLOUD_TEST_CONFIG, "--output", str(output), "--dry-run",
        ]) == 0
        captured = capsys.readouterr()
        assert "output status:       ALREADY EXISTS (16 bytes)" in captured.out
        assert "already exists" in captured.err
        assert output.read_bytes() == b"existing product"
    finally:
        output.unlink()


def test_dry_run_creates_no_output(capsys):
    output = ROOT / "validation" / "generated" / "test_dry_run_not_created.nc"
    try:
        assert main(["--config", CLOUD_TEST_CONFIG, "--output", str(output), "--dry-run"]) == 0
        captured = capsys.readouterr().out
        assert "output status:       does not exist yet" in captured
        assert not output.exists()
    finally:
        output.unlink(missing_ok=True)


def test_lut_grid_is_never_selected_implicitly(capsys):
    # Without --lut or --test there is no default: a production-sized grid must
    # never be chosen merely because --test is absent.
    with pytest.raises(SystemExit) as raised:
        main([
            "--platform", "meteosat-10", "--instrument", "seviri",
            "--microphysics", "liquid-water_stg.mm", "--channels", "1", "--dry-run",
        ])
    assert raised.value.code == 2
    message = capsys.readouterr().err
    assert "--lut liquid-water-cloud.lut" in message
    assert "--test" in message and "liquid-water-cloud_test.lut" in message


def test_test_flag_conflicts_with_an_explicit_grid(capsys):
    with pytest.raises(SystemExit) as raised:
        main([
            "--platform", "meteosat-10", "--instrument", "seviri",
            "--microphysics", "liquid-water_stg.mm", "--lut", "liquid-water-cloud.lut",
            "--test", "--dry-run",
        ])
    assert raised.value.code == 2
    assert "--test cannot be combined with an explicit LUT grid" in capsys.readouterr().err


def test_explicit_full_grid_is_resolved_without_running_it(capsys):
    assert main([
        "--platform", "meteosat-10", "--instrument", "seviri",
        "--microphysics", "liquid-water_stg.mm", "--lut", "liquid-water-cloud.lut",
        "--dry-run",
    ]) == 0
    captured = capsys.readouterr().out
    assert "liquid-water-cloud.lut  [full grid]" in captured
    assert "grid points/channel  374,000" in captured
    assert "estimated workload:  153,340 DISORT states" in captured


def test_cli_derives_formulation_and_legacy_output_name_like_makerunfile_v2(capsys):
    # makerunfile_v2: material 'aerosol' -> FM 'aerosol', lutfile 'aerosol[_test].lut'
    assert main([
        "--platform", "meteosat-10", "--instrument", "seviri",
        "--microphysics", "aerosol_a79.mm", "--test", "--gas", "--channels", "1", "--dry-run",
    ]) == 0
    captured = capsys.readouterr().out
    assert "forward model:       aerosol" in captured
    assert "instrument file:     create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst" in captured
    assert "LUT grid definition: create_orac_lut/input_files/lut/aerosol_test.lut" in captured
    assert "output:              create_orac_lut/luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc" in captured


CLOUD_HEADER = "input_files\nmeteosat-10_seviri_v1.inst\nliquid-water_stg.mm\nliquid-water-cloud_test.lut\n2\n"


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        # Unknown instrument, misspelled model, missing grid: never substituted.
        ("input_files\nmeteosat-10_sevri_v1.inst\nliquid-water_stg.mm\nliquid-water-cloud_test.lut\n2\n",
         "instrument file not found"),
        ("input_files\nmeteosat-10_seviri_v1.inst\nliquid-water_stgg.mm\nliquid-water-cloud_test.lut\n2\n",
         "microphysics file not found"),
        ("input_files\nmeteosat-10_seviri_v1.inst\nliquid-water_stg.mm\nliquid-water-cloud_tiny.lut\n2\n",
         "LUT file not found"),
        ("input_files\nmeteosat-10_seviri_v1.inst\nliquid-water_stg.mm\ndefault.lut\n2\n",
         "No LUT grid selected"),
        (CLOUD_HEADER + "channelid=[1,99]\n", "are not available on meteosat-10/seviri"),
        (CLOUD_HEADER + "channelid=[1,two]\n", "channels must be integers"),
        (CLOUD_HEADER + "forward_model=clod\n", "must be one of ('cloud', 'aerosol')"),
        (CLOUD_HEADER + "forward_model=aerosol\n", "requires a LUT grid with a sixth surface-pressure block"),
        (CLOUD_HEADER + "gas=1\n",
         "gas absorption is currently implemented for the aerosol formulation only"),
        (CLOUD_HEADER + "srf_quad=3\n", "srf_quad must be 1 (monochromatic) or 2"),
        (CLOUD_HEADER + "srf_quad=mono\n", "expects an integer"),
        (CLOUD_HEADER + "gas=maybe\n", "expects 1/0"),
        (CLOUD_HEADER + "phase_order=50\n", "must be at least twice the DISORT stream count"),
        (CLOUD_HEADER + "nstreams=8\n", "Unrecognised configuration key 'nstreams'"),
        (CLOUD_HEADER + "tmatrix_path=/somewhere\n", "not supported by the Python generator"),
        ("inputs\nmeteosat-10_seviri_v1.inst\nliquid-water_stg.mm\nliquid-water-cloud_test.lut\n2\n",
         "input root must be the repository input hierarchy"),
        ("input_files\nmeteosat-10_seviri_v1.inst\nliquid-water_stg.mm\nliquid-water-cloud_test.lut\n9\n",
         "atmosphere must be a MODTRAN code 0-6"),
    ],
)
def test_invalid_configurations_fail_early_with_a_useful_message(
    text, expected, temporary_configuration, capsys
):
    config = temporary_configuration(text)
    with pytest.raises(SystemExit) as raised:
        main(["--config", config, "--dry-run"])
    assert raised.value.code == 2
    assert expected in capsys.readouterr().err


def test_ignored_legacy_workflow_options_do_not_block_a_run(temporary_configuration, capsys):
    config = temporary_configuration(CLOUD_HEADER + "reuse_scat=1\n")
    assert main(["--config", config, "--dry-run"]) == 0
    captured = capsys.readouterr()
    assert "has no effect in the Python generator" in captured.err
    assert "dry-run: no Mie or DISORT calculation was executed" in captured.out


def test_legacy_lut_filename_matches_archive_convention():
    assert legacy_lut_filename(
        platform="meteosat-10", instrument="seviri", srf_quad=1, substance="liquid-water",
        shortname="stg", atmosphere=2, gas=False, rayleigh=True, revision=21,
    ) == "meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"
    assert legacy_lut_filename(
        platform="meteosat-10", instrument="seviri", srf_quad=1, substance="aerosol",
        shortname="a79", atmosphere=2, gas=True, rayleigh=True, revision=21,
    ) == "meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc"
    assert legacy_lut_filename(
        platform="aqua", instrument="modis", srf_quad=2, substance="liquid-water",
        shortname="stg", atmosphere=2, gas=False, rayleigh=False, revision=21,
    ).endswith("_b_liquid-water_a00_pstg_v21.nc")


def test_production_lut_definitions_follow_load_lutstr_grid_rules():
    aerosol = read_lut_grid(INPUTS / "lut" / "aerosol.lut")
    assert aerosol.surface_pressure.tolist() == [1013.0]
    assert aerosol.optical_depth.size == 20 and aerosol.spacings[0] == "uneven_logarithmic"
    cloud = read_lut_grid(INPUTS / "lut" / "liquid-water-cloud.lut")
    assert cloud.surface_pressure is None
    assert cloud.optical_depth.size == 17 and cloud.effective_radius.tolist()[:2] == [1.0, 3.0]


def test_cli_refuses_existing_output_before_scientific_execution(tmp_path):
    output = tmp_path / "existing.nc"
    # CLI output paths are intentionally repository-local.
    output = ROOT / "validation" / "generated" / "test_existing_output.nc"
    output.write_bytes(b"do not replace")
    try:
        with pytest.raises(SystemExit) as raised:
            main([
                "--platform", "meteosat-10", "--instrument", "seviri",
                "--microphysics", "liquid-water_stg.mm", "--test",
                "--channels", "1", "--output", str(output),
            ])
        assert raised.value.code == 2
        assert output.read_bytes() == b"do not replace"
    finally:
        output.unlink()


def test_aerosol_inputs_preserve_components_and_pressure_grid():
    model = read_microphysics(INPUTS / "microphysics" / "aerosol_a79.mm")
    grid = read_lut_grid(INPUTS / "lut" / "aerosol_test.lut")
    assert [component.component for component in model.components] == ["waf", "saf"]
    assert [component.mixing_ratio for component in model.components] == [0.625, 0.375]
    assert grid.surface_pressure.tolist() == [950.0, 1013.0, 1050.0]


@pytest.mark.slow
def test_cloud_pipeline_and_v2_writer(tmp_path):
    result = generate_cloud(
        ROOT,
        lut_path=INPUTS / "lut" / "liquid-water-cloud_test.lut",
        microphysics_path=INPUTS / "microphysics" / "liquid-water_stg.mm",
        instrument_path=INSTRUMENT,
        channels=(1,),
        phase_order=120,
    )
    output = tmp_path / "cloud.nc"
    write_generation(output, result)
    lut = read_lut(output)
    assert lut.variable_dimensions["T_dv"] == ("satellite_zenith", "optical_depth", "effective_radius", "channels")
    assert lut.variables["T_dv"].shape == (2, 2, 3, 1)


@pytest.mark.slow
def test_aerosol_pipeline_and_pressure_aware_v2_writer(tmp_path):
    result = generate_aerosol(
        ROOT,
        lut_path=INPUTS / "lut" / "aerosol_test.lut",
        microphysics_path=INPUTS / "microphysics" / "aerosol_a79.mm",
        instrument_path=INSTRUMENT,
        channels=(1,),
        phase_order=120,
    )
    output = tmp_path / "aerosol.nc"
    write_generation(output, result)
    lut = read_lut(output)
    assert lut.variable_dimensions["T_dv"] == (
        "satellite_zenith", "optical_depth", "effective_radius", "surface_pressure", "channels"
    )
    assert lut.variables["R_0v"].shape == (2, 2, 2, 2, 2, 3, 1)
    assert lut.variables["surface_pressure"].tolist() == [950.0, 1013.0, 1050.0]
