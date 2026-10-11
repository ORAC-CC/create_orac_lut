from pathlib import Path
from types import SimpleNamespace

import numpy as np

import pytest

from oraclut.master_run import expand_master, print_preflight, read_master, write_manifest, write_manifest_task
from create_orac_luts import _print_execution_configuration


ROOT = Path(__file__).parents[1]
MASTER = ROOT / "makerunfile"


def test_master_lives_at_repository_root_only():
    assert MASTER.is_file()
    assert not (ROOT / "makerunfile.run").exists()
    assert not (ROOT / "runs" / "makerunfile.run").exists()


def _master_text(**changes):
    values = {
        "platform": "noaa-20",
        "instrument": "viirs",
        "forward_model": "cloud",
        "in_path": "create_orac_lut/input_files",
        "instfile": "noaa-20_viirs_v1.inst",
        "mmfiles": ["liquid-water_stg.mm"],
        "lutfile": "auto",
        "atmospheres": 2,
        "channelid": [1, 2, 3],
        "gas": 0,
        "no_rayleigh": 0,
        "reuse_scat": 0,
        "scat_only": 0,
        "srf_quad": 1,
        "nstreams": 60,
        "nmom": 1000,
        "version": 22,
        "out_path": "luts",
        "tmatrix_path": "/network/aopp/matin/eodg/shared/dubovik_tmatrix",
    }
    values.update(changes)
    return "\n".join(f"{key} = {value!r}" for key, value in values.items()) + "\n"


def test_current_master_expands_validation_experiment():
    expanded = expand_master(MASTER)
    assert len(expanded) == 1
    assert expanded[0]["mmfile"] == "liquid-water_old.mm"
    assert expanded[0]["lutfile"] == "liquid-water-cloud_test.lut"
    assert expanded[0]["channelid"] == [1, 10]
    assert expanded[0]["output"].name == "meteosat-10_seviri_m_liquid-water_a01_pold_v22.nc"


def test_master_preflight_reports_resolved_science_and_grid(capsys):
    print_preflight(MASTER)
    output = capsys.readouterr().out
    for expected in (
        "Platform:             meteosat-10",
        "Instrument:           seviri",
        "Forward model:        cloud",
        "Python LUT version:   22",
        "Instrument file:      meteosat-10_seviri_v1.inst",
        "Microphysics:         liquid-water_old.mm",
        "LUT definition:       liquid-water-cloud_test.lut",
        "Channels:             [1, 10]",
        "Atmosphere:           2",
        "Gas absorption:       off",
        "Rayleigh scattering:  on",
        "Scattering only:      no",
        "SRF quadrature:       2 (segmented SRF integration)",
        "DISORT streams:       60",
        "Legendre moments:     1000",
        "Optical Depth         2",
        "Effective Radius      3",
        "Solar Zenith          2",
        "Satellite Zenith      2",
        "Relative Azimuth      2",
        "Output:            luts/meteosat-10_seviri_m_liquid-water_a01_pold_v22.nc",
    ):
        assert expected in output


def test_worker_reporting_uses_actual_srf_point_counts(capsys):
    inst = SimpleNamespace(platform="meteosat-10", instrument="seviri", channelid=np.array([1, 10]))
    mm = SimpleNamespace(substance="liquid-water", shortname="old")
    lut = SimpleNamespace(opd_n=2, efr_n=3, soz_n=2, saz_n=2, raa_n=2)
    qm2_srf = [SimpleNamespace(nwvl=20), SimpleNamespace(nwvl=23)]
    _print_execution_configuration(
        forward_model="cloud", inststr=inst, mmstr=mm, mmfile="liquid-water_old.mm",
        lutfile="liquid-water-cloud_test.lut", lutstr=lut, qm=2, atmospheres=2,
        gas_flag=0, rayleigh_flag=1, scat_only=0, nstreams=60, nmom=1000,
        version=22, srfstrarr=qm2_srf, output="luts/example.nc",
    )
    output = capsys.readouterr().out
    assert "Channel 1: SRF integration points: 20" in output
    assert "Channel 10: SRF integration points: 23" in output
    assert "LUT definition:       liquid-water-cloud_test.lut" in output

    capsys.readouterr()
    _print_execution_configuration(
        forward_model="cloud", inststr=inst, mmstr=mm, mmfile="liquid-water_old.mm",
        lutfile="liquid-water-cloud_test.lut", lutstr=lut, qm=1, atmospheres=2,
        gas_flag=0, rayleigh_flag=1, scat_only=0, nstreams=60, nmom=1000,
        version=22, srfstrarr=[SimpleNamespace(nwvl=1), SimpleNamespace(nwvl=1)],
        output="luts/example.nc",
    )
    assert "Channel 1: SRF integration points: 1" in capsys.readouterr().out


def test_auto_preflight_associates_each_task_with_its_resolved_lut(tmp_path, capsys):
    path = tmp_path / "multi.run"
    path.write_text(_master_text(mmfiles=["liquid-water_stg.mm", "water-ice_sph.mm"]))
    print_preflight(path)
    output = capsys.readouterr().out
    assert "task 0:" in output and "liquid-water_stg.mm" in output
    assert "LUT definition:   liquid-water-cloud.lut" in output
    assert "task 1:" in output and "water-ice_sph.mm" in output
    assert "LUT definition:   ice-cloud.lut" in output


def test_manifest_freezes_deterministic_tasks_and_materialises_one_task(tmp_path):
    manifest = tmp_path / "manifest.json"
    write_manifest(MASTER, manifest)
    payload = manifest.read_text()
    assert '"index": 0' in payload
    concrete = tmp_path / "task.run"
    write_manifest_task(manifest, 0, concrete)
    assert "mmfile = 'liquid-water_old.mm'" in concrete.read_text()


@pytest.mark.parametrize("models", [["liquid-water_stg.mm"], ["liquid-water_stg.mm", "water-ice_sph.mm"]])
def test_one_and_two_composition_expansion(tmp_path, models):
    path = tmp_path / "master.run"
    path.write_text(_master_text(mmfiles=models))
    expanded = expand_master(path)
    assert [item["mmfile"] for item in expanded] == models


def test_five_compositions_preserve_explicit_channels(tmp_path):
    path = tmp_path / "master.run"
    models = [f"aerosol_a{value}.mm" for value in (70, 75, 76, 77, 79)]
    path.write_text(_master_text(forward_model="aerosol", mmfiles=models, channelid=[1, 4, 9], gas=1))
    expanded = expand_master(path)
    assert len(expanded) == 5
    assert all(item["channelid"] == [1, 4, 9] for item in expanded)
    assert all("_a12_" in item["output"].name for item in expanded)


@pytest.mark.parametrize("channels", [[], [1, 1], [999]])
def test_invalid_channel_selection_is_rejected(tmp_path, channels):
    path = tmp_path / "master.run"
    path.write_text(_master_text(channelid=channels))
    with pytest.raises(ValueError):
        expand_master(path)


def test_missing_microphysics_fails_before_expansion(tmp_path):
    path = tmp_path / "master.run"
    path.write_text(_master_text(mmfiles=["missing.mm"]))
    with pytest.raises(FileNotFoundError, match="microphysics"):
        expand_master(path)


def test_duplicate_models_are_rejected_as_duplicate_outputs(tmp_path):
    path = tmp_path / "master.run"
    path.write_text(_master_text(mmfiles=["liquid-water_stg.mm", "liquid-water_stg.mm"]))
    with pytest.raises(ValueError, match="not unique"):
        expand_master(path)


def test_master_rejects_version_21_and_empty_mmfiles(tmp_path):
    version_path = tmp_path / "version.run"
    version_path.write_text(_master_text(version=21))
    with pytest.raises(ValueError, match="version"):
        read_master(version_path)
    empty_path = tmp_path / "empty.run"
    empty_path.write_text(_master_text(mmfiles=[]))
    with pytest.raises(ValueError, match="mmfiles"):
        read_master(empty_path)


def test_single_composition_accepts_explicit_validation_lut(tmp_path):
    path = tmp_path / "validation.run"
    path.write_text(_master_text(lutfile="liquid-water-cloud_test.lut", srf_quad=2))
    expanded = expand_master(path)
    assert len(expanded) == 1
    assert expanded[0]["lutfile"] == "liquid-water-cloud_test.lut"


def test_multiple_compositions_reject_explicit_lut(tmp_path):
    path = tmp_path / "unsafe.run"
    path.write_text(_master_text(
        mmfiles=["liquid-water_stg.mm", "water-ice_sph.mm"],
        lutfile="liquid-water-cloud_test.lut",
    ))
    with pytest.raises(ValueError, match="exactly one mmfile"):
        read_master(path)
