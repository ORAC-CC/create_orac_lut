"""Code-version identification, LUT-version consistency and provenance records.

The authoritative definition is the CODE_VERSION file; Git supplies the
commit and working-tree state.  These tests check the parsing and validation
of that file, the banner, the detection of a modified tree and of missing Git
metadata (in temporary repositories, never in this checkout), the rejection
of a LUT version the source release does not produce by both generator entry
points before any calculation, the provenance sidecar and its verification.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

import create_orac_luts
from oraclut import provenance, version
from oraclut.generate import main as generate_main

ROOT = Path(__file__).resolve().parents[1]
GIT = version._git_executable()


def git(root: Path, *args: str) -> str:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update({"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@x", "GIT_COMMITTER_NAME": "t",
                "GIT_COMMITTER_EMAIL": "t@x", "LC_ALL": "C"})
    return subprocess.run([GIT, "-C", str(root), *args], check=True, capture_output=True, text=True, env=env).stdout


def write_code_version(root: Path, **overrides) -> Path:
    fields = {"code_release": "lut-code-v25", "lut_version": "25", "compatible_lut_versions": "24",
              "scientific_config": "v25", "scientific_config_description": "test description"}
    fields.update(overrides)
    path = root / "CODE_VERSION"
    path.write_text("# test\n" + "".join(f"{k} = {v}\n" for k, v in fields.items()))
    return path


@pytest.fixture
def repo(tmp_path):
    """A small Git repository with a CODE_VERSION file, tagged lut-code-v25."""

    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-q")
    write_code_version(root)
    (root / "src").mkdir()
    (root / "src" / "module.py").write_text("x = 1\n")
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "release")
    git(root, "tag", "-a", "lut-code-v25", "-m", "lut-code-v25: test release")
    return root


# ---------------------------------------------------------------------------
# CODE_VERSION
# ---------------------------------------------------------------------------

def test_repository_code_version_is_valid_and_matches_a_release_tag_name():
    code = version.read_code_version(ROOT / "CODE_VERSION")
    assert version.RELEASE_TAG_PATTERN.match(code.code_release)
    assert code.lut_version >= 25
    assert code.lut_version not in code.compatible_lut_versions
    assert code.scientific_config and code.scientific_config_description


@pytest.mark.parametrize("overrides, message", [
    ({"code_release": "v25"}, "annotated tag name"),
    ({"lut_version": "twenty-five"}, "must be integers"),
    ({"compatible_lut_versions": "25"}, "must not repeat"),
    ({"scientific_config": ""}, "must not be empty"),
])
def test_malformed_code_version_is_rejected(tmp_path, overrides, message):
    path = write_code_version(tmp_path, **overrides)
    with pytest.raises(version.VersionError, match=message):
        version.read_code_version(path)


def test_missing_or_unknown_settings_are_rejected(tmp_path):
    path = tmp_path / "CODE_VERSION"
    path.write_text("code_release = lut-code-v1\n")
    with pytest.raises(version.VersionError, match="missing settings"):
        version.read_code_version(path)
    write_code_version(tmp_path)
    path.write_text(path.read_text() + "extra = 1\n")
    with pytest.raises(version.VersionError, match="unknown settings"):
        version.read_code_version(path)
    with pytest.raises(version.VersionError, match="not found"):
        version.read_code_version(tmp_path / "absent")


# ---------------------------------------------------------------------------
# Git state and banner
# ---------------------------------------------------------------------------

def test_clean_tagged_release_is_identified(repo):
    summary = version.summarise(repo)
    assert summary.git.available and summary.git.working_tree == "CLEAN"
    assert summary.git.commit == git(repo, "rev-parse", "HEAD").strip()
    assert summary.git.release_relation("lut-code-v25") == "HEAD is the tagged release lut-code-v25"
    banner = summary.banner(lut_version=25)
    for line in ("ORAC LUT Generator", "LUT version:       25 (requested 25)", "Code release:      lut-code-v25",
                 f"Git commit:        {summary.git.commit[:7]}", "Working tree:      CLEAN",
                 "Scientific config: v25"):
        assert line in banner, banner


def test_modified_tree_is_never_reported_clean(repo):
    (repo / "src" / "module.py").write_text("x = 2\n")
    state = version.git_state(repo, "lut-code-v25")
    assert state.working_tree == "MODIFIED" and state.tracked_changes == ("src/module.py",)
    banner = version.summarise(repo).banner()
    assert "Working tree:      MODIFIED" in banner and "modified:      src/module.py" in banner
    git(repo, "checkout", "--", "src/module.py")
    (repo / "src" / "new_production_file.py").write_text("y = 1\n")
    state = version.git_state(repo, "lut-code-v25")
    assert state.working_tree == "MODIFIED" and state.untracked_production == ("src/new_production_file.py",)
    (repo / "src" / "new_production_file.py").unlink()
    (repo / "notes.txt").write_text("not production content\n")
    assert version.git_state(repo, "lut-code-v25").working_tree == "CLEAN"


def test_commits_after_the_release_tag_are_reported(repo):
    (repo / "src" / "module.py").write_text("x = 3\n")
    git(repo, "commit", "-q", "-am", "later change")
    state = version.git_state(repo, "lut-code-v25")
    assert state.working_tree == "CLEAN"
    assert state.release_relation("lut-code-v25") == "HEAD is 1 commit(s) after tag lut-code-v25"
    assert state.describe.startswith("lut-code-v25-1-g")
    assert "not present" in version.git_state(repo, "lut-code-v26").release_relation()


def test_missing_git_metadata_is_reported_not_invented(tmp_path):
    write_code_version(tmp_path)
    state = version.git_state(tmp_path, "lut-code-v25")
    assert not state.available and state.commit is None
    assert state.working_tree.startswith("UNKNOWN (Git metadata unavailable")
    banner = version.summarise(tmp_path).banner()
    assert "Git commit:        unavailable" in banner and "UNKNOWN" in banner


def test_release_listing_reads_annotated_tags(repo):
    releases = version.list_releases(repo)
    assert [r["tag"] for r in releases] == ["lut-code-v25"]
    assert releases[0]["annotated"] and releases[0]["commit"] == git(repo, "rev-parse", "HEAD").strip()
    assert releases[0]["subject"] == "lut-code-v25: test release"


def test_version_command_line(repo, capsys):
    assert version.main(["--root", str(repo)]) == 0
    assert "Code release:      lut-code-v25" in capsys.readouterr().out
    assert version.main(["--root", str(repo), "--json"]) == 0
    record = json.loads(capsys.readouterr().out)
    assert record["git"]["working_tree"] == "CLEAN" and record["software"]["python"]
    assert version.main(["--root", str(repo), "--list-releases"]) == 0
    assert "lut-code-v25" in capsys.readouterr().out
    assert version.main(["--root", str(repo), "--check-lut-version", "24"]) == 0
    assert version.main(["--root", str(repo), "--check-lut-version", "22"]) == 1
    assert "requests LUT version 22" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# LUT-version consistency
# ---------------------------------------------------------------------------

def test_check_lut_version_accepts_declared_and_compatible_versions_only():
    code = version.read_code_version(ROOT / "CODE_VERSION")
    assert version.check_lut_version(code.lut_version, code) == code.lut_version
    for compatible in code.compatible_lut_versions:
        assert version.check_lut_version(compatible, code) == compatible
    with pytest.raises(version.LutVersionMismatch, match="requests LUT version 22"):
        version.check_lut_version(22, code, source="run file x.run")
    with pytest.raises(version.LutVersionMismatch, match="does not state the LUT version"):
        version.check_lut_version(None, code)


def test_run_file_with_a_foreign_lut_version_is_refused_before_any_calculation(tmp_path, capsys, monkeypatch):
    text = (ROOT / "runs" / "template.run").read_text()
    code = version.read_code_version(ROOT / "CODE_VERSION")
    foreign = 22 if not code.accepts(22) else 1
    text = text.replace(f"version       = {code.lut_version}", f"version       = {foreign}")
    text = text.replace("out_path      = 'luts'", f"out_path      = '{tmp_path / 'out'}'")
    run = tmp_path / "foreign.run"
    run.write_text(text)
    called = []
    monkeypatch.setattr(create_orac_luts, "create_orac_cloud_lut", lambda *a, **k: called.append(1) or 0)
    assert create_orac_luts.main([str(run)]) == 2
    captured = capsys.readouterr()
    assert "ORAC LUT Generator" in captured.out and f"(requested {foreign})" in captured.out
    assert f"requests LUT version {foreign}" in captured.err
    assert not called and not (tmp_path / "out").exists()


def test_run_file_with_the_declared_lut_version_passes_the_check(tmp_path, capsys, monkeypatch):
    text = (ROOT / "runs" / "template.run").read_text().replace("out_path      = 'luts'",
                                                               f"out_path      = '{tmp_path / 'out'}'")
    run = tmp_path / "declared.run"
    run.write_text(text)
    monkeypatch.setattr(create_orac_luts, "create_orac_cloud_lut", lambda *a, **k: 0)
    assert create_orac_luts.main([str(run)]) == 0
    assert "making ..." in capsys.readouterr().out


def test_configuration_entry_point_refuses_a_foreign_lut_version(tmp_path, capsys):
    output = ROOT / "validation" / "generated" / "test_version_mismatch_output.nc"
    output.unlink(missing_ok=True)
    try:
        with pytest.raises(SystemExit) as raised:
            generate_main(["--platform", "meteosat-10", "--instrument", "seviri", "--microphysics",
                           "liquid-water_stg.mm", "--test", "--channels", "1", "--version", "2",
                           "--output", str(output)])
        assert raised.value.code == 2
        captured = capsys.readouterr()
        assert "ORAC LUT Generator" in captured.out and "requests LUT version 2" in captured.err
        assert not output.exists()
    finally:
        output.unlink(missing_ok=True)


def test_configuration_entry_point_dry_run_prints_the_banner_without_rejecting(capsys):
    assert generate_main(["--platform", "meteosat-10", "--instrument", "seviri", "--microphysics",
                          "liquid-water_stg.mm", "--test", "--channels", "1", "--version", "2", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "ORAC LUT Generator" in out and "(requested 2)" in out and "dry-run" in out


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def _small_lut(path: Path) -> None:
    with Dataset(path, "w") as ds:
        ds.createDimension("optical_depth", 2)
        v = ds.createVariable("optical_depth", "f4", ("optical_depth",))
        v[:] = [1.0, 2.0]


def test_provenance_record_and_sidecar(tmp_path):
    lut = tmp_path / "x_m_liquid-water_a01_pstg_v25.nc"
    _small_lut(lut)
    config = tmp_path / "x.run"
    config.write_text("version = 25\n")
    record = provenance.build_record(
        lut, generator="test", lut_version=25,
        configuration={"platform": "x", "channels": [1, 2], "array": np.array([1.5], dtype=np.float32)},
        configuration_file=config, input_files={"instrument": config, "missing": tmp_path / "absent"},
        grid={"opd": {"n": 2, "values": np.array([1.0, 2.0])}}, numerics={"nmom": np.int32(1000)},
        radiative_transfer={"nstreams": 60})
    sidecar = provenance.write_provenance(lut, record)
    assert sidecar == tmp_path / "x_m_liquid-water_a01_pstg_v25.provenance.json"
    loaded = provenance.read_provenance(sidecar)
    assert loaded["lut_version"] == 25 and loaded["generator"] == "test"
    assert loaded["output"]["file"] == lut.name and loaded["output"]["sha256"] == version.sha256_file(lut)
    assert loaded["source"]["code_release"] == version.read_code_version(ROOT / "CODE_VERSION").code_release
    assert loaded["source"]["git"]["commit"] and loaded["source"]["git"]["working_tree"] in ("CLEAN", "MODIFIED")
    assert loaded["generated_utc"].endswith("Z")
    assert loaded["configuration_file"]["sha256"] == version.sha256_file(config)
    assert loaded["configuration_file"]["text"] == "version = 25\n"
    assert loaded["input_files"]["missing"] == {"path": str(tmp_path / "absent"), "exists": False}
    assert loaded["configuration"]["array"] == [1.5] and loaded["numerics"]["nmom"] == 1000
    assert loaded["source"]["software"]["numpy"] == np.__version__


def test_provenance_verification_detects_a_changed_product(tmp_path, capsys):
    lut = tmp_path / "y_v25.nc"
    _small_lut(lut)
    record = provenance.build_record(lut, generator="test", lut_version=25, configuration={}, configuration_file=None,
                                     input_files={}, grid={}, numerics={}, radiative_transfer={})
    sidecar = provenance.write_provenance(lut, record)
    result = version.verify_provenance(sidecar, root=ROOT)
    assert result["commit_in_repository"] is True and result["lut_sha256_matches"] is True
    assert result["consistent"] == (result["working_tree"] == "CLEAN")
    with Dataset(lut, "r+") as ds:
        ds.variables["optical_depth"][0] = 9.0
    result = version.verify_provenance(sidecar, root=ROOT)
    assert result["lut_sha256_matches"] is False and not result["consistent"]
    assert version.main(["--verify", str(sidecar), "--root", str(ROOT)]) == 1
    assert '"lut_sha256_matches": false' in capsys.readouterr().out
