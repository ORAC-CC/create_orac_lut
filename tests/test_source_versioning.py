"""Per-file source versioning, the sealed V22 manifest and release launcher.

The fast tests check that each working source has an identical latest numbered
copy, that every sealed *_v01 copy equals the lut-code-v22.1 blob, the manifest
rules, the launcher's refusals and module binding.  The slow test runs a compact
LUT through the sealed launcher and checks provenance, the existing LUT reader
and the inspection tool.  It deliberately does not compare a historical
release with a newer unversioned working implementation.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]
# The sealed V22 release is driven by the current launcher (v03 selects among
# several V22 manifests and serves the sealed CODE_VERSION from CODE_VERSION_v01).
LAUNCHER = ROOT / "create_orac_luts_release_v03.py"
MANIFEST = ROOT / "releases" / "lut_v22.release"
INSPECT = ROOT / "validation" / "source_versioning" / "inspect_lut_provenance.py"
PYTHON = sys.executable
GIT = shutil.which("git") or "/usr/bin/git"

CLOSURE = {
    "create_orac_luts": "create_orac_luts.py",
    "oraclut": "src/oraclut/__init__.py",
    "oraclut.version": "src/oraclut/version.py",
    "oraclut.provenance": "src/oraclut/provenance.py",
    "oraclut.idl_mirror": "src/oraclut/idl_mirror/__init__.py",
    "oraclut.io": "src/oraclut/io/__init__.py",
    "oraclut.io.lut": "src/oraclut/io/lut.py",
    "oraclut.io.v2": "src/oraclut/io/v2.py",
    "oraclut.io.string_attributes": "src/oraclut/io/string_attributes.py",
    "oraclut.io.netcdf_c": "src/oraclut/io/netcdf_c.py",
    "oraclut.optics": "src/oraclut/optics/__init__.py",
    "oraclut.optics.legacy": "src/oraclut/optics/legacy.py",
    "oraclut.optics.legacy_mie": "src/oraclut/optics/legacy_mie.py",
    "oraclut.radiative_transfer": "src/oraclut/radiative_transfer/__init__.py",
    "oraclut.radiative_transfer.legacy_disort": "src/oraclut/radiative_transfer/legacy_disort.py",
}
for _name in ("load_inststr", "load_lutstr", "load_srfstrarr", "segment", "load_mmdat", "load_atmstr", "load_gasstr",
              "interpol", "generate_scattering_properties", "create_bwgp", "baum", "dubovik", "setup_disort",
              "call_disort", "write_v2_lut"):
    CLOSURE[f"oraclut.idl_mirror.{_name}"] = f"src/oraclut/idl_mirror/{_name}.py"


def v01(path):
    path = Path(path)
    return path.with_name(path.stem + "_v01" + path.suffix)


def latest_version(path):
    path = Path(path)
    versions = sorted(path.parent.glob(path.stem + "_v[0-9][0-9]" + path.suffix))
    assert versions, f"no SourceVault version exists for {path}"
    return versions[-1]


def load_launcher():
    spec = importlib.util.spec_from_file_location("release_launcher_under_test", LAUNCHER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(args, cwd=ROOT, timeout=900, manifest="lut_v22.release"):
    environment = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}
    environment.pop("ORAC_LUT_RELEASE_MANIFEST", None)
    if manifest:
        environment["ORAC_LUT_RELEASE_MANIFEST"] = manifest
    return subprocess.run([PYTHON, *args], cwd=cwd, capture_output=True, text=True, timeout=timeout, env=environment)


def write_run_file(path, *, version=22, out_path, extra=""):
    text = (ROOT / "validation" / "source_versioning" / "runs" / "case034_launcher.run").read_text()
    text = text.replace("version       = 22", f"version       = {version}")
    text = "\n".join(line for line in text.splitlines() if not line.startswith("out_path"))
    path.write_text(text + f"\nout_path      = '{out_path}'\n{extra}")
    return path


# ---------------------------------------------------------------------------
# Transition guard and manifest
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("relative", sorted(set(CLOSURE.values()) | {
    "src/oraclut/idl_total.py", "src/oraclut/idl_mirror/total.py", "src/oraclut/pipeline.py",
}))
def test_latest_sourcevault_copy_is_byte_identical_to_the_working_file(relative):
    original = ROOT / relative
    archived = latest_version(original)
    assert archived.read_bytes() == original.read_bytes(), (
        f"{relative} differs from its latest SourceVault copy {archived.name}")


def test_v01_copies_equal_the_lut_code_v22_1_files():
    for logical, relative in CLOSURE.items():
        blob = subprocess.run([GIT, "-C", str(ROOT), "show", f"lut-code-v22.1:{relative}"], capture_output=True)
        if blob.returncode != 0:
            pytest.skip("Git or the lut-code-v22.1 tag is unavailable")
        assert v01(ROOT / relative).read_bytes() == blob.stdout, relative


def test_manifest_lists_exactly_the_production_closure_and_matches_the_files():
    launcher = load_launcher()
    manifest = launcher.read_manifest(MANIFEST)
    assert manifest["header"]["lut_release"] == 22
    python = {entry["logical"]: entry for entry in manifest["entries"] if entry["role"] in launcher.PYTHON_ROLES}
    assert set(python) == set(CLOSURE)
    for logical, entry in python.items():
        assert entry["path"] == str(v01(CLOSURE[logical]))
        assert entry["file_version"] == "v01"
    roles = sorted({entry["role"] for entry in manifest["entries"]})
    assert roles == ["adapter", "config", "driver", "kernel", "library", "module", "package"]
    launcher.verify_manifest_files(manifest, launcher.ROLES)     # raises on any mismatch


@pytest.mark.parametrize("bad_line, message", [
    ("module   oraclut.extra   v01  " + "0" * 64 + "  src/oraclut/extra.py", "must be a versioned"),
    ("module   oraclut.extra   v02  " + "0" * 64 + "  src/oraclut/extra_v01.py", "does not match the file name"),
    ("kernel   k.f             pinned  " + "0" * 64 + "  k_v01.f", "has a version suffix"),
    ("thing    x               pinned  " + "0" * 64 + "  x.f", "unknown role"),
    ("module   oraclut.version v01  " + "0" * 64 + "  src/oraclut/other_v01.py", "duplicate logical"),
    ("module   oraclut.extra   v01  " + "z" * 64 + "  src/oraclut/extra_v01.py", "malformed SHA-256"),
])
def test_malformed_manifest_is_refused(tmp_path, bad_line, message):
    launcher = load_launcher()
    path = tmp_path / "lut_v22.release"
    path.write_text(MANIFEST.read_text() + bad_line + "\n")
    with pytest.raises(launcher.ReleaseRefusal, match=message):
        launcher.read_manifest(path)


def test_run_file_version_is_read_with_the_driver_comment_rules(tmp_path):
    launcher = load_launcher()
    path = tmp_path / "x.run"
    path.write_text("; version = 99 (comment)\nversion = 22 ; trailing comment\n# version = 98\n")
    assert launcher.read_run_file_version(path) == 22
    path.write_text("version = '22'\n")
    with pytest.raises(launcher.ReleaseRefusal, match="integer version"):
        launcher.read_run_file_version(path)
    path.write_text("platform = 'x'\n")
    with pytest.raises(launcher.ReleaseRefusal, match="integer version"):
        launcher.read_run_file_version(path)


# ---------------------------------------------------------------------------
# Launcher refusals and binding (subprocesses; no LUT calculation)
# ---------------------------------------------------------------------------

def test_a_run_file_without_a_release_manifest_is_refused(tmp_path):
    run_file = write_run_file(tmp_path / "v23.run", version=23, out_path=str(tmp_path / "out"))
    result = run([str(LAUNCHER), str(run_file)], manifest=None)
    assert result.returncode == 2
    assert "no release manifest defines LUT release V23" in result.stderr
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("option", ["reuse_scat", "scat_only"])
def test_cached_optics_and_scattering_only_runs_are_refused(tmp_path, option):
    out = ROOT / "validation" / "tmp" / f"pytest-source-versioning-{uuid.uuid4().hex[:8]}"
    run_file = write_run_file(tmp_path / "r.run", out_path=str(out))
    import re
    text, count = re.subn(rf"^{option}\s*=\s*0", f"{option} = 1", run_file.read_text(), flags=re.MULTILINE)
    assert count == 1
    run_file.write_text(text)
    result = run([str(LAUNCHER), str(run_file)])
    assert result.returncode == 2, result.stdout + result.stderr
    assert f"{option} = 1 is not allowed in a release run" in result.stderr
    assert not out.exists()


def test_a_modified_sealed_file_is_refused_before_anything_is_imported(tmp_path):
    launcher = load_launcher()
    manifest = launcher.read_manifest(MANIFEST)
    copy = tmp_path / "repo"
    for relative in [LAUNCHER.name, "releases/lut_v22.release", "CODE_VERSION_v01"] + [entry["path"] for entry in manifest["entries"]]:
        (copy / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, copy / relative)
    target = copy / "src" / "oraclut" / "idl_mirror" / "create_bwgp_v01.py"
    target.write_text(target.read_text() + "\n# an edit\n")
    run_file = write_run_file(tmp_path / "r.run", out_path=str(copy / "validation" / "out"))
    result = run([str(copy / LAUNCHER.name), str(run_file)], cwd=copy)
    assert result.returncode == 2
    assert "MODIFIED src/oraclut/idl_mirror/create_bwgp_v01.py" in result.stderr
    assert "new file version" in result.stderr


def test_the_finder_binds_every_module_to_its_v01_file_and_refuses_others():
    probe = f"""
import importlib, importlib.util, json, sys
spec = importlib.util.spec_from_file_location("launcher", {str(LAUNCHER)!r})
launcher = importlib.util.module_from_spec(spec); spec.loader.exec_module(launcher)
manifest = launcher.read_manifest({str(MANIFEST)!r})
finder = launcher.ReleaseFinder(manifest)
sys.meta_path.insert(0, finder)
driver = importlib.import_module("create_orac_luts")
launcher.verify_bindings(finder)
files = {{name: module.__file__ for name, module in sys.modules.items() if name.split(".")[0] in ("oraclut", "create_orac_luts")}}
refused = {{}}
for name in ("oraclut.generate", "oraclut.io.v2_v01", "create_orac_luts_v01"):
    try:
        importlib.import_module(name)
        refused[name] = False
    except ImportError as exc:
        refused[name] = str(exc)
print(json.dumps({{"files": files, "refused": refused, "loaded": sorted(finder.loaded)}}))
"""
    result = run(["-c", probe])
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout.splitlines()[-1])
    assert set(report["files"]) == set(CLOSURE)
    for name, file in report["files"].items():
        assert Path(file).resolve() == v01(ROOT / CLOSURE[name]).resolve(), name
    assert sorted(report["loaded"]) == sorted(CLOSURE)
    assert "not part of LUT release V22" in report["refused"]["oraclut.generate"]
    assert "versioned files are loaded only under their logical names" in report["refused"]["oraclut.io.v2_v01"]
    assert "versioned files are loaded only under their logical names" in report["refused"]["create_orac_luts_v01"]


def test_a_process_that_already_imported_oraclut_is_refused(tmp_path):
    run_file = write_run_file(tmp_path / "r.run", out_path="validation/tmp/never-written")
    probe = f"""
import sys, importlib.util
sys.path.insert(0, {str(ROOT / 'src')!r})
import oraclut.version
spec = importlib.util.spec_from_file_location("launcher", {str(LAUNCHER)!r})
launcher = importlib.util.module_from_spec(spec); spec.loader.exec_module(launcher)
sys.exit(launcher.main([{str(run_file)!r}]))
"""
    result = run(["-c", probe])
    assert result.returncode == 2
    assert "was imported before the release finder was installed" in result.stderr


def test_products_outside_validation_need_committed_release_files(tmp_path):
    launcher = load_launcher()
    state = launcher.git_state(launcher.read_manifest(MANIFEST), [LAUNCHER.name, "releases/lut_v22.release"])
    if not state.get("available") or not state["release_files_not_in_commit"]:
        pytest.skip("the release files are committed; this run would start a real calculation")
    run_file = write_run_file(tmp_path / "r.run", out_path=str(tmp_path / "production"))
    result = run([str(LAUNCHER), str(run_file)])
    assert result.returncode == 2
    assert "products outside validation/ require every release file" in result.stderr
    assert not (tmp_path / "production").exists()


# ---------------------------------------------------------------------------
# Inspection tool
# ---------------------------------------------------------------------------

def make_product(path, source_root):
    (source_root / "src").mkdir(parents=True)
    (source_root / "src" / "a_v01.py").write_text("x = 1\n")
    (source_root / "k.f").write_text("      END\n")
    (source_root / "releases").mkdir()
    sha = lambda p: __import__("hashlib").sha256(Path(p).read_bytes()).hexdigest()
    manifest = source_root / "releases" / "lut_v22.release"
    manifest.write_text(f"module a v01 {sha(source_root / 'src' / 'a_v01.py')} src/a_v01.py\n"
                        f"kernel k.f pinned {sha(source_root / 'k.f')} k.f\n")
    (source_root / "launch_v01.py").write_text("pass\n")
    record = {
        "provenance_schema": "oraclut-lut-provenance-1", "generated_utc": "2026-10-10T00:00:00Z", "started_utc": "2026-10-10T00:00:00Z",
        "driver": {"path": "src/a_v01.py", "file_version": "v01"},
        "release_manifest": {"path": "releases/lut_v22.release", "sha256": sha(manifest)},
        "launcher": {"path": "launch_v01.py", "sha256": sha(source_root / "launch_v01.py")},
        "run_file": {"path": "r.run", "sha256": "0" * 64}, "generation_command": "python launch_v01.py r.run",
        "git_at_start": {"available": False, "detail": "test"}, "repository": {"public": "https://example"},
        "source_files": [{"role": "module", "logical": "a", "file_version": "v01", "sha256": sha(source_root / "src" / "a_v01.py"), "path": "src/a_v01.py"},
                         {"role": "kernel", "logical": "k.f", "file_version": "pinned", "sha256": sha(source_root / "k.f"), "path": "k.f"}],
        "inputs_read": [],
    }
    with Dataset(path, "w") as dataset:
        dataset.setncattr_string("lut_release", "V22")
        dataset.setncattr_string("provenance_json", json.dumps(record))


def test_inspection_tool_reports_match_mismatch_and_missing(tmp_path):
    product = tmp_path / "p.nc"
    source_root = tmp_path / "sources"
    make_product(product, source_root)
    result = run([str(INSPECT), str(product), "--source-root", str(source_root)])
    assert result.returncode == 0, result.stdout + result.stderr
    assert "MATCH" in result.stdout and "the embedded source list equals releases/lut_v22.release" in result.stdout
    (source_root / "k.f").write_text("      STOP\n      END\n")
    result = run([str(INSPECT), str(product), "--source-root", str(source_root)])
    assert result.returncode == 1 and "MISMATCH" in result.stdout
    (source_root / "src" / "a_v01.py").unlink()
    result = run([str(INSPECT), str(product), "--source-root", str(source_root)])
    assert result.returncode == 1 and "MISSING" in result.stdout


def test_inspection_tool_infers_nothing_for_a_product_without_provenance(tmp_path):
    product = tmp_path / "old.nc"
    with Dataset(product, "w") as dataset:
        dataset.createDimension("x", 1)
    result = run([str(INSPECT), str(product)])
    assert result.returncode == 3
    assert "no embedded release provenance" in result.stdout


# ---------------------------------------------------------------------------
# End to end (compact LUT; minutes)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_launcher_product_carries_provenance(tmp_path):
    from oraclut.io import netcdf_c
    from oraclut.io.lut import read_lut
    from oraclut.io.string_attributes import audit_file

    base = ROOT / "validation" / "tmp" / f"pytest-source-versioning-{uuid.uuid4().hex[:8]}"
    try:
        launcher_run = write_run_file(tmp_path / "launcher.run", out_path=str(base / "launcher"))
        result = run([str(LAUNCHER), str(launcher_run)], timeout=3600)
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
        name = "meteosat-10_seviri_m_liquid-water_a01_pstg_v22.nc"
        released = base / "launcher" / name

        with Dataset(released) as a:
            assert a.getncattr("lut_release") == "V22"
            record = json.loads(a.getncattr("provenance_json"))
        assert audit_file(released).compliant
        types = {item.name: item.type_name for item in netcdf_c.list_attributes(released) if item.variable is None}
        assert types["lut_release_version"] == "NC_INT"
        assert all(t == "NC_STRING" for n, t in types.items() if n != "lut_release_version")
        assert all(item["loaded"] for item in record["source_files"] if item["role"] in ("driver", "package", "module"))
        assert {"create_orac_lut/input_files/sun/Gueymard2018.sssi", "create_orac_lut/input_files/atm/mls.atm"} <= {
            item["path"] for item in record["inputs_read"]}
        lut = read_lut(released)                       # the existing reader is unaffected
        assert lut.variables["R_0v"].ndim == 6
        result = run([str(INSPECT), str(released), "--verify-inputs"])
        assert result.returncode == 0, result.stdout
    finally:
        shutil.rmtree(base, ignore_errors=True)
