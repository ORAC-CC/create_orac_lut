"""Release launcher v04: the corrected single-precision V22 release, revision 2 (original Intel DISORT kernel)."""

from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = ROOT / "create_orac_luts_release_v04.py"
MANIFEST = ROOT / "releases" / "lut_v22_single-r2.release"
SMOKE = ROOT / "validation/v21_v22_terra_modis/v21_v22_single/launcher_smoke.run"

pytestmark = pytest.mark.skipif(not MANIFEST.is_file(), reason="lut_v22_single.release not present")


def _launch(*args, **environment):
   env = {"HOME": str(Path.home()), "PATH": "/usr/bin:/bin", "ORAC_LUT_RELEASE_MANIFEST": "lut_v22_single-r2.release", **environment}
   return subprocess.run([sys.executable, str(LAUNCHER), *args], cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)


def _load():
   import importlib.util
   spec = importlib.util.spec_from_file_location("launcher_v04", LAUNCHER)
   module = importlib.util.module_from_spec(spec)
   spec.loader.exec_module(module)
   return module


def test_manifest_declares_single_precision_and_the_intel_kernel():
   launcher = _load()
   manifest = launcher.read_manifest(MANIFEST)
   assert manifest["header"]["disort_precision"] == "single" and manifest["header"]["disort_sp_kernel"] == "intel_dlm"
   logical = {entry["logical"]: entry for entry in manifest["entries"]}
   assert logical["disort-intel-library"]["path"] == "create_orac_lut/disort2/src/DISORTtoIDL.so"
   assert "disort-library" not in logical and "disort-double-library" not in logical
   runtime = [entry for entry in manifest["entries"] if entry["logical"].startswith("intel-runtime-")]
   assert len(runtime) == 7 and all(entry["role"] == "library" and Path(entry["path"]).is_absolute() for entry in runtime)


def test_launcher_v03_still_parses_every_manifest():
   # v03 scans releases/ when no manifest is named; the new manifest must not stop it.
   import importlib.util
   spec = importlib.util.spec_from_file_location("launcher_v03", ROOT / "create_orac_luts_release_v03.py")
   v03 = importlib.util.module_from_spec(spec)
   spec.loader.exec_module(v03)
   for path in sorted((ROOT / "releases").glob("lut_v*.release")):
      v03.read_manifest(path)


def test_runtime_entries_must_be_absolute_and_pinned(tmp_path):
   launcher = _load()
   text = MANIFEST.read_text()
   bad = tmp_path / "lut_v22_bad.release"
   bad.write_text(text.replace(" /network/software/", " network/software/", 1))
   with pytest.raises(launcher.ReleaseRefusal, match="absolute path"):
      launcher.read_manifest(bad)


def test_check_verifies_every_file():
   result = _launch("--check", "22")
   assert result.returncode == 0, result.stderr
   assert "all 58 files match the manifest" in result.stdout


@pytest.mark.parametrize(("environment", "message"), [
   ({}, "single-precision release, but ORAC_LUT_DISORT_PRECISION selects double"),
   ({"ORAC_LUT_DISORT_PRECISION": "single", "ORAC_LUT_DISORT_SP_KERNEL": "gfortran"}, "ORAC_LUT_DISORT_SP_KERNEL selects gfortran"),
])
def test_release_refuses_another_precision_or_kernel(environment, message):
   result = _launch(str(SMOKE), **environment)
   assert result.returncode == 2 and message in result.stderr
