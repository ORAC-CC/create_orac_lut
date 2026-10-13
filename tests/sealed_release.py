"""Run test code with the scientific code of a sealed LUT release (not the working tree).

The V21/V22 reproduction tests check the lut-code-v22.4 scientific code, which
the working tree no longer holds once a later release (V24) changes it.  They
run that code exactly as the release launcher does: a fresh Python process in
which create_orac_luts_release_v04.ReleaseFinder maps every oraclut module and
the driver to the SHA-256-verified files of the release manifest
(releases/lut_v23.release = the lut-code-v22.4 code with Grid B; the grid is
irrelevant to these tests, which supply their own definitions).
"""

import os
from pathlib import Path
import subprocess
import sys
import textwrap

ROOT = Path(__file__).resolve().parents[1]
V22_SCIENCE_MANIFEST = ROOT / "releases" / "lut_v23.release"

PRELUDE = """
import importlib, importlib.util, sys
spec = importlib.util.spec_from_file_location("release_launcher_v04", {launcher!r})
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)
manifest = launcher.read_manifest({manifest!r})
launcher.verify_manifest_files(manifest, ("driver", "package", "module", "config", "adapter", "kernel"))
sys.meta_path.insert(0, launcher.ReleaseFinder(manifest))
create_orac_luts = importlib.import_module("create_orac_luts")
launcher.substitute_config_files()
"""


def run_with_release(code, manifest=V22_SCIENCE_MANIFEST, env=None, timeout=3600):
   """Run ``code`` (Python source) after installing the release's import finder; raise on failure."""

   source = PRELUDE.format(launcher=str(ROOT / "create_orac_luts_release_v04.py"), manifest=str(manifest))
   source += textwrap.dedent(code)
   environment = {key: value for key, value in os.environ.items() if not key.startswith("PYTEST")}
   environment.pop("PYTHONPATH", None)
   environment.update(env or {})
   result = subprocess.run([sys.executable, "-c", source], cwd=ROOT, env=environment, capture_output=True, text=True,
                           timeout=timeout)
   if result.returncode != 0:
      raise AssertionError(f"sealed-release run failed ({result.returncode}):\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}")
   return result
