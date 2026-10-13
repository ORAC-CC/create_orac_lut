"""Select the scientific code of a validation process: the sealed V23 release or the working (V24) tree.

   import codebase; codebase.select("v23")    # before anything imports oraclut or create_orac_luts

"v23" installs the release launcher's import finder for releases/lut_v23.release
(create_orac_luts_release_v04.ReleaseFinder), so every oraclut module and the
driver are the exact, SHA-256-verified files of LUT release V23 (the
lut-code-v22.4 scientific code).  "v24" uses the working tree (src/ and
create_orac_luts.py).  The two cannot be mixed in one process.  Read-only.
"""

import importlib
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
LAUNCHER = ROOT / "create_orac_luts_release_v04.py"
SELECTED = None
# The files that differ between V23 and V24 (their digests are printed when "v24" is selected)
V24_FILES = ("create_orac_luts.py", "src/oraclut/idl_mirror/create_bwgp.py",
             "src/oraclut/idl_mirror/generate_scattering_properties.py", "src/oraclut/idl_mirror/legendre_expansion.py")


def select(code):
   """Make 'v23' (sealed release modules) or 'v24' (working tree) importable; returns the driver module."""

   global SELECTED
   if SELECTED is not None:
      raise RuntimeError(f"the code base is already selected ({SELECTED})")
   if any(name.split(".", 1)[0] in ("oraclut", "create_orac_luts") for name in sys.modules):
      raise RuntimeError("oraclut was imported before the code base was selected")
   if code == "v23":
      spec = importlib.util.spec_from_file_location("oraclut_release_launcher_v04", LAUNCHER)
      launcher = importlib.util.module_from_spec(spec)
      spec.loader.exec_module(launcher)
      manifest = launcher.read_manifest(ROOT / "releases" / "lut_v23.release")
      launcher.verify_manifest_files(manifest, ("driver", "package", "module", "config", "adapter", "kernel"))
      sys.meta_path.insert(0, launcher.ReleaseFinder(manifest))
      driver = importlib.import_module("create_orac_luts")
      launcher.substitute_config_files()            # the sealed CODE_VERSION once the working one is V24's
   elif code == "v24":
      sys.path.insert(0, str(ROOT / "src"))
      sys.path.insert(0, str(ROOT))
      driver = importlib.import_module("create_orac_luts")
      # Record the exact working-tree revision of the V24 scattering code used
      import hashlib
      for name in V24_FILES:
         print(f"V24 source {name}: sha256 {hashlib.sha256((ROOT / name).read_bytes()).hexdigest()}", flush=True)
   else:
      raise ValueError(f"unknown code base {code!r}: 'v23' or 'v24'")
   SELECTED = code
   return driver


def module(name):
   """An oraclut module of the selected code base (the package attributes can shadow submodules)."""

   return importlib.import_module(name)
