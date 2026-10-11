"""Write releases/lut_v22_single.release: the corrected single-precision V22 release manifest.

   python scripts/make_lut_v22_single_manifest.py

LUT release V22, single precision: the Python generator with IDL-compatible
arithmetic that reproduces the V21 IDL lookup tables, solving DISORT with the
original V21 Intel kernel create_orac_lut/disort2/src/DISORTtoIDL.so (called
unmodified, with IDL's array memory layout).  The source set is the closure of
the active driver create_orac_luts.py (every oraclut module it imports, each
listed as its highest SourceVault version *_vNN.py, which must be
byte-identical to the working file); CODE_VERSION; the Python-owned C
adapters; the IDL-shared Fortran/C kernel sources and the DLM wrapper; the
gfortran Mie library (the Mie kernel build is numerically immaterial: the
original IDL Mie DLM gives bitwise-identical scattering properties on the full
Terra MODIS grids); the original Intel DISORT DLM; and the Intel 2022 Fortran
run-time libraries it loads (pinned by absolute path).

Products are made with launcher v04:
   ORAC_LUT_RELEASE_MANIFEST=lut_v22_single.release ORAC_LUT_DISORT_PRECISION=single \\
   ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>.run

The script refuses to overwrite an existing manifest (sealed once written).
"""

import datetime
import hashlib
import importlib
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "releases" / "lut_v22_single.release"


def sha256(path):
   return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def latest_version(working):
   working = Path(working)
   stem, suffix = working.stem, working.suffix
   versions = sorted(working.parent.glob(f"{stem}_v[0-9][0-9]{suffix}"))
   if not versions:
      sys.exit(f"{working} has no SourceVault version; run: version {working}")
   if sha256(versions[-1]) != sha256(working):
      sys.exit(f"{working} differs from its latest version {versions[-1].name}; run: version {working}")
   return versions[-1]


def main():
   if OUT.exists():
      sys.exit(f"refusing to overwrite sealed manifest {OUT}")
   os.environ["ORAC_LUT_DISORT_PRECISION"] = "single"
   os.environ["ORAC_LUT_DISORT_SP_KERNEL"] = "intel_dlm"
   sys.path.insert(0, str(ROOT / "src"))
   sys.path.insert(0, str(ROOT))
   importlib.import_module("create_orac_luts")
   legacy_disort = sys.modules["oraclut.radiative_transfer.legacy_disort"]
   modules = sorted(name for name in sys.modules if name == "oraclut" or name.startswith("oraclut."))
   rows = []

   def add(role, logical, path, version=None):
      path = Path(path)
      shown = path if path.is_absolute() and not str(path).startswith(str(ROOT) + os.sep) else path.relative_to(ROOT)
      rows.append(f"{role:<8s} {logical:<52s} {version or 'pinned':<7s} {sha256(path)}  {shown}")

   driver = latest_version(ROOT / "create_orac_luts.py")
   add("driver", "create_orac_luts", driver, driver.stem[-3:])
   for name in modules:
      working = Path(sys.modules[name].__file__)
      role = "package" if working.name == "__init__.py" else "module"
      version = latest_version(working)
      add(role, name, version, version.stem[-3:])
   add("config", "CODE_VERSION", ROOT / "CODE_VERSION")
   add("adapter", "mie-adapter", ROOT / "src/oraclut/optics/legacy_mie.c")
   add("adapter", "disort-adapter", ROOT / "src/oraclut/radiative_transfer/legacy_disort.c")
   add("kernel", "mie-mieint.f", ROOT / "mie/dlm-code/mieint.f")
   for name in ("DISORTfunctions.f", "BDREF.f", "PRTFIN.f", "ambrals-fortran.c", "ambralsfor.c", "ambralsfor.h",
                "disort/ErrPack.f", "disort/LINPAK.f", "disort/RDI1MACH.f", "DISORTwrapper.c", "DISORTtoIDL.c", "DISORTtoIDL.dlm"):
      add("kernel", "disort-" + Path(name).name, ROOT / "create_orac_lut/disort2/src" / name)
   add("library", "mie-library", ROOT / "build/oraclut/liboraclut_mie_gfortran.so")
   add("library", "disort-intel-library", legacy_disort.INTEL_DLM_LIBRARY)
   runtime = legacy_disort.intel_runtime_directory()
   for name in legacy_disort.INTEL_RUNTIME_LIBRARIES:
      add("runtime", "intel-runtime-" + name, runtime / name)

   head = subprocess.run(["/usr/bin/git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
   header = f"""\
# ORAC LUT release manifest
#
# LUT release V22, corrected single precision: the Python generator with
# IDL-compatible arithmetic that reproduces the IDL V21 lookup tables, with
# DISORT solved by the original V21 Intel kernel DISORTtoIDL.so (unmodified,
# IDL array memory layout).  NSTR 60, NMOM 1000, original LUT grids, float32
# NetCDF storage.  Validation: validation/v21_v22_terra_modis/
# REPORT_single_precision_equivalence.md.  A product is made with
#
#    ORAC_LUT_RELEASE_MANIFEST=lut_v22_single.release ORAC_LUT_DISORT_PRECISION=single \\
#    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>.run
#
# which refuses to run unless every file below has exactly this content and the
# environment selects single precision and the Intel kernel.  The historical
# single-precision V22 release (gfortran DISORT, before the IDL-compatible
# corrections) is releases/lut_v22.release; the double-precision V22 release is
# releases/lut_v22_double.release.  Both remain sealed.  This file is sealed:
# once a product has been made from it, it is never edited.
#
# LUT release numbers (V21, V22, ...) and source-file versions (_v01, _v02, ...)
# are independent: each file's version advances only when that file changes.
# "pinned" files have no versioned name; their content is fixed by SHA-256.
# runtime entries are outside the repository (absolute paths).
#
# Columns: role  logical-name  file-version  sha256  path

manifest_schema = oraclut-release-manifest-1
lut_release = 22
disort_precision = single
disort_sp_kernel = intel_dlm
library_checkout = {ROOT}
created = {datetime.date.today().isoformat()}
git_head_at_creation = {head}

"""
   OUT.write_text(header + "\n".join(rows) + "\n")
   print(f"wrote {OUT} ({len(rows)} entries)")


if __name__ == "__main__":
   main()
