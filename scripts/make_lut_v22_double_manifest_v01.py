"""Write releases/lut_v22_double.release: the source manifest of LUT release V22 with double-precision DISORT.

   python scripts/make_lut_v22_double_manifest.py

The source set is the closure of the active driver create_orac_luts.py: every
oraclut module it imports (determined by importing it), each listed as the
highest SourceVault version *_vNN.py, which must be byte-identical to the
working file; CODE_VERSION; the Python-owned C adapters; the IDL-shared
Fortran/C kernels; the double-precision DISORT build sources; and the three
compiled libraries (Mie, single-precision DISORT for GETMOM/PLKAVG, and the
double-precision DISORT kernel create_orac_lut/disort2_dp/libdisort_dp.so).

Products are made with
   ORAC_LUT_RELEASE_MANIFEST=lut_v22_double.release python create_orac_luts_release_v03.py runs/<product>.run

The script refuses to overwrite an existing manifest (sealed once written).
"""

import hashlib
import importlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "releases" / "lut_v22_double.release"


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
   sys.path.insert(0, str(ROOT / "src"))
   sys.path.insert(0, str(ROOT))
   importlib.import_module("create_orac_luts")
   modules = sorted(name for name in sys.modules if name == "oraclut" or name.startswith("oraclut."))
   rows = []

   def add(role, logical, path, version=None):
      path = Path(path)
      if version is None:
         version = "pinned"
      rows.append(f"{role:<8s} {logical:<52s} {version:<7s} {sha256(path)}  {path.relative_to(ROOT)}")

   driver = latest_version(ROOT / "create_orac_luts.py")
   add("driver", "create_orac_luts", driver, driver.stem[-3:])
   for name in modules:
      module = sys.modules[name]
      working = Path(module.__file__)
      role = "package" if working.name == "__init__.py" else "module"
      version = latest_version(working)
      add(role, name, version, version.stem[-3:])
   add("config", "CODE_VERSION", ROOT / "CODE_VERSION")
   add("adapter", "mie-adapter", ROOT / "src/oraclut/optics/legacy_mie.c")
   add("adapter", "disort-adapter", ROOT / "src/oraclut/radiative_transfer/legacy_disort.c")
   add("kernel", "mie-mieint.f", ROOT / "mie/dlm-code/mieint.f")
   for name in ("DISORTfunctions.f", "BDREF.f", "PRTFIN.f", "ambrals-fortran.c", "ambralsfor.c", "ambralsfor.h",
                "disort/ErrPack.f", "disort/LINPAK.f", "disort/RDI1MACH.f"):
      add("kernel", "disort-" + Path(name).name, ROOT / "create_orac_lut/disort2/src" / name)
   for name in ("dp_support.f", "build_dlm.sh", "BUILD_INFO.txt"):
      add("kernel", "disort-double-" + name, ROOT / "create_orac_lut/disort2_dp" / name)
   add("library", "mie-library", ROOT / "build/oraclut/liboraclut_mie_gfortran.so")
   add("library", "disort-library", ROOT / "build/oraclut/liboraclut_disort_gfortran.so")
   add("library", "disort-double-library", ROOT / "create_orac_lut/disort2_dp/libdisort_dp.so")

   head = subprocess.run(["/usr/bin/git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
   header = f"""\
# ORAC LUT release manifest
#
# LUT release V22 with double-precision DISORT: the Python generator that
# reproduces the IDL V21 lookup tables (active revisions with IDL-compatible
# TOTAL reductions), with DISORT solved by the double-precision kernel shared
# with the IDL DLM DISORTtoIDL_DP (create_orac_lut/disort2_dp).  NSTR 60,
# NMOM 1000, original LUT grids, float32 NetCDF storage.  A product is made with
#
#    ORAC_LUT_RELEASE_MANIFEST=lut_v22_double.release python create_orac_luts_release_v03.py runs/<product>.run
#
# which refuses to run unless every file below has exactly this content (and
# ORAC_LUT_DISORT_PRECISION selects double, the default).  The historical
# single-precision V22 release is releases/lut_v22.release (sealed; its
# products are reproduced from Git tag lut-code-v22.1).  This file is sealed:
# once a product has been made from it, it is never edited.
#
# LUT release numbers (V21, V22, ...) and source-file versions (_v01, _v02, ...)
# are independent: each file's version advances only when that file changes.
# "pinned" files have no versioned name; their content is fixed by SHA-256.
# The library entries are built in the checkout named by library_checkout.
#
# Columns: role  logical-name  file-version  sha256  path (relative to the repository)

manifest_schema = oraclut-release-manifest-1
lut_release = 22
disort_precision = double
library_checkout = {ROOT}
created = {__import__('datetime').date.today().isoformat()}
git_head_at_creation = {head}

"""
   OUT.write_text(header + "\n".join(rows) + "\n")
   print(f"wrote {OUT} ({len(rows)} entries)")


if __name__ == "__main__":
   main()
