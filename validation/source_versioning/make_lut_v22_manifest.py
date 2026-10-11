"""Write releases/lut_v22.release, the source manifest of LUT release V22 (record of how it was made).

   python validation/source_versioning/make_lut_v22_manifest.py

The V22 source set is the production closure traced on 2026-10-10
(validation/source_versioning/REPORT_source_versioning_v01.md): the driver,
the 29 oraclut modules and package initialisers it loads, the CODE_VERSION
file it reads, the Python-owned C adapters and IDL-shared Fortran/C kernels
compiled into the two kernel libraries, and the two libraries themselves as
built from those sources in this checkout.  Every Python file is listed as
its *_v01 copy; the compiled sources and libraries are listed under their
own names with file version "pinned" (SHA-256 fixed, no versioned name yet).

The script refuses to overwrite an existing manifest: a manifest is sealed
once written.
"""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "releases" / "lut_v22.release"

ENTRIES = [
   ("driver", "create_orac_luts", "create_orac_luts_v01.py"),
   ("package", "oraclut", "src/oraclut/__init___v01.py"),
   ("module", "oraclut.version", "src/oraclut/version_v01.py"),
   ("module", "oraclut.provenance", "src/oraclut/provenance_v01.py"),
   ("package", "oraclut.idl_mirror", "src/oraclut/idl_mirror/__init___v01.py"),
]
for name in ("load_inststr", "load_lutstr", "load_srfstrarr", "segment", "load_mmdat", "load_atmstr", "load_gasstr",
             "interpol", "generate_scattering_properties", "create_bwgp", "baum", "dubovik", "setup_disort",
             "call_disort", "write_v2_lut"):
   ENTRIES.append(("module", f"oraclut.idl_mirror.{name}", f"src/oraclut/idl_mirror/{name}_v01.py"))
ENTRIES += [
   ("package", "oraclut.io", "src/oraclut/io/__init___v01.py"),
   ("module", "oraclut.io.lut", "src/oraclut/io/lut_v01.py"),
   ("module", "oraclut.io.v2", "src/oraclut/io/v2_v01.py"),
   ("module", "oraclut.io.string_attributes", "src/oraclut/io/string_attributes_v01.py"),
   ("module", "oraclut.io.netcdf_c", "src/oraclut/io/netcdf_c_v01.py"),
   ("package", "oraclut.optics", "src/oraclut/optics/__init___v01.py"),
   ("module", "oraclut.optics.legacy", "src/oraclut/optics/legacy_v01.py"),
   ("module", "oraclut.optics.legacy_mie", "src/oraclut/optics/legacy_mie_v01.py"),
   ("package", "oraclut.radiative_transfer", "src/oraclut/radiative_transfer/__init___v01.py"),
   ("module", "oraclut.radiative_transfer.legacy_disort", "src/oraclut/radiative_transfer/legacy_disort_v01.py"),
   ("config", "CODE_VERSION", "CODE_VERSION"),
   ("adapter", "mie-adapter", "src/oraclut/optics/legacy_mie.c"),
   ("adapter", "disort-adapter", "src/oraclut/radiative_transfer/legacy_disort.c"),
   ("kernel", "mie-mieint.f", "mie/dlm-code/mieint.f"),
]
for name in ("DISORTfunctions.f", "BDREF.f", "PRTFIN.f", "ambrals-fortran.c", "ambralsfor.c", "ambralsfor.h",
             "disort/ErrPack.f", "disort/LINPAK.f", "disort/RDI1MACH.f"):
   ENTRIES.append(("kernel", "disort-" + Path(name).name, f"create_orac_lut/disort2/src/{name}"))
ENTRIES += [
   ("library", "mie-library", "build/oraclut/liboraclut_mie_gfortran.so"),
   ("library", "disort-library", "build/oraclut/liboraclut_disort_gfortran.so"),
]

HEADER = """\
# ORAC LUT release manifest
#
# LUT release V22: the Python generator that reproduces the IDL V21 lookup
# tables.  This file lists, with SHA-256, every source file that makes up the
# V22 generator.  A V22 product is made with
#
#    python create_orac_luts_release_v01.py runs/<product>_v22.run
#
# which refuses to run unless every file below has exactly this content, loads
# the Python modules from exactly these files and records this manifest in the
# product.  This file is sealed: once a product has been made from it, it is
# never edited.  A new LUT release gets a new manifest.
#
# LUT release numbers (V21, V22, ...) and source-file versions (_v01, _v02, ...)
# are independent: each file's version advances only when that file changes.
# "pinned" files have no versioned name yet; their content is fixed by SHA-256.
# The two library entries are the kernel libraries built from the pinned
# sources in the checkout named by library_checkout (the DISORT library embeds
# its absolute source paths).
#
# Columns: role  logical-name  file-version  sha256  path (relative to the repository)
# Roles: driver, package, module (Python, loaded by logical name), config (read
#   at run time), adapter (Python-owned C source), kernel (IDL-shared Fortran/C
#   source), library (compiled kernel library).

manifest_schema = oraclut-release-manifest-1
lut_release = 22
created = 2026-10-10
library_checkout = /home/g/grainger/project-oraclut

"""


def file_version(path):
   stem = Path(path).stem
   if len(stem) > 4 and stem[-4:-2] == "_v" and stem[-2:].isdigit():
      return stem[-3:]
   return "pinned"


def main():
   if OUT.exists():
      print(f"{OUT} exists; a release manifest is sealed and is never rewritten", file=sys.stderr)
      return 1
   lines = [HEADER]
   for role, logical, path in ENTRIES:
      data = (ROOT / path).read_bytes()
      lines.append(f"{role:<8} {logical:<51} {file_version(path):<7} {hashlib.sha256(data).hexdigest()}  {path}\n")
   OUT.parent.mkdir(exist_ok=True)
   OUT.write_text("".join(lines))
   print(f"wrote {OUT}: {len(ENTRIES)} entries")
   return 0


if __name__ == "__main__":
   sys.exit(main())
