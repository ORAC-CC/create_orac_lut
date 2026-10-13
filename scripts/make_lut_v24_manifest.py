"""Write releases/lut_v24.release: LUT release V24 = the V23 release with the three V24 scattering refinements.

   python scripts/make_lut_v24_manifest.py

V24 (V24_RELEASE.md) changes the scientific code of V23 (releases/lut_v23.release,
lut-code-v22.4) in exactly four files, and nothing else:

* the driver create_orac_luts (v04 -> v05: max(nmom, L) DISORT moments per
  phase function, provenance description, scattering-cache scheme);
* oraclut.idl_mirror.create_bwgp (v03 -> v05: 3.5 r_e lattice limit, nested
  radius refinement, blocked Mie sums, double-precision phase function on
  request; v04 differs from v03 only in a docstring);
* oraclut.idl_mirror.generate_scattering_properties (v03 -> v04: the adaptive
  path for Mie classes);
* oraclut.idl_mirror.legendre_expansion (new, v01).

Every other entry of the V23 manifest (modules, adapters, DISORT and Mie
kernel sources, the Intel DISORT library and run time, and the pinned input
files: Grid B liquid-water-cloud.lut / ice-cloud.lut, instrument, microphysics,
refractive indices, atmosphere, SRFs, solar spectrum) is copied unchanged,
verified against the files on disk.  CODE_VERSION is pinned with its V24
content (code_release lut-code-v24.1, lut_version 24).

Products are made with launcher v04:
   ORAC_LUT_RELEASE_MANIFEST=lut_v24.release ORAC_LUT_DISORT_PRECISION=single \\
   ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>_v24.run

The script refuses to overwrite an existing manifest (sealed once written).
"""

import datetime
import hashlib
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
V23 = ROOT / "releases" / "lut_v23.release"
OUT = ROOT / "releases" / "lut_v24.release"
# logical name -> (role, file version, path) of the V24 files
CHANGED = {
   "create_orac_luts": ("driver", "v05", "create_orac_luts_v05.py"),
   "oraclut.idl_mirror.create_bwgp": ("module", "v05", "src/oraclut/idl_mirror/create_bwgp_v05.py"),
   "oraclut.idl_mirror.generate_scattering_properties": ("module", "v04", "src/oraclut/idl_mirror/generate_scattering_properties_v04.py"),
}
ADDED = {"oraclut.idl_mirror.legendre_expansion": ("module", "v01", "src/oraclut/idl_mirror/legendre_expansion_v01.py")}
# The working files must equal these SourceVault copies (no unversioned change can enter the release)
WORKING = {"create_orac_luts_v05.py": "create_orac_luts.py",
           "src/oraclut/idl_mirror/create_bwgp_v05.py": "src/oraclut/idl_mirror/create_bwgp.py",
           "src/oraclut/idl_mirror/generate_scattering_properties_v04.py": "src/oraclut/idl_mirror/generate_scattering_properties.py",
           "src/oraclut/idl_mirror/legendre_expansion_v01.py": "src/oraclut/idl_mirror/legendre_expansion.py"}


def sha256(path):
   return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def row(role, logical, version, digest, path):
   return f"{role:<8s} {logical:<52s} {version:<7s} {digest}  {path}"


def main():
   if OUT.exists():
      sys.exit(f"refusing to overwrite sealed manifest {OUT}")
   for copy, working in WORKING.items():
      if sha256(ROOT / copy) != sha256(ROOT / working):
         sys.exit(f"{working} differs from its SourceVault copy {copy}")
   rows, replaced = [], set()
   for line in V23.read_text().splitlines():
      fields = line.split("#", 1)[0].split()
      if len(fields) != 5 or fields[0] == "config":
         continue
      role, logical, version, digest, path = fields
      full = Path(path) if Path(path).is_absolute() else ROOT / path
      if sha256(full) != digest:
         sys.exit(f"{path} no longer matches the V23 manifest")
      if logical in CHANGED:
         new_role, new_version, new_path = CHANGED[logical]
         if new_role != role:
            sys.exit(f"{logical}: role {role} in V23, {new_role} expected")
         rows.append(row(role, logical, new_version, sha256(ROOT / new_path), new_path))
         replaced.add(logical)
         if logical == "oraclut.idl_mirror.generate_scattering_properties":
            role, version, path = ADDED["oraclut.idl_mirror.legendre_expansion"]
            rows.append(row(role, "oraclut.idl_mirror.legendre_expansion", version, sha256(ROOT / path), path))
         continue
      rows.append(row(role, logical, version, digest, path))
   if replaced != set(CHANGED):
      sys.exit(f"V23 manifest lacks {sorted(set(CHANGED) - replaced)}")
   rows.append(row("config", "CODE_VERSION", "pinned", sha256(ROOT / "CODE_VERSION"), "CODE_VERSION"))
   head = subprocess.run(["/usr/bin/git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
   header = f"""\
# ORAC LUT release manifest
#
# LUT release V24: the V23 release with three refinements of the scattering
# properties of Mie classes (V24_RELEASE.md).
#   LUT release      V24
#   scientific code  lut-code-v24.1 = lut-code-v22.4 (V22/V23) with the driver
#                    create_orac_luts_v05, create_bwgp_v05,
#                    generate_scattering_properties_v04 and the new
#                    legendre_expansion_v01; every other entry is copied from
#                    releases/lut_v23.release
#   refinements      liquid-water size integration to 3.5 r_e on the legacy
#                    radius lattice; nested radius refinement (liquid water,
#                    ice spheres); adaptive Legendre expansion, DISORT
#                    NMOM = max(1000, L); Baum and T-matrix classes as V23
#   grid             Grid B, as V23 (liquid-water-cloud.lut, ice-cloud.lut:
#                    SourceVault _v02), pinned below
#   DISORT           the original single-precision Intel kernel DISORTtoIDL.so, NSTR 60
# A product is made with
#
#    ORAC_LUT_RELEASE_MANIFEST=lut_v24.release ORAC_LUT_DISORT_PRECISION=single \\
#    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>_v24.run
#
# which refuses to run unless every file below has exactly this content.  The
# git-ignored Baum table and the external run time are listed in
# releases/lut_v24.dependencies once the products exist.  This file is sealed:
# once a product has been made from it, it is never edited.
#
# Columns: role  logical-name  file-version  sha256  path

manifest_schema = oraclut-release-manifest-1
lut_release = 24
disort_precision = single
disort_sp_kernel = intel_dlm
scientific_code = lut-code-v24.1
scattering = v24-scattering-refinements
grid = Grid-B
library_checkout = {ROOT}
created = {datetime.date.today().isoformat()}
git_head_at_creation = {head}

"""
   OUT.write_text(header + "\n".join(rows) + "\n")
   print(f"wrote {OUT} ({len(rows)} entries)")


if __name__ == "__main__":
   main()
