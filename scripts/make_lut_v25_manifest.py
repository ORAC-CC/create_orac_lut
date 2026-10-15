"""Write releases/lut_v25.release: LUT release V25 = the V24 release with vertically varying cloud temperature.

   python scripts/make_lut_v25_manifest.py

V25 (V25_RELEASE.md) changes the scientific code of V24
(releases/lut_v24.release) in exactly three files, and adds one input:

* the driver create_orac_luts (v05 -> v07: run-file setting
  cloud_vertical_profile, the profile set up once per LUT, the emission call
  on the sub-layers of cloud_temperature.emission_layers, E_md valid range,
  provenance text; v06 was the 2026-10-15 candidate);
* oraclut.idl_mirror.write_v2_lut (v03 -> v04: optional E_md valid range);
* oraclut.idl_mirror.cloud_temperature (new, v02: the profiles and the
  emission layering; v01 was the candidate's 40 equal sub-layers);
* the input references/data/ocalut_cloudprofile_Cirrostratus.dat (pinned).

Every other entry of the V24 manifest (modules, adapters, DISORT and Mie
kernel sources, the Intel DISORT library and run time, and the pinned input
files: Grid B liquid-water-cloud.lut / ice-cloud.lut, instrument, microphysics,
refractive indices, atmosphere, SRFs, solar spectrum) is copied unchanged,
verified against the files on disk.  CODE_VERSION is pinned with its V25
content (code_release lut-code-v25.2, lut_version 25).

Products are made with launcher v04:
   ORAC_LUT_RELEASE_MANIFEST=lut_v25.release ORAC_LUT_DISORT_PRECISION=single \\
   ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>_v25.run

The script refuses to overwrite an existing manifest (sealed once written).
"""

import datetime
import hashlib
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
V24 = ROOT / "releases" / "lut_v24.release"
OUT = ROOT / "releases" / "lut_v25.release"
# logical name -> (role, file version, path) of the V25 files
CHANGED = {
   "create_orac_luts": ("driver", "v07", "create_orac_luts_v07.py"),
   "oraclut.idl_mirror.write_v2_lut": ("module", "v04", "src/oraclut/idl_mirror/write_v2_lut_v04.py"),
}
ADDED_AFTER = {   # new entry -> the V24 entry it follows
   "oraclut.idl_mirror.cloud_temperature": ("oraclut.idl_mirror.call_disort", "module", "v02",
                                            "src/oraclut/idl_mirror/cloud_temperature_v02.py"),
   "input-cloudprofile/ocalut_cloudprofile_Cirrostratus.dat": ("input-bbconstants.pro", "kernel", "pinned",
                                                               "references/data/ocalut_cloudprofile_Cirrostratus.dat"),
}
# The working files must equal these SourceVault copies (no unversioned change can enter the release)
WORKING = {"create_orac_luts_v07.py": "create_orac_luts.py",
           "src/oraclut/idl_mirror/write_v2_lut_v04.py": "src/oraclut/idl_mirror/write_v2_lut.py",
           "src/oraclut/idl_mirror/cloud_temperature_v02.py": "src/oraclut/idl_mirror/cloud_temperature.py",
           "CODE_VERSION_v08": "CODE_VERSION"}


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
   rows, replaced, placed = [], set(), set()
   for line in V24.read_text().splitlines():
      fields = line.split("#", 1)[0].split()
      if len(fields) != 5 or fields[0] == "config":
         continue
      role, logical, version, digest, path = fields
      full = Path(path) if Path(path).is_absolute() else ROOT / path
      if sha256(full) != digest:
         sys.exit(f"{path} no longer matches the V24 manifest")
      if logical in CHANGED:
         new_role, new_version, new_path = CHANGED[logical]
         if new_role != role:
            sys.exit(f"{logical}: role {role} in V24, {new_role} expected")
         rows.append(row(role, logical, new_version, sha256(ROOT / new_path), new_path))
         replaced.add(logical)
      else:
         rows.append(row(role, logical, version, digest, path))
      for name, (after, new_role, new_version, new_path) in ADDED_AFTER.items():
         if after == logical:
            rows.append(row(new_role, name, new_version, sha256(ROOT / new_path), new_path))
            placed.add(name)
   if replaced != set(CHANGED) or placed != set(ADDED_AFTER):
      sys.exit(f"V24 manifest lacks {sorted(set(CHANGED) - replaced)} or the anchors of {sorted(set(ADDED_AFTER) - placed)}")
   rows.append(row("config", "CODE_VERSION", "pinned", sha256(ROOT / "CODE_VERSION"), "CODE_VERSION"))
   head = subprocess.run(["/usr/bin/git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
   header = f"""\
# ORAC LUT release manifest
#
# LUT release V25: the V24 release with a vertically varying cloud
# temperature in the thermal-emission calculation (V25_RELEASE.md).
#   LUT release      V25
#   scientific code  lut-code-v25.2 = the V24 code (releases/lut_v24.release)
#                    with the driver create_orac_luts_v07, write_v2_lut_v04 and
#                    the new cloud_temperature_v02; every other entry is copied
#                    from releases/lut_v24.release
#   cloud emission   cloud_vertical_profile: water ice 'cirrostratus' (the
#                    supplied profile, pinned below), liquid water
#                    'wet_adiabat'; cloud top 250 K, E_md / B(250 K); geometric
#                    emission sub-layers with a projected Planck source, or
#                    isothermal cells for optically thin clouds
#   scattering       the V24 refinements, unchanged (Mie: 3.5 r_e, nested
#                    radius refinement, adaptive Legendre, DISORT
#                    NMOM = max(1000, L); Baum and T-matrix as V23)
#   grid             Grid B (liquid-water-cloud.lut, ice-cloud.lut: SourceVault
#                    _v02), pinned below
#   DISORT           the original single-precision Intel kernel DISORTtoIDL.so, NSTR 60
# A product is made with
#
#    ORAC_LUT_RELEASE_MANIFEST=lut_v25.release ORAC_LUT_DISORT_PRECISION=single \\
#    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>_v25.run
#
# which refuses to run unless every file below has exactly this content.  The
# git-ignored Baum table and the external run time are listed in
# releases/lut_v25.dependencies once the products exist.  This file is sealed:
# once a product has been made from it, it is never edited.
#
# Columns: role  logical-name  file-version  sha256  path

manifest_schema = oraclut-release-manifest-1
lut_release = 25
disort_precision = single
disort_sp_kernel = intel_dlm
scientific_code = lut-code-v25.2
scattering = v24-scattering-refinements
cloud_emission = v25-cloud-temperature-profiles
grid = Grid-B
library_checkout = {ROOT}
created = {datetime.date.today().isoformat()}
git_head_at_creation = {head}

"""
   OUT.write_text(header + "\n".join(rows) + "\n")
   print(f"wrote {OUT} ({len(rows)} entries)")


if __name__ == "__main__":
   main()
