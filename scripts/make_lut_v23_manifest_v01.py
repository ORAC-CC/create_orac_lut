"""Write releases/lut_v23.release: LUT release V23 = the lut-code-v22.4 scientific code with the Grid B definitions.

   python scripts/make_lut_v23_manifest.py

V23 is a grid-only successor of the definitive V22 (releases/lut_v22_single-r2.release):

* every code, adapter, kernel-source, library and Intel run-time entry of the
  V22 manifest is copied unchanged (same file versions and SHA-256, verified
  against the files on disk), so the release runs exactly the V22 scientific
  code (driver create_orac_luts_v04.py, legacy_disort_v04, idl_math_v02, ...);
* CODE_VERSION is pinned with its V23 content (code_release lut-code-v22.4,
  lut_version 23);
* every Git-tracked input file read by the V22 products is pinned (role
  "kernel": the launcher verifies it by SHA-256 at the start and end of a run
  and refuses any difference), with the LUT definitions liquid-water-cloud.lut
  and ice-cloud.lut holding Grid B (SourceVault _v02).  The git-ignored Baum
  table is recorded in releases/lut_v23.dependencies instead.

Products are made with launcher v04:
   ORAC_LUT_RELEASE_MANIFEST=lut_v23.release ORAC_LUT_DISORT_PRECISION=single \\
   ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>_v23.run

The script refuses to overwrite an existing manifest (sealed once written).
"""

import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
V22 = ROOT / "releases" / "lut_v22_single-r2.release"
OUT = ROOT / "releases" / "lut_v23.release"
V22_PROVENANCE = ROOT / "releases" / "products" / "lut_v22"
GRIDS = {"create_orac_lut/input_files/lut/liquid-water-cloud.lut": "liquid-water-cloud_v02.lut",
         "create_orac_lut/input_files/lut/ice-cloud.lut": "ice-cloud_v02.lut"}


def sha256(path):
   return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
   if OUT.exists():
      sys.exit(f"refusing to overwrite sealed manifest {OUT}")
   tracked = set(subprocess.run(["/usr/bin/git", "-C", str(ROOT), "ls-files"], capture_output=True, text=True).stdout.split("\n"))
   rows = []
   for line in V22.read_text().splitlines():
      fields = line.split("#", 1)[0].split()
      if len(fields) != 5 or fields[0] == "config":
         continue
      role, logical, version, digest, path = fields
      full = Path(path) if Path(path).is_absolute() else ROOT / path
      if sha256(full) != digest:
         sys.exit(f"{path} no longer matches the V22 manifest")
      rows.append(f"{role:<8s} {logical:<52s} {version:<7s} {digest}  {path}")
   rows.append(f"{'config':<8s} {'CODE_VERSION':<52s} {'pinned':<7s} {sha256(ROOT / 'CODE_VERSION')}  CODE_VERSION")
   inputs = set()
   for record in sorted(V22_PROVENANCE.glob("*.provenance_embedded.json")):
      for item in json.loads(record.read_text())["inputs_read"]:
         inputs.add(item["path"])
   for path in sorted(inputs):
      if path == "CODE_VERSION" or path.startswith("runs/"):
         continue
      if path not in tracked:
         continue                                      # the git-ignored Baum table: lut_v23.dependencies
      note = GRIDS.get(path)
      if note and sha256(ROOT / path) != sha256(ROOT / Path(path).parent / note):
         sys.exit(f"{path} does not hold the Grid B content of {note}")
      rows.append(f"{'kernel':<8s} {'input-' + path.split('input_files/', 1)[1]:<52s} {'pinned':<7s} {sha256(ROOT / path)}  {path}")
   head = subprocess.run(["/usr/bin/git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
   header = f"""\
# ORAC LUT release manifest
#
# LUT release V23: a grid-only successor of the definitive V22.
#   LUT release      V23
#   scientific code  lut-code-v22.4, unchanged: every code, library and run-time
#                    entry below is copied from releases/lut_v22_single-r2.release
#   grid             Grid B (commit 4bb0ac7, validation/v23 value-and-Jacobian
#                    study), installed as liquid-water-cloud.lut and ice-cloud.lut
#                    (SourceVault liquid-water-cloud_v02.lut, ice-cloud_v02.lut)
# A product is made with
#
#    ORAC_LUT_RELEASE_MANIFEST=lut_v23.release ORAC_LUT_DISORT_PRECISION=single \\
#    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/<product>_v23.run
#
# which refuses to run unless every file below has exactly this content.  The
# "kernel" entries named input-* are the input files (LUT definitions,
# instrument, microphysics, refractive indices, atmosphere, SRFs, solar
# spectrum): pinned, so a run can never silently use another grid.  The
# git-ignored Baum table is listed in releases/lut_v23.dependencies.  This file
# is sealed: once a product has been made from it, it is never edited.
#
# Columns: role  logical-name  file-version  sha256  path

manifest_schema = oraclut-release-manifest-1
lut_release = 23
disort_precision = single
disort_sp_kernel = intel_dlm
scientific_code = lut-code-v22.4
grid = Grid-B
library_checkout = {ROOT}
created = {datetime.date.today().isoformat()}
git_head_at_creation = {head}

"""
   OUT.write_text(header + "\n".join(rows) + "\n")
   print(f"wrote {OUT} ({len(rows)} entries)")


if __name__ == "__main__":
   main()
