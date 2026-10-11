"""Report and verify the provenance embedded in an ORAC LUT (no Git needed).

   python validation/source_versioning/inspect_lut_provenance.py <lut.nc> [--source-root DIR] [--verify-inputs] [--json]

Answers "exactly which source files, at which file versions, generated this
LUT?" from the product's NetCDF global attributes (written by
create_orac_luts_release_vNN.py), and checks each listed file's SHA-256
against a source tree (default: this repository; any copy or archive
extraction of it works):

   MATCH      the file exists with the recorded SHA-256
   MISMATCH   the file exists with different content
   MISSING    the file does not exist under the source root

The kernel libraries are build products; a MISMATCH there means the library
in the source tree is not the one that ran (for example a different checkout
path), not that the sources differ.  --verify-inputs checks the input files
the calculation read in the same way.

Exit status: 0 every source file (and, if requested, input) matches;
1 something does not match or is missing; 3 the product carries no release
provenance (made before source-file versioning, or outside the release
launcher), in which case nothing is inferred.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[2]


def sha256_file(path):
   digest = hashlib.sha256()
   with open(path, "rb") as handle:
      for block in iter(lambda: handle.read(1 << 20), b""):
         digest.update(block)
   return digest.hexdigest()


def check(root, relative_path, sha256):
   path = Path(root) / relative_path
   if not path.is_file():
      return "MISSING"
   return "MATCH" if sha256_file(path) == sha256 else "MISMATCH"


def read_provenance(lut_file):
   with Dataset(lut_file, "r") as dataset:
      attributes = {name: dataset.getncattr(name) for name in dataset.ncattrs()}
   if "lut_release" not in attributes or "provenance_json" not in attributes:
      return attributes, None
   return attributes, json.loads(attributes["provenance_json"])


def inspect(lut_file, source_root, verify_inputs):
   attributes, record = read_provenance(lut_file)
   if record is None:
      print(f"{lut_file}: no embedded release provenance (a product made before source-file versioning, or outside the "
            "release launcher); its source cannot be identified from the file")
      return 3, None
   source_root = Path(source_root)
   report = {"lut_file": str(lut_file), "lut_release": attributes["lut_release"], "source_root": str(source_root),
             "sources": [], "inputs": [], "problems": 0}

   print(f"LUT file:          {lut_file}")
   print(f"LUT release:       {attributes['lut_release']}   (provenance schema {record['provenance_schema']})")
   print(f"Generated (UTC):   {record['generated_utc']}  (started {record['started_utc']})")
   print(f"Driver:            {record['driver']['path']}  file version {record['driver']['file_version']}")
   print(f"Release manifest:  {record['release_manifest']['path']}  sha256 {record['release_manifest']['sha256']}")
   print(f"Launcher:          {record['launcher']['path']}  sha256 {record['launcher']['sha256']}")
   print(f"Run file:          {record['run_file']['path']}  sha256 {record['run_file']['sha256']}")
   print(f"Command:           {record['generation_command']}")
   git = record["git_at_start"]
   if git.get("available"):
      outside = git.get("release_files_not_in_commit", [])
      print(f"Git commit:        {git['commit']}  ({git['working_tree']}; "
            + ("every release file in this commit" if not outside else f"{len(outside)} release file(s) not in this commit")
            + ")")
   else:
      print(f"Git commit:        unavailable ({git.get('detail')})")
   print(f"Repository:        {record['repository']['public']}")
   print()
   print(f"Source files (checked under {source_root}):")
   rows = [("launcher", "-", record["launcher"]["path"], record["launcher"]["sha256"]),
           ("manifest", "-", record["release_manifest"]["path"], record["release_manifest"]["sha256"])]
   rows += [(item["role"], item["file_version"], item["path"], item["sha256"]) for item in record["source_files"]]
   # A sealed CODE_VERSION served from its SourceVault copy by the launcher
   # (create_orac_luts_release_v03.py) is checked at the copy that was read.
   served = record.get("config_substitutions", {})
   for role, file_version, path, sha256 in rows:
      status = check(source_root, served.get(path, path), sha256)
      if path in served:
         path = f"{path} (served from {served[path]})"
      if status != "MATCH" and role != "library":
         report["problems"] += 1
      report["sources"].append({"role": role, "file_version": file_version, "path": path, "sha256": sha256, "status": status})
      note = "   (build product)" if role == "library" and status != "MATCH" else ""
      print(f"   {status:<8}  {role:<8} {file_version:<7} {sha256}  {path}{note}")

   manifest_path = source_root / record["release_manifest"]["path"]
   if manifest_path.is_file() and sha256_file(manifest_path) == record["release_manifest"]["sha256"]:
      listed = []
      for raw in manifest_path.read_text().splitlines():
         fields = raw.split("#", 1)[0].split()
         if len(fields) == 5:
            listed.append((fields[0], fields[1], fields[2], fields[3], fields[4]))
      embedded = [(i["role"], i["logical"], i["file_version"], i["sha256"], i["path"]) for i in record["source_files"]]
      if listed != embedded:
         report["problems"] += 1
         print("   PROBLEM   the embedded source list differs from the release manifest it names")
      else:
         print(f"   the embedded source list equals {record['release_manifest']['path']} ({len(listed)} entries)")

   print()
   print(f"Inputs read: {len(record['inputs_read'])} files" + ("" if verify_inputs else " (not verified; use --verify-inputs)"))
   for item in record["inputs_read"]:
      status = check(source_root, item["path"], item["sha256"]) if verify_inputs else "-"
      if verify_inputs and status == "MISSING" and Path(item["realpath"]).is_file():
         status = "MATCH" if sha256_file(item["realpath"]) == item["sha256"] else "MISMATCH"
      if verify_inputs and status != "MATCH":
         report["problems"] += 1
      report["inputs"].append({**item, "status": status})
      print(f"   {status:<8}  {item['sha256']}  {item['path']}  ({item['via']})")
   print()
   if report["problems"]:
      print(f"RESULT: {report['problems']} file(s) do not match the provenance recorded in the LUT")
      return 1, report
   print("RESULT: every source file" + (" and input" if verify_inputs else "") + " matches the provenance recorded in the LUT")
   return 0, report


def main(argv=None):
   parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
   parser.add_argument("lut_file", type=Path)
   parser.add_argument("--source-root", type=Path, default=ROOT, help="tree holding the source files (default: this repository)")
   parser.add_argument("--verify-inputs", action="store_true", help="also verify the input files the calculation read")
   parser.add_argument("--json", type=Path, help="write the verification result as JSON")
   args = parser.parse_args(argv)
   status, report = inspect(args.lut_file, args.source_root, args.verify_inputs)
   if args.json and report is not None:
      args.json.write_text(json.dumps(report, indent=1))
   return status


if __name__ == "__main__":
   sys.exit(main())
