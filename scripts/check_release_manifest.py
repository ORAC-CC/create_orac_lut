"""Verify every file of an ORAC release manifest (Python or IDL) against its SHA-256.

   python scripts/check_release_manifest.py releases/idl_lut_v21_single.release [more manifests]
   python scripts/check_release_manifest.py --products releases/idl_lut_v21_double.products

Reads the entry lines of releases/*.release (oraclut-release-manifest-1 or
oraclut-idl-release-manifest-1) or a .products file, hashes each listed path
(relative to the repository, or absolute) and reports any missing or changed
file.  Exit status 0 only if every file matches.  Nothing is modified.
"""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
   digest = hashlib.sha256()
   with open(path, "rb") as handle:
      for block in iter(lambda: handle.read(1 << 20), b""):
         digest.update(block)
   return digest.hexdigest()


def entries(manifest):
   for line in Path(manifest).read_text().splitlines():
      fields = line.split()
      if line.lstrip().startswith("#") or len(fields) < 5:
         continue
      if len(fields[3]) != 64 or any(c not in "0123456789abcdef" for c in fields[3]):
         continue
      yield fields[0], fields[1], fields[3], " ".join(fields[4:])


def check(manifest):
   bad = 0
   count = 0
   for role, logical, digest, path in entries(manifest):
      count += 1
      full = Path(path) if Path(path).is_absolute() else ROOT / path
      if not full.is_file():
         print(f"MISSING  {role:<11s} {logical}  {path}")
         bad += 1
      elif sha256(full) != digest:
         print(f"CHANGED  {role:<11s} {logical}  {path}")
         bad += 1
   print(f"{manifest}: {count} entries, {bad} problem(s)")
   return bad == 0 and count > 0


def main(argv):
   paths = [a for a in argv if a != "--products"]
   if not paths:
      sys.exit(__doc__)
   ok = all([check(p) for p in paths])
   return 0 if ok else 1


if __name__ == "__main__":
   raise SystemExit(main(sys.argv[1:]))
