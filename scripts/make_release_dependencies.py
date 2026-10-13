"""Write releases/<release>.dependencies from the provenance records of the release's products.

   python scripts/make_release_dependencies.py <release manifest> <output> <product>.nc ...

The dependencies file lists what a release's products depended on beyond the
manifest itself, as recorded by the release launcher in each product's
provenance (the provenance_json attribute of each NetCDF product:
inputs_read, shared_objects_mapped), in the manifest columns
(role  logical-name  version  sha256  path), so that
scripts/check_release_manifest.py can verify it:

   input      every input file read (each path must have one digest in all products);
   external   every shared object mapped from outside the repository (Intel
              run time, system glibc, the conda environment).

The format follows releases/lut_v22_single-r2.dependencies.  The script refuses
to overwrite an existing file.
"""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main(argv):
   if len(argv) < 3:
      sys.exit(__doc__)
   manifest, output, records = Path(argv[0]).resolve(), Path(argv[1]).resolve(), [Path(p) for p in argv[2:]]
   if output.exists():
      sys.exit(f"refusing to overwrite {output}")
   inputs, external, jobs = {}, {}, []
   from netCDF4 import Dataset
   for record in records:
      with Dataset(record) as product:
         provenance = json.loads(product.getncattr("provenance_json"))
      if provenance["release_manifest"]["path"] != str(manifest.relative_to(ROOT)):
         sys.exit(f"{record} was made by {provenance['release_manifest']['path']}, not {manifest}")
      jobs.append(str(provenance["runtime"]["environment"].get("SLURM_JOB_ID", "interactive")))
      for item in provenance["inputs_read"]:
         if inputs.setdefault(item["path"], item["sha256"]) != item["sha256"]:
            sys.exit(f"{item['path']} has different digests in the products")
      for path, digest in provenance["shared_objects_mapped"].items():
         if not Path(path).resolve().is_relative_to(ROOT) and external.setdefault(path, digest) != digest:
            sys.exit(f"{path} has different digests in the products")
   lines = [f"# Dependencies of the products of {manifest.relative_to(ROOT)},",
            "# extracted from the provenance embedded in the products (inputs read and",
            f"# shared objects mapped by the generating processes, jobs {', '.join(sorted(jobs))}).",
            "# Same columns as the release manifests: role  logical-name  version  sha256  path;",
            f"# verify with: python scripts/check_release_manifest.py {output.relative_to(ROOT)}",
            "#",
            "# input        every input file read (all in Git except the third-party Baum table,",
            "#              which is git-ignored and must be kept in place)",
            "# external     shared objects outside the repository: Intel 2022.1.0 run time,",
            "#              system glibc and the conda environment 'science'",
            ""]
   for path in sorted(inputs):
      lines.append(f"{'input':<8s} {path:<90s} {'pinned':<7s} {inputs[path]}  {path}")
   for path in sorted(external):
      lines.append(f"{'external':<8s} {Path(path).name:<90s} {'pinned':<7s} {external[path]}  {path}")
   output.write_text("\n".join(lines) + "\n")
   print(f"wrote {output}: {len(inputs)} inputs, {len(external)} external shared objects")


if __name__ == "__main__":
   main(sys.argv[1:])
