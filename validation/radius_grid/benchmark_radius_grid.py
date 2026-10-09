"""Cost of the refined radius integration (V24 radius-grid benchmark).

Validation-only code.  Times generate_scattering_properties (single thread:
OMP/OPENBLAS/MKL_NUM_THREADS = 1, one process per variant) on the V23 Grid B
grids at representative channels, for

   v23        source revision 9d663e9 (the V23 products: 0.001-100 um, nmom 1000)
   dbcc42c    source revision dbcc42c (3.5 r_e lattice limit, adaptive Legendre,
              legacy radius spacing)
   adopted    the current source (refined radius grid)
   reference  the current source with every modified-gamma grid forced to
              refinement level 6 (the convergence reference)

Each revision runs from an export made with git archive under
validation/tmp/radius_grid/rev_<hash> (links to build/, mie/, create_orac_lut/).
The Mie kernel is wrapped to count radius nodes and node x angle
evaluations and to time the kernel separately from the rest.

   python validation/radius_grid/benchmark_radius_grid.py <case> <variant>
   python validation/radius_grid/benchmark_radius_grid.py summary

Cases: liquid (liquid-water_stg, liquid-water Grid B, 24 radii) and ice_sph
(water-ice_sph, ice Grid B, 29 radii), each at Terra MODIS channels
3, 1, 6, 7, 20, 31 (0.469-11.03 um) in one call including the 0.55 um reference.
"""

import csv
import importlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
   os.environ[variable] = "1"

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TMP = ROOT / "validation" / "tmp" / "radius_grid"
RESULTS = HERE / "results" / "benchmark"
CHANNELS = [3, 1, 6, 7, 20, 31]
CASES = {"liquid": ("liquid-water_stg.mm", "liquid-water-cloud-grid-b.lut"),
         "ice_sph": ("water-ice_sph.mm", "ice-cloud-grid-b.lut")}
REVISIONS = {"v23": "9d663e9", "dbcc42c": "dbcc42c"}
VARIANTS = ["v23", "dbcc42c", "adopted", "reference"]


def export(revision):
   tree = TMP / f"rev_{revision}"
   if not (tree / "create_orac_luts.py").exists():
      tree.mkdir(parents=True, exist_ok=True)
      archive = subprocess.run(["git", "-C", str(ROOT), "archive", revision, "create_orac_luts.py", "src"],
                               check=True, capture_output=True).stdout
      subprocess.run(["tar", "-x", "-C", str(tree)], input=archive, check=True)
      for name in ("build", "create_orac_lut", "mie"):
         (tree / name).symlink_to(ROOT / name)
   return tree


def cpu_model():
   try:
      for line in Path("/proc/cpuinfo").read_text().splitlines():
         if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
   except OSError:
      pass
   return platform.processor()


def run(case, variant):
   root = export(REVISIONS[variant]) if variant in REVISIONS else ROOT
   sys.path.insert(0, str(root / "src"))
   from oraclut.idl_mirror import generate_scattering_properties, load_inststr, load_lutstr, load_mmdat, load_srfstrarr

   bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
   if variant == "reference":
      bwgp.radius_refinement_level = lambda *args, **kwargs: 6
   counts = {"mie_calls": 0, "radius_nodes": 0, "node_angle_evaluations": 0, "mie_cpu_s": 0.0}
   kernel = bwgp.mie_single_batch

   def counted(size_parameters, refractive_index, scattering_cosines):
      started = time.process_time()
      result = kernel(size_parameters, refractive_index, scattering_cosines)
      counts["mie_cpu_s"] += time.process_time() - started
      counts["mie_calls"] += 1
      counts["radius_nodes"] += int(size_parameters.size)
      counts["node_angle_evaluations"] += int(size_parameters.size) * int(scattering_cosines.size)
      return result

   bwgp.mie_single_batch = counted
   mmfile, grid = CASES[case]
   inputs = ROOT / "create_orac_lut" / "input_files"
   inststr = load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=CHANNELS)
   lutstr = load_lutstr(inputs / "lut" / grid, inststr.max_sat_zenith)
   srfstrarr, nwvl_max = load_srfstrarr(inststr, inputs / "sun" / "Gueymard2018.sssi", 1, inputs)
   mmstr = load_mmdat(inputs / "microphysics" / mmfile, inputs)
   started_wall, started_cpu = time.time(), time.process_time()
   generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, 1000 if variant == "v23" else None)
   row = {"case": case, "variant": variant, "revision": REVISIONS.get(variant, "working tree"),
          "host": socket.gethostname(), "cpu": cpu_model(), "threads": 1,
          "wall_s": time.time() - started_wall, "cpu_s": time.process_time() - started_cpu, **counts}
   row["other_cpu_s"] = row["cpu_s"] - row["mie_cpu_s"]
   RESULTS.mkdir(parents=True, exist_ok=True)
   (RESULTS / f"{case}_{variant}.json").write_text(json.dumps(row, indent=2) + "\n")
   print(row, flush=True)


def summary():
   rows = [json.loads(path.read_text()) for path in sorted(RESULTS.glob("*.json"))]
   with open(RESULTS / "benchmark_radius_grid.csv", "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)
   for row in rows:
      print(f"{row['case']:8s} {row['variant']:10s} wall {row['wall_s']:8.0f} s  cpu {row['cpu_s']:8.0f} s  "
            f"Mie {row['mie_cpu_s']:8.0f} s  nodes {row['radius_nodes']:>11d}  node x angle {row['node_angle_evaluations']:.3e}")


def main():
   if sys.argv[1] == "summary":
      summary()
   else:
      run(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
   main()
