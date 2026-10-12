"""V23 pre-production smoke test and benchmark (compact, through the release launchers).

Runs liquid-water p240 and water-ice agg for Terra MODIS channels 1 and 20 at
two optical depths (1, 8) and two effective radii (5, 9 um), nodes common to the
V22 coarse grids and to Grid B:

* V23: releases/lut_v23.release, with the full Grid B angles
  (21 solar zenith, 21 satellite zenith, 37 relative azimuth);
* V22: releases/lut_v22_single-r2.release, with the V22 coarse angles.

The compact grid definitions are written to a private input directory (the
pinned working definitions are not touched; the launcher still verifies them).
Reports the wall time per (channel, optical depth, radius) block of each run,
the dimensions and release metadata of the V23 products, and the V23/V22
differences at the exactly shared nodes (solar zenith 0/30/60, satellite
zenith 0, relative azimuth 0/90/180).

Usage: python validation/v23_terra_modis/smoke.py   -> validation/v23_terra_modis/smoke/ and smoke_results.json
"""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "smoke"
PYTHON = sys.executable
GRIDB_ANGLES = """21 uneven_linear # solar zenith
0.0 3.75 7.5 11.25 15.0 18.75 22.5 26.25 30.0 33.75 37.5 41.25 45.0 48.75 52.5 56.25 60.0 63.75 67.5 71.25 75.0
21 uneven_linear # satellite zenith
0.0 3.75 7.5 11.25 15.0 18.75 22.5 26.25 30.0 33.75 37.5 41.25 45.0 48.75 52.5 56.25 60.0 63.75 67.5 71.25 75.0
37 linear # relative azimuth
0 180
"""
V22_ANGLES = """10 linear # solar zenith
0 90
10 linear # satellite zenith
0 90
11 linear # relative azimuth
0 180
"""
COMMON = "2 uneven_logarithmic # optical depth\n1.0 8.0\n2 uneven_linear # effective radius\n5.0 9.0\n"
PRODUCTS = {"p240": ("liquid-water_240.mm", "liquid-water-cloud.lut", "terra_modis_m_liquid-water_a01_p240_v{v}.nc"),
            "pagg": ("water-ice_agg.mm", "ice-cloud.lut", "terra_modis_m_water-ice_a01_pagg_v{v}.nc")}
RELEASES = {23: ("lut_v23.release", GRIDB_ANGLES), 22: ("lut_v22_single-r2.release", V22_ANGLES)}


def stage(version):
   inputs = WORK / f"inputs_v{version}"
   (inputs / "lut").mkdir(parents=True)
   for entry in (ROOT / "create_orac_lut" / "input_files").iterdir():
      if entry.name != "lut":
         (inputs / entry.name).symlink_to(entry)
   for _, lutfile, _ in PRODUCTS.values():
      (inputs / "lut" / lutfile).write_text("# V23 smoke-test grid (compact)\n" + COMMON + RELEASES[version][1])
   return inputs


def run(version, product):
   mmfile, lutfile, _ = PRODUCTS[product]
   inputs = WORK / f"inputs_v{version}"
   out = WORK / f"out_v{version}"
   run_file = WORK / f"{product}_v{version}.run"
   run_file.write_text(f"""platform = 'terra'
instrument = 'modis'
forward_model = 'cloud'
in_path = '{inputs.relative_to(ROOT)}'
instfile = 'terra_modis_v1.inst'
mmfile = '{mmfile}'
lutfile = '{lutfile}'
atmospheres = 2
channelid = [1, 20]
gas = 0
no_rayleigh = 0
srf_quad = 1
nstreams = 60
nmom = 1000
version = {version}
out_path = '{out.relative_to(ROOT)}'
""")
   env = {"HOME": os.environ["HOME"], "PATH": "/usr/bin:/bin", "ORAC_LUT_RELEASE_MANIFEST": RELEASES[version][0],
          "ORAC_LUT_DISORT_PRECISION": "single", "ORAC_LUT_DISORT_SP_KERNEL": "intel_dlm"}
   start = time.time()
   result = subprocess.run([PYTHON, "create_orac_luts_release_v04.py", str(run_file.relative_to(ROOT))], cwd=ROOT, env=env,
                           capture_output=True, text=True)
   elapsed = time.time() - start
   (WORK / f"{product}_v{version}.log").write_text(result.stdout + result.stderr)
   if result.returncode:
      sys.exit(f"{product} V{version} failed:\n{result.stderr[-2000:]}")
   return elapsed, out / PRODUCTS[product][2].format(v=version)


def shared(reference, candidate):
   rows = {}
   with Dataset(reference) as r, Dataset(candidate) as c:
      r.set_auto_maskandscale(False)
      c.set_auto_maskandscale(False)
      for name in ("R_0v", "R_0d", "T_0d", "T_00", "R_dv", "T_dv", "R_dd", "T_dd", "E_md"):
         if name not in r.variables:
            continue
         a, b = r[name], c[name]
         index_a, index_b = [], []
         for dim in a.dimensions:
            x = np.asarray(r[dim][:]) if dim in r.variables else np.arange(len(r.dimensions[dim]))
            y = np.asarray(c[dim][:]) if dim in c.variables else np.arange(len(c.dimensions[dim]))
            common = sorted(set(np.round(x.astype(float), 4)) & set(np.round(y.astype(float), 4)))
            index_a.append([int(np.argmin(abs(x - v))) for v in common])
            index_b.append([int(np.argmin(abs(y - v))) for v in common])
         va = np.asarray(a[:])[np.ix_(*index_a)]
         vb = np.asarray(b[:])[np.ix_(*index_b)]
         differing = va.view(np.uint32) != vb.view(np.uint32)
         rows[name] = {"shared_elements": int(va.size), "differing": int(differing.sum()),
                       "max_abs": float(np.max(np.abs(va.astype(float) - vb))) if va.size else 0.0}
   return rows


def main():
   if WORK.exists():
      shutil.rmtree(WORK)
   WORK.mkdir(parents=True)
   for version in RELEASES:
      stage(version)
   results = {}
   for product in PRODUCTS:
      t23, p23 = run(23, product)
      t22, p22 = run(22, product)
      with Dataset(p23) as d:
         dims = {k: len(v) for k, v in d.dimensions.items()}
         meta = {"lut_release": d.getncattr("lut_release"), "manifest": d.getncattr("provenance_release_manifest")[:60],
                 "nan_or_fill": {n: int(np.count_nonzero(~np.isfinite(np.asarray(d[n][:], float)))) for n in d.variables
                                 if np.issubdtype(d[n].dtype, np.floating)}}
      blocks = 2 * 2 * 2
      results[product] = {"seconds_v23": t23, "seconds_v22": t22, "seconds_per_block_v23": t23 / blocks,
                          "seconds_per_block_v22": t22 / blocks, "dimensions_v23": dims, "metadata_v23": meta,
                          "shared_nodes": shared(p22, p23)}
      print(product, json.dumps({k: v for k, v in results[product].items() if k != "metadata_v23"}), flush=True)
   (HERE / "smoke_results.json").write_text(json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
   main()
