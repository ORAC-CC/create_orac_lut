"""Compact LUT comparison of the adaptive Legendre expansion (Phase 3).

Validation-only code.  Runs the production generator on compact grids and
compares three Legendre treatments of the same optics:

  adaptive   current source: Gauss-Legendre order above each wavelength's
             degree bound, King-criterion expansion length L per phase function
  long       current source with twice the quadrature order and every exact
             coefficient up to each phase function's degree bound D_r kept
             (a deliberately long reference expansion, testing the stopping
             criterion itself)
  fixed      the fixed nmom = 1000 expansion of source revision c6ad545 (the
             integration-limit change without the Legendre change), run from
             an export of that revision under validation/tmp/rev_c6ad545

Cases: MODIS and dual-view SLSTR liquid-water cloud on
validation/size_distribution_limits/compact_liquid_grid.lut, and SEVIRI
aerosol classes a79 (two-mode), a76 (the class of the production PA76 Grid B
aerosol LUTs) and a78 (three Mie modes) on compact_aerosol_grid.lut (a75 has
a T-matrix dust component and therefore keeps the fixed expansion).  For MODIS/SLSTR "fixed" reuses the products of
validation/size_distribution_limits/compact_lut_comparison.py "adopted",
which were generated from the same c6ad545 source.

   python compact_lut_legendre.py generate <case> <adaptive|long|fixed>
   python compact_lut_legendre.py compare
"""

import argparse
import csv
import importlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
MIRROR = ROOT / "validation" / "tmp" / "compact_inputs"
PRODUCTS = ROOT / "validation" / "tmp" / "compact_luts_legendre"
OLD_TREE = ROOT / "validation" / "tmp" / "rev_c6ad545"
PHASE2_PRODUCTS = ROOT / "validation" / "tmp" / "compact_luts"
RESULTS = HERE / "results" / "compact_lut"
CASES = {
   # name: (forward model, instrument file, microphysics, grid, channels, gas)
   "modis": ("cloud", "terra_modis_v1.inst", "liquid-water_stg.mm", "compact_liquid_grid.lut", [1, 2, 6, 7, 20, 31], 0),
   "slstr": ("cloud", "sentinel-3a_slstr_v1.inst", "liquid-water_stg.mm", "compact_liquid_grid.lut", [1, 3, 5, 6, 7, 8], 0),
   "seviri_a79": ("aerosol", "meteosat-10_seviri_v1.inst", "aerosol_a79.mm", "compact_aerosol_grid.lut", [1, 9], 1),
   "seviri_a76": ("aerosol", "meteosat-10_seviri_v1.inst", "aerosol_a76.mm", "compact_aerosol_grid.lut", [1, 9], 1),
   "seviri_a78": ("aerosol", "meteosat-10_seviri_v1.inst", "aerosol_a78.mm", "compact_aerosol_grid.lut", [1, 9], 1),
}


def make_mirror():
   (MIRROR / "lut").mkdir(parents=True, exist_ok=True)
   for item in INPUTS.iterdir():
      link = MIRROR / item.name
      if item.name != "lut" and not link.exists():
         link.symlink_to(item)
   for grid in (ROOT / "validation" / "size_distribution_limits" / "compact_liquid_grid.lut", HERE / "compact_aerosol_grid.lut"):
      (MIRROR / "lut" / grid.name).write_text(grid.read_text())


def make_old_tree():
   """Export source revision c6ad545 (generator and src/) with links to the kernels and inputs."""

   if (OLD_TREE / "create_orac_luts.py").exists():
      return
   OLD_TREE.mkdir(parents=True, exist_ok=True)
   archive = subprocess.run(["git", "-C", str(ROOT), "archive", "c6ad545", "create_orac_luts.py", "src"],
                            check=True, capture_output=True).stdout
   subprocess.run(["tar", "-x", "-C", str(OLD_TREE)], input=archive, check=True)
   for name in ("build", "create_orac_lut", "mie"):
      (OLD_TREE / name).symlink_to(ROOT / name)


def generate(case, variant):
   forward, instfile, mmfile, grid, channels, gas = CASES[case]
   root = OLD_TREE if variant == "fixed" else ROOT
   sys.path.insert(0, str(root))
   sys.path.insert(0, str(root / "src"))
   import create_orac_luts

   if variant == "long":
      gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
      legendre = importlib.import_module("oraclut.idl_mirror.legendre_expansion")
      original_order = gsp.initial_quadrature_order
      gsp.initial_quadrature_order = lambda degree: 2 * original_order(degree)

      def every_exact_coefficient(omega, degree, check_mu, check_phase):
         length, noise, error = legendre.expansion_length(omega, degree, check_mu, check_phase)
         length = int(degree) + 1
         reconstructed = legendre.reconstruct_phase_function(omega, length, check_mu)
         return length, noise, float(max(abs(reconstructed / check_phase - 1.0)))

      gsp.expansion_length = every_exact_coefficient

   timing = {}
   scattering = create_orac_luts.generate_scattering_properties

   def timed(*args, **kwargs):
      started = time.process_time()
      result = scattering(*args, **kwargs)
      timing["scattering_cpu_s"] = time.process_time() - started
      return result

   create_orac_luts.generate_scattering_properties = timed
   out = PRODUCTS / case / variant
   out.mkdir(parents=True, exist_ok=True)
   nmom = 1000 if variant == "fixed" else None
   started = time.process_time()
   if forward == "cloud":
      status = create_orac_luts.create_orac_cloud_lut(MIRROR, instfile, mmfile, grid, out, 2, channelid=channels, srf_quad=1,
                                                      version=None, nstreams=60, nmom=nmom, work_path=out / "work")
   else:
      status = create_orac_luts.create_orac_aerosol_lut(MIRROR, instfile, mmfile, grid, out, 2, channelid=channels, gas=gas,
                                                        srf_quad=1, version=None, nstreams=60, nmom=nmom, work_path=out / "work")
   timing["total_cpu_s"] = time.process_time() - started
   timing["status"] = status
   (out / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
   print(case, variant, timing, flush=True)


def product(case, variant):
   if variant == "fixed" and case in ("modis", "slstr"):
      return PHASE2_PRODUCTS / case / "adopted"
   return PRODUCTS / case / variant


def compare():
   import numpy as np
   sys.path.insert(0, str(ROOT / "src"))
   from oraclut.io.lut import read_lut

   RESULTS.mkdir(parents=True, exist_ok=True)
   rows, moments = [], []
   for case in CASES:
      luts = {v: read_lut(next(product(case, v).glob("*.nc"))) for v in ("adaptive", "long", "fixed")}
      caches = {v: np.load(product(case, v) / "work" / "scatfile.npz") for v in ("adaptive", "long", "fixed")}
      for name in luts["long"].variable_names:
         reference = np.asarray(luts["long"].variables[name])
         if reference.dtype.kind != "f":
            continue
         reference = reference.astype(np.float64)
         for variant in ("adaptive", "fixed"):
            other = np.asarray(luts[variant].variables[name], dtype=np.float64)
            finite = np.isfinite(reference) & np.isfinite(other) & (np.abs(reference) < 1e30)
            if not np.any(finite):
               continue
            d = np.abs(other[finite] - reference[finite])
            scale = np.abs(reference[finite])
            rel = d[scale > 1e-6 * scale.max()] / scale[scale > 1e-6 * scale.max()]
            rows.append({"case": case, "variant": variant, "variable": name, "max_abs_vs_long": float(d.max()),
                         "max_rel_vs_long": float(rel.max()) if rel.size else 0.0, "max_value": float(scale.max())})
      # moments actually passed to DISORT (single precision), against the long expansion
      long_m = caches["long"]["amom"].astype(np.float64)
      long_l = caches["long"]["lmom"]
      for variant in ("adaptive", "fixed"):
         m = caches[variant]["amom"].astype(np.float64)
         lengths = caches[variant]["lmom"] if "lmom" in caches[variant] else None
         n = min(m.shape[0], long_m.shape[0])
         d = np.abs(m[:n] - long_m[:n])
         moments.append({"case": case, "variant": variant,
                         "max_L": int(np.max(lengths)) if lengths is not None else m.shape[0],
                         "min_L": int(np.min(lengths[lengths > 0])) if lengths is not None else m.shape[0],
                         "long_max_L": int(np.max(long_l)),
                         "max_dchi_0_60": float(d[:61].max()), "max_dchi_61_up": float(d[61:].max()) if n > 61 else 0.0,
                         "max_long_chi_beyond_compared": float(np.abs(long_m[n:]).max()) if long_m.shape[0] > n else 0.0,
                         "max_dg": float(np.max(np.abs(caches[variant]["g"].astype(np.float64) - caches["long"]["g"]))),
                         "max_dbext_rel": float(np.max(np.abs(caches[variant]["bext"].astype(np.float64) / np.where(caches["long"]["bext"] == 0, 1, caches["long"]["bext"]) - np.where(caches["long"]["bext"] == 0, 0, 1))))})
   timing = []
   for case in CASES:
      for variant in ("adaptive", "long", "fixed"):
         timing.append(dict(case=case, variant=variant, **json.loads((product(case, variant) / "timing.json").read_text())))
   for filename, table in (("lut_vs_long.csv", rows), ("moments_vs_long.csv", moments), ("timing.csv", timing)):
      with open(RESULTS / filename, "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(table[0]))
         writer.writeheader()
         writer.writerows(table)
   print("results in", RESULTS, "; all finite:",
         all(np.isfinite(v) for table in (rows, moments) for row in table for v in row.values() if isinstance(v, float)))


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   sub = parser.add_subparsers(dest="command", required=True)
   g = sub.add_parser("generate")
   g.add_argument("case", choices=sorted(CASES))
   g.add_argument("variant", choices=["adaptive", "long", "fixed"])
   sub.add_parser("compare")
   sub.add_parser("all")
   args = parser.parse_args()
   make_mirror()
   if args.command == "generate":
      if args.variant == "fixed":
         make_old_tree()
      generate(args.case, args.variant)
   elif args.command == "compare":
      compare()
   else:
      make_old_tree()
      jobs = []
      for case in CASES:
         for variant in ("adaptive", "long", "fixed"):
            if variant == "fixed" and case in ("modis", "slstr"):
               continue
            jobs.append(subprocess.Popen([sys.executable, __file__, "generate", case, variant]))
      if any(job.wait() for job in jobs):
         raise SystemExit("a generation run failed")
      compare()


if __name__ == "__main__":
   main()
