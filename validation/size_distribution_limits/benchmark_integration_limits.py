"""Computational cost of the liquid-cloud upper radius limit (Phase 2D benchmark).

Validation-only code.  Times the production create_bwgp (Mie size integration
at the production 1000 phase-function points) for the V23 Grid B liquid
effective radii, 1-39 um (create_orac_lut/input_files/lut/liquid-water-cloud-grid-b.lut),
at representative MODIS wavelengths, for

   reference   0.001-100 um (create_bwgp radius_upper_factor=None, the IDL limits)
   adopted     the production liquid-water rule: legacy lattice up to the first
               node >= LIQUID_UPPER_RADIUS_FACTOR x re

and records single-thread CPU time and radius-point counts.

   PYTHONPATH=src python validation/size_distribution_limits/benchmark_integration_limits.py
"""

import csv
import importlib
import os
import time
from pathlib import Path

# Timings are single-thread CPU: numpy's thread pool otherwise inflates CPU time
# (process_time counts every thread) without speeding up the Mie kernel.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
   os.environ[variable] = "1"

import numpy as np   # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

create_bwgp_module = importlib.import_module("oraclut.idl_mirror.create_bwgp")
from oraclut.idl_mirror.generate_scattering_properties import _interpol_complex   # noqa: E402
from oraclut.idl_mirror.load_mmdat import load_mmdat                             # noqa: E402

INPUTS = ROOT / "create_orac_lut" / "input_files"
RADII = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0,
         29.0, 31.0, 33.0, 35.0, 37.0, 39.0]
WAVELENGTHS = (0.6462791, 1.6291, 2.1142, 3.7850, 11.0262)
UPPER_FACTOR = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties").LIQUID_UPPER_RADIUS_FACTOR


def points(upper, wavelength):
   return max(int(2.0 * np.pi * (upper - 0.001) / wavelength / 0.4), 200)


def main():
   table = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS).comp[0]
   abscissas, _ = create_bwgp_module.quadrature("g", 1000)
   qv = -abscissas
   rows = []
   for wavelength in WAVELENGTHS:
      ri = _interpol_complex(table.cm, table.wl, [wavelength])
      for variant, factor in (("reference", None), ("adopted", UPPER_FACTOR)):
         started = time.process_time()
         for re in RADII:
            create_bwgp_module.create_bwgp("modified_gamma", re, 0.1111111, ri, [wavelength], qv, radius_upper_factor=factor)
         cpu = time.process_time() - started
         count = sum(points(100.0, wavelength) if factor is None else
                     create_bwgp_module.lattice_upper_radius(re, factor, 1.0 / wavelength)[1] for re in RADII)
         rows.append({"wavelength_um": wavelength, "variant": variant, "cpu_s": cpu, "radius_points": count})
         print(wavelength, variant, f"{cpu:.1f} s", count, flush=True)
   out = HERE / "results" / "integration_limit_benchmark.csv"
   out.parent.mkdir(parents=True, exist_ok=True)
   with open(out, "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)
   total = {v: sum(row["cpu_s"] for row in rows if row["variant"] == v) for v in ("reference", "adopted")}
   print("total CPU s:", total, "speed-up:", total["reference"] / total["adopted"])


if __name__ == "__main__":
   main()
