"""Compact Terra MODIS cloud LUT with the V23 or the V24 code, identical settings (read-only for the code).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v24/scattering_refinements/compact_lut.py v23|v24 <model> <out_dir> <channel> ...

The production create_orac_cloud_lut of the selected code (codebase.py) is
called with the settings of the standard Terra MODIS run files (atmosphere 2,
no gas, Rayleigh on, SRF quadrature 1, NSTR 60, nmom 1000) and the compact
Grid B subset compact_liquid_grid_b.lut (liquid water) or
compact_ice_grid_b.lut (ice).  The product is written to <out_dir> (under
validation/v24/scattering_refinements/work/, git-ignored for NetCDF); the wall
time is appended to <out_dir>/timing.json.
"""

import json
import sys
import time
from pathlib import Path

import codebase

code, model, out_dir, channels = sys.argv[1], sys.argv[2], Path(sys.argv[3]), [int(c) for c in sys.argv[4:]]
driver = codebase.select(code)
legacy_disort = codebase.module("oraclut.radiative_transfer.legacy_disort")
if legacy_disort.disort_precision() != "single" or legacy_disort.single_kernel() != "intel_dlm":
   raise SystemExit("set ORAC_LUT_DISORT_PRECISION=single and ORAC_LUT_DISORT_SP_KERNEL=intel_dlm (the V23/V24 kernel)")
here = Path(__file__).resolve().parent
grid = here / ("compact_ice_grid_b.lut" if model.startswith("water-ice") else "compact_liquid_grid_b.lut")
out_dir.mkdir(parents=True, exist_ok=True)
start = time.perf_counter()
status = driver.create_orac_cloud_lut(codebase.ROOT / "create_orac_lut" / "input_files", "terra_modis_v1.inst",
                                      f"{model}.mm", str(grid), out_dir, 2, channelid=channels, gas=0, no_rayleigh=0,
                                      srf_quad=1, version=int(code[1:]), nstreams=60, nmom=1000)
elapsed = time.perf_counter() - start
if status != 0:
   raise SystemExit(f"create_orac_cloud_lut returned {status}")
timing = out_dir / "timing.json"
records = json.loads(timing.read_text()) if timing.is_file() else []
records.append({"code": code, "model": model, "channels": channels, "seconds": elapsed})
timing.write_text(json.dumps(records, indent=1) + "\n")
print(f"{code} {model}: {elapsed:.0f} s -> {out_dir}")
