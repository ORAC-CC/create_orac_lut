"""Stage 7: full-grid scattering properties of the corrected V22 generator against the production V21 LUTs.

The production V21 single-precision LUTs (releases/idl_lut_v21_single.release)
store, for all 36 Terra MODIS channels and every effective radius of the
original V2.1 grid, the scattering properties BextOUT, BextRatOUT, SSAOUT and
GOUT of generate_scattering_properties.pro.  This script runs the corrected
Python generator with scat_only=1 (same run configuration as the V22
production runs: atmosphere 2, channels 1-36, srf_quad 1, nmom 1000) and
compares those four arrays bitwise.

The comparison of the two Mie kernels (results/full_grid/*_intel_dlm.json,
*_gfortran.json and mie_kernels_scatfile.json) was made on 2026-10-11 with a
development revision of oraclut.optics.legacy_mie that could call the
original IDL Mie DLM (ORAC_LUT_MIE_KERNEL=intel_dlm).  Both kernels gave
bitwise-identical scattering properties on the full pold and p240 grids, so
that backend was not retained; the production generator, and this script,
use the gfortran build.  The kernel label below records the environment
setting of the run.

Usage: python full_grid_scattering.py <pold|p240|pagg> <work_dir>
   -> results/full_grid/<product>_<kernel>.json
"""

import importlib
import json
import os
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
V21 = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
PRODUCTS = {"pold": ("liquid-water_old.mm", "liquid-water-cloud.lut", "terra_modis_m_liquid-water_a01_pold_v21.nc"),
            "p240": ("liquid-water_240.mm", "liquid-water-cloud.lut", "terra_modis_m_liquid-water_a01_p240_v21.nc"),
            "pagg": ("water-ice_agg.mm", "ice-cloud.lut", "terra_modis_m_water-ice_a01_pagg_v21.nc")}
VARIABLES = {"extinction_coefficient": "bext", "extinction_coefficient_ratio": "bextrat", "single_scatter_albedo": "w",
             "asymmetry_parameter": "g"}


def ulps(a, b):
   a = np.asarray(a, np.float32).ravel().view(np.int32).astype(np.int64)
   b = np.asarray(b, np.float32).ravel().view(np.int32).astype(np.int64)
   a = np.where(a < 0, -(a & 0x7FFFFFFF), a)
   b = np.where(b < 0, -(b & 0x7FFFFFFF), b)
   return np.abs(a - b)


def main():
   product, work = sys.argv[1], Path(sys.argv[2])
   mmfile, lutfile, v21name = PRODUCTS[product]
   kernel = os.environ.get("ORAC_LUT_MIE_KERNEL", "").strip().lower() or "gfortran"
   production = importlib.import_module("create_orac_luts")
   production.create_orac_cloud_lut(ROOT / "create_orac_lut" / "input_files", "terra_modis_v1.inst", mmfile, lutfile, work / "luts", "2",
                                    channelid=list(range(1, 37)), gas=0, no_rayleigh=0, srf_quad=1, reuse_scat=0, scat_only=1,
                                    tmatrix_path=None, version=22, driver=None, nstreams=60, nmom=1000, work_path=work / "work")
   scat = np.load(work / "work" / "scatfile.npz")
   m = scat["bext"].shape[0] // 2                       # IDL: m = Nwvl_Max / 2
   result = {"product": product, "mie_kernel": kernel, "v21": str(V21 / v21name), "variables": {}}
   with Dataset(V21 / v21name) as d:
      d.set_auto_maskandscale(False)
      for name, key in VARIABLES.items():
         ref = np.asarray(d[name][:], np.float32)               # (effective_radius, channels)
         mine = np.asarray(scat[key][m], np.float32).T          # [l, r] -> (r, l)
         u = ulps(ref, mine)
         differing = ref.view(np.uint32) != mine.view(np.uint32)
         result["variables"][name] = {"shape": list(ref.shape), "differing": int(differing.sum()), "max_ulp": int(u.max()),
                                      "max_abs": float(np.max(np.abs(ref.astype(np.float64) - mine))),
                                      "differing_at": [[int(i), int(j)] for i, j in zip(*np.nonzero(differing))][:20]}
   out = HERE / "results" / "full_grid"
   out.mkdir(parents=True, exist_ok=True)
   (out / f"{product}_{kernel}.json").write_text(json.dumps(result, indent=1) + "\n")
   print(product, kernel, {k: (v["differing"], v["max_ulp"]) for k, v in result["variables"].items()})


if __name__ == "__main__":
   main()
