"""Post-production checks of the V25 cloud LUT products (validation only).

For every expected V25 product (the 40 of scripts/submit_v25_cloud_luts.sh):
existence and size, the V2 structure (read_lut validation), the Grid B
dimensions, finite values, the V25 global attributes, and a comparison with
the V24 product: every reflection / transmission operator within the
single-precision reproducibility floor of the production calculation
(validation/v24/results/REPORT_disort_options.md: up to about 4e-4 relative in
the diffuse operators between compute nodes), E_md changed smoothly (the
ratio V25 / V24 varies without jumps along the optical-depth grid).  Writes
validation/v25/results/products.csv and prints a summary.

   python validation/v25/check_products.py [PRODUCT.nc ...]
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.io.lut import read_lut   # noqa: E402

LUT_DIR = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
RESULTS = HERE / "results"
LIQUID = ("old", "stg", "240", "253", "263", "273")
ICE = ("sph", "agg", "ghm", "src")
GRID = {"liquid": (24, 24, 21, 21, 37), "ice": (24, 29, 21, 21, 37)}
RT_TOLERANCE = 2e-3                                   # absolute, on operators of order 1 (reproducibility floor x 5)


def expected_products():
   for instrument, platforms in (("modis", ("aqua", "terra")), ("slstr", ("sentinel-3a", "sentinel-3b"))):
      for platform in platforms:
         for family, names, substance in (("liquid", LIQUID, "liquid-water"), ("ice", ICE, "water-ice")):
            for name in names:
               yield family, LUT_DIR / f"{platform}_{instrument}_m_{substance}_a01_p{name}_v25.nc"


def check(path, family):
   row = {"product": path.name, "exists": path.is_file(), "size_MB": None, "structure": None, "grid": None,
          "finite": None, "v25_attributes": None, "rt_max_abs_diff_vs_v24": None, "rt_within_floor": None,
          "e_md_min": None, "e_md_max": None, "e_md_ratio_min": None, "e_md_ratio_max": None,
          "e_md_ratio_max_jump": None, "e_md_smooth": None, "status": "MISSING"}
   if not path.is_file():
      return row
   row["size_MB"] = round(path.stat().st_size / 1e6, 1)
   try:
      lut = read_lut(path)
      row["structure"] = True
   except Exception as exc:                                  # noqa: BLE001
      row["structure"] = f"error: {exc}"
      row["status"] = "BAD STRUCTURE"
      return row
   d = lut.dimensions
   dims = (d["optical_depth"], d["effective_radius"], d["solar_zenith"], d["satellite_zenith"], d["relative_azimuth"])
   row["grid"] = dims == GRID[family]
   finite = all(np.all(np.isfinite(np.asarray(lut.variables[n]))) for n in lut.variable_names
                if np.asarray(lut.variables[n]).dtype.kind == "f" and n not in ("Bext", "BextRat", "SSA", "G"))
   row["finite"] = bool(finite)
   g = lut.global_attributes
   row["v25_attributes"] = g.get("cloud_temperature_profile") == "adiabatic" and float(g.get("cloud_top_temperature_K", 0)) == 240.0
   v24 = path.with_name(path.name.replace("_v25.nc", "_v24.nc"))
   if v24.is_file():
      old = read_lut(v24)
      worst = 0.0
      for name in ("R_0v", "R_0d", "R_dd", "R_dv", "T_00", "T_0d", "T_dd", "T_dv"):
         if name in lut.variables:
            x = np.asarray(lut.variables[name], dtype=np.float64)
            y = np.asarray(old.variables[name], dtype=np.float64)
            worst = max(worst, float(np.nanmax(np.abs(x - y))))
      row["rt_max_abs_diff_vs_v24"] = worst
      row["rt_within_floor"] = worst <= RT_TOLERANCE
      e_new = np.asarray(lut.variables["E_md"], dtype=np.float64)      # (saz, opd, efr, thermal)
      e_old = np.asarray(old.variables["E_md"], dtype=np.float64)
      row["e_md_min"], row["e_md_max"] = float(e_new.min()), float(e_new.max())
      ok = e_old > 1e-3
      ratio = np.where(ok, e_new / np.where(ok, e_old, 1.0), np.nan)
      row["e_md_ratio_min"], row["e_md_ratio_max"] = float(np.nanmin(ratio)), float(np.nanmax(ratio))
      jump = np.nanmax(np.abs(np.diff(ratio, axis=1)))                 # along the optical-depth grid
      row["e_md_ratio_max_jump"] = float(jump)
      row["e_md_smooth"] = bool(jump < 0.5)
   checks = [row["structure"] is True, row["grid"], row["finite"], row["v25_attributes"]]
   if row["rt_within_floor"] is not None:
      checks += [row["rt_within_floor"], row["e_md_smooth"]]
   row["status"] = "OK" if all(checks) else "CHECK"
   return row


def main():
   RESULTS.mkdir(parents=True, exist_ok=True)
   rows = []
   if len(sys.argv) > 1:
      targets = [("liquid" if "liquid-water" in p else "ice", Path(p)) for p in sys.argv[1:]]
   else:
      targets = list(expected_products())
   for family, path in targets:
      rows.append(check(path, family))
      print(f"{rows[-1]['status']:14s} {path.name}  rt_diff={rows[-1]['rt_max_abs_diff_vs_v24']}  "
            f"E_md={rows[-1]['e_md_min']}..{rows[-1]['e_md_max']}  ratio_jump={rows[-1]['e_md_ratio_max_jump']}", flush=True)
   with open(RESULTS / "products.csv", "w", newline="") as handle:
      w = csv.DictWriter(handle, fieldnames=list(rows[0]))
      w.writeheader()
      w.writerows(rows)
   counts = {s: sum(r["status"] == s for r in rows) for s in ("OK", "CHECK", "MISSING", "BAD STRUCTURE")}
   print(counts)


if __name__ == "__main__":
   main()
