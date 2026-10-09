"""V24 production smoke tests (validation-only code).

For representative V24 run files (MODIS and dual-view SLSTR; liquid water,
ice spheres and a Baum ice class) it writes a derived run file that keeps
every setting except: the compact validation grid instead of Grid B, three
channels (visible, 3.7 um, thermal), and an output directory under
validation/tmp/v24_smoke.  The derived run file is executed through
create_orac_luts.run (the production run-file path), from a source tree
given on the command line (default: this working tree), with peak memory
recorded.

   python validation/v24/smoke_v24.py generate <case> [SOURCE_TREE]
   python validation/v24/smoke_v24.py check

check verifies, for every case:
   * V24 file name, written to the smoke directory only;
   * every LUT variable finite where defined;
   * for Mie classes: the scattering cache's radius-grid marker, liquid
     upper limit at the lattice node >= 3.5 r_e, refined node counts and
     adaptive expansion lengths (> 1000 where expected);
   * for the Baum class: bitwise identity with the same run from source
     revision dbcc42c (the radius-grid change must not touch Baum);
and that the Grid B coordinates expanded by load_lutstr are those of the
V23 products (effective radius, optical depth, geometry), so that the full
V24 products use the V23 sampling.
"""

import importlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SMOKE = ROOT / "validation" / "tmp" / "v24_smoke"
MIRROR = ROOT / "validation" / "tmp" / "radius_grid" / "compact_inputs"
LUT_DIR = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
CASES = {
   "modis_liquid": ("terra_modis_cloud_liquid-water_stg_v24.run", "compact_radius_liquid_grid.lut", [1, 20, 31]),
   "modis_ice_sph": ("aqua_modis_cloud_water-ice_sph_v24.run", "compact_radius_ice_grid.lut", [1, 20, 31]),
   "modis_ice_agg": ("terra_modis_cloud_water-ice_agg_v24.run", "compact_radius_ice_grid.lut", [1, 20, 31]),
   "slstr_liquid": ("sentinel-3a_slstr_cloud_liquid-water_263_v24.run", "compact_radius_liquid_grid.lut", [1, 7, 8]),
   "slstr_ice_sph": ("sentinel-3b_slstr_cloud_water-ice_sph_v24.run", "compact_radius_ice_grid.lut", [1, 7, 8]),
}


def derived_run(case, out):
   runfile, grid, channels = CASES[case]
   text = (ROOT / "runs" / runfile).read_text()
   lines = []
   for line in text.splitlines():
      key = line.split("=")[0].strip()
      if key == "in_path":
         line = f"in_path       = '{MIRROR}'"
      elif key == "lutfile":
         line = f"lutfile       = '{grid}'"
      elif key == "channelid":
         line = f"channelid     = {channels}"
      elif key == "out_path":
         line = f"out_path      = '{out}'"
      lines.append(line)
   path = out / runfile
   path.write_text("\n".join(lines) + "\n")
   return path


def generate(case, tree):
   tree = Path(tree)
   label = "dbcc42c" if tree != ROOT else "v24"
   out = SMOKE / case / label
   out.mkdir(parents=True, exist_ok=True)
   run = derived_run(case, out)
   sys.path.insert(0, str(tree))
   sys.path.insert(0, str(tree / "src"))
   import create_orac_luts
   started = time.time()
   status = create_orac_luts.run(str(run))
   record = {"case": case, "source": str(tree), "revision": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
             capture_output=True, text=True).stdout.strip() if tree == ROOT else "dbcc42c (export)",
             "status": status, "wall_s": time.time() - started,
             "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0}
   (out / "smoke.json").write_text(json.dumps(record, indent=2) + "\n")
   print(record, flush=True)


def check():
   sys.path.insert(0, str(ROOT))
   sys.path.insert(0, str(ROOT / "src"))
   from oraclut.io.lut import read_lut
   from oraclut.idl_mirror import load_lutstr, load_mmdat
   import create_orac_luts
   bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
   gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
   inputs = ROOT / "create_orac_lut" / "input_files"
   ok = True

   def report(flag, message):
      nonlocal ok
      ok &= bool(flag)
      print(("PASS  " if flag else "FAIL  ") + message)

   for case, (runfile, grid, channels) in CASES.items():
      out = SMOKE / case / "v24"
      products = sorted(out.glob("*.nc"))
      record = json.loads((out / "smoke.json").read_text())
      report(len(products) == 1 and products[0].name.endswith("_v24.nc"),
             f"{case}: one product {[p.name for p in products]} (status {record['status']}, "
             f"{record['wall_s']:.0f} s, peak RSS {record['peak_rss_mb']:.0f} MB, revision {record['revision'][:7]})")
      lut = read_lut(products[0])
      bad = [n for n in lut.variable_names if np.asarray(lut.variables[n]).dtype.kind == "f"
             and not np.all(np.isfinite(np.asarray(lut.variables[n])))]
      report(not bad, f"{case}: all floating-point variables finite {bad if bad else ''}")
      settings = create_orac_luts.read_runfile(str(ROOT / "runs" / runfile))
      mmstr = load_mmdat(inputs / "microphysics" / settings["mmfile"], inputs)
      cache = np.load(out / "work" / "scatfile.npz") if (out / "work" / "scatfile.npz").exists() else None
      if mmstr.distname[0] == "modified_gamma":
         xres = gsp.refined_xres(mmstr, 0)
         report(xres == (gsp.LIQUID_REFINED_XRES if "liquid" in case else gsp.ICE_SPHERE_REFINED_XRES),
                f"{case}: refined size-parameter step {xres}")
         factor = gsp.radius_upper_factor(mmstr, 0)
         re = 10.0
         wn = 1.0 / 0.6462791
         params, npts = bwgp.mie_integration_limits("modified_gamma", re, mmstr.s[0], wn, factor, xres)
         legacy_params, legacy_npts = bwgp.mie_integration_limits("modified_gamma", re, mmstr.s[0], wn, factor)
         if legacy_npts is None:
            legacy_npts = max(int(2.0 * np.pi * (params[3] - params[2]) * wn / bwgp.MIE_XRES), 200)
         upper_ok = (params[3] >= 3.5 * re and params[3] - bwgp.legacy_radius_spacing(wn) < 3.5 * re) if factor \
            else params[3] == 100.0
         report(upper_ok and params == legacy_params,
                f"{case}: r_e 10 um at 0.646 um integrates {params[2]}-{params[3]:.3f} um "
                f"({'lattice node >= 3.5 r_e' if factor else 'legacy 0.001-100 um'})")
         report((npts - 1) % (legacy_npts - 1) == 0 and (npts - 1) // (legacy_npts - 1) == (16 if "liquid" in case else 8),
                f"{case}: {npts} nodes = ({legacy_npts} - 1) x {(npts - 1) // (legacy_npts - 1)} + 1")
         if cache is not None:
            report(str(cache["radius_grid"]) == create_orac_luts.RADIUS_GRID and np.max(cache["lmom"]) > 1000,
                   f"{case}: cache marker {cache['radius_grid']}, adaptive L {int(np.min(cache['lmom'][cache['lmom'] > 0]))}"
                   f"-{int(np.max(cache['lmom']))}")
         else:
            print(f"INFO  {case}: private scattering cache not kept (work area removed); cache checks skipped")
      else:
         reference = SMOKE / case / "dbcc42c"
         other = read_lut(next(reference.glob("*.nc")))
         equal = lambda a, b: np.array_equal(a, b, equal_nan=a.dtype.kind == "f")
         same = all(equal(np.asarray(lut.variables[n]), np.asarray(other.variables[n])) for n in lut.variable_names)
         report(same and set(lut.variable_names) == set(other.variable_names),
                f"{case}: Baum product bitwise identical to dbcc42c")
   # Full Grid B coordinates against the V23 products
   for grid, product in (("liquid-water-cloud-grid-b.lut", "terra_modis_m_liquid-water_a01_pstg_v23.nc"),
                         ("ice-cloud-grid-b.lut", "sentinel-3a_slstr_m_water-ice_a01_psph_v23.nc")):
      lutstr = load_lutstr(inputs / "lut" / grid, 75.0)
      v23 = read_lut(LUT_DIR / product)
      pairs = (("optical_depth", lutstr.opd), ("effective_radius", lutstr.efr), ("solar_zenith", lutstr.soz),
               ("satellite_zenith", lutstr.saz), ("relative_azimuth", lutstr.raa))
      same = all(np.allclose(np.asarray(v23.variables[n], dtype=np.float64), np.asarray(v, dtype=np.float64),
                             rtol=1e-6, atol=1e-6) for n, v in pairs)
      report(same, f"{grid}: expanded coordinates equal those of {product}")
   print("ALL PASSED" if ok else "FAILURES")
   return ok


def main():
   if sys.argv[1] == "check":
      sys.exit(0 if check() else 1)
   if sys.argv[1] != "generate":
      sys.exit("usage: smoke_v24.py generate <case> [SOURCE_TREE] | check")
   generate(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else ROOT)


if __name__ == "__main__":
   main()
