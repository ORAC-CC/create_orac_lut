"""Predicted Mie cost of full-channel V24 cloud runs per radius-refinement level.

Validation-only code.  The Mie kernel's work for one radius node is
proportional to NStop(x) x (number of angles).  This script sums
NStop x angles over every node of every (wavelength, effective radius) Mie
evaluation of generate_scattering_properties (production limits, adaptive
quadrature order Nq + check angles, the 0.55 um reference included; the
quadrature retries, which are rare, are not), for a forced refinement level k.
The constant (CPU seconds per unit) is calibrated on the measured six-channel
benchmark (results/benchmark/*_dbcc42c.json, which also counts nodes), so
the result is the predicted single-thread Mie CPU time on the benchmark host.

   python validation/radius_grid/cost_model.py
"""

import importlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
le = importlib.import_module("oraclut.idl_mirror.legendre_expansion")
from oraclut.idl_mirror import load_inststr, load_lutstr, load_mmdat, load_srfstrarr   # noqa: E402

CASES = {"liquid": ("liquid-water_stg.mm", "liquid-water-cloud-grid-b.lut"),
         "ice_sph": ("water-ice_sph.mm", "ice-cloud-grid-b.lut")}
INSTRUMENTS = {"modis": ("terra_modis_v1.inst", list(range(1, 37))), "slstr": ("sentinel-3a_slstr_v1.inst", list(range(1, 10))),
               "benchmark": ("terra_modis_v1.inst", [3, 1, 6, 7, 20, 31])}


def nstop(x):
   return np.where(x < 0.02, 2.0, x + 4.05 * np.cbrt(x) + 2.0)


def units(mmstr, radii, wavelengths, level):
   factor = gsp.radius_upper_factor(mmstr, 0)
   variance = float(mmstr.s[0])
   angles_check = le.CHECK_THETA.size
   total = 0.0
   nodes = 0
   for wl in wavelengths:
      wn = 1.0 / wl
      for re in radii:
         params, npts = bwgp.mie_integration_limits("modified_gamma", re, variance, wn, factor)
         if npts is None:
            npts = max(int(2.0 * np.pi * (params[3] - params[2]) * wn / bwgp.MIE_XRES), 200)
         if level is None:
            level_r = bwgp.radius_refinement_level((params[3] - params[2]) / (npts - 1), wn, re, variance,
                                                   gsp.refined_xres(mmstr, 0))
         else:
            level_r = level
         n = (npts - 1) * 2**level_r + 1
         r = np.linspace(params[2], params[3], n)
         x = 2.0 * np.pi * r * wn
         nq = le.initial_quadrature_order(le.mie_phase_degree(2.0 * np.pi * params[3] * wn))
         total += float(np.sum(nstop(x))) * (nq + angles_check)
         nodes += n
   return total, nodes


def main():
   rows = []
   calibration = {}
   for case, (mmfile, grid) in CASES.items():
      mmstr = load_mmdat(INPUTS / "microphysics" / mmfile, INPUTS)
      for instrument, (instfile, channels) in INSTRUMENTS.items():
         inststr = load_inststr(INPUTS / "inst" / instfile, requestedchannelid=channels)
         lutstr = load_lutstr(INPUTS / "lut" / grid, inststr.max_sat_zenith)
         srfstrarr, _ = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
         wavelengths = [0.55] + [float(s.wvl_centre) for s in srfstrarr]
         radii = np.asarray(lutstr.efr, dtype=float)
         for level in (0, 2, 3, 4, None):
            u, n = units(mmstr, radii, wavelengths, level)
            rows.append({"case": case, "instrument": instrument, "level": "rule" if level is None else level, "units": u, "nodes": n})
      measured = json.loads((HERE / "results" / "benchmark" / f"{case}_dbcc42c.json").read_text())
      bench0 = next(r for r in rows if r["case"] == case and r["instrument"] == "benchmark" and r["level"] == 0)
      calibration[case] = measured["mie_cpu_s"] / bench0["units"]
      print(f"{case}: benchmark nodes model {bench0['nodes']} measured {measured['radius_nodes']} "
            f"(the 0.55 um reference and repeated radii make these differ slightly); {calibration[case]:.3e} s per unit")
   for row in rows:
      row["predicted_mie_cpu_h"] = row["units"] * calibration[row["case"]] / 3600.0
      print(f"{row['case']:8s} {row['instrument']:9s} level {str(row['level']):4s} nodes {row['nodes']:>10d}  "
            f"Mie CPU {row['predicted_mie_cpu_h']:7.2f} h")
   out = HERE / "results" / "cost_model.json"
   out.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
   main()
