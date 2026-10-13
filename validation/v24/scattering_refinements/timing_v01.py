"""Cost of the three V24 refinements in the scattering step, per Grid B liquid-water radius (read-only).

   OMP_NUM_THREADS=1 python validation/v24/scattering_refinements/timing.py <model> <wavelength um> <out.json>

For each Grid B effective radius at one wavelength, the wall time of the
production routines (working tree) for:

   v23        create_bwgp with the IDL limits and grid on the 1000-point Gauss
              rule (the V23 Mie step per wavelength and radius)
   m1         + liquid-water limit 3.5 r_e on the legacy lattice, 1000 points
   m1_m3      + nested radius refinement, 1000 points
   m1_m2      + adaptive Legendre expansion, legacy-step lattice (refinement
              disabled in this process only, by replacing refined_xres)
   v24        all three: generate_scattering_properties._adaptive_mie_wavelength

The 1000-point rule is that of the IDL quadrature.pro.  One CPU thread.
"""

import json
import sys
import time

import numpy as np

import codebase

codebase.select("v24")
mirror = codebase.module("oraclut.idl_mirror")
gsp = codebase.module("oraclut.idl_mirror.generate_scattering_properties")
cb = codebase.module("oraclut.idl_mirror.create_bwgp")

model, wavelength, out_file = sys.argv[1], float(sys.argv[2]), sys.argv[3]
inputs = codebase.ROOT / "create_orac_lut" / "input_files"
mm = mirror.load_mmdat(inputs / "microphysics" / f"{model}.mm", inputs)
grid = "ice-cloud.lut" if model.startswith("water-ice") else "liquid-water-cloud.lut"
inst = mirror.load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=[1])
radii = mirror.load_lutstr(inputs / "lut" / grid, inst.max_sat_zenith).efr.astype(np.float64)
ri = gsp._interpol_complex(mm.comp[0].cm, mm.comp[0].wl, [wavelength])
qv = -cb.quadrature("G", 1000)[0]
factor, refine = gsp.radius_upper_factor(mm, 0), gsp.refined_xres(mm, 0)


def timed(function):
   start = time.perf_counter()
   value = function()
   return time.perf_counter() - start, value


rows = []
for radius in radii:
   row = {"effective_radius": float(radius)}
   row["v23"] = timed(lambda: cb.create_bwgp(mm.distname[0], radius, mm.s[0], ri, [wavelength], qv))[0]
   row["m1"] = timed(lambda: cb.create_bwgp(mm.distname[0], radius, mm.s[0], ri, [wavelength], qv,
                                            radius_upper_factor=factor))[0]
   row["m1_m3"] = timed(lambda: cb.create_bwgp(mm.distname[0], radius, mm.s[0], ri, [wavelength], qv,
                                               radius_upper_factor=factor, refined_xres=refine))[0]
   seconds, result = timed(lambda: gsp._adaptive_mie_wavelength(mm, np.ones((1, 1)), np.asarray([[radius]]), ri,
                                                                wavelength, "timing"))
   row["v24"], row["L"], row["Nq"] = seconds, int(result[6][0]), int(result[7]["order"][0])
   rows.append(row)
   print(json.dumps(row), flush=True)

# the adaptive expansion without the refinement (this process only)
gsp.refined_xres = lambda mmstr, c: None
for row in rows:
   radius = row["effective_radius"]
   row["m1_m2"] = timed(lambda: gsp._adaptive_mie_wavelength(mm, np.ones((1, 1)), np.asarray([[radius]]), ri,
                                                             wavelength, "timing"))[0]
totals = {key: float(sum(row[key] for row in rows)) for key in ("v23", "m1", "m1_m3", "m1_m2", "v24")}
with open(out_file, "w") as handle:
   json.dump({"model": model, "wavelength": wavelength, "totals_seconds": totals, "rows": rows}, handle, indent=1)
print("totals", json.dumps(totals))
