"""Class scattering properties of a microphysical model with the V23 or the V24 code (read-only).

   python validation/v24/scattering_refinements/optics.py v23|v24 <model> <out.npz> --channels 1 20 --efr 5 21 39

<model> is a Terra MODIS cloud microphysics file stem (liquid-water_240,
liquid-water_old, water-ice_agg, water-ice_sph).  The run uses the production
generate_scattering_properties of the selected code (codebase.py) with the
Terra MODIS instrument, the Gueymard 2018 solar spectrum, SRF quadrature 1,
nmom = 1000 and the given channels and effective radii (nodes of the active
Grid B definition).  Saved: bext550, w550, g550, bextrat, bext, w, g, vavg,
the Legendre moments amom, the number of moments of each phase function lmom,
and the wall time.
"""

import argparse
import time

import numpy as np

import codebase

parser = argparse.ArgumentParser()
parser.add_argument("code", choices=("v23", "v24"))
parser.add_argument("model")
parser.add_argument("out")
parser.add_argument("--channels", type=int, nargs="+", required=True)
parser.add_argument("--efr", type=float, nargs="+", required=True)
args = parser.parse_args()

codebase.select(args.code)
mirror = codebase.module("oraclut.idl_mirror")
gsp = codebase.module("oraclut.idl_mirror.generate_scattering_properties")

inputs = codebase.ROOT / "create_orac_lut" / "input_files"
inst = mirror.load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=args.channels)
grid = "ice-cloud.lut" if args.model.startswith("water-ice") else "liquid-water-cloud.lut"
lut = mirror.load_lutstr(inputs / "lut" / grid, inst.max_sat_zenith)
efr = np.asarray(args.efr, dtype=lut.efr.dtype)
if not np.all(np.isin(efr, lut.efr)):
   raise SystemExit(f"effective radii {args.efr} are not all nodes of {grid}")
lut.efr, lut.efr_n = efr, efr.size
srf, nwvl = mirror.load_srfstrarr(inst, inputs / "sun" / "Gueymard2018.sssi", 1, inputs)
mm = mirror.load_mmdat(inputs / "microphysics" / f"{args.model}.mm", inputs)

start = time.perf_counter()
result = gsp.generate_scattering_properties(srf, nwvl, inst, mm, lut, 1000)
elapsed = time.perf_counter() - start
first, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom = result
lmom = np.full(bext.shape, first, dtype=np.int64) if np.isscalar(first) else np.asarray(first)
np.savez(args.out, code=args.code, model=args.model, channels=np.asarray(args.channels), efr=efr,
         wavelengths=np.asarray([s.wvl_centre for s in srf]), bext550=bext550, w550=w550, g550=g550, bextrat=bextrat,
         bext=bext, w=w, g=g, vavg=vavg, amom=amom, lmom=lmom, seconds=elapsed)
print(f"{args.code} {args.model}: {elapsed:.1f} s, lmom {lmom.min()}-{lmom.max()}, wrote {args.out}")
