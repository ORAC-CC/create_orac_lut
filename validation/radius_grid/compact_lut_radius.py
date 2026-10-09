"""Compact LUT comparison of radius-integration refinement levels (V24 radius grid).

Validation-only code.  Runs the production generator
(create_orac_luts.create_orac_cloud_lut) on compact grids with the Mie radius
integration forced to a nested dyadic refinement level k of the legacy
lattice (spacing h_0 / 2^k on the same endpoints; k = 0 is the V23 grid as
modified by c6ad545), or with the production rule as it stands ("adopted").
Everything else (limits, SRF treatment, Mie kernel, adaptive Legendre,
DISORT, writer) is the production path at the current source revision.

Cases: MODIS and dual-view SLSTR liquid-water cloud (stg) on
compact_radius_liquid_grid.lut, and MODIS ice spheres (sph) on
compact_radius_ice_grid.lut.

   python compact_lut_radius.py generate <case> <k0..k6|adopted>
   python compact_lut_radius.py compare [REFERENCE]       (default k6)

The level is imposed at run time: on the production rule
create_bwgp.radius_refinement_level when it exists, otherwise by refining the
node count that create_bwgp.mie_integration_limits returns (modified-gamma
classes only; log-normal classes are untouched).  Products are written under
validation/tmp/radius_grid/compact; summaries to results/compact_lut.  Forced
levels use the memory-bounded Mie summation of batched_mie.py (same sums,
blockwise; rounding-level differences only).
"""

import argparse
import csv
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

INPUTS = ROOT / "create_orac_lut" / "input_files"
MIRROR = ROOT / "validation" / "tmp" / "radius_grid" / "compact_inputs"
PRODUCTS = ROOT / "validation" / "tmp" / "radius_grid" / "compact"
RESULTS = HERE / "results" / "compact_lut"
CASES = {
   # name: (instrument file, microphysics, grid, channels)
   # MODIS 0.469, 0.646, 0.858, 1.63, 2.13, 3.79 (solar + thermal) and 11.03 um
   "modis_liquid": ("terra_modis_v1.inst", "liquid-water_stg.mm", "compact_radius_liquid_grid.lut", [3, 1, 2, 6, 7, 20, 31]),
   # SLSTR (dual view) 0.555, 0.865, 1.61, 2.25, 3.74 and 10.85 um
   "slstr_liquid": ("sentinel-3a_slstr_v1.inst", "liquid-water_stg.mm", "compact_radius_liquid_grid.lut", [1, 3, 5, 6, 7, 8]),
   "modis_ice_sph": ("terra_modis_v1.inst", "water-ice_sph.mm", "compact_radius_ice_grid.lut", [1, 2, 6, 7, 20, 31]),
}
VARIANTS = [f"k{k}" for k in range(7)] + ["adopted"]


def make_mirror():
   (MIRROR / "lut").mkdir(parents=True, exist_ok=True)
   for item in INPUTS.iterdir():
      link = MIRROR / item.name
      if item.name != "lut" and not link.exists():
         link.symlink_to(item)
   for grid in sorted({case[2] for case in CASES.values()}):
      (MIRROR / "lut" / grid).write_text((HERE / grid).read_text())


def force_level(level):
   bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
   if hasattr(bwgp, "radius_refinement_level"):
      bwgp.radius_refinement_level = lambda *args, **kwargs: level
      return
   production = bwgp.mie_integration_limits

   def refined(distname, rm, s, wavenumber, radius_upper_factor=None):
      params, npts = production(distname, rm, s, wavenumber, radius_upper_factor)
      if distname == "modified_gamma":
         if npts is None:
            npts = max(int(2.0 * np.pi * (params[3] - params[2]) * wavenumber / bwgp.MIE_XRES), 200)
         npts = (npts - 1) * 2**level + 1
      return params, npts

   bwgp.mie_integration_limits = refined


def generate(case, variant):
   import create_orac_luts

   if variant != "adopted":
      sys.path.insert(0, str(HERE))
      import batched_mie
      batched_mie.install(importlib.import_module("oraclut.idl_mirror.create_bwgp"))
      force_level(int(variant[1:]))
   timing = {}
   scattering = create_orac_luts.generate_scattering_properties

   def timed(*args, **kwargs):
      started = time.process_time()
      result = scattering(*args, **kwargs)
      timing["scattering_cpu_s"] = time.process_time() - started
      return result

   create_orac_luts.generate_scattering_properties = timed
   instfile, mmfile, grid, channels = CASES[case]
   out = PRODUCTS / case / variant
   out.mkdir(parents=True, exist_ok=True)
   started = time.process_time()
   status = create_orac_luts.create_orac_cloud_lut(MIRROR, instfile, mmfile, grid, out, 2, channelid=channels, srf_quad=1,
                                                   version=None, nstreams=60, nmom=None, work_path=out / "work")
   timing["total_cpu_s"] = time.process_time() - started
   timing["status"] = status
   (out / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
   print(case, variant, timing, flush=True)


def _rel(a, b, floor):
   a = np.asarray(a, dtype=np.float64)
   b = np.asarray(b, dtype=np.float64)
   mask = np.abs(a) > floor
   return float(np.max(np.abs(b[mask] / a[mask] - 1.0))) if np.any(mask) else 0.0


def compare(reference="k6"):
   from oraclut.io.lut import read_lut
   from oraclut.idl_mirror import load_inststr
   from oraclut.validation.quadrature import cell_midpoints, interpolate_2d, measurement_uncertainty

   RESULTS.mkdir(parents=True, exist_ok=True)
   optics_rows, lut_rows, retrieval_rows, timing_rows = [], [], [], []
   for case, (instfile, mmfile, grid, channels) in CASES.items():
      base = PRODUCTS / case
      if not (base / reference / "timing.json").exists():
         continue
      inst = load_inststr(INPUTS / "inst" / instfile, requestedchannelid=channels)
      ref_cache = np.load(base / reference / "work" / "scatfile.npz")
      ref_lut = read_lut(next((base / reference).glob("*.nc")))
      radius = np.asarray(ref_lut.variables["effective_radius"], dtype=np.float64)
      tau = np.asarray(ref_lut.variables["optical_depth"], dtype=np.float64)
      m = ref_cache["bext"].shape[0] // 2
      for variant in VARIANTS:
         if variant == reference or not (base / variant / "timing.json").exists():
            continue
         timing_rows.append(dict(case=case, variant=variant, **json.loads((base / variant / "timing.json").read_text())))
         cache = np.load(base / variant / "work" / "scatfile.npz")
         for l, channel in enumerate(channels):
            for r, re in enumerate(radius):
               n = min(int(cache["lmom"][m, l, r]), int(ref_cache["lmom"][m, l, r]))
               chi_r = ref_cache["amom"][:, m, l, r].astype(np.float64)
               chi_c = cache["amom"][:, m, l, r].astype(np.float64)
               k = min(chi_r.size, chi_c.size)
               optics_rows.append({
                  "case": case, "variant": variant, "channel": channel, "effective_radius_um": re,
                  "rel_bext": _rel(ref_cache["bext"][m, l, r], cache["bext"][m, l, r], 0.0),
                  "rel_bextrat": _rel(ref_cache["bextrat"][m, l, r], cache["bextrat"][m, l, r], 0.0),
                  "abs_ssa": float(abs(cache["w"][m, l, r] - ref_cache["w"][m, l, r])),
                  "abs_g": float(abs(cache["g"][m, l, r] - ref_cache["g"][m, l, r])),
                  "L": int(cache["lmom"][m, l, r]), "L_reference": int(ref_cache["lmom"][m, l, r]),
                  "chi_0_60": float(np.max(np.abs(chi_c[:61] - chi_r[:61]))),
                  "chi_61_up": float(np.max(np.abs(chi_c[61:k] - chi_r[61:k]))) if k > 61 else 0.0,
                  "common_L": n})
         lut = read_lut(next((base / variant).glob("*.nc")))
         for name in ref_lut.variable_names:
            a = np.asarray(ref_lut.variables[name])
            if a.dtype.kind != "f":
               continue
            b = np.asarray(lut.variables[name])
            dims = ref_lut.variable_dimensions[name]
            channel_axis = next((i for i, d in enumerate(dims) if "channel" in d), None)
            if channel_axis is None or a.ndim < 2:
               continue
            for c in range(a.shape[channel_axis]):
               x = np.take(a, c, axis=channel_axis).astype(np.float64)
               y = np.take(b, c, axis=channel_axis).astype(np.float64)
               finite = np.isfinite(x) & np.isfinite(y)
               if not np.any(finite):
                  continue
               lut_rows.append({"case": case, "variant": variant, "variable": name, "dimension": dims[channel_axis], "index": c,
                                "max_abs": float(np.max(np.abs(y[finite] - x[finite]))),
                                "max_rel": _rel(x[finite], y[finite], 1e-6 * max(1e-30, float(np.max(np.abs(x[finite]))))),
                                "max_value": float(np.max(np.abs(x[finite])))})
         # retrieval level: R_0v at cell midpoints in (log10 tau, r_e) and the Jacobians of the interpolant
         dims = ref_lut.variable_dimensions["R_0v"]
         a = np.asarray(ref_lut.variables["R_0v"], dtype=np.float64)
         b = np.asarray(lut.variables["R_0v"], dtype=np.float64)
         order = [dims.index("optical_depth"), dims.index("effective_radius")]
         rest = [i for i in range(a.ndim) if i not in order]
         a = np.transpose(a, order + rest)
         b = np.transpose(b, order + rest)
         channel_dim = dims[[i for i in rest if "channel" in dims[i]][0]]
         solar = [int(ch) for ch, flag in zip(inst.channelid, inst.solar_channel_flag) if flag]
         tq = cell_midpoints(tau, "log10")
         rq = cell_midpoints(radius, "linear")
         pad = [None] * (a.ndim - 2)
         va = interpolate_2d(tau, radius, a, tq, rq, "log10", "linear")
         vb = interpolate_2d(tau, radius, b, tq, rq, "log10", "linear")
         ja_t = np.diff(interpolate_2d(tau, radius, a, tau, rq, "log10", "linear"), axis=0) / np.diff(np.log10(tau))[:, None, *pad]
         jb_t = np.diff(interpolate_2d(tau, radius, b, tau, rq, "log10", "linear"), axis=0) / np.diff(np.log10(tau))[:, None, *pad]
         ja_r = np.diff(interpolate_2d(tau, radius, a, tq, radius, "log10", "linear"), axis=1) / np.diff(radius)[None, :, *pad]
         jb_r = np.diff(interpolate_2d(tau, radius, b, tq, radius, "log10", "linear"), axis=1) / np.diff(radius)[None, :, *pad]
         channel_axis = 2 + rest.index(dims.index(channel_dim))
         for c, channel in enumerate(solar):
            i = list(inst.channelid).index(channel)
            u = measurement_uncertainty(channel, wavelength_um=1.0, solar=True, thermal=False, oldnefr=float(inst.oldnefr[i]),
                                        snr=float(inst.snr[i]) or None, nedt=None, refbt=None)
            sigma = u.solar_reflectance_sigma
            take = lambda array: np.take(array, c, axis=channel_axis)
            dv = np.abs(take(vb) - take(va))
            djt = np.abs(take(jb_t) - take(ja_t))
            djr = np.abs(take(jb_r) - take(ja_r))
            retrieval_rows.append({
               "case": case, "variant": variant, "solar_channel": channel, "sigma_reflectance": sigma,
               "max_abs_dR_midpoint": float(np.max(dv)), "max_dR_over_sigma": float(np.max(dv) / sigma) if sigma else None,
               "max_rel_dR_midpoint": _rel(take(va), take(vb), 1e-6),
               "max_rel_dJ_logtau": _rel(take(ja_t), take(jb_t), 1e-3 * float(np.max(np.abs(take(ja_t))))),
               "max_rel_dJ_re": _rel(take(ja_r), take(jb_r), 1e-3 * float(np.max(np.abs(take(ja_r))))),
               "max_dJ_logtau_over_sigma": float(np.max(djt) / sigma) if sigma else None,
               "max_dJ_re_over_sigma_per_um": float(np.max(djr) / sigma) if sigma else None})
      timing_rows.append(dict(case=case, variant=reference, **json.loads((base / reference / "timing.json").read_text())))
   for filename, rows in (("optics.csv", optics_rows), ("lut_variables.csv", lut_rows),
                          ("retrieval_level.csv", retrieval_rows), ("timing.csv", timing_rows)):
      with open(RESULTS / f"{reference}_{filename}", "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
         writer.writeheader()
         writer.writerows(rows)
   finite = all(np.isfinite(v) for rows in (optics_rows, lut_rows) for row in rows for v in row.values() if isinstance(v, float))
   print("results in", RESULTS, "; all finite:", finite)


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   sub = parser.add_subparsers(dest="command", required=True)
   g = sub.add_parser("generate")
   g.add_argument("case", choices=sorted(CASES))
   g.add_argument("variant", choices=VARIANTS)
   c = sub.add_parser("compare")
   c.add_argument("reference", nargs="?", default="k6", choices=VARIANTS)
   args = parser.parse_args()
   if args.command == "generate":
      make_mirror()
      generate(args.case, args.variant)
   else:
      compare(args.reference)


if __name__ == "__main__":
   main()
