"""Compact LUT comparison of liquid-cloud radius-integration limits (Phase 2B/2C).

Validation-only code.  It runs the production generator
(create_orac_luts.create_orac_cloud_lut) on the compact grid
compact_liquid_grid.lut for one instrument, twice:

  reference   the production size integration, 0.001-100 um
  candidate   the same with the modified-gamma upper limit 3 x effective radius
              (substituted at run time into create_bwgp's call of
              mie_size_dist_new; production source is not edited)
  lattice     the 3 x effective-radius limit on the reference node lattice
              (same spacing as the 0.001-100 um grid, ending at the first node
              >= 3 re, beyond 100 um when needed)
  adopted     the production code as it stands (after the change it
              implements the lattice limit itself; reference and the run-time
              variants then re-impose the legacy limits first)

Everything else (SRF treatment, Mie kernel, Legendre expansion, DISORT,
NetCDF writer) is the production path at the current source revision.  The
scattering cache of each run is kept so that the optical properties and the
1000 production Legendre moments can be compared as well as every LUT
variable.

   generate:  python compact_lut_comparison.py generate <modis|slstr> <reference|candidate>
   compare:   python compact_lut_comparison.py compare [candidate|lattice]

Inputs are read through a mirror of create_orac_lut/input_files in
validation/tmp/compact_inputs (symbolic links plus the compact grid), and all
products are written under validation/tmp/compact_luts.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

INPUTS = ROOT / "create_orac_lut" / "input_files"
MIRROR = ROOT / "validation" / "tmp" / "compact_inputs"
PRODUCTS = ROOT / "validation" / "tmp" / "compact_luts"
RESULTS_BASE = HERE / "results" / "compact_lut"
GRID = "compact_liquid_grid.lut"
MICROPHYSICS = "liquid-water_stg.mm"
UPPER_FACTOR = 3.5
INSTRUMENTS = {
   # MODIS 0.646, 0.858, 1.63, 2.13, 3.79 (solar + thermal) and 11.03 um
   "modis": ("terra_modis_v1.inst", [1, 2, 6, 7, 20, 31]),
   # SLSTR (dual view) 0.555, 0.865, 1.61, 2.25, 3.74 (solar + thermal) and 10.85 um
   "slstr": ("sentinel-3a_slstr_v1.inst", [1, 3, 5, 6, 7, 8]),
}


def make_mirror():
   (MIRROR / "lut").mkdir(parents=True, exist_ok=True)
   for item in INPUTS.iterdir():
      link = MIRROR / item.name
      if item.name != "lut" and not link.exists():
         link.symlink_to(item)
   (MIRROR / "lut" / GRID).write_text((HERE / GRID).read_text())


def generate(instrument, variant):
   import importlib
   import create_orac_luts
   # the module (oraclut.idl_mirror exports a function of the same name)
   create_bwgp_module = importlib.import_module("oraclut.idl_mirror.create_bwgp")

   gsp_module = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
   if variant != "adopted" and hasattr(gsp_module, "radius_upper_factor"):
      # From the production change onwards: impose the legacy 0.001-100 um
      # limits so that "reference" reproduces the V23 integration and the
      # run-time candidates start from it.
      gsp_module.radius_upper_factor = lambda mmstr, c: None
   if variant == "candidate":
      production = create_bwgp_module.mie_size_dist_new

      def adaptive_upper(distname, nd, params, wavenumber, cm, dqv, xres=0.1):
         if distname == "modified_gamma":
            params = [params[0], params[1], params[2], UPPER_FACTOR * params[0]]
         return production(distname, nd, params, wavenumber, cm, dqv, xres=xres)

      create_bwgp_module.mie_size_dist_new = adaptive_upper
   elif variant == "lattice":
      # The same upper limit, but on the reference node lattice: the spacing of
      # the legacy 0.001-100 um grid is kept and the grid stops at the first
      # lattice node >= 3 re (beyond 100 um when needed), so that the nodes up
      # to 100 um are those of the reference.
      production = create_bwgp_module.mie_size_dist_new
      quadrature = create_bwgp_module.quadrature
      override = {}

      def lattice_quadrature(quadtype, npts):
         return quadrature(quadtype, override.pop("npts", npts) if quadtype.upper() == "T" else npts)

      def lattice_upper(distname, nd, params, wavenumber, cm, dqv, xres=0.1):
         if distname == "modified_gamma":
            lower = params[2]
            reference_points = max(int(2.0 * np.pi * (params[3] - lower) * wavenumber / xres), 200)
            spacing = (params[3] - lower) / (reference_points - 1)
            intervals = int(np.ceil((UPPER_FACTOR * params[0] - lower) / spacing))
            params = [params[0], params[1], lower, lower + intervals * spacing]
            override["npts"] = intervals + 1
         return production(distname, nd, params, wavenumber, cm, dqv, xres=xres)

      create_bwgp_module.quadrature = lattice_quadrature
      create_bwgp_module.mie_size_dist_new = lattice_upper
   elif variant not in ("reference", "adopted"):
      raise ValueError(variant)

   timing = {}
   scattering = create_orac_luts.generate_scattering_properties

   def timed(*args, **kwargs):
      started = time.process_time()
      result = scattering(*args, **kwargs)
      timing["scattering_cpu_s"] = time.process_time() - started
      return result

   create_orac_luts.generate_scattering_properties = timed
   instfile, channels = INSTRUMENTS[instrument]
   out = PRODUCTS / instrument / variant
   out.mkdir(parents=True, exist_ok=True)
   started = time.process_time()
   status = create_orac_luts.create_orac_cloud_lut(
      MIRROR, instfile, MICROPHYSICS, GRID, out, 2, channelid=channels, srf_quad=1,
      version=None, nstreams=60, nmom=1000, work_path=out / "work")
   timing["total_cpu_s"] = time.process_time() - started
   timing["status"] = status
   (out / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
   print(instrument, variant, timing)


def _rel(a, b, floor):
   a = np.asarray(a, dtype=np.float64)
   b = np.asarray(b, dtype=np.float64)
   mask = np.abs(a) > floor
   return float(np.max(np.abs(b[mask] / a[mask] - 1.0))) if np.any(mask) else 0.0


def compare(variant="candidate"):
   from oraclut.io.lut import read_lut
   from oraclut.idl_mirror import load_inststr
   from oraclut.validation.quadrature import cell_midpoints, interpolate_2d, measurement_uncertainty

   RESULTS = RESULTS_BASE / variant
   RESULTS.mkdir(parents=True, exist_ok=True)
   optics_rows, lut_rows, retrieval_rows, timing_rows = [], [], [], []
   for instrument, (instfile, channels) in INSTRUMENTS.items():
      base = PRODUCTS / instrument
      ref_cache = np.load(base / "reference" / "work" / "scatfile.npz")
      cand_cache = np.load(base / variant / "work" / "scatfile.npz")
      inst = load_inststr(INPUTS / "inst" / instfile, requestedchannelid=channels)
      efr = np.array([5.0, 10.0, 20.0, 30.0, 40.0])
      m = ref_cache["bext"].shape[0] // 2
      for l, channel in enumerate(channels):
         for r, re in enumerate(efr):
            chi_r = ref_cache["amom"][:, m, l, r].astype(np.float64)
            chi_c = cand_cache["amom"][:, m, l, r].astype(np.float64)
            optics_rows.append({
               "instrument": instrument, "channel": channel, "effective_radius_um": re,
               "rel_bext": _rel(ref_cache["bext"][m, l, r], cand_cache["bext"][m, l, r], 0.0),
               "rel_bextrat": _rel(ref_cache["bextrat"][m, l, r], cand_cache["bextrat"][m, l, r], 0.0),
               "abs_ssa": float(abs(cand_cache["w"][m, l, r] - ref_cache["w"][m, l, r])),
               "abs_g": float(abs(cand_cache["g"][m, l, r] - ref_cache["g"][m, l, r])),
               "chi_0_127": float(np.max(np.abs(chi_c[:128] - chi_r[:128]))),
               "chi_128_999": float(np.max(np.abs(chi_c[128:] - chi_r[128:]))),
            })
      ref_lut = read_lut(next((base / "reference").glob("*.nc")))
      cand_lut = read_lut(next((base / variant).glob("*.nc")))
      for name in ref_lut.variable_names:
         a = np.asarray(ref_lut.variables[name])
         if a.dtype.kind != "f":
            continue
         b = np.asarray(cand_lut.variables[name])
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
            lut_rows.append({"instrument": instrument, "variable": name, "dimension": dims[channel_axis], "index": c,
                             "max_abs": float(np.max(np.abs(y[finite] - x[finite]))),
                             "max_rel": _rel(x[finite], y[finite], 1e-6 * max(1e-30, float(np.max(np.abs(x[finite]))))),
                             "max_value": float(np.max(np.abs(x[finite])))})
      # retrieval-level: R_0v at cell midpoints in (log10 tau, r_e) and its Jacobians
      dims = ref_lut.variable_dimensions["R_0v"]
      tau = np.asarray(ref_lut.variables["optical_depth"], dtype=np.float64)
      radius = np.asarray(ref_lut.variables["effective_radius"], dtype=np.float64)
      a = np.asarray(ref_lut.variables["R_0v"], dtype=np.float64)
      b = np.asarray(cand_lut.variables["R_0v"], dtype=np.float64)
      order = [dims.index("optical_depth"), dims.index("effective_radius")]
      rest = [i for i in range(a.ndim) if i not in order]
      a = np.transpose(a, order + rest)
      b = np.transpose(b, order + rest)
      channel_dim = dims[[i for i in rest if "channel" in dims[i]][0]]
      solar = [int(ch) for ch, flag in zip(inst.channelid, inst.solar_channel_flag) if flag]
      tq = cell_midpoints(tau, "log10")
      rq = cell_midpoints(radius, "linear")
      va = interpolate_2d(tau, radius, a, tq, rq, "log10", "linear")
      vb = interpolate_2d(tau, radius, b, tq, rq, "log10", "linear")
      # Jacobians of the piecewise-bilinear interpolant at the cell midpoints
      ja_t = np.diff(interpolate_2d(tau, radius, a, tau, rq, "log10", "linear"), axis=0) / np.diff(np.log10(tau))[:, None, *([None] * (a.ndim - 2))]
      jb_t = np.diff(interpolate_2d(tau, radius, b, tau, rq, "log10", "linear"), axis=0) / np.diff(np.log10(tau))[:, None, *([None] * (a.ndim - 2))]
      ja_r = np.diff(interpolate_2d(tau, radius, a, tq, radius, "log10", "linear"), axis=1) / np.diff(radius)[None, :, *([None] * (a.ndim - 2))]
      jb_r = np.diff(interpolate_2d(tau, radius, b, tq, radius, "log10", "linear"), axis=1) / np.diff(radius)[None, :, *([None] * (a.ndim - 2))]
      channel_axis = 2 + rest.index(dims.index(channel_dim))
      for c, channel in enumerate(solar):
         u = measurement_uncertainty(channel, wavelength_um=1.0, solar=True, thermal=False,
                                     oldnefr=float(inst.oldnefr[list(inst.channelid).index(channel)]),
                                     snr=float(inst.snr[list(inst.channelid).index(channel)]) or None, nedt=None, refbt=None)
         sigma = u.solar_reflectance_sigma
         take = lambda array: np.take(array, c, axis=channel_axis)
         dv = np.abs(take(vb) - take(va))
         djt = np.abs(take(jb_t) - take(ja_t))
         djr = np.abs(take(jb_r) - take(ja_r))
         retrieval_rows.append({
            "instrument": instrument, "solar_channel": channel, "sigma_reflectance": sigma,
            "max_abs_dR_midpoint": float(np.max(dv)), "max_dR_over_sigma": float(np.max(dv) / sigma) if sigma else None,
            "max_rel_dR_midpoint": _rel(take(va), take(vb), 1e-6),
            "max_rel_dJ_logtau": _rel(take(ja_t), take(jb_t), 1e-3 * float(np.max(np.abs(take(ja_t))))),
            "max_rel_dJ_re": _rel(take(ja_r), take(jb_r), 1e-3 * float(np.max(np.abs(take(ja_r))))),
            "max_dJ_logtau_over_sigma": float(np.max(djt) / sigma) if sigma else None,
            "max_dJ_re_over_sigma_per_um": float(np.max(djr) / sigma) if sigma else None,
         })
      for name in ("reference", variant):
         timing = json.loads((base / name / "timing.json").read_text())
         timing_rows.append(dict(instrument=instrument, variant=name, **timing))
   for filename, rows in (("optics.csv", optics_rows), ("lut_variables.csv", lut_rows),
                          ("retrieval_level.csv", retrieval_rows), ("timing.csv", timing_rows)):
      with open(RESULTS / filename, "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
         writer.writeheader()
         writer.writerows(rows)
   finite = all(np.isfinite(v) for rows in (optics_rows, lut_rows) for row in rows for v in row.values() if isinstance(v, float))
   print("results in", RESULTS, "; all finite:", finite)


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   sub = parser.add_subparsers(dest="command", required=True)
   g = sub.add_parser("generate")
   g.add_argument("instrument", choices=sorted(INSTRUMENTS))
   g.add_argument("variant", choices=["reference", "candidate", "lattice", "adopted"])
   c = sub.add_parser("compare")
   c.add_argument("variant", nargs="?", default="candidate", choices=["candidate", "lattice", "adopted"])
   sub.add_parser("all")
   args = parser.parse_args()
   if args.command == "generate":
      make_mirror()
      generate(args.instrument, args.variant)
   elif args.command == "compare":
      compare(args.variant)
   else:
      make_mirror()
      jobs = [subprocess.Popen([sys.executable, __file__, "generate", instrument, variant])
              for instrument in INSTRUMENTS for variant in ("reference", "candidate")]
      if any(job.wait() for job in jobs):
         raise SystemExit("a generation run failed")
      compare()


if __name__ == "__main__":
   main()
