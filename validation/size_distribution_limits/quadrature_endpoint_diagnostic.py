"""Radius-quadrature endpoint diagnostic for modified-gamma liquid clouds (Phase 2A).

Validation-only code.  It separates the two effects of changing the upper
radius limit of the production size integration
(src/oraclut/idl_mirror/create_bwgp.py: mie_size_dist_new):

* the physical tail removed beyond the new limit, and
* the numerical change caused by moving every node of the linear-radius
  trapezoidal rule, whose spacing and node positions both depend on the
  endpoint (npts = max(int(2 pi (ru - rl) / (lambda xres)), 200), xres = 0.4).

For each wavelength and effective radius it computes, with the preserved Mie
kernel and the production distribution weights:

  truth       dense trapezoid, xres = 0.05, 0.001-150 um (converged tail and nodes)
  reference   production grid, 0.001-100 um
  candidate   production grid, 0.001-3 re
  jitter      production grids ending at 100 um x (1 +- 0.0025, 1 +- 0.005);
              the tail beyond 99.5 um is negligible here, so their spread is
              pure node-placement noise
  lattice     the reference node lattice (same spacing) up to the first node
              >= 3 re, extended beyond 100 um when needed; up to 100 um its
              nodes are those of the reference, so for 3 re <= 100 um its
              difference from the reference is the removed tail alone

Bulk extinction, scattering, single-scattering albedo, asymmetry parameter,
the phase function and its Legendre coefficients are compared.  Coefficients
are computed on common Gauss-Legendre nodes, Newton-refined so that their
weights are accurate (numpy's leggauss weights lose accuracy at these orders),
with an order above the Mie degree bound of the widest (150 um) integration so
that no coefficient is aliased.

Usage, from the repository root (about 15 minutes with nice):

   PYTHONPATH=src python validation/size_distribution_limits/quadrature_endpoint_diagnostic.py
"""

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

from oraclut.idl_mirror.create_bwgp import create_bwgp, quadrature, shift_quadrature   # noqa: E402
from oraclut.idl_mirror.generate_scattering_properties import _interpol_complex       # noqa: E402
from oraclut.idl_mirror.load_mmdat import load_mmdat                                  # noqa: E402
from oraclut.optics.legacy_mie import mie_single_batch                                # noqa: E402

INPUTS = ROOT / "create_orac_lut" / "input_files"
EFFECTIVE_VARIANCE = 0.1111111
RADII = (5.0, 10.0, 20.0, 30.0, 40.0)
# MODIS 0.646 / 1.628 / 2.114 / 3.785 / 11.026 um and SLSTR 0.868 um centre wavelengths (as float32 in the .inst files)
WAVELENGTHS = (0.6462791, 0.8680, 1.6291, 2.1142, 3.7850, 11.0262)
XRES = 0.4
TRUTH_XRES = 0.05
TRUTH_UPPER = 150.0
LOWER = 0.001
UPPER = 100.0
UPPER_FACTOR = 3.5
JITTER = (-0.005, -0.0025, 0.0025, 0.005)
BATCH = 1500


def gauss_legendre(npts):
   """numpy leggauss nodes refined by Newton iteration; weights from the derivative."""

   x = np.polynomial.legendre.leggauss(npts)[0]
   for iteration in range(3):
      pnm1, pn = np.ones_like(x), x.copy()
      for n in range(2, npts + 1):
         pnm1, pn = pn, ((2.0 * n - 1.0) * x * pn - (n - 1.0) * pnm1) / n
      dpn = npts * (pnm1 - x * pn) / ((1.0 - x) * (1.0 + x))
      x = x - pn / dpn
   pnm1, pn = np.ones_like(x), x.copy()
   for n in range(2, npts + 1):
      pnm1, pn = pn, ((2.0 * n - 1.0) * x * pn - (n - 1.0) * pnm1) / n
   dpn = npts * (pnm1 - x * pn) / ((1.0 - x) * (1.0 + x))
   return x, 2.0 / ((1.0 - x) * (1.0 + x) * dpn**2)


def nstop(x):
   return 2 if x < 0.02 else int(x + 4.05 * x**(1.0 / 3.0) + 2.0) + 1


def production_grid(lower, upper, wavelength, xres):
   """The mie_size_dist_new radius nodes and trapezoid weights (create_bwgp.py)."""

   wavenumber = 1.0 / wavelength
   npts = max(int(2.0 * np.pi * (upper - lower) * wavenumber / xres), 200)
   absc, wght = quadrature("T", npts)
   return shift_quadrature(absc, wght, lower, upper)


def modified_gamma_weights(r, w1, effective_radius):
   """Production number weights W1P (create_bwgp.py, mie_size_dist_new, modified gamma)."""

   alpha = (1.0 - 3.0 * EFFECTIVE_VARIANCE) / EFFECTIVE_VARIANCE
   b = 1.0 / (effective_radius * EFFECTIVE_VARIANCE)
   n = b ** (-alpha - 1.0) * math.gamma(alpha + 1.0)
   return w1 / n * r**alpha * np.exp(-b * r)


def integrate(r, w1, wavelength, cm, mu, radii):
   """Bulk properties and phase functions for every effective radius on one radius grid.

   The Mie calculation is shared by all effective radii (the grid does not
   depend on them); the sums are those of mie_size_dist_new.
   """

   nre = len(radii)
   ext = np.zeros(nre)
   sca = np.zeros(nre)
   gsca = np.zeros(nre)
   phase = np.zeros((nre, mu.size))
   weights = np.stack([modified_gamma_weights(r, w1, re) for re in radii])     # (nre, npts)
   area = weights * np.pi * r**2
   for start in range(0, r.size, BATCH):
      sl = slice(start, start + BATCH)
      mie = mie_single_batch(2.0 * np.pi * r[sl] / wavelength, cm, mu)
      ext += area[:, sl] @ mie["qext"]
      sca += area[:, sl] @ mie["qsca"]
      gsca += area[:, sl] @ (mie["g"] * mie["qsca"])
      phase += (area[:, sl] * mie["qsca"][None, :]) @ mie["f11"]
   return {"ext": ext, "sca": sca, "ssa": sca / ext, "g": gsca / sca, "phase": phase / sca[:, None]}


def coefficients(phase, mu, weights):
   """Legendre coefficients omega_l (Grainger 1990 eq. 4.4.2) for each row of phase."""

   npts = mu.size
   weighted = phase * weights[None, :]
   omega = np.zeros((phase.shape[0], npts))
   pnm1, pn = np.ones_like(mu), mu.copy()
   omega[:, 0] = weighted @ pnm1 / 2.0
   omega[:, 1] = 3.0 * (weighted @ pn) / 2.0
   for n in range(2, npts):
      pnm1, pn = pn, ((2.0 * n - 1.0) / n) * mu * pn - ((n - 1.0) / n) * pnm1
      omega[:, n] = (2.0 * n + 1.0) * (weighted @ pn) / 2.0
   return omega


def compare(result, truth, omega, omega_truth, index, mu):
   """Differences of one variant from the truth for effective-radius index ``index``."""

   chi = omega[index] / (2.0 * np.arange(omega.shape[1]) + 1.0)
   chi_truth = omega_truth[index] / (2.0 * np.arange(omega.shape[1]) + 1.0)
   dchi = np.abs(chi - chi_truth)
   theta = np.degrees(np.arccos(mu))
   relp = np.abs(result["phase"][index] / truth["phase"][index] - 1.0)
   return {
      "rel_ext": abs(result["ext"][index] / truth["ext"][index] - 1.0),
      "rel_sca": abs(result["sca"][index] / truth["sca"][index] - 1.0),
      "abs_ssa": abs(result["ssa"][index] - truth["ssa"][index]),
      "abs_g": abs(result["g"][index] - truth["g"][index]),
      "phase_rel_30_180": float(np.max(relp[theta >= 30.0])),
      "phase_rel_all": float(np.max(relp)),
      "chi_0_127": float(np.max(dchi[:128])),
      "chi_128_999": float(np.max(dchi[128:1000])) if dchi.size > 128 else 0.0,
      "chi_1000_up": float(np.max(dchi[1000:])) if dchi.size > 1000 else 0.0,
   }


def run(results, wavelengths, radii):
   ri_table = load_mmdat(INPUTS / "microphysics" / "liquid-water_stg.mm", INPUTS).comp[0]
   rows = []
   check = None
   for wavelength in wavelengths:
      started = time.time()
      cm = complex(_interpol_complex(ri_table.cm, ri_table.wl, [wavelength])[0])
      degree = 2 * nstop(2.0 * np.pi * TRUTH_UPPER / wavelength)
      npts = degree + 1 + 32 + degree // 32
      x, w = gauss_legendre(npts)
      mu = -x
      variants = {}
      r, w1 = production_grid(LOWER, TRUTH_UPPER, wavelength, TRUTH_XRES)
      variants["truth"] = (integrate(r, w1, wavelength, cm, mu, radii), r.size)
      r_ref, w1_ref = production_grid(LOWER, UPPER, wavelength, XRES)
      variants["reference"] = (integrate(r_ref, w1_ref, wavelength, cm, mu, radii), r_ref.size)
      for jitter in JITTER:
         r, w1 = production_grid(LOWER, UPPER * (1.0 + jitter), wavelength, XRES)
         variants[f"jitter{jitter:+.4f}"] = (integrate(r, w1, wavelength, cm, mu, radii), r.size)
      # candidate (own production grid) and lattice (reference-node subset), one radius at a time
      for j, re in enumerate(radii):
         r, w1 = production_grid(LOWER, UPPER_FACTOR * re, wavelength, XRES)
         variants[f"candidate{j}"] = (integrate(r, w1, wavelength, cm, mu, [re]), r.size)
         h = (UPPER - LOWER) / (r_ref.size - 1)
         intervals = int(np.ceil((UPPER_FACTOR * re - LOWER) / h))
         r_lat = LOWER + h * np.arange(intervals + 1)        # reference nodes, extended beyond 100 um when needed
         w_lat = np.full(r_lat.size, h)
         w_lat[0] = w_lat[-1] = h / 2.0
         variants[f"lattice{j}"] = (integrate(r_lat, w_lat, wavelength, cm, mu, [re]), r_lat.size)
      # one consistency check against the production routine itself
      if check is None:
         bext, ssa, g, phi, _ = create_bwgp("modified_gamma", radii[0], EFFECTIVE_VARIANCE, [cm], [wavelength], mu[:64])
         ref = integrate(r_ref, w1_ref, wavelength, cm, mu[:64], radii[:1])
         check = {"wavelength": wavelength, "re": radii[0],
                  "rel_ext": float(abs(ref["ext"][0] / bext[0] - 1.0)), "abs_ssa": float(abs(ref["ssa"][0] - ssa[0])),
                  "abs_g": float(abs(ref["g"][0] - g[0])),
                  "rel_phase": float(np.max(np.abs(ref["phase"][0] / phi[:, 0] - 1.0)))}
      omega = {name: coefficients(result["phase"], mu, w) for name, (result, count) in variants.items()}
      truth = variants["truth"][0]
      for j, re in enumerate(radii):
         base = {"wavelength_um": wavelength, "effective_radius_um": re, "legendre_points": npts, "degree_bound_150um": degree,
                 "truth_radius_points": variants["truth"][1]}
         for name in ["reference", "candidate", "lattice"] + [f"jitter{jitter:+.4f}" for jitter in JITTER]:
            key = name + (str(j) if name in ("candidate", "lattice") else "")
            result, count = variants[key]
            index = 0 if name in ("candidate", "lattice") else j
            row = dict(base, variant=name, radius_points=count,
                       upper_um=float(result_upper(name, re, wavelength)))
            row.update(compare(result, truth, omega[key], omega["truth"], index, mu) if name not in ("candidate", "lattice")
                       else compare(result, slice_truth(truth, j), omega[key], omega["truth"][j:j + 1], 0, mu))
            rows.append(row)
         # candidate and lattice relative to the reference (what a production change would show)
         for name in ("candidate", "lattice"):
            result, count = variants[name + str(j)]
            row = dict(base, variant=name + "_minus_reference", radius_points=count, upper_um=UPPER_FACTOR * re)
            row.update(compare(result, slice_truth(variants["reference"][0], j), omega[name + str(j)],
                               omega["reference"][j:j + 1], 0, mu))
            rows.append(row)
      print(f"{wavelength:.4f} um: Nq = {npts}, truth points = {variants['truth'][1]}, "
            f"reference points = {variants['reference'][1]}, {time.time() - started:.0f} s", flush=True)
   results.mkdir(parents=True, exist_ok=True)
   with open(results / "quadrature_endpoint_diagnostic.csv", "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)
   (results / "quadrature_endpoint_check.json").write_text(json.dumps(check, indent=2) + "\n")
   return rows, check


def slice_truth(truth, j):
   return {"ext": truth["ext"][j:j + 1], "sca": truth["sca"][j:j + 1], "ssa": truth["ssa"][j:j + 1],
           "g": truth["g"][j:j + 1], "phase": truth["phase"][j:j + 1]}


def result_upper(name, re, wavelength):
   if name == "reference":
      return UPPER
   if name in ("candidate", "lattice"):
      return UPPER_FACTOR * re
   return UPPER * (1.0 + float(name[len("jitter"):]))


def main():
   global UPPER_FACTOR
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("--results", type=Path, default=HERE / "results" / "quadrature_endpoint")
   parser.add_argument("--wavelengths", type=float, nargs="*", default=list(WAVELENGTHS))
   parser.add_argument("--radii", type=float, nargs="*", default=list(RADII))
   parser.add_argument("--factor", type=float, default=UPPER_FACTOR, help="upper limit as a multiple of effective radius")
   args = parser.parse_args()
   UPPER_FACTOR = args.factor
   rows, check = run(args.results, args.wavelengths, args.radii)
   print("production-routine consistency:", check)
   finite = all(np.isfinite(v) for row in rows for v in row.values() if isinstance(v, float))
   print("all values finite:", finite)


if __name__ == "__main__":
   main()
