"""Radius-integration convergence study (V24 radius-grid investigation).

Validation-only code.  For one size distribution, wavelength and effective
radius it evaluates the preserved Mie kernel once on the finest of a family
of nested dyadic refinements of the legacy radius lattice

   r_j = r_l + j h_0 / 2^K,   j = 0 ... J,   r_l = 0.001 um

where h_0 is the spacing of the legacy 0.001-100 um grid
(npts = max(int(2 pi (100 - 0.001) / (lambda 0.4)), 200)), and the upper
limit r_u is exactly the production one (the liquid-water lattice node at or
beyond 3.5 r_e, create_bwgp.lattice_upper_radius, or 100 um otherwise).  Every
coarser level k < K is a subset of the finest, so one Mie evaluation gives the
trapezoid rule at every level k = 0 (the legacy production grid) ... K, plus
the Richardson (Simpson) combination (4 T_k - T_{k-1}) / 3.  Only the
numerical radius integration differs between levels: size distribution,
refractive index, Mie kernel, limits and phase-function definition are those
of production (mie_size_dist_new, verified bitwise at k = 0).

Quantities, each per level: extinction, scattering, single-scattering
albedo, asymmetry parameter, the size-averaged phase function F11 at
(a) the production Gauss-Legendre nodes of this phase function (adaptive
order, legendre_expansion) and (b) ~1400 dense independent angles, and the
Legendre coefficients chi_l = omega_l / (2l+1) of each level's phase function.

   PYTHONPATH=src python validation/radius_grid/radius_convergence.py \
       --microphysics liquid-water_stg.mm --radius 10 --wavelength 0.6462791 --levels 4 --out FILE.npz
"""

import argparse
import importlib
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"

bwgp = importlib.import_module("oraclut.idl_mirror.create_bwgp")
gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
le = importlib.import_module("oraclut.idl_mirror.legendre_expansion")
from oraclut.idl_mirror.load_mmdat import load_mmdat     # noqa: E402
from oraclut.optics.legacy_mie import mie_single_batch    # noqa: E402

DENSE_THETA = np.unique(np.concatenate((np.arange(0.0, 1.0, 0.005), np.arange(1.0, 30.0, 0.05),
                                        np.arange(30.0, 180.0001, 0.25))))
BATCH = 800


def legacy_spacing(wavelength):
   npts = max(int(2.0 * np.pi * (bwgp.MIE_RADIUS_UPPER - bwgp.MIE_RADIUS_LOWER) / wavelength / bwgp.MIE_XRES), 200)
   return (bwgp.MIE_RADIUS_UPPER - bwgp.MIE_RADIUS_LOWER) / (npts - 1)


def number_weights(r, effective_radius, variance):
   """mie_size_dist_new modified-gamma number density (without the quadrature weight)."""

   alpha = (1.0 - 3.0 * variance) / variance
   b = 1.0 / (effective_radius * variance)
   n = b ** (-alpha - 1.0) * math.gamma(alpha + 1.0)
   return r**alpha * np.exp(-b * r) / n


def run(mmfile, effective_radius, wavelength, levels):
   mmstr = load_mmdat(INPUTS / "microphysics" / mmfile, INPUTS)
   if mmstr.distname[0] != "modified_gamma" or mmstr.ncomp != 1:
      raise SystemExit("this study covers single-component modified-gamma (liquid water, ice spheres)")
   variance = float(mmstr.s[0])
   cm = complex(gsp._interpol_complex(mmstr.comp[0].cm, mmstr.comp[0].wl, [wavelength])[0])
   factor = gsp.radius_upper_factor(mmstr, 0)
   params, npts0 = bwgp.mie_integration_limits("modified_gamma", effective_radius, variance, 1.0 / wavelength, factor)
   upper = params[3]
   h0 = legacy_spacing(wavelength)
   intervals0 = int(round((upper - bwgp.MIE_RADIUS_LOWER) / h0))
   if npts0 is not None and intervals0 + 1 != npts0:
      raise SystemExit("lattice mismatch")
   if npts0 is None:                                   # legacy 0.001-100 um grid
      intervals0 = max(int(2.0 * np.pi * (upper - bwgp.MIE_RADIUS_LOWER) / wavelength / bwgp.MIE_XRES), 200) - 1
      h0 = (upper - bwgp.MIE_RADIUS_LOWER) / intervals0
   finest = intervals0 * 2**levels
   r = bwgp.MIE_RADIUS_LOWER + (upper - bwgp.MIE_RADIUS_LOWER) * np.arange(finest + 1) / finest
   # trapezoid weights of every level on the finest nodes
   weights = np.zeros((levels + 1, r.size))
   for k in range(levels + 1):
      step = 2 ** (levels - k)
      h = (upper - bwgp.MIE_RADIUS_LOWER) / (intervals0 * 2**k)
      weights[k, ::step] = h
      weights[k, 0] = weights[k, -1] = h / 2.0
   # angular nodes: the production Gauss-Legendre order of this phase function, plus dense angles
   degree = le.mie_phase_degree(2.0 * np.pi * upper / wavelength)
   nq = le.initial_quadrature_order(degree)
   x_gauss, w_gauss = le.gauss_legendre(nq)
   mu = np.concatenate((-x_gauss, np.cos(np.deg2rad(DENSE_THETA))))
   if mu.size > le.MIE_MAX_ANGLES:
      raise SystemExit("too many angles for the Mie kernel")
   area = number_weights(r, effective_radius, variance) * np.pi * r**2
   ext = np.zeros(levels + 1)
   sca = np.zeros(levels + 1)
   gsca = np.zeros(levels + 1)
   phase = np.zeros((levels + 1, mu.size))
   started = time.time()
   for start in range(0, r.size, BATCH):
      sl = slice(start, start + BATCH)
      mie = mie_single_batch(2.0 * np.pi * r[sl] / wavelength, cm, mu)
      wa = weights[:, sl] * area[None, sl]
      ext += wa @ mie["qext"]
      sca += wa @ mie["qsca"]
      gsca += wa @ (mie["g"] * mie["qsca"])
      phase += (wa * mie["qsca"][None, :]) @ mie["f11"]
   mie_seconds = time.time() - started
   phase /= sca[:, None]
   omega = le.legendre_coefficients(-x_gauss, w_gauss, phase[:, :nq].T)          # (nq, levels+1)
   # production consistency at k = 0
   b, s, g, phi, _ = bwgp.create_bwgp("modified_gamma", effective_radius, variance, [cm], [wavelength], mu[:16],
                                      radius_upper_factor=factor)
   check = {"rel_ext": abs(ext[0] / b[0] - 1.0), "abs_ssa": abs(sca[0] / ext[0] - s[0]), "abs_g": abs(gsca[0] / sca[0] - g[0]),
            "rel_phase": float(np.max(np.abs(phase[0, :16] / phi[:, 0] - 1.0)))}
   return dict(microphysics=mmfile, effective_radius=effective_radius, wavelength=wavelength, levels=levels,
               upper=upper, h0=h0, intervals0=intervals0, nodes=r.size, degree=degree, nq=nq,
               dense_theta=DENSE_THETA, ext=ext, sca=sca, gsca=gsca, phase=phase, omega=omega,
               mie_seconds=mie_seconds, check=np.array([check[k] for k in ("rel_ext", "abs_ssa", "abs_g", "rel_phase")]))


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("--microphysics", required=True)
   parser.add_argument("--radius", type=float, required=True)
   parser.add_argument("--wavelength", type=float, required=True)
   parser.add_argument("--levels", type=int, default=4)
   parser.add_argument("--out", type=Path, required=True)
   args = parser.parse_args()
   result = run(args.microphysics, args.radius, args.wavelength, args.levels)
   args.out.parent.mkdir(parents=True, exist_ok=True)
   np.savez(args.out, **result)
   print(f"{args.microphysics} re={args.radius} wl={args.wavelength}: nodes {result['nodes']}, Nq {result['nq']}, "
         f"Mie {result['mie_seconds']:.0f} s, production check {result['check']}", flush=True)


if __name__ == "__main__":
   main()
