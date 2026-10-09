"""Figures for the V24 radius-grid report.

Validation-only code.  Reads results/convergence_matrix.csv,
results/compact_lut/k6_r0v_by_geometry.csv and one convergence case
(validation/tmp/radius_grid/matrix/liquid-water_stg_re5_wl0.6462791.npz, if
present) and writes PNG figures to results/figures.

   MPLCONFIGDIR=validation/generated/matplotlib-cache python validation/radius_grid/make_figures.py
"""

import collections
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
FIGURES = RESULTS / "figures"
CASE = ROOT / "validation" / "tmp" / "radius_grid" / "matrix" / "liquid-water_stg_re5_wl0.6462791.npz"


def convergence():
   rows = [r for r in csv.DictReader(open(RESULTS / "convergence_matrix.csv")) if r["rule"] == "trapezoid"]
   by = collections.defaultdict(list)
   for r in rows:
      by[int(r["level"])].append(r)
   levels = sorted(by)[:-1]
   metrics = {"extinction (rel.)": "rel_ext", "g (abs.)": "abs_g", "phase 5-150° (max rel.)": "phase_max_rel_side_5_150",
              "phase 150-180° (max rel.)": "phase_max_rel_back_150_180", "phase RMS (rel.)": "phase_rms_rel"}
   fig, ax = plt.subplots(figsize=(6.5, 4.5))
   for label, key in metrics.items():
      ax.semilogy(levels, [max(float(r[key]) for r in by[k]) for k in levels], "o-", label=label)
   ax.set_xlabel("refinement level k (size-parameter step 0.4 / 2^k)")
   ax.set_ylabel("worst case over 54 cases, vs level 6")
   ax.set_title("Radius-quadrature error of the size-averaged optics")
   ax.axvline(3, color="0.6", ls=":")
   ax.axvline(4, color="0.6", ls="--")
   ax.text(3.05, 2e-6, "ice spheres", rotation=90, color="0.4")
   ax.text(4.05, 2e-6, "liquid water", rotation=90, color="0.4")
   ax.grid(True, which="both", alpha=0.3)
   ax.legend(fontsize=8)
   fig.tight_layout()
   fig.savefig(FIGURES / "convergence_by_level.png", dpi=150)


def compact():
   rows = list(csv.DictReader(open(RESULTS / "compact_lut" / "k6_r0v_by_geometry.csv")))
   fig, ax = plt.subplots(figsize=(6.5, 4.5))
   for case, colour in (("modis_liquid", "C0"), ("slstr_liquid", "C1"), ("modis_ice_sph", "C2")):
      for region, style in (("glory_ge170", "o-"), ("other_lt170", "s--")):
         levels = list(range(6))
         values = [max(float(r[f"max_abs_dR_{region}"]) for r in rows if r["case"] == case and r["variant"] == f"k{k}")
                   for k in levels]
         ax.semilogy(levels, values, style, color=colour,
                     label=f"{case} {'Θ ≥ 170°' if region.startswith('glory') else 'Θ < 170°'}")
   ax.axhline(1.5e-3, color="k", lw=0.8)
   ax.text(0.1, 1.7e-3, "tolerance 1.5e-3 (0.3 σ SLSTR NeFR)", fontsize=8)
   ax.set_xlabel("refinement level k")
   ax.set_ylabel("max |ΔR_0v| over compact LUT nodes, vs level 6")
   ax.set_title("Compact LUT: bidirectional reflectance")
   ax.grid(True, which="both", alpha=0.3)
   ax.legend(fontsize=7)
   fig.tight_layout()
   fig.savefig(FIGURES / "compact_lut_r0v_by_level.png", dpi=150)


def phase_case():
   if not CASE.exists():
      return
   d = np.load(CASE)
   nq = int(d["nq"])
   theta = d["dense_theta"]
   ref = d["phase"][-1][nq:]
   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 6), sharex=True)
   ax1.semilogy(theta, ref, "k", lw=0.8)
   ax1.set_ylabel("F11 (level 6)")
   ax1.set_title("Liquid water, r_e = 5 µm, 0.646 µm")
   for k, colour in ((0, "C3"), (3, "C1"), (4, "C0")):
      ax2.semilogy(theta, np.abs(d["phase"][k][nq:] / ref - 1.0) + 1e-9, color=colour, lw=0.7,
                   label=f"level {k}" + (" (V23)" if k == 0 else ""))
   ax2.set_xlabel("scattering angle (deg)")
   ax2.set_ylabel("|relative error|")
   ax2.legend(fontsize=8)
   ax2.grid(True, which="both", alpha=0.3)
   fig.tight_layout()
   fig.savefig(FIGURES / "phase_error_re5_0.646um.png", dpi=150)


def main():
   FIGURES.mkdir(parents=True, exist_ok=True)
   convergence()
   compact()
   phase_case()
   print("figures in", FIGURES)


if __name__ == "__main__":
   main()
