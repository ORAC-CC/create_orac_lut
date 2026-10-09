"""Figures for REPORT_disort_options.md (validation only).

Reads the saved forward-model and timing tables in
validation/v24/results/disort_options (analyse_fm.py, propagate.py, timing.py)
and writes PNGs to validation/v24/results/disort_options/figures:

   fm_stream_convergence.png  max r = |dy|/sigma (production ORAC Sy) of the complete
                              ORAC measurement against the reference (double, NSTR 236),
                              by NSTR, for solar, mixed (3.7 um) and thermal channels;
                              away from backscatter (< 170 deg) and at exact backscatter
   fm_options.png             max r of each NSTR 60 option experiment against the
                              production baseline, away from backscatter and overall
   runtime_vs_streams.png     DISORT CPU time relative to NSTR 60 (repeated timings)

   MPLCONFIGDIR=validation/generated/matplotlib-cache python make_figures.py
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "validation" / "v24" / "results" / "disort_options"
FIGURES = RESULTS / "figures"
CASES = ("modis_liquid", "modis_ice_agg", "slstr_liquid", "slstr_ice_sph")
SURFACES = ("dark", "vegetation", "snow")
OPTIONS = ("accur0.01", "accur0", "nocorint", "nodeltam_nocorint", "moments_nstr1", "double", "nstr056", "nstr068")


def table(name):
   return list(csv.DictReader(open(RESULTS / name)))


def worst(rows, key):
   values = [float(r[key]) for r in rows if r[key] not in ("", "nan")]
   return max(values) if values else np.nan


def convergence():
   ladder = table("fm_ladder.csv")
   fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
   for ax, kind in zip(axes, ("solar", "mixed", "thermal")):
      for n_case, case in enumerate(CASES):
         colour = f"C{n_case}"
         for angles, style in (("off", "-"), ("exact", "--"), ("all", "-")):
            if (kind == "thermal") != (angles == "all"):
               continue
            sel = [r for r in ladder if r["case"] == case and r["kind"] == kind and r["angles"] == angles and r["surface"] in SURFACES]
            single = sorted({int(r["nstreams"]) for r in sel if r["precision"] == "single" and r["nstreams"]})
            if not single:
               continue
            y = [worst([r for r in sel if r["precision"] == "single" and r["nstreams"] == str(n)], "max_r") for n in single]
            label = f"{case}" + ("" if angles != "exact" else " (exact backscatter)")
            ax.semilogy(single, y, style, marker="o", ms=3, color=colour, label=label)
            for n in (60, 136):
               d = [r for r in sel if r["precision"] == "double" and r["nstreams"] == str(n)]
               if d:
                  ax.semilogy([n], [worst(d, "max_r")], "s", mfc="none", color=colour, ms=7)
      for level, lw in ((2.0, 0.8), (0.5, 0.6), (0.1, 0.5)):
         ax.axhline(level, color="k", lw=lw, ls=":" if level != 2.0 else "-")
      ax.axvline(60, color="0.5", ls="--")
      ax.set_xlabel("NSTR (open squares: double precision)")
      ax.set_title({"solar": "solar reflectance (dark, vegetation, snow)", "mixed": "mixed 3.7 µm BT",
                    "thermal": "thermal BT"}[kind], fontsize=10)
      ax.grid(True, which="both", alpha=0.3)
   axes[0].set_ylabel("max |Δy| / σ  vs reference (double, NSTR 236)")
   axes[0].legend(fontsize=7)
   fig.tight_layout()
   fig.savefig(FIGURES / "fm_stream_convergence.png", dpi=140)


def options():
   meas = [r for r in table("fm_measurement.csv") if r["against"] == "baseline" and r["surface"] in SURFACES]
   fig, axes = plt.subplots(1, 2, figsize=(14, 4.6), sharey=True)
   x = np.arange(len(OPTIONS))
   for ax, key, title in ((axes[0], "max_r_off_glory", "away from backscatter (< 170°)"), (axes[1], "max_r", "all geometries")):
      for n_case, case in enumerate(CASES):
         for kind, marker in (("solar", "o"), ("mixed", "^"), ("thermal", "s")):
            y = [worst([r for r in meas if r["case"] == case and r["kind"] == kind and r["experiment"] == e], key) for e in OPTIONS]
            ax.semilogy(x + 0.08 * (n_case - 1.5), np.maximum(y, 1e-4), marker, color=f"C{n_case}", ls="none",
                        label=f"{case} {kind}" if ax is axes[0] else None)
      for level, lw in ((2.0, 0.8), (0.5, 0.6), (0.1, 0.5)):
         ax.axhline(level, color="k", lw=lw, ls=":" if level != 2.0 else "-")
      ax.set_xticks(x, OPTIONS, rotation=30, ha="right")
      ax.set_title(title, fontsize=10)
      ax.grid(True, which="both", alpha=0.3)
   axes[0].set_ylabel("max |Δy| / σ  vs production baseline (floor 1e-4)")
   axes[0].legend(fontsize=6, ncol=2)
   fig.tight_layout()
   fig.savefig(FIGURES / "fm_options.png", dpi=140)


def runtime():
   rows = table("runtime.csv")
   fig, ax = plt.subplots(figsize=(6.4, 4.6))
   for n_case, case in enumerate(CASES):
      pts = sorted((int(r["nstreams"]), float(r["ratio_to_baseline"])) for r in rows if r["case"] == case and
                   (r["experiment"] == "baseline" or r["experiment"].startswith("nstr") and "_" not in r["experiment"]))
      ax.loglog([p[0] for p in pts], [p[1] for p in pts], "o-", color=f"C{n_case}", label=case)
      d = [float(r["ratio_to_baseline"]) for r in rows if r["case"] == case and r["experiment"] == "double"]
      if d:
         ax.loglog([60], d, "s", mfc="none", color=f"C{n_case}", ms=8)
   n = np.array([30, 140])
   ax.loglog(n, (n / 60.0) ** 3, "k:", lw=0.8, label="(NSTR/60)³")
   ax.axvline(60, color="0.5", ls="--")
   ax.set_xlabel("NSTR (open squares: double precision at 60)")
   ax.set_ylabel("DISORT CPU time / production NSTR 60 (median of 3)")
   ax.grid(True, which="both", alpha=0.3)
   ax.legend(fontsize=8)
   fig.tight_layout()
   fig.savefig(FIGURES / "runtime_vs_streams.png", dpi=140)


def main():
   FIGURES.mkdir(parents=True, exist_ok=True)
   convergence()
   options()
   runtime()
   print("figures in", FIGURES)


if __name__ == "__main__":
   main()
