"""Nakajima-King diagrams of the complete ORAC forward model (validation only).

Reads results/disort_options/fm_arrays_<case>_<surface>.npz (propagate.py) and
draws, in the measurement space ORAC compares with observations, the
production baseline (NSTR 60, solid) and the numerical reference
(double precision, NSTR 236, dashed) for a visible / non-absorbing channel
(x) against an absorbing channel (y), at one geometry.  Optical-depth
isolines are labelled by tau and effective-radius isolines by re (um); the
production measurement uncertainty (ORAC Sy: noise + homogeneity +
co-registration) is drawn as error bars at the reference nodes.  A second
panel shows the NSTR 60 - reference displacement in units of sigma.

   python nk_fm.py CASE SURFACE XCH YCH SZA VZA RAA
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import orac_fm as fm                                   # noqa: E402
import propagate as pg                                 # noqa: E402

FIGURES = pg.RESULTS / "figures"


def main():
   case, surf, xch, ych = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
   sza, vza, raa = (float(v) for v in sys.argv[5:8])
   d = np.load(pg.RESULTS / f"fm_arrays_{case}_{surf}.npz")
   node = d["state_kind"] == "node"
   tau = np.unique(d["state_tau"][node])
   re = np.unique(d["state_re"][node])
   i, j, k = (int(np.argmin(np.abs(d[n] - v))) for n, v in (("sza", sza), ("vza", vza), ("raa", raa)))
   const = pg.channel_constants(pg.lut(case, pg.REFERENCE), pg.PLATFORM[case][1])

   def grid(experiment, ch):
      y = d[f"{experiment}__{ch}__y"][node][:, i, j, k]
      return y.reshape(tau.size, re.size)

   fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
   ax = axes[0]
   for experiment, style, label in (("baseline", "-", "NSTR 60 (production)"), (pg.REFERENCE, "--", "reference: NSTR 236, double")):
      X, Y = grid(experiment, xch), grid(experiment, ych)
      for a in range(tau.size):
         ax.plot(X[a], Y[a], style, color="C0", lw=1)
      for b in range(re.size):
         ax.plot(X[:, b], Y[:, b], style, color="C3", lw=1)
      if experiment == pg.REFERENCE:
         for a in range(tau.size):
            ax.annotate(f"τ={tau[a]:g}", (X[a, -1], Y[a, -1]), fontsize=7, color="C0")
         for b in range(re.size):
            ax.annotate(f"rₑ={re[b]:g}", (X[-1, b], Y[-1, b]), fontsize=7, color="C3")
      ax.plot([], [], style, color="k", label=label)
   X, Y = grid(pg.REFERENCE, xch), grid(pg.REFERENCE, ych)
   sx = fm.sigma_solar(X, const[xch]["noise_rel"], True)
   sy = fm.sigma_solar(Y, const[ych]["noise_rel"], True)
   ax.errorbar(X.ravel(), Y.ravel(), xerr=sx.ravel(), yerr=sy.ravel(), fmt="none", ecolor="0.6", lw=0.6,
               label="production σ (ORAC Sy)")
   ax.set_xlabel(f"channel {xch} reflectance (ORAC measurement)")
   ax.set_ylabel(f"channel {ych} reflectance (ORAC measurement)")
   ax.set_title(f"{case}, {surf} surface, SZA {d['sza'][i]:g}°, VZA {d['vza'][j]:g}°, RAA {d['raa'][k]:g}°", fontsize=9)
   ax.legend(fontsize=7)
   ax.grid(alpha=0.3)
   ax = axes[1]
   B = {ch: grid("baseline", ch) for ch in (xch, ych)}
   rx = (B[xch] - X) / sx
   ry = (B[ych] - Y) / sy
   q = ax.quiver(np.log10(tau)[:, None] * np.ones_like(rx), np.ones_like(rx) * re[None, :], rx, ry,
                 np.hypot(rx, ry), angles="xy", scale_units="xy", scale=2.0, cmap="viridis")
   fig.colorbar(q, ax=ax, label="|NSTR 60 − reference| / σ (2-channel)")
   ax.set_xlabel("log10 τ")
   ax.set_ylabel("rₑ (µm)")
   ax.set_title("displacement of the production baseline, in units of σ (arrow: x = ch %d, y = ch %d)" % (xch, ych), fontsize=9)
   ax.grid(alpha=0.3)
   fig.tight_layout()
   FIGURES.mkdir(parents=True, exist_ok=True)
   out = FIGURES / f"nk_fm_{case}_{surf}_ch{xch}_ch{ych}_sza{d['sza'][i]:g}_vza{d['vza'][j]:g}_raa{d['raa'][k]:g}.png"
   fig.savefig(out, dpi=140)
   print("written", out)


if __name__ == "__main__":
   main()
