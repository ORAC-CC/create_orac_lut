"""Angular extent of the backscatter error in the complete ORAC measurement (validation only).

Reads validation/tmp/disort_options/fm_full/<case>_bs_<surface>.npz, written by

   python propagate.py --cases modis_liquid_bs,modis_ice_agg_bs,slstr_liquid_bs,slstr_ice_sph_bs \
                       --out validation/v24/results/disort_options/near_backscatter

for the near-backscatter geometry sets (SZA, VZA = 30, 33.75 deg; RAA 0-180 every
5 deg; Grid B nodes), and writes to results/disort_options/near_backscatter/:

   near_backscatter.csv  max and 99th-percentile r = |dy|/sigma (production ORAC Sy)
                         against the reference (double, NSTR 236) per case, channel
                         kind, surface, experiment and scattering-angle bin
   extent.csv            per case, kind and experiment: the smallest scattering
                         angle at which r exceeds 0.5, 1 and 2 (worst over
                         channels, surfaces and states)
   ../figures/near_backscatter.png
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyse_fm as af                                # noqa: E402
import propagate as pg                                 # noqa: E402

OUT = pg.RESULTS / "near_backscatter"
CASES = ("modis_liquid_bs", "modis_ice_agg_bs", "slstr_liquid_bs", "slstr_ice_sph_bs")
EXPERIMENTS = ("baseline", "nstr100", "nstr136", "nstr200", "nstr236", "double_nstr136")
NSTR = {"baseline": 60, "double_nstr136": 136}
BINS = [(150, 160), (160, 165), (165, 170), (170, 172), (172, 174), (174, 176), (176, 177), (177, 178), (178, 179),
        (179, 179.9), (179.9, 180.01)]
SURFACES = ("dark", "vegetation", "snow")


def main():
   OUT.mkdir(parents=True, exist_ok=True)
   rows, extent, curves = [], [], {}
   for case in CASES:
      for surf in SURFACES:
         path = af.FULL / f"{case}_{surf}.npz"
         if not path.exists():
            print("missing", path)
            continue
         d = np.load(path)
         kinds = dict(k.split(":") for k in d["kinds"])
         theta = d["theta"]
         for ch_s, kind in kinds.items():
            if kind == "thermal":
               continue
            ch = int(ch_s)
            y_ref = d[f"{pg.REFERENCE}__{ch}__y"].astype(np.float64)
            sig = af.sigma(case[:-3], ch, kind, y_ref)
            t = np.broadcast_to(theta, y_ref.shape)
            for e in EXPERIMENTS:
               key = f"{e}__{ch}__y"
               if key not in d.files:
                  continue
               r = np.abs(d[key] - y_ref) / sig
               ok = np.isfinite(r)
               for lo, hi in BINS:
                  m = ok & (t >= lo) & (t < hi)
                  if not m.any():
                     continue
                  rows.append({"case": case, "surface": surf, "channel": ch, "kind": kind, "experiment": e,
                               "nstreams": NSTR.get(e, int(e[4:]) if e.startswith("nstr") else None),
                               "precision": "double" if e.startswith("double") else "single",
                               "theta_lo": lo, "theta_hi": hi, "points": int(m.sum()),
                               "max_r": float(r[m].max()), "p99_r": float(np.percentile(r[m], 99))})
               # worst r as a function of scattering angle, per distinct angle
               for ang in np.unique(np.round(theta, 2)):
                  m = ok & (np.round(t, 2) == ang)
                  k = (case, kind, e)
                  curves.setdefault(k, {})
                  curves[k][float(ang)] = max(curves[k].get(float(ang), 0.0), float(r[m].max()))
   for (case, kind, e), c in sorted(curves.items()):
      ang = np.array(sorted(c))
      val = np.array([c[a] for a in ang])
      row = {"case": case, "kind": kind, "experiment": e}
      for level in (0.5, 1.0, 2.0):
         above = ang[val > level]
         row[f"theta_first_r_gt_{level:g}"] = float(above.min()) if above.size else None
      row["max_r_theta_lt_170"] = float(val[ang < 170].max()) if (ang < 170).any() else None
      row["max_r_theta_170_178"] = float(val[(ang >= 170) & (ang < 178)].max()) if ((ang >= 170) & (ang < 178)).any() else None
      row["max_r_theta_ge_178"] = float(val[ang >= 178].max())
      extent.append(row)
   for name, table in (("near_backscatter.csv", rows), ("extent.csv", extent)):
      with open(OUT / name, "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(table[0]))
         writer.writeheader()
         writer.writerows(table)

   fig, axes = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)
   for ax, kind in zip(axes, ("solar", "mixed")):
      for n_case, case in enumerate(CASES):
         for e, style in (("baseline", "-"), ("nstr136", "--"), ("nstr236", ":")):
            c = curves.get((case, kind, e))
            if not c:
               continue
            ang = np.array(sorted(c))
            ax.semilogy(ang, [c[a] for a in ang], style, color=f"C{n_case}", marker=".", ms=3,
                        label=f"{case[:-3]} NSTR {NSTR.get(e, e[4:])}")
      for level, lw in ((2.0, 0.8), (0.5, 0.6), (0.1, 0.5)):
         ax.axhline(level, color="k", lw=lw, ls=":" if level != 2.0 else "-")
      ax.set_xlim(150, 180.5)
      ax.set_xlabel("scattering angle Θ (deg)")
      ax.set_title({"solar": "solar reflectance (dark, vegetation, snow)", "mixed": "mixed 3.7 µm BT"}[kind], fontsize=10)
      ax.grid(True, which="both", alpha=0.3)
   axes[0].set_ylabel("max |Δy| / σ  vs reference (double, NSTR 236)")
   axes[0].legend(fontsize=6, ncol=2)
   fig.tight_layout()
   fig.savefig(pg.RESULTS / "figures" / "near_backscatter.png", dpi=140)
   print("written", OUT)


if __name__ == "__main__":
   main()
