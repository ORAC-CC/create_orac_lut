"""Summarise radius_convergence.py results: error of every level (and its Richardson value) against the finest.

   python validation/radius_grid/analyse_convergence.py RESULTS_DIR [--csv OUT.csv]

For each case and level k the metrics are, against the reference level
(the finest K, after checking that K-1 and K agree):

   extinction, scattering   relative difference
   ssa, g                   absolute difference
   phase                    max relative difference in the forward (0-5 deg),
                            side (5-150 deg) and back (150-180 deg) regions of
                            the dense angles, and the RMS relative difference
                            over 0-180 deg
   chi_l                    max |delta chi_l| over l <= 60, 61-999 and >= 1000
                            (Gauss-Legendre coefficients up to the degree bound)
"""

import argparse
import csv
from pathlib import Path

import numpy as np


def metrics(data, a, b):
   """Differences of quantity set a from set b (each a dict of per-level arrays selected)."""

   theta = data["dense_theta"]
   nq = int(data["nq"])
   rel_phase = np.abs(a["phase"][nq:] / b["phase"][nq:] - 1.0)
   regions = {"forward_0_5": theta < 5.0, "side_5_150": (theta >= 5.0) & (theta < 150.0), "back_150_180": theta >= 150.0}
   degree = int(data["degree"])
   l = np.arange(min(degree + 1, a["omega"].size))
   dchi = np.abs((a["omega"][l] - b["omega"][l]) / (2.0 * l + 1.0))
   out = {"rel_ext": abs(a["ext"] / b["ext"] - 1.0), "rel_sca": abs(a["sca"] / b["sca"] - 1.0),
          "abs_ssa": abs(a["sca"] / a["ext"] - b["sca"] / b["ext"]), "abs_g": abs(a["gsca"] / a["sca"] - b["gsca"] / b["sca"]),
          "phase_rms_rel": float(np.sqrt(np.mean(rel_phase**2)))}
   for name, mask in regions.items():
      out[f"phase_max_rel_{name}"] = float(rel_phase[mask].max())
   out["chi_0_60"] = float(dchi[:61].max())
   out["chi_61_999"] = float(dchi[61:1000].max()) if dchi.size > 61 else 0.0
   out["chi_1000_up"] = float(dchi[1000:].max()) if dchi.size > 1000 else 0.0
   return out


def level(data, k):
   return {"ext": data["ext"][k], "sca": data["sca"][k], "gsca": data["gsca"][k], "phase": data["phase"][k],
           "omega": data["omega"][:, k]}


def richardson(data, k):
   return {key: (4.0 * level(data, k)[key] - level(data, k - 1)[key]) / 3.0 for key in ("ext", "sca", "gsca", "phase", "omega")}


def summarise(path):
   data = np.load(path)
   kmax = int(data["levels"])
   reference = level(data, kmax)
   rows = []
   for k in range(kmax + 1):
      for label, values in (("trapezoid", level(data, k)), ("richardson", richardson(data, k) if k > 0 else None)):
         if values is None:
            continue
         row = {"microphysics": str(data["microphysics"]), "effective_radius_um": float(data["effective_radius"]),
                "wavelength_um": float(data["wavelength"]), "upper_um": float(data["upper"]), "level": k,
                "rule": label, "delta_x": float(2.0 * np.pi * data["h0"] / data["wavelength"] / 2**k),
                "delta_r_um": float(data["h0"] / 2**k), "nodes": int(data["intervals0"]) * 2**k + 1}
         row.update(metrics(data, values, reference))
         rows.append(row)
   # convergence of the reference itself: K-1 against K
   successive = metrics(data, level(data, kmax - 1), reference)
   return rows, successive


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("results", type=Path)
   parser.add_argument("--csv", type=Path)
   args = parser.parse_args()
   rows = []
   for path in sorted(args.results.glob("*.npz")):
      case_rows, successive = summarise(path)
      rows += case_rows
      print(f"{path.stem}: K-1 vs K  ext {successive['rel_ext']:.1e} ssa {successive['abs_ssa']:.1e} g {successive['abs_g']:.1e} "
            f"phase side {successive['phase_max_rel_side_5_150']:.1e} back {successive['phase_max_rel_back_150_180']:.1e} "
            f"chi {max(successive['chi_0_60'], successive['chi_61_999'], successive['chi_1000_up']):.1e}")
      for row in case_rows:
         if row["rule"] == "trapezoid":
            print(f"   k={row['level']} dx={row['delta_x']:.3f} dr={row['delta_r_um']:.4f}: ext {row['rel_ext']:.1e} "
                  f"ssa {row['abs_ssa']:.1e} g {row['abs_g']:.1e} phase fwd {row['phase_max_rel_forward_0_5']:.1e} "
                  f"side {row['phase_max_rel_side_5_150']:.1e} back {row['phase_max_rel_back_150_180']:.1e} "
                  f"rms {row['phase_rms_rel']:.1e} chi {max(row['chi_0_60'], row['chi_61_999'], row['chi_1000_up']):.1e}")
   if args.csv:
      args.csv.parent.mkdir(parents=True, exist_ok=True)
      with open(args.csv, "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
         writer.writeheader()
         writer.writerows(rows)


if __name__ == "__main__":
   main()
