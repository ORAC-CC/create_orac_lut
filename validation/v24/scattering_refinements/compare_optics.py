"""Compare the V23 and V24 class scattering properties written by optics.py (read-only).

   python validation/v24/scattering_refinements/compare_optics.py <model> ... -> results/optics_<model>.json

For each channel (SRF centre wavelength) and effective radius: relative change
of the extinction coefficient and of the extinction ratio (the quantity the
LUT scales optical depth with), absolute change of the single-scattering
albedo and asymmetry parameter, the number of Legendre moments (1000 in V23,
L in V24) and the largest change of the first 1000 moments.  Physical
consistency of V24: 0 <= w <= 1, -1 < g < 1, chi_0 = 1, |chi_l| <= 1 and
chi_1 = g (the first moment of the expanded phase function and the Mie
asymmetry parameter are independent calculations).
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main(models):
   for model in models:
      a = np.load(HERE / "work" / f"optics_v23_{model}.npz")
      b = np.load(HERE / "work" / f"optics_v24_{model}.npz")
      if not (np.array_equal(a["channels"], b["channels"]) and np.array_equal(a["efr"], b["efr"])):
         raise SystemExit(f"{model}: V23 and V24 runs used different channels or radii")
      m = 0
      rows = []
      for l, channel in enumerate(a["channels"]):
         for r, radius in enumerate(a["efr"]):
            chi23, chi24 = a["amom"][:, m, l, r], b["amom"][:, m, l, r]
            n = min(1000, chi24.size)
            rows.append({"channel": int(channel), "wavelength": float(a["wavelengths"][l]), "effective_radius": float(radius),
                         "extinction_rel": float(b["bext"][m, l, r] / a["bext"][m, l, r] - 1.0),
                         "extinction_ratio_rel": float(b["bextrat"][m, l, r] / a["bextrat"][m, l, r] - 1.0),
                         "w_abs": float(b["w"][m, l, r] - a["w"][m, l, r]), "g_abs": float(b["g"][m, l, r] - a["g"][m, l, r]),
                         "L_v24": int(b["lmom"][m, l, r]), "max_chi_change_l_lt_1000": float(np.max(np.abs(chi24[:n] - chi23[:n]))),
                         "v24_chi1_minus_g": float(chi24[1] - b["g"][m, l, r])})
      consistency = {"w_in_0_1": bool(np.all((b["w"] >= 0) & (b["w"] <= 1))), "g_in_-1_1": bool(np.all(np.abs(b["g"]) < 1)),
                     "chi0_is_1": bool(np.all(b["amom"][0] == 1.0)), "max_abs_chi": float(np.max(np.abs(b["amom"]))),
                     "max_abs_chi1_minus_g": float(max(abs(row["v24_chi1_minus_g"]) for row in rows)),
                     "reference_550nm": {"extinction_rel": (b["bext550"] / a["bext550"] - 1.0).tolist(),
                                         "w_abs": (b["w550"] - a["w550"]).tolist(), "g_abs": (b["g550"] - a["g550"]).tolist()},
                     "seconds": {"v23": float(a["seconds"]), "v24": float(b["seconds"])}}
      out = {"model": model, "consistency": consistency, "rows": rows}
      (HERE / "results" / f"optics_{model}.json").write_text(json.dumps(out, indent=1) + "\n")
      print(f"== {model}: V23 {consistency['seconds']['v23']:.0f} s, V24 {consistency['seconds']['v24']:.0f} s; "
            f"w in [0,1] {consistency['w_in_0_1']}, chi_0 = 1 {consistency['chi0_is_1']}, max|chi| {consistency['max_abs_chi']:.3g}, "
            f"max|chi_1 - g| {consistency['max_abs_chi1_minus_g']:.1e}")
      print(f"   {'ch':>3s} {'um':>6s} {'re':>4s} {'dext':>9s} {'dratio':>9s} {'dw':>9s} {'dg':>9s} {'L':>5s} {'dchi<1000':>9s}")
      for row in rows:
         print(f"   {row['channel']:3d} {row['wavelength']:6.3f} {row['effective_radius']:4.0f} {row['extinction_rel']:9.1e} "
               f"{row['extinction_ratio_rel']:9.1e} {row['w_abs']:9.1e} {row['g_abs']:9.1e} {row['L_v24']:5d} {row['max_chi_change_l_lt_1000']:9.1e}")


if __name__ == "__main__":
   main(sys.argv[1:])
