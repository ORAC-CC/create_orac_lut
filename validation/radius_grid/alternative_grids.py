"""Alternative radius quadratures at equal node count (V24 radius-grid investigation).

Validation-only code.  For a liquid-water case it compares, at the same number
of radius nodes N on the production interval [0.001 um, r_u]:

   uniform     trapezoid on the nested refinement of the legacy lattice (the
               production candidate family)
   log         trapezoid in ln r (nodes geometrically spaced)
   gauss       one-panel Gauss-Legendre in r

against the converged level-6 uniform reference written by
radius_convergence.py (same inputs and dense angles).

   PYTHONPATH=src python validation/radius_grid/alternative_grids.py REFERENCE.npz OUT.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import radius_convergence as rc   # noqa: E402
from oraclut.idl_mirror.load_mmdat import load_mmdat     # noqa: E402
from oraclut.optics.legacy_mie import mie_single_batch    # noqa: E402


def integrate(r, w, effective_radius, variance, wavelength, cm, mu):
   area = rc.number_weights(r, effective_radius, variance) * np.pi * r**2 * w
   ext = sca = gsca = 0.0
   phase = np.zeros(mu.size)
   for start in range(0, r.size, rc.BATCH):
      sl = slice(start, start + rc.BATCH)
      mie = mie_single_batch(2.0 * np.pi * r[sl] / wavelength, cm, mu)
      ext += area[sl] @ mie["qext"]
      sca += area[sl] @ mie["qsca"]
      gsca += area[sl] @ (mie["g"] * mie["qsca"])
      phase += (area[sl] * mie["qsca"]) @ mie["f11"]
   return ext, sca, gsca / sca, phase / sca


def main():
   reference = np.load(sys.argv[1])
   mmstr = load_mmdat(rc.INPUTS / "microphysics" / str(reference["microphysics"]), rc.INPUTS)
   re, wl, variance = float(reference["effective_radius"]), float(reference["wavelength"]), float(mmstr.s[0])
   cm = complex(rc.gsp._interpol_complex(mmstr.comp[0].cm, mmstr.comp[0].wl, [wl])[0])
   lower, upper = rc.bwgp.MIE_RADIUS_LOWER, float(reference["upper"])
   nq = int(reference["nq"])
   theta = reference["dense_theta"]
   mu = np.cos(np.deg2rad(theta))
   K = int(reference["levels"])
   ref_ext, ref_sca = reference["ext"][K], reference["sca"][K]
   ref_g = reference["gsca"][K] / ref_sca
   ref_phase = reference["phase"][K][nq:]
   rows = []
   for k in range(0, 5):
      n = int(reference["intervals0"]) * 2**k + 1
      grids = {}
      r = lower + (upper - lower) * np.arange(n) / (n - 1)
      w = np.full(n, (upper - lower) / (n - 1)); w[0] = w[-1] = w[1] / 2.0
      grids["uniform"] = (r, w)
      t = np.linspace(np.log(lower), np.log(upper), n)
      r = np.exp(t)
      w = np.full(n, t[1] - t[0]) * r; w[0] /= 2.0; w[-1] /= 2.0
      grids["log"] = (r, w)
      x, gw = np.polynomial.legendre.leggauss(n)
      grids["gauss"] = (lower + (upper - lower) * (x + 1.0) / 2.0, gw * (upper - lower) / 2.0)
      for name, (r, w) in grids.items():
         ext, sca, g, phase = integrate(r, w, re, variance, wl, cm, mu)
         rel = np.abs(phase / ref_phase - 1.0)
         rows.append({"case": Path(sys.argv[1]).stem, "nodes": n, "grid": name, "rel_ext": abs(ext / ref_ext - 1.0),
                      "abs_g": abs(g - ref_g), "phase_max_rel_side_5_150": float(rel[(theta >= 5) & (theta < 150)].max()),
                      "phase_max_rel_back_150_180": float(rel[theta >= 150].max()), "phase_rms_rel": float(np.sqrt(np.mean(rel**2)))})
         print(rows[-1], flush=True)
   with open(sys.argv[2], "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)


if __name__ == "__main__":
   main()
