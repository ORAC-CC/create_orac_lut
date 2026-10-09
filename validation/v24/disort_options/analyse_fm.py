"""Detailed analysis of the complete ORAC forward-model results (validation only).

Reads validation/tmp/disort_options/fm_full/<case>_<surface>.npz (propagate.py)
and writes compact tables to validation/v24/results/disort_options/:

   fm_ladder.csv        stream-count ladder (single precision, production kernel
                        semantics) and the double-precision points, against the
                        reference (double, NSTR 236): r = |dy|/sigma (production
                        ORAC Sy) by channel kind and scattering-angle class
   fm_breakdown.csv     NSTR 60 (production) against the reference: share of
                        points with r > 0.1, 0.5, 1, 2 by channel and angle class,
                        and the same with the reference's own uncertainty
                        (|double NSTR 136 - double NSTR 236|) taken into account
   fm_jacobian_ladder.csv  the same ladder for dy/dlog10(tau) and dy/dre,
                        |dJ| x (half the grid step) / sigma, i.e. the measurement
                        error the Jacobian error implies over half a LUT cell

Scattering-angle classes: "off" (< 170 deg), "near" (170-179 deg), "exact"
(180 deg: solar zenith = view zenith, relative azimuth 0).
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import orac_fm as fm                                  # noqa: E402
import propagate as pg                                # noqa: E402

FULL = pg.TMP / "fm_full"
LADDER = ["nstr008", "nstr016", "nstr024", "nstr032", "nstr040", "nstr048", "nstr056", "baseline", "nstr068", "nstr072",
          "nstr080", "nstr088", "nstr096", "nstr100", "nstr112", "nstr136", "nstr152", "nstr168", "nstr200", "nstr236",
          "double", "double_nstr136"]
NSTR = {"baseline": 60, "double": 60, "double_nstr136": 136}
SURFACES = ("dark", "vegetation", "snow")


def nstr(e):
   return NSTR.get(e, int(e[4:]) if e.startswith("nstr") and "_" not in e else None)


def load(case, surf):
   return np.load(FULL / f"{case}_{surf}.npz")


def sigma(case, ch, kind, y):
   const = pg.channel_constants(pg.lut(case, pg.REFERENCE), pg.PLATFORM[case][1])[ch]
   if kind == "solar":
      return np.where(y > 0, fm.sigma_solar(y, const["noise_rel"], True), 1e3)
   return fm.sigma_bt(y, const["nedt"], const["refbt"], const["planck"], mixed=kind == "mixed", production=True)


def angle_classes(theta, shape):
   t = np.broadcast_to(theta, shape)
   return {"off": t < 170.0, "near": (t >= 170.0) & (t < 179.5), "exact": t >= 179.5, "all": np.ones(shape, bool)}


def main():
   ladder_rows, breakdown_rows, jac_rows = [], [], []
   for case in [c for c in pg.CASES if not (c.endswith("_bs") or c.startswith("defect"))]:
      for surf in SURFACES:
         path = FULL / f"{case}_{surf}.npz"
         if not path.exists():
            continue
         d = load(case, surf)
         kinds = dict(k.split(":") for k in d["kinds"])
         axes_re = np.unique(d["state_re"][d["state_kind"] == "node"])
         axes_lt = np.log10(np.unique(d["state_tau"][d["state_kind"] == "node"]))
         half_step = {0: 0.5 * np.min(np.diff(axes_lt)), 1: 0.5 * np.min(np.diff(axes_re))}
         for ch_s, kind in kinds.items():
            ch = int(ch_s)
            y_ref = d[f"{pg.REFERENCE}__{ch}__y"].astype(np.float64)
            j_ref = d[f"{pg.REFERENCE}__{ch}__J"].astype(np.float64)
            sig = sigma(case, ch, kind, y_ref)
            classes = angle_classes(d["theta"], y_ref.shape) if kind != "thermal" else {"all": np.ones(y_ref.shape, bool)}
            ref_unc = np.abs(d[f"double_nstr136__{ch}__y"] - y_ref) / sig
            for e in LADDER:
               key = f"{e}__{ch}__y"
               if key not in d.files:
                  continue
               r = np.abs(d[key] - y_ref) / sig
               for cls, mask in classes.items():
                  m = mask & np.isfinite(r)
                  if not m.any():
                     continue
                  ladder_rows.append({"case": case, "surface": surf, "channel": ch, "kind": kind, "experiment": e,
                                      "nstreams": nstr(e), "precision": "double" if e.startswith("double") else "single",
                                      "angles": cls, "max_r": float(r[m].max()), "p99_r": float(np.percentile(r[m], 99)),
                                      "rms_r": float(np.sqrt(np.mean(r[m] ** 2)))})
               jkey = f"{e}__{ch}__J"
               for ix, name in ((0, "dlog10tau"), (1, "dre")):
                  dj = np.abs(d[jkey][:, ix].astype(np.float64) - j_ref[:, ix]) * half_step[ix] / sig
                  m = np.isfinite(dj)
                  jac_rows.append({"case": case, "surface": surf, "channel": ch, "kind": kind, "experiment": e,
                                   "nstreams": nstr(e), "derivative": name, "max": float(dj[m].max()),
                                   "p99": float(np.percentile(dj[m], 99)), "max_off_glory":
                                   float(dj[m & classes.get("off", classes["all"])].max())})
               if e == "baseline":
                  for cls, mask in classes.items():
                     m = mask & np.isfinite(r) & np.isfinite(ref_unc)
                     if not m.any():
                        continue
                     resolved = m & (r > 2.0 * ref_unc)        # baseline error clearly above the reference uncertainty
                     row = {"case": case, "surface": surf, "channel": ch, "kind": kind, "angles": cls, "points": int(m.sum())}
                     for level in (0.1, 0.5, 1.0, 2.0):
                        row[f"share_r_gt_{level:g}"] = float((r[m] > level).mean())
                        row[f"share_r_gt_{level:g}_resolved"] = float(((r > level) & resolved)[m].mean())
                     row["max_r"] = float(r[m].max())
                     row["max_reference_uncertainty"] = float(ref_unc[m].max())
                     breakdown_rows.append(row)
   for name, rows in (("fm_ladder.csv", ladder_rows), ("fm_breakdown.csv", breakdown_rows), ("fm_jacobian_ladder.csv", jac_rows)):
      with open(pg.RESULTS / name, "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
         writer.writeheader()
         writer.writerows(rows)
   print("written", pg.RESULTS)


if __name__ == "__main__":
   main()
