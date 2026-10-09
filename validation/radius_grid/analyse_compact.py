"""Geometry-resolved comparison of the compact radius-level LUTs (V24 radius grid).

Validation-only code, run after compact_lut_radius.py.  For each case and
level against a reference level it reports, per solar channel, the largest
absolute and relative change of the bidirectional reflectance R_0v at the LUT
nodes, separately for scattering angles >= 170 deg (the glory / exact
backscatter, where the radius-quadrature error of the phase function is
largest) and < 170 deg, and in units of the channel's measurement
uncertainty (the instrument file's SAD noise-equivalent reflectance) where
one is defined.  The scattering angle uses the ORAC convention (relative
azimuth 0 = satellite on the sun's side):
   cos(Theta) = -cos(sza) cos(vza) - sin(sza) sin(vza) cos(raa).

   python validation/radius_grid/analyse_compact.py REFERENCE
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
from compact_lut_radius import CASES, INPUTS, PRODUCTS, RESULTS, VARIANTS   # noqa: E402
from oraclut.io.lut import read_lut                                          # noqa: E402
from oraclut.idl_mirror import load_inststr                                  # noqa: E402


def main():
   reference = sys.argv[1]
   rows = []
   for case, (instfile, mmfile, grid, channels) in CASES.items():
      base = PRODUCTS / case
      if not (base / reference / "timing.json").exists():
         continue
      inst = load_inststr(INPUTS / "inst" / instfile, requestedchannelid=channels)
      ref = read_lut(next((base / reference).glob("*.nc")))
      dims = ref.variable_dimensions["R_0v"]
      a = np.asarray(ref.variables["R_0v"], dtype=np.float64)
      order = [dims.index(d) for d in ("solar_zenith", "satellite_zenith", "relative_azimuth")]
      sza, vza, raa = (np.deg2rad(np.asarray(ref.variables[d], dtype=np.float64))
                       for d in ("solar_zenith", "satellite_zenith", "relative_azimuth"))
      cos_theta = (-np.cos(sza)[:, None, None] * np.cos(vza)[None, :, None]
                   - np.sin(sza)[:, None, None] * np.sin(vza)[None, :, None] * np.cos(raa)[None, None, :])
      theta = np.rad2deg(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
      # (sza, vza, raa) moved onto the axes of R_0v, broadcast over the others
      theta_full = np.broadcast_to(np.moveaxis(theta.reshape(theta.shape + (1,) * (a.ndim - 3)), [0, 1, 2], order), a.shape)
      channel_axis = dims.index("solar_channels")
      solar = [int(ch) for ch, flag in zip(inst.channelid, inst.solar_channel_flag) if flag]
      nviews = a.shape[channel_axis] // len(solar)
      for variant in VARIANTS:
         if variant == reference or not (base / variant / "timing.json").exists():
            continue
         b = np.asarray(read_lut(next((base / variant).glob("*.nc"))).variables["R_0v"], dtype=np.float64)
         for c in range(a.shape[channel_axis]):
            channel = solar[c % len(solar)]
            i = list(inst.channelid).index(channel)
            sigma = float(inst.oldnefr[i]) if float(inst.oldnefr[i]) > 0 else None
            x = np.take(a, c, axis=channel_axis)
            y = np.take(b, c, axis=channel_axis)
            t = np.take(theta_full, c, axis=channel_axis)
            d = np.abs(y - x)
            rel = d / np.maximum(np.abs(x), 1e-3)
            row = {"case": case, "reference": reference, "variant": variant, "channel": channel,
                   "view": c // len(solar) if nviews > 1 else 0, "sigma_reflectance": sigma}
            for label, mask in (("glory_ge170", t >= 170.0), ("other_lt170", t < 170.0)):
               row[f"max_abs_dR_{label}"] = float(d[mask].max()) if np.any(mask) else 0.0
               row[f"max_rel_dR_{label}"] = float(rel[mask].max()) if np.any(mask) else 0.0
               row[f"max_dR_over_sigma_{label}"] = float(d[mask].max() / sigma) if sigma and np.any(mask) else None
            row["p99_rel_dR_all"] = float(np.percentile(rel, 99))
            rows.append(row)
   out = RESULTS / f"{reference}_r0v_by_geometry.csv"
   with open(out, "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)
   for case in CASES:
      for variant in VARIANTS:
         sel = [r for r in rows if r["case"] == case and r["variant"] == variant]
         if not sel:
            continue
         worst = lambda key: max((r[key] for r in sel if r[key] is not None), default=float("nan"))
         print(f"{case:14s} {variant:7s} vs {reference}: glory max |dR| {worst('max_abs_dR_glory_ge170'):.1e} "
               f"(rel {worst('max_rel_dR_glory_ge170'):.1e}, {worst('max_dR_over_sigma_glory_ge170'):.2f} sigma); "
               f"other max |dR| {worst('max_abs_dR_other_lt170'):.1e} (rel {worst('max_rel_dR_other_lt170'):.1e}, "
               f"{worst('max_dR_over_sigma_other_lt170'):.2f} sigma); p99 rel {worst('p99_rel_dR_all'):.1e}")
   print("written", out)


if __name__ == "__main__":
   main()
