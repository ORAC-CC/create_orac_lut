"""Census of negative TOA bidirectional reflectance (R_0v) in ORAC LUT products (validation only).

A TOA radiance cannot be negative.  For every product given (default: all
*_v23.nc and *_v24.nc cloud products in the permanent LUT directory) this
counts R_0v < 0 and < -1e-4, the minimum, the scattering angles (ORAC
convention: relative azimuth 0 = backscatter side) and the channels and
effective radii where R_0v < -1e-4 occurs, read-only.

   python negative_radiance.py [PRODUCT.nc ...]   -> results/disort_options/negative_r0v.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.io.lut import read_lut   # noqa: E402

LUT_DIR = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
OUT = ROOT / "validation" / "v24" / "results" / "disort_options" / "negative_r0v.csv"


def census(path):
   lut = read_lut(path)
   dims = lut.variable_dimensions["R_0v"]
   x = np.asarray(lut.variables["R_0v"])
   sza, vza, raa = (np.deg2rad(np.asarray(lut.variables[d], dtype=np.float64))
                    for d in ("solar_zenith", "satellite_zenith", "relative_azimuth"))
   cos_theta = (-np.cos(sza)[:, None, None] * np.cos(vza)[None, :, None]
                - np.sin(sza)[:, None, None] * np.sin(vza)[None, :, None] * np.cos(raa)[None, None, :])
   theta = np.rad2deg(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
   order = [dims.index(d) for d in ("solar_zenith", "satellite_zenith", "relative_azimuth")]
   theta = np.broadcast_to(np.moveaxis(theta.reshape(theta.shape + (1,) * (x.ndim - 3)), [0, 1, 2], order), x.shape)
   bad = x < -1e-4
   idx = np.argwhere(bad)
   channels = np.asarray(lut.variables["solar_channel_id"])
   radii = np.asarray(lut.variables["effective_radius"])
   taus = np.asarray(lut.variables["optical_depth"])
   return {"product": path.name, "values": int(x.size), "negative": int((x < 0).sum()), "below_minus_1e-4": int(bad.sum()),
           "min": float(x.min()),
           "share_theta_ge_175": float((theta[bad] >= 175.0).mean()) if bad.any() else None,
           "theta_min": float(theta[bad].min()) if bad.any() else None,
           "channels": " ".join(str(int(c)) for c in np.unique(channels[idx[:, dims.index("solar_channels")]])) if bad.any() else "",
           "effective_radii": " ".join(f"{r:g}" for r in np.unique(radii[idx[:, dims.index("effective_radius")]])) if bad.any() else "",
           "optical_depths": " ".join(f"{t:g}" for t in np.unique(taus[idx[:, dims.index("optical_depth")]])) if bad.any() else ""}


def main():
   paths = [Path(p) for p in sys.argv[1:]] or sorted(p for p in LUT_DIR.glob("*_v2[34].nc") if "liquid-water" in p.name or "water-ice" in p.name)
   rows = []
   for path in paths:
      row = census(path)
      rows.append(row)
      print(row, flush=True)
   OUT.parent.mkdir(parents=True, exist_ok=True)
   with open(OUT, "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)


if __name__ == "__main__":
   main()
