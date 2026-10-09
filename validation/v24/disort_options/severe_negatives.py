"""List the severe negative R_0v values (< -0.05) of ORAC cloud LUT products (validation only).

For each product the distinct (channel, tau, re, sza, vza) combinations where
R_0v < -0.05 at any relative azimuth, with the minimum and the number of
azimuths affected.  Read-only.

   python severe_negatives.py PRODUCT.nc [...]   -> results/disort_options/severe_negative_r0v.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.io.lut import read_lut   # noqa: E402

OUT = ROOT / "validation" / "v24" / "results" / "disort_options" / "severe_negative_r0v.csv"


def main():
   rows = []
   for p in sys.argv[1:]:
      lut = read_lut(p)
      dims = lut.variable_dimensions["R_0v"]
      x = np.asarray(lut.variables["R_0v"])
      ax = {d: np.asarray(lut.variables[d]) for d in ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith",
                                                      "relative_azimuth")}
      ch = np.asarray(lut.variables["solar_channel_id"])
      raa_axis = dims.index("relative_azimuth")
      m = np.moveaxis(x, raa_axis, -1)
      rest = [d for d in dims if d != "relative_azimuth"]
      mins = m.min(axis=-1)
      count = (m < -0.05).sum(axis=-1)
      for idx in np.argwhere(mins < -0.05):
         key = dict(zip(rest, idx))
         rows.append({"product": Path(p).name, "channel": int(ch[key["solar_channels"]]),
                      "tau": float(ax["optical_depth"][key["optical_depth"]]), "re": float(ax["effective_radius"][key["effective_radius"]]),
                      "sza": float(ax["solar_zenith"][key["solar_zenith"]]), "vza": float(ax["satellite_zenith"][key["satellite_zenith"]]),
                      "min_r0v": float(mins[tuple(idx)]), "azimuths_below": int(count[tuple(idx)])})
      print(Path(p).name, sum(1 for r in rows if r["product"] == Path(p).name), flush=True)
   OUT.parent.mkdir(parents=True, exist_ok=True)
   with open(OUT, "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=["product", "channel", "tau", "re", "sza", "vza", "min_r0v", "azimuths_below"])
      writer.writeheader()
      writer.writerows(rows)


if __name__ == "__main__":
   main()
