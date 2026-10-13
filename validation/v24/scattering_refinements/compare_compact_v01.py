"""Compare the compact V23 and V24 Terra MODIS LUTs of compact_lut.py (read-only).

   python validation/v24/scattering_refinements/compare_compact.py <model> ... -> results/compact_<model>.json

Both products come from the same compact Grid B subset and settings.  Checks:
identical variable list, dimensions and coordinates; finite values; for each
RT operator the median, 99th percentile and maximum |V24 - V23| with its
location, the relative maximum where |V23| >= 0.01, and the largest change
at each effective radius; the bulk optics written to the LUT.  Isolated
DISORT failures: for R_0v each direct-beam call (channel, optical depth,
radius, solar zenith) is summarised by its median |V24 - V23| over the view
geometries, and calls more than ten times the median of the other solar
zeniths of the same (channel, optical depth, radius) are listed.
"""

import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
RT = ("R_0v", "R_0d", "T_0d", "T_00", "R_dv", "T_dv", "R_dd", "T_dd", "E_md")
OPTICS = ("extinction_coefficient", "extinction_coefficient_ratio", "single_scatter_albedo", "asymmetry_parameter")
SHORT = {"liquid-water_240": "liquid-water_a01_p240", "liquid-water_old": "liquid-water_a01_pold", "water-ice_agg": "water-ice_a01_pagg"}


def location(dataset, dims, index):
   out = {}
   for dim, i in zip(dims, index):
      if dim in dataset.variables:
         out[dim] = float(dataset[dim][i])
      elif dim.endswith("channels"):
         ids = {"channels": "channel_id", "solar_channels": "solar_channel_id", "thermal_channels": "thermal_channel_id"}[dim]
         out[dim] = int(dataset[ids][i])
   return out


def main(models):
   for model in models:
      a = Dataset(HERE / "work" / f"compact_v23_{model}" / f"terra_modis_m_{SHORT[model]}_v23.nc")
      b = Dataset(HERE / "work" / f"compact_v24_{model}" / f"terra_modis_m_{SHORT[model]}_v24.nc")
      for d in (a, b):
         d.set_auto_maskandscale(False)
      out = {"model": model, "same_variables": list(a.variables) == list(b.variables),
             "same_dimensions": {k: len(v) for k, v in a.dimensions.items()} == {k: len(v) for k, v in b.dimensions.items()},
             "same_coordinates": all(np.array_equal(a[c][:], b[c][:]) for c in
                                     ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")),
             "v24_finite": all(bool(np.all(np.isfinite(np.asarray(b[n][:])))) for n in RT + OPTICS),
             "bitwise_identical_variables": [n for n in a.variables if np.asarray(a[n][:]).tobytes() == np.asarray(b[n][:]).tobytes()],
             "operators": {}, "optics": {}}
      radii = np.asarray(a["effective_radius"][:])
      for name in RT + OPTICS:
         x, y = np.asarray(a[name][:], np.float64), np.asarray(b[name][:], np.float64)
         if name in OPTICS:
            x, y = x[:, :8], y[:, :8]                             # the calculated channels (the rest is fill)
         d = np.abs(y - x)
         imax = np.unravel_index(np.argmax(d), d.shape)
         big = np.abs(x) >= 0.01
         rel = np.where(big, d / np.where(big, np.abs(x), 1.0), 0.0)
         irel = np.unravel_index(np.argmax(rel), rel.shape)
         axis = a[name].dimensions.index("effective_radius")
         by_radius = {f"{radii[r]:g}": float(np.max(np.take(d, r, axis=axis))) for r in range(radii.size)}
         entry = {"median": float(np.median(d)), "p99": float(np.percentile(d, 99)), "max": float(d.max()),
                  "max_at": location(a, a[name].dimensions, imax), "values_at_max": [float(x[imax]), float(y[imax])],
                  "max_rel": float(rel.max()), "max_rel_at": location(a, a[name].dimensions, irel), "max_by_radius": by_radius}
         (out["optics"] if name in OPTICS else out["operators"])[name] = entry
      # isolated direct-beam call outliers in R_0v (dims: raa, saz, soz, opd, efr, solar channel)
      d = np.abs(np.asarray(b["R_0v"][:], np.float64) - np.asarray(a["R_0v"][:], np.float64))
      call = np.median(d.reshape(-1, *d.shape[2:]), axis=0)          # (soz, opd, efr, channel)
      outliers = []
      for s in range(call.shape[0]):
         others = np.median(np.delete(call, s, axis=0), axis=0)
         ratio = call[s] / np.maximum(others, 1e-7)
         for idx in np.argwhere((ratio > 10.0) & (call[s] > 1e-4)):
            o, r, c = idx
            outliers.append({"solar_zenith": float(a["solar_zenith"][s]), "optical_depth": float(a["optical_depth"][o]),
                             "effective_radius": float(radii[r]), "channel": int(a["solar_channel_id"][c]),
                             "median_view_change": float(call[s, o, r, c]), "others": float(others[o, r, c])})
      out["r0v_call_outliers"] = outliers
      (HERE / "results" / f"compact_{model}.json").write_text(json.dumps(out, indent=1) + "\n")
      print(f"== {model}: variables {out['same_variables']}, dims {out['same_dimensions']}, coords {out['same_coordinates']}, "
            f"V24 finite {out['v24_finite']}; bitwise identical: {len(out['bitwise_identical_variables'])} of {len(a.variables)}; "
            f"R_0v call outliers {len(outliers)}")
      for group in ("optics", "operators"):
         for name, e in out[group].items():
            print(f"   {name:30s} med {e['median']:.1e} p99 {e['p99']:.1e} max {e['max']:.2e} rel {e['max_rel']:.1e} at {e['max_at']} "
                  f"by r_e {[f'{v:.0e}' for v in e['max_by_radius'].values()]}")


if __name__ == "__main__":
   main(sys.argv[1:])
