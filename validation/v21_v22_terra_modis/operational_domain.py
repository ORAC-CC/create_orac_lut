"""Reanalyse the Terra MODIS V21/V22 LUT pairs inside ORAC's operational angular domain.

ORAC (ORAC-CC/orac eb9a1233) classifies a pixel as day only when
0 <= solar zenith < Ctrl%MaxSolZen (75 by default; src/determine_illum.F90,
src/read_driver.F90:614); only day pixels use solar LUT terms.  Satellite
zenith is rejected above Ctrl%MaxSatZen (default 90, src/get_geometry.F90);
the stated operational limit of 75 degrees is applied here as well.  LUT
angles are interpolated linearly between the bracketing nodes iX0 = locate(x)
and iX0 + 1 (src/set_gzero.F90, src/locate.F90, src/int_lut_*.F90), and the
coefficients are interpolated before the forward model is evaluated.

The operational solar-zenith set is therefore every node below the limit plus
interpolated slices (1 - w) F(70) + w F(80) at w = 0.25 and w -> 0.5 (solar
zenith 72.5 and 75-), the supremum of the weight ORAC can give the 80 degree
support node.  The forward-model and uncertainty definitions mirror
``radiance_impact.py``; with ``--validate`` the full-grid results are checked
against that utility's stored JSON (the historical V21 / lut-code-v22.1 analysis
results, removed from the working tree on 2026-10-12: restore
validation/v21_v22_terra_modis/results/ from commit f9b0b3b first).

Usage:
   python operational_domain.py [--validate]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import radiance_impact as impact                          # noqa: E402
from radiance_impact import fm, atmosphere_terms          # noqa: E402

LUT_DIRECTORY = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
CASES = {
   "pold": ("terra_modis_m_liquid-water_a01_pold", "liquid"),
   "p240": ("terra_modis_m_liquid-water_a01_p240", "liquid"),
   "pagg": ("terra_modis_m_water-ice_a01_pagg", "ice"),
}
MAX_SOL_ZEN = 75.0          # ORAC default Ctrl%MaxSolZen; day iff SZA < limit
MAX_SAT_ZEN = 75.0          # stated operational limit (ORAC default is 90); inclusive
RESULTS = HERE / "results" / "operational_domain"
PLOTS = ROOT / "plots"
SZA_VARIABLES = ("R_0v", "R_0d", "T_0d", "T_00")
RT_VARIABLES = ("R_0v", "R_0d", "T_0d", "T_00", "T_dv", "R_dv", "T_dd", "R_dd", "E_md")
THRESHOLDS = (1.0e-4, 1.0e-3, 1.0e-2)
MINIMUM_SIGNAL = 0.01       # as radiance_impact.py for relative/sigma statistics


# ---------------------------------------------------------------------------
# Angular domain
# ---------------------------------------------------------------------------

def orac_bracket(grid, x):
   """ORAC locate(bounded=.true.) lower index and weight for an ascending grid."""
   low = int(np.searchsorted(grid, x, side="right")) - 1
   low = min(max(low, 0), grid.size - 2)
   return low, (x - grid[low]) / (grid[low + 1] - grid[low])


def classify_nodes(grid, limit, inclusive):
   """Classify nodes as in-domain, interpolation support or unused.

   Every operational angle in [0, limit) (or [0, limit]) is sampled finely; a
   node is support if it lies outside the domain but receives non-zero weight.
   """
   upper = limit if inclusive else np.nextafter(limit, -np.inf)
   samples = np.unique(np.concatenate([np.linspace(0.0, upper, 200001), grid[grid <= upper]]))
   maximum_weight = np.zeros(grid.size)
   for x in samples:
      low, w = orac_bracket(grid, x)
      maximum_weight[low] = max(maximum_weight[low], 1.0 - w)
      maximum_weight[low + 1] = max(maximum_weight[low + 1], w)
   inside = grid <= limit if inclusive else grid < limit
   classes = np.where(inside, "in_domain", np.where(maximum_weight > 0.0, "support", "unused"))
   return [{"node": float(g), "class": str(c), "maximum_weight": float(m)}
           for g, c, m in zip(grid, classes, maximum_weight)]


def operational_sza_weights(grid, limit):
   """Rows of weights over the solar-zenith nodes for each operational geometry."""
   geometries, rows = [], []
   for index, node in enumerate(grid):
      if node < limit:
         row = np.zeros(grid.size)
         row[index] = 1.0
         geometries.append(float(node))
         rows.append(row)
   low, _ = orac_bracket(grid, np.nextafter(limit, -np.inf))
   if grid[low + 1] >= limit:
      for w in (0.25, 0.5):
         row = np.zeros(grid.size)
         row[low], row[low + 1] = 1.0 - w, w
         geometries.append(float(grid[low] + w * (grid[low + 1] - grid[low])))
         rows.append(row)
   return np.asarray(geometries), np.asarray(rows)


def apply_sza(values, axis, weights):
   """Linear angular interpolation of coefficients, as ORAC does."""
   return np.moveaxis(np.tensordot(weights, np.moveaxis(values, axis, 0), axes=(1, 0)), 0, axis)


# ---------------------------------------------------------------------------
# Coefficient statistics
# ---------------------------------------------------------------------------

def coordinates(dataset, dimensions, shape, sza_values=None):
   result = []
   for dimension, length in zip(dimensions, shape):
      if dimension == "solar_zenith" and sza_values is not None:
         result.append(np.asarray(sza_values))
      else:
         result.append(impact.comparison._coordinate_values(dataset, dimension, length))
   return result


def statistics(reference, candidate, coords, dimensions):
   difference = candidate - reference
   finite = np.isfinite(difference)
   d = np.where(finite, np.abs(difference), -1.0)
   index = np.unravel_index(int(np.argmax(d)), d.shape)
   meaningful = finite & (np.abs(reference) >= MINIMUM_SIGNAL)
   relative = np.where(meaningful, np.abs(difference) / np.where(meaningful, np.abs(reference), 1.0), -1.0)
   rel_index = np.unravel_index(int(np.argmax(relative)), relative.shape)
   n = int(finite.sum())
   record = {
      "elements": n,
      "bitwise_equal_fraction": float(np.mean(reference[finite] == candidate[finite])),
      "maximum_absolute_difference": float(d[index]),
      "location_of_maximum": {dim: float(c[i]) for dim, c, i in zip(dimensions, coords, index)},
      "reference_at_maximum": float(reference[index]),
      "candidate_at_maximum": float(candidate[index]),
      "rms_difference": float(np.sqrt(np.mean(difference[finite] ** 2))),
      "maximum_relative_difference_reference_ge_0.01": float(relative[rel_index]) if meaningful.any() else None,
      "location_of_maximum_relative": {dim: float(c[i]) for dim, c, i in zip(dimensions, coords, rel_index)},
   }
   for threshold in THRESHOLDS:
      count = int(np.sum(np.abs(difference[finite]) > threshold))
      record[f"count_over_{threshold:g}"] = count
      record[f"fraction_over_{threshold:g}"] = count / n if n else None
   return record


def coefficient_analysis(reference, candidate, sza_grid, sza_geometries, sza_weights, sza_classes):
   rows, dependence = [], []
   support = [i for i, c in enumerate(sza_classes) if c["class"] == "support"]
   unused = [i for i, c in enumerate(sza_classes) if c["class"] == "unused"]
   for name in RT_VARIABLES:
      dimensions = tuple(reference.variables[name].dimensions)
      ref = impact._values(reference, name, dimensions)
      cand = impact._values(candidate, name, dimensions)
      subsets = {}
      if "solar_zenith" in dimensions:
         axis = dimensions.index("solar_zenith")
         subsets["operational"] = (apply_sza(ref, axis, sza_weights), apply_sza(cand, axis, sza_weights), sza_geometries)
         subsets["support_node_raw"] = (np.take(ref, support, axis), np.take(cand, support, axis), sza_grid[support])
         subsets["unused_nodes"] = (np.take(ref, unused, axis), np.take(cand, unused, axis), sza_grid[unused])
         subsets["full_lut"] = (ref, cand, sza_grid)
      else:
         subsets["operational"] = (ref, cand, None)
      for subset, (r, c, sza_values) in subsets.items():
         coords = coordinates(reference, dimensions, r.shape, sza_values)
         record = statistics(r, c, coords, dimensions)
         record.update({"variable": name, "subset": subset})
         rows.append(record)
         if subset != "operational":
            continue
         difference = np.abs(c - r)
         for axis, dimension in enumerate(dimensions):
            other = tuple(a for a in range(difference.ndim) if a != axis)
            maximum = np.nanmax(difference, axis=other)
            rms = np.sqrt(np.nanmean((c - r) ** 2, axis=other))
            for value, m, s in zip(coords[axis], maximum, rms):
               dependence.append({"variable": name, "dimension": dimension, "coordinate": float(value),
                                  "maximum_absolute_difference": float(m), "rms_difference": float(s)})
   return rows, dependence


# ---------------------------------------------------------------------------
# Direct-beam calls
# ---------------------------------------------------------------------------

CALL_ORDER = ("solar_channels", "optical_depth", "effective_radius", "solar_zenith", "satellite_zenith",
              "relative_azimuth")


def call_analysis(reference, candidate, sza_classes):
   """Per DISORT direct-beam call (channel, tau, r_e, SZA node) statistics of R_0v."""
   ref = impact._values(reference, "R_0v", CALL_ORDER)
   cand = impact._values(candidate, "R_0v", CALL_ORDER)
   difference = np.abs(cand - ref)
   snr = impact._values(reference, "snr")
   # Mixed channels carry no SNR in the LUT; for them only the production
   # calibration terms (2.2 % of R) are used as a reflectance-level proxy.  Their
   # operational significance is assessed in brightness temperature below.
   noise = np.where(snr > 0.0, 1.0 / np.where(snr > 0.0, snr, 1.0), 0.0)
   sigma = fm.sigma_solar(ref, noise[:, None, None, None, None, None], production=True)
   ratio = np.where(np.abs(ref) >= MINIMUM_SIGNAL, difference / sigma, np.nan)
   call_max = np.nanmax(difference, axis=(4, 5))
   call_sigma = np.nanmax(np.where(np.isfinite(ratio), ratio, -1.0), axis=(4, 5))
   call_median = np.nanmedian(difference, axis=(4, 5))
   shape = call_max.shape
   # Neighbour context: same channel, adjacent tau / r_e / SZA nodes.
   neighbour = np.full(shape, np.nan)
   for index in np.ndindex(shape):
      values = []
      for axis in (1, 2, 3):
         for step in (-1, 1):
            j = list(index)
            j[axis] += step
            if 0 <= j[axis] < shape[axis]:
               values.append(call_max[tuple(j)])
      neighbour[index] = np.median(values)
   channels = impact._values(reference, "solar_channel_id").astype(int)
   mixed = set(impact._values(reference, "mixed_channel_id").astype(int))
   tau = impact._values(reference, "optical_depth")
   re = impact._values(reference, "effective_radius")
   sza = impact._values(reference, "solar_zenith")
   rows = []
   for index in np.ndindex(shape):
      c, t, r, s = index
      rows.append({
         "channel": int(channels[c]), "channel_type": "mixed" if channels[c] in mixed else "solar",
         "optical_depth": float(tau[t]), "effective_radius": float(re[r]), "solar_zenith": float(sza[s]),
         "sza_class": sza_classes[s]["class"],
         "maximum_absolute_R_0v_difference": float(call_max[index]),
         "median_absolute_R_0v_difference": float(call_median[index]),
         "neighbour_median_of_call_maxima": float(neighbour[index]),
         "isolation_ratio": float(call_max[index] / neighbour[index]) if neighbour[index] > 0 else None,
         "maximum_black_surface_sigma": float(call_sigma[index]) if call_sigma[index] >= 0 else None,
         "reference_R_0v_minimum": float(np.nanmin(ref[index])), "candidate_R_0v_minimum": float(np.nanmin(cand[index])),
      })
   return rows


def classify_calls(rows):
   """Isolated failure: >= 1 sigma and >= 10x the neighbouring calls' maxima."""
   for row in rows:
      ratio = row["isolation_ratio"]
      sig = row["maximum_black_surface_sigma"] or 0.0
      if sig >= 1.0 and ratio is not None and ratio >= 10.0:
         row["call_class"] = "isolated_significant"
      elif sig >= 1.0:
         row["call_class"] = "significant_not_isolated"
      else:
         row["call_class"] = "below_1_sigma"
   return rows


# ---------------------------------------------------------------------------
# Measurement impact (mirrors radiance_impact.py, with SZA geometry weights)
# ---------------------------------------------------------------------------

def summary(difference, sigma, signal, coords, dimensions, minimum_signal=None):
   record = impact._ratio_summary(difference, sigma, signal, minimum_signal)
   if record.get("elements"):
      valid = np.isfinite(difference) & np.isfinite(sigma) & (sigma > 0.0)
      if minimum_signal is not None:
         valid &= np.abs(signal) >= minimum_signal
      ratio = np.abs(difference[valid]) / sigma[valid]
      record["count_over_1_sigma"] = int(np.sum(ratio > 1.0))
      record["count_over_0.1_sigma"] = int(np.sum(ratio > 0.1))
      record["location_of_maximum"] = {d: float(c[i]) for d, c, i in zip(dimensions, coords, record["maximum_index"])}
   return record


def solar_black(reference, candidate, sza_geometries, sza_weights):
   dims = ("solar_channels", "optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")
   ref = apply_sza(impact._values(reference, "R_0v", dims), 3, sza_weights)
   cand = apply_sza(impact._values(candidate, "R_0v", dims), 3, sza_weights)
   channels = impact._values(reference, "solar_channel_id").astype(int)
   mixed = set(impact._values(reference, "mixed_channel_id").astype(int))
   snr = impact._values(reference, "snr")
   coords = coordinates(reference, dims[1:], ref.shape[1:], sza_geometries)
   rows = []
   for index, channel in enumerate(channels):
      noise = 1.0 / snr[index] if snr[index] > 0.0 else np.nan
      signal = ref[index]
      for uncertainty, production in (("production", True), ("instrument_only", False)):
         sigma = fm.sigma_solar(signal, noise, production=production)
         record = summary(cand[index] - signal, sigma, signal, coords, dims[1:], MINIMUM_SIGNAL)
         record.update({"measurement": "solar_reflectance_black_surface", "channel": int(channel),
                        "channel_type": "mixed" if channel in mixed else "solar",
                        "uncertainty": uncertainty, "units": "reflectance"})
         rows.append(record)
   return rows


def lambertian(dataset, albedo, sza_geometries, sza_weights):
   _, channel_indices = impact._channel_indices(dataset, "solar_channel_id")
   r0v = impact._values(dataset, "R_0v", ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth",
                                          "effective_radius", "solar_channels"))[0, 0]
   r0v = apply_sza(r0v, 0, sza_weights)
   t00_all = impact._values(dataset, "T_00", ("solar_zenith", "optical_depth", "effective_radius", "solar_channels"))
   t0d = apply_sza(impact._values(dataset, "T_0d", ("solar_zenith", "optical_depth", "effective_radius",
                                                     "solar_channels")), 0, sza_weights)
   t00 = apply_sza(t00_all, 0, sza_weights)
   tdv = np.take(impact._values(dataset, "T_dv", ("satellite_zenith", "optical_depth", "effective_radius",
                                                  "channels"))[0], channel_indices, axis=-1)
   rdd = np.take(impact._values(dataset, "R_dd", ("optical_depth", "effective_radius", "channels")),
                 channel_indices, axis=-1)
   tvv = t00_all[0]                    # nadir view: T_vv is the exact SZA = 0 node
   mu0 = np.cos(np.deg2rad(sza_geometries))[:, None, None, None]
   a = 1.0 - albedo * rdd[None, ...]
   b = albedo * (t00 + t0d)
   c = tdv[None, ...] + albedo * rdd[None, ...] * tvv[None, ...]
   return r0v + mu0 * (albedo * t00 * tvv[None, ...] + albedo * t0d * tvv[None, ...] + b * c / a)


def solar_lambertian(reference, candidate, sza_geometries, sza_weights):
   channels = impact._values(reference, "solar_channel_id").astype(int)
   mixed = set(impact._values(reference, "mixed_channel_id").astype(int))
   snr = impact._values(reference, "snr")
   dims = ("solar_zenith", "optical_depth", "effective_radius")
   rows = []
   for albedo in (0.05, 0.2, 0.8):
      ref = lambertian(reference, albedo, sza_geometries, sza_weights)
      cand = lambertian(candidate, albedo, sza_geometries, sza_weights)
      coords = coordinates(reference, dims, ref.shape[:3], sza_geometries)
      for index, channel in enumerate(channels):
         signal = ref[..., index]
         noise = 1.0 / snr[index] if snr[index] > 0.0 else np.nan
         sigma = fm.sigma_solar(signal, noise, production=True)
         record = summary(cand[..., index] - signal, sigma, signal, coords, dims, MINIMUM_SIGNAL)
         record.update({"measurement": "solar_reflectance_lambertian_nadir", "channel": int(channel),
                        "channel_type": "mixed" if channel in mixed else "solar",
                        "uncertainty": "production", "units": "reflectance", "surface_albedo": albedo})
         rows.append(record)
   return rows


def thermal_and_mixed(reference, candidate, cloud, sza_geometries, sza_weights):
   v = impact._values
   all_channels = v(reference, "channel_id").astype(int)
   solar_channels, _ = impact._channel_indices(reference, "solar_channel_id")
   thermal_channels, _ = impact._channel_indices(reference, "thermal_channel_id")
   mixed = set(v(reference, "mixed_channel_id").astype(int))
   solar_lookup = {channel: index for index, channel in enumerate(solar_channels)}
   thermal_lookup = {channel: index for index, channel in enumerate(thermal_channels)}
   planck = {int(channel): tuple(float(v(reference, name)[index]) for name in ("B1", "B2", "T1", "T2"))
             for index, channel in enumerate(thermal_channels)}
   nedt, refbt = v(reference, "nedt"), v(reference, "refbt")
   f0_reference, f0_candidate = v(reference, "F0"), v(candidate, "F0")
   vza = v(reference, "satellite_zenith")
   tdv_dims = ("satellite_zenith", "optical_depth", "effective_radius", "channels")
   em_dims = ("satellite_zenith", "optical_depth", "effective_radius", "thermal_channels")
   r0v_dims = ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth", "effective_radius",
               "solar_channels")
   tdv = [v(d, "T_dv", tdv_dims) for d in (reference, candidate)]
   rdv = [v(d, "R_dv", tdv_dims) for d in (reference, candidate)]
   em = [v(d, "E_md", em_dims) for d in (reference, candidate)]
   r0v = [apply_sza(v(d, "R_0v", r0v_dims), 2, sza_weights) for d in (reference, candidate)]
   rows = []
   for vza_index, angle in enumerate(vza):
      if angle > MAX_SAT_ZEN:
         continue
      atmosphere = atmosphere_terms("terra", "modis", list(all_channels), planck, cloud, float(angle))
      for channel in thermal_channels:
         ti = thermal_lookup[channel]
         ai = int(np.flatnonzero(all_channels == channel)[0])
         entry = atmosphere[int(channel)]
         radiance = [((entry["Rbc_up"] * tdv[k][vza_index, :, :, ai] + entry["B_c"] * em[k][vza_index, :, :, ti] +
                       entry["Rac_dwn"] * rdv[k][vza_index, :, :, ai]) * entry["Tac_lw"] + entry["Rac_up"])
                     for k in (0, 1)]
         is_mixed = int(channel) in mixed
         if is_mixed:
            si = solar_lookup[channel]
            secv = 1.0 / np.cos(np.deg2rad(angle))
            transmission = entry["Tac"] ** (1.0 / np.cos(np.deg2rad(sza_geometries))) * entry["Tac"] ** secv
            totals = []
            for k, f0 in ((0, f0_reference), (1, f0_candidate)):
               reflectance = np.transpose(r0v[k][:, vza_index, :, :, :, si], (1, 2, 3, 0))
               totals.append(radiance[k][None, :, :, None] + f0[si] * transmission[:, None, None, None] * reflectance)
            dims = ("solar_zenith", "optical_depth", "effective_radius", "relative_azimuth")
         else:
            totals = radiance
            dims = ("optical_depth", "effective_radius")
         bt = [fm.r2t(total, planck[int(channel)])[0] for total in totals]
         sigma = fm.sigma_bt(bt[0], nedt[ti], refbt[ti], planck[int(channel)], mixed=is_mixed, production=True)
         coords = coordinates(reference, dims, bt[0].shape, sza_geometries)
         record = summary(bt[1] - bt[0], sigma, bt[0], coords, dims)
         record.update({"measurement": "mixed_brightness_temperature" if is_mixed else "thermal_brightness_temperature",
                        "channel": int(channel), "channel_type": "mixed" if is_mixed else "thermal",
                        "uncertainty": "production", "units": "K", "satellite_zenith": float(angle),
                        "cloud_profile": cloud})
         rows.append(record)
   return rows


def measurement_rows(reference, candidate, cloud, sza_geometries, sza_weights):
   rows = solar_black(reference, candidate, sza_geometries, sza_weights)
   rows += solar_lambertian(reference, candidate, sza_geometries, sza_weights)
   rows += thermal_and_mixed(reference, candidate, cloud, sza_geometries, sza_weights)
   return rows


# ---------------------------------------------------------------------------
# Validation against radiance_impact.py
# ---------------------------------------------------------------------------

def validate(case, reference, candidate, cloud):
   """With the full SZA grid as identity weights, reproduce the stored utility output."""
   sza = impact._values(reference, "solar_zenith")
   rows = measurement_rows(reference, candidate, cloud, sza, np.eye(sza.size))
   stored = json.loads((HERE / "results" / f"terra_modis_{case}_v21_v22_radiance_impact.json").read_text())
   keys = ("maximum_sigma", "maximum_absolute_measurement_difference", "rms_sigma", "share_over_1_sigma", "elements")
   if len(rows) != len(stored):
      raise AssertionError(f"{case}: {len(rows)} rows versus {len(stored)} stored")
   worst = 0.0
   for mine, theirs in zip(rows, stored):
      for key in ("measurement", "channel"):
         if mine[key] != theirs[key]:
            raise AssertionError(f"{case}: row order differs at {mine[key]} / {theirs[key]}")
      for key in keys:
         a, b = mine.get(key), theirs.get(key)
         if a is None and b is None:
            continue
         scale = max(abs(b), 1e-300)
         worst = max(worst, abs(a - b) / scale)
   if worst > 1e-12:
      raise AssertionError(f"{case}: full-grid reproduction relative error {worst:.3e}")
   return {"case": case, "rows": len(rows), "maximum_relative_mismatch": worst}


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_csv(path, rows):
   fieldnames = []
   for row in rows:
      for name in row:
         if name not in fieldnames:
            fieldnames.append(name)
   with path.open("w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=fieldnames)
      writer.writeheader()
      for row in rows:
         writer.writerow({k: (json.dumps(v) if isinstance(v, (dict, list)) else v) for k, v in row.items()})


def plots(all_calls, all_measurements, full_measurements):
   os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "validation" / "generated" / "matplotlib-cache"))
   Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
   import matplotlib
   matplotlib.use("Agg")
   import matplotlib.pyplot as plt

   paths = []
   # 1. Call maxima of |dR_0v| against SZA node, coloured by domain class.
   fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
   colours = {"in_domain": "tab:blue", "support": "tab:orange", "unused": "0.6"}
   for axis, (case, rows) in zip(axes, all_calls.items()):
      for cls, colour in colours.items():
         sel = [r for r in rows if r["sza_class"] == cls]
         x = np.asarray([r["solar_zenith"] for r in sel]) + np.random.default_rng(0).uniform(-2, 2, len(sel))
         y = np.maximum([r["maximum_absolute_R_0v_difference"] for r in sel], 1e-9)
         axis.scatter(x, y, s=4, color=colour, alpha=0.5, label=cls.replace("_", " "))
      axis.axvline(MAX_SOL_ZEN, color="red", linestyle="--", linewidth=1)
      axis.set_yscale("log")
      axis.set_title(f"Terra MODIS {case}")
      axis.set_xlabel("solar zenith node (deg)")
      axis.grid(True, which="both", alpha=0.25)
   axes[0].set_ylabel("per-call max |R_0v(V22) - R_0v(V21)|")
   axes[0].legend(loc="lower left", markerscale=3, fontsize=8)
   fig.suptitle("V21/V22 direct-beam call differences; ORAC day limit SZA < 75 deg (red)")
   fig.tight_layout()
   path = PLOTS / "terra_modis_v21_v22_operational_call_differences.png"
   fig.savefig(path, dpi=160)
   plt.close(fig)
   paths.append(path)

   # 2. Measurement impact per channel, full LUT versus operational domain.
   fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
   for axis, case in zip(axes, all_measurements):
      def best(rows):
         out = {}
         for r in rows:
            if r.get("uncertainty") != "production" or "maximum_sigma" not in r:
               continue
            if r["measurement"] == "solar_reflectance_black_surface" and r["channel_type"] == "mixed":
               continue
            key = r["channel"]
            out[key] = max(out.get(key, 0.0), r["maximum_sigma"])
         return out
      full = best(full_measurements[case])
      oper = best(all_measurements[case])
      channels = sorted(oper)
      y = np.arange(len(channels))
      axis.barh(y + 0.2, [full[c] for c in channels], height=0.4, color="0.7", label="full LUT grid")
      axis.barh(y - 0.2, [oper[c] for c in channels], height=0.4, color="tab:blue", label="operational domain")
      axis.set_yticks(y, [f"ch{c}" for c in channels])
      axis.axvline(1.0, color="red", linestyle="--", linewidth=1)
      axis.set_xscale("log")
      axis.set_title(f"Terra MODIS {case}")
      axis.set_xlabel("max |V22 - V21| / production uncertainty")
      axis.grid(True, which="both", axis="x", alpha=0.25)
   axes[0].legend(loc="lower right", fontsize=8)
   fig.suptitle("Representative measurement impact (max over surfaces, atmosphere and geometry)")
   fig.tight_layout()
   path = PLOTS / "terra_modis_v21_v22_operational_measurement_impact.png"
   fig.savefig(path, dpi=160)
   plt.close(fig)
   paths.append(path)
   return paths


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("--validate", action="store_true", help="check full-grid reproduction of radiance_impact.py")
   args = parser.parse_args()
   RESULTS.mkdir(parents=True, exist_ok=True)
   PLOTS.mkdir(parents=True, exist_ok=True)
   overall = {"orac_source": "ORAC-CC/orac eb9a1233c1136ef95fa0c37f650fe519e9134ed4",
              "max_sol_zen": MAX_SOL_ZEN, "max_sat_zen": MAX_SAT_ZEN, "cases": {}, "validation": []}
   all_calls, all_measurements, full_measurements = {}, {}, {}
   for case, (stem, cloud) in CASES.items():
      reference_path = LUT_DIRECTORY / f"{stem}_v21.nc"
      candidate_path = LUT_DIRECTORY / f"{stem}_v22.nc"
      with Dataset(reference_path) as reference, Dataset(candidate_path) as candidate:
         sza = impact._values(reference, "solar_zenith")
         vza = impact._values(reference, "satellite_zenith")
         sza_classes = classify_nodes(sza, MAX_SOL_ZEN, inclusive=False)
         vza_classes = classify_nodes(vza, MAX_SAT_ZEN, inclusive=True)
         if any(c["class"] != "in_domain" for c in vza_classes):
            raise NotImplementedError("satellite-zenith nodes outside the domain are not handled")
         geometries, weights = operational_sza_weights(sza, MAX_SOL_ZEN)
         if args.validate:
            overall["validation"].append(validate(case, reference, candidate, cloud))
         coefficient_rows, dependence = coefficient_analysis(reference, candidate, sza, geometries, weights, sza_classes)
         calls = classify_calls(call_analysis(reference, candidate, sza_classes))
         measurements = measurement_rows(reference, candidate, cloud, geometries, weights)
         full_measurements[case] = measurement_rows(reference, candidate, cloud, sza, np.eye(sza.size))
      all_calls[case] = calls
      all_measurements[case] = measurements
      write_csv(RESULTS / f"{case}_coefficients.csv", coefficient_rows)
      write_csv(RESULTS / f"{case}_coordinate_dependence.csv", dependence)
      ranked = sorted(calls, key=lambda r: -r["maximum_absolute_R_0v_difference"])
      write_csv(RESULTS / f"{case}_direct_beam_calls.csv", ranked[:1000])
      write_csv(RESULTS / f"{case}_measurement_impact.csv", measurements)
      operational_calls = [r for r in calls if r["sza_class"] == "in_domain"]
      overall["cases"][case] = {
         "reference": str(reference_path), "candidate": str(candidate_path),
         "solar_zenith_nodes": sza_classes, "satellite_zenith_nodes": vza_classes,
         "operational_solar_zenith_geometries": geometries.tolist(),
         "coefficients": coefficient_rows,
         "call_classes": {cls: {sza_cls: sum(1 for r in calls if r["call_class"] == cls and r["sza_class"] == sza_cls)
                                for sza_cls in ("in_domain", "support", "unused")}
                          for cls in ("isolated_significant", "significant_not_isolated", "below_1_sigma")},
         "significant_calls": [r for r in ranked if r["call_class"] != "below_1_sigma"],
         "in_domain_call_maximum_percentiles": {
            str(p): float(np.percentile([r["maximum_absolute_R_0v_difference"] for r in operational_calls], p))
            for p in (50, 90, 99, 99.9, 100)},
         "measurements": measurements,
      }
   overall["plots"] = [str(p) for p in plots(all_calls, all_measurements, full_measurements)]
   (RESULTS / "operational_domain_summary.json").write_text(json.dumps(overall, indent=1, allow_nan=False) + "\n")
   print(json.dumps({"validation": overall["validation"], "plots": overall["plots"]}, indent=1))
   return 0


if __name__ == "__main__":
   raise SystemExit(main())
