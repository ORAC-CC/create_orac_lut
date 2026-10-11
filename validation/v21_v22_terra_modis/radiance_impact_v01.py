"""Propagate a V21/V22 LUT pair through representative operational ORAC measurements.

The forward equations and uncertainty definitions are the validated mirror in
``validation/v24/disort_options/orac_fm.py``.  Solar black-surface results use
all LUT geometries.  Lambertian results use nadir view, where T_vv is an exact
LUT node, and the McGarragh et al. (2018) equation-form-3 surface coupling.
Thermal and mixed-channel results use the repository mid-latitude-summer
atmosphere and gas profiles through ``atmosphere.py``.
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
FM_DIRECTORY = ROOT / "validation" / "v24" / "disort_options"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(FM_DIRECTORY))
import compare_luts as comparison                         # noqa: E402
import orac_fm as fm                                      # noqa: E402
from atmosphere import terms as atmosphere_terms          # noqa: E402


def _values(dataset: Dataset, name: str, dimensions: tuple[str, ...] | None = None) -> np.ndarray:
   variable = dataset.variables[name]
   values = comparison._float_with_nan(comparison._decoded(variable))
   if dimensions is None:
      return values
   if set(dimensions) != set(variable.dimensions) or len(dimensions) != variable.ndim:
      raise ValueError(f"{name}: requested dimensions {dimensions} do not match {variable.dimensions}")
   return np.transpose(values, [variable.dimensions.index(dimension) for dimension in dimensions])


def _channel_indices(dataset: Dataset, subset_name: str) -> tuple[np.ndarray, np.ndarray]:
   all_channels = _values(dataset, "channel_id").astype(int)
   subset = _values(dataset, subset_name).astype(int)
   indices = np.asarray([int(np.flatnonzero(all_channels == channel)[0]) for channel in subset])
   return subset, indices


def _ratio_summary(difference: np.ndarray, sigma: np.ndarray, signal: np.ndarray, minimum_signal: float | None = None) -> dict:
   valid = np.isfinite(difference) & np.isfinite(sigma) & (sigma > 0.0)
   if minimum_signal is not None:
      valid &= np.abs(signal) >= minimum_signal
   ratio = np.full(difference.shape, np.nan, dtype=np.float64)
   ratio[valid] = np.abs(difference[valid]) / sigma[valid]
   values = ratio[valid]
   if not values.size:
      return {"elements": 0}
   flat_index = int(np.nanargmax(ratio))
   index = tuple(int(item) for item in np.unravel_index(flat_index, ratio.shape))
   return {
      "elements": int(values.size),
      "maximum_absolute_measurement_difference": float(np.max(np.abs(difference[valid]))),
      "rms_measurement_difference": float(np.sqrt(np.mean(difference[valid] ** 2))),
      "maximum_sigma": float(np.max(values)),
      "p99_sigma": float(np.percentile(values, 99.0)),
      "p999_sigma": float(np.percentile(values, 99.9)),
      "rms_sigma": float(np.sqrt(np.mean(values ** 2))),
      "share_over_0.1_sigma": float(np.mean(values > 0.1)),
      "share_over_0.5_sigma": float(np.mean(values > 0.5)),
      "share_over_1_sigma": float(np.mean(values > 1.0)),
      "share_over_2_sigma": float(np.mean(values > 2.0)),
      "maximum_index": list(index),
      "reference_measurement_at_maximum": float(signal[index]),
      "candidate_minus_reference_at_maximum": float(difference[index]),
   }


def _solar_black(reference: Dataset, candidate: Dataset) -> list[dict]:
   dimensions = ("solar_channels", "optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")
   reference_r0v = _values(reference, "R_0v", dimensions)
   candidate_r0v = _values(candidate, "R_0v", dimensions)
   channels = _values(reference, "solar_channel_id").astype(int)
   snr = _values(reference, "snr")
   rows = []
   for index, channel in enumerate(channels):
      noise = 1.0 / snr[index] if snr[index] > 0.0 else np.nan
      signal = reference_r0v[index]
      difference = candidate_r0v[index] - signal
      production_sigma = fm.sigma_solar(signal, noise, production=True)
      instrument_sigma = fm.sigma_solar(signal, noise, production=False)
      for uncertainty, sigma in (("production", production_sigma), ("instrument_only", instrument_sigma)):
         record = _ratio_summary(difference, sigma, signal, minimum_signal=0.01)
         record.update({"measurement": "solar_reflectance_black_surface", "channel": int(channel),
                        "uncertainty": uncertainty, "units": "reflectance", "minimum_abs_reference_signal": 0.01})
         rows.append(record)
   return rows


def _lambertian_reflectance(dataset: Dataset, albedo: float) -> tuple[np.ndarray, tuple[str, ...]]:
   solar_channels, channel_indices = _channel_indices(dataset, "solar_channel_id")
   r0v = _values(dataset, "R_0v", ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth",
                                      "effective_radius", "solar_channels"))[0, 0]
   t00 = _values(dataset, "T_00", ("solar_zenith", "optical_depth", "effective_radius", "solar_channels"))
   t0d = _values(dataset, "T_0d", ("solar_zenith", "optical_depth", "effective_radius", "solar_channels"))
   tdv = np.take(_values(dataset, "T_dv", ("satellite_zenith", "optical_depth", "effective_radius", "channels"))[0],
                 channel_indices, axis=-1)
   rdd = np.take(_values(dataset, "R_dd", ("optical_depth", "effective_radius", "channels")), channel_indices, axis=-1)
   tvv = t00[0]
   mu0 = np.cos(np.deg2rad(_values(dataset, "solar_zenith")))[:, None, None, None]
   a = 1.0 - albedo * rdd[None, ...]
   b = albedo * (t00 + t0d)
   c = tdv[None, ...] + albedo * rdd[None, ...] * tvv[None, ...]
   reflectance = r0v + mu0 * (albedo * t00 * tvv[None, ...] + albedo * t0d * tvv[None, ...] + b * c / a)
   return reflectance, ("solar_zenith", "optical_depth", "effective_radius", "solar_channels")


def _solar_lambertian(reference: Dataset, candidate: Dataset) -> list[dict]:
   channels = _values(reference, "solar_channel_id").astype(int)
   snr = _values(reference, "snr")
   rows = []
   for albedo in (0.05, 0.2, 0.8):
      reference_measurement, _ = _lambertian_reflectance(reference, albedo)
      candidate_measurement, _ = _lambertian_reflectance(candidate, albedo)
      for channel_index, channel in enumerate(channels):
         signal = reference_measurement[..., channel_index]
         difference = candidate_measurement[..., channel_index] - signal
         noise = 1.0 / snr[channel_index] if snr[channel_index] > 0.0 else np.nan
         sigma = fm.sigma_solar(signal, noise, production=True)
         record = _ratio_summary(difference, sigma, signal, minimum_signal=0.01)
         record.update({"measurement": "solar_reflectance_lambertian_nadir", "channel": int(channel),
                        "uncertainty": "production", "units": "reflectance", "surface_albedo": albedo,
                        "minimum_abs_reference_signal": 0.01})
         rows.append(record)
   return rows


def _thermal_and_mixed(reference: Dataset, candidate: Dataset, cloud: str) -> list[dict]:
   all_channels = _values(reference, "channel_id").astype(int)
   solar_channels, _ = _channel_indices(reference, "solar_channel_id")
   thermal_channels, thermal_all_indices = _channel_indices(reference, "thermal_channel_id")
   mixed = set(_values(reference, "mixed_channel_id").astype(int))
   solar_lookup = {channel: index for index, channel in enumerate(solar_channels)}
   thermal_lookup = {channel: index for index, channel in enumerate(thermal_channels)}
   planck = {int(channel): tuple(float(_values(reference, name)[index]) for name in ("B1", "B2", "T1", "T2"))
             for index, channel in enumerate(thermal_channels)}
   nedt = _values(reference, "nedt")
   refbt = _values(reference, "refbt")
   f0_reference = _values(reference, "F0")
   f0_candidate = _values(candidate, "F0")
   vza = _values(reference, "satellite_zenith")
   sza = _values(reference, "solar_zenith")
   reference_tdv = _values(reference, "T_dv", ("satellite_zenith", "optical_depth", "effective_radius", "channels"))
   candidate_tdv = _values(candidate, "T_dv", ("satellite_zenith", "optical_depth", "effective_radius", "channels"))
   reference_rdv = _values(reference, "R_dv", ("satellite_zenith", "optical_depth", "effective_radius", "channels"))
   candidate_rdv = _values(candidate, "R_dv", ("satellite_zenith", "optical_depth", "effective_radius", "channels"))
   reference_em = _values(reference, "E_md", ("satellite_zenith", "optical_depth", "effective_radius", "thermal_channels"))
   candidate_em = _values(candidate, "E_md", ("satellite_zenith", "optical_depth", "effective_radius", "thermal_channels"))
   reference_r0v = _values(reference, "R_0v", ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth",
                                                        "effective_radius", "solar_channels"))
   candidate_r0v = _values(candidate, "R_0v", ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth",
                                                        "effective_radius", "solar_channels"))
   rows = []
   for vza_index, angle in enumerate(vza):
      atmosphere = atmosphere_terms("terra", "modis", list(all_channels), planck, cloud, float(angle))
      for channel in thermal_channels:
         thermal_index = thermal_lookup[channel]
         all_index = int(np.flatnonzero(all_channels == channel)[0])
         entry = atmosphere[int(channel)]
         def radiance(tdv, em, rdv):
            return ((entry["Rbc_up"] * tdv + entry["B_c"] * em + entry["Rac_dwn"] * rdv) *
                    entry["Tac_lw"] + entry["Rac_up"])
         reference_radiance = radiance(reference_tdv[vza_index, :, :, all_index], reference_em[vza_index, :, :, thermal_index],
                                        reference_rdv[vza_index, :, :, all_index])
         candidate_radiance = radiance(candidate_tdv[vza_index, :, :, all_index], candidate_em[vza_index, :, :, thermal_index],
                                        candidate_rdv[vza_index, :, :, all_index])
         is_mixed = int(channel) in mixed
         if is_mixed:
            solar_index = solar_lookup[channel]
            secv = 1.0 / np.cos(np.deg2rad(angle))
            solar_transmission = entry["Tac"] ** (1.0 / np.cos(np.deg2rad(sza))) * entry["Tac"] ** secv
            reference_reflectance = np.transpose(reference_r0v[:, vza_index, :, :, :, solar_index], (1, 2, 3, 0))
            candidate_reflectance = np.transpose(candidate_r0v[:, vza_index, :, :, :, solar_index], (1, 2, 3, 0))
            reference_total = reference_radiance[None, :, :, None] + f0_reference[solar_index] * solar_transmission[:, None, None, None] * reference_reflectance
            candidate_total = candidate_radiance[None, :, :, None] + f0_candidate[solar_index] * solar_transmission[:, None, None, None] * candidate_reflectance
         else:
            reference_total = reference_radiance
            candidate_total = candidate_radiance
         reference_bt = fm.r2t(reference_total, planck[int(channel)])[0]
         candidate_bt = fm.r2t(candidate_total, planck[int(channel)])[0]
         sigma = fm.sigma_bt(reference_bt, nedt[thermal_index], refbt[thermal_index], planck[int(channel)],
                             mixed=is_mixed, production=True)
         record = _ratio_summary(candidate_bt - reference_bt, sigma, reference_bt)
         record.update({"measurement": "mixed_brightness_temperature" if is_mixed else "thermal_brightness_temperature",
                        "channel": int(channel), "uncertainty": "production", "units": "K",
                        "satellite_zenith": float(angle), "cloud_profile": cloud,
                        "atmosphere": "mid-latitude summer with repository MODTRAN channel profiles"})
         rows.append(record)
   return rows


def _write_csv(path: Path, rows: list[dict]) -> None:
   fieldnames = []
   for row in rows:
      for name in row:
         if name not in fieldnames:
            fieldnames.append(name)
   with path.open("w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=fieldnames)
      writer.writeheader()
      writer.writerows(rows)


def _plot(rows: list[dict], path: Path, label: str) -> None:
   os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "validation" / "generated" / "matplotlib-cache"))
   Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
   import matplotlib
   matplotlib.use("Agg")
   import matplotlib.pyplot as plt

   selected = [row for row in rows if row.get("uncertainty") == "production" and
               (row["measurement"] == "solar_reflectance_black_surface" or row["measurement"].endswith("brightness_temperature"))]
   grouped = {}
   for row in selected:
      key = (row["measurement"], row["channel"])
      grouped[key] = max(grouped.get(key, 0.0), row.get("maximum_sigma", 0.0))
   labels = [f"{measurement.split('_')[0]} ch{channel}" for measurement, channel in grouped]
   values = np.asarray(list(grouped.values()))
   order = np.argsort(values)
   fig, axis = plt.subplots(figsize=(11, max(5, 0.27 * len(values))))
   axis.barh(np.arange(len(values)), values[order])
   axis.set_yticks(np.arange(len(values)), [labels[index] for index in order])
   axis.axvline(0.1, color="0.5", linestyle="--", linewidth=1)
   axis.axvline(1.0, color="red", linestyle="--", linewidth=1)
   axis.set_xscale("log")
   axis.set_xlabel("maximum |V22 - V21| / production ORAC measurement uncertainty")
   axis.set_title(label + ": representative operational measurement impact")
   axis.grid(True, which="both", axis="x", alpha=0.25)
   fig.tight_layout()
   fig.savefig(path, dpi=180)
   plt.close(fig)


def main() -> int:
   parser = argparse.ArgumentParser(description=__doc__)
   parser.add_argument("reference", type=Path)
   parser.add_argument("candidate", type=Path)
   parser.add_argument("--cloud", choices=("liquid", "ice"), required=True)
   parser.add_argument("--output-directory", type=Path, required=True)
   parser.add_argument("--plot-directory", type=Path, required=True)
   parser.add_argument("--label")
   args = parser.parse_args()
   label = args.label or args.candidate.stem
   args.output_directory.mkdir(parents=True, exist_ok=True)
   args.plot_directory.mkdir(parents=True, exist_ok=True)
   with Dataset(args.reference, "r") as reference, Dataset(args.candidate, "r") as candidate:
      rows = _solar_black(reference, candidate)
      rows.extend(_solar_lambertian(reference, candidate))
      rows.extend(_thermal_and_mixed(reference, candidate, args.cloud))
   csv_path = args.output_directory / (label + "_radiance_impact.csv")
   json_path = args.output_directory / (label + "_radiance_impact.json")
   plot_path = args.plot_directory / (label + "_radiance_impact.png")
   _write_csv(csv_path, rows)
   json_path.write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
   _plot(rows, plot_path, label)
   print(json.dumps({"rows": len(rows), "csv": str(csv_path), "json": str(json_path), "plot": str(plot_path)}, indent=2))
   return 0


if __name__ == "__main__":
   raise SystemExit(main())
