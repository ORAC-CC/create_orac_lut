"""Compare two ORAC NetCDF LUTs by decoded physical value.

The reference and candidate paths are explicit.  Variables are aligned by
dimension name before comparison, so a harmless storage-order difference is
not reported as a physical difference.  CF scale/offset decoding and
fill/missing masking are applied explicitly, while packing attributes and mask
differences are recorded separately.  `valid_range` does not hide diagnostic
outliers.

Outputs are JSON plus flat CSV tables suitable for independent checking.  The
coordinate table contains a residual profile along every dimension of every
differing numerical variable.  For R_0v, complete direct-beam DISORT calls are
also ranked after collapsing only the view-angle dimensions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
from netCDF4 import Dataset, default_fillvals


RT_OPERATORS = ("T_dv", "T_dd", "R_dv", "R_dd", "R_0v", "R_0d", "T_0d", "T_00", "E_md")
PACKING_ATTRIBUTES = ("_FillValue", "missing_value", "scale_factor", "add_offset", "_Unsigned")
CALL_DIMENSIONS = ("solar_channels", "optical_depth", "effective_radius", "solar_zenith")


def _json_value(value: Any) -> Any:
   if isinstance(value, np.ma.MaskedArray):
      return _json_value(value.filled(np.nan))
   if isinstance(value, np.ndarray):
      return [_json_value(item) for item in value.tolist()]
   if isinstance(value, np.generic):
      return _json_value(value.item())
   if isinstance(value, bytes):
      return value.decode("utf-8", errors="replace")
   if isinstance(value, (list, tuple)):
      return [_json_value(item) for item in value]
   if isinstance(value, float) and not np.isfinite(value):
      return str(value)
   return value


def _attributes(owner: Any) -> dict[str, Any]:
   return {name: _json_value(owner.getncattr(name)) for name in owner.ncattrs()}


def _sha256(path: Path) -> str:
   digest = hashlib.sha256()
   with path.open("rb") as handle:
      for block in iter(lambda: handle.read(1024 * 1024), b""):
         digest.update(block)
   return digest.hexdigest()


def _file_record(path: Path, dataset: Dataset) -> dict[str, Any]:
   stat = path.stat()
   return {
      "path": str(path.resolve()),
      "size_bytes": stat.st_size,
      "mtime_ns": stat.st_mtime_ns,
      "mtime_utc": np.datetime64(stat.st_mtime_ns, "ns").astype("datetime64[ms]").astype(str) + "Z",
      "sha256": _sha256(path),
      "data_model": dataset.data_model,
      "disk_format": dataset.disk_format,
   }


def _dimension_record(dataset: Dataset) -> dict[str, int]:
   return {name: len(dimension) for name, dimension in dataset.dimensions.items()}


def _variable_schema(variable: Any) -> dict[str, Any]:
   return {
      "dimensions": list(variable.dimensions),
      "shape": list(variable.shape),
      "dtype": str(variable.dtype),
      "datatype": str(variable.datatype),
      "packing": {name: _json_value(variable.getncattr(name)) for name in PACKING_ATTRIBUTES if name in variable.ncattrs()},
      "attributes": _attributes(variable),
   }


def _attribute_differences(reference: dict[str, Any], candidate: dict[str, Any], scope: str) -> list[dict[str, Any]]:
   rows = []
   for name in sorted(set(reference) | set(candidate)):
      if reference.get(name) != candidate.get(name):
         rows.append({"scope": scope, "attribute": name, "reference": reference.get(name), "candidate": candidate.get(name)})
   return rows


def _align_candidate(reference_variable: Any, candidate_variable: Any, candidate_values: Any) -> tuple[Any, list[int]]:
   reference_dimensions = tuple(reference_variable.dimensions)
   candidate_dimensions = tuple(candidate_variable.dimensions)
   if reference_dimensions == candidate_dimensions:
      return candidate_values, list(range(len(candidate_dimensions)))
   if len(reference_dimensions) != len(candidate_dimensions) or set(reference_dimensions) != set(candidate_dimensions):
      raise ValueError("dimension names differ")
   order = [candidate_dimensions.index(name) for name in reference_dimensions]
   return np.transpose(candidate_values, order), order


def _coordinate_values(dataset: Dataset, dimension: str, length: int) -> np.ndarray:
   aliases = {"channels": "channel_id", "solar_channels": "solar_channel_id", "thermal_channels": "thermal_channel_id",
              "mixed_channels": "mixed_channel_id"}
   variable_name = dimension if dimension in dataset.variables else aliases.get(dimension)
   if variable_name is None or variable_name not in dataset.variables:
      return np.arange(length)
   variable = dataset.variables[variable_name]
   if variable.ndim != 1 or variable.shape[0] != length:
      return np.arange(length)
   values = _decoded(variable)
   if np.issubdtype(values.dtype, np.number):
      return _float_with_nan(values)
   return np.asarray(values)


def _location(dataset: Dataset, dimensions: tuple[str, ...], shape: tuple[int, ...], index: tuple[int, ...]) -> dict[str, Any]:
   location = {}
   for axis, dimension in enumerate(dimensions):
      coordinates = _coordinate_values(dataset, dimension, shape[axis])
      location[dimension] = {"index": int(index[axis]), "coordinate": _json_value(coordinates[index[axis]])}
   return location


def _decoded(variable: Any) -> np.ma.MaskedArray:
   """Decode packing while masking fill/missing values, not valid_range.

   netCDF4's combined automatic mask/scale also masks values outside
   ``valid_range``.  Such values are scientifically important diagnostics in
   these LUTs, so the CF packing operations are applied explicitly here.
   """

   variable.set_auto_maskandscale(False)
   raw = np.asarray(variable[:])
   if raw.dtype.kind not in "biufc":
      return np.ma.asarray(raw)
   mask = np.zeros(raw.shape, dtype=bool)
   markers = []
   for name in ("_FillValue", "missing_value"):
      if name in variable.ncattrs():
         markers.extend(np.atleast_1d(variable.getncattr(name)).tolist())
   fill_key = raw.dtype.str[1:]
   if fill_key in default_fillvals:
      markers.append(default_fillvals[fill_key])
   for marker in markers:
      if isinstance(marker, (float, np.floating)) and np.isnan(marker):
         mask |= np.isnan(raw)
      else:
         mask |= raw == np.asarray(marker, dtype=raw.dtype)
   values = raw
   if raw.dtype.kind == "i" and str(getattr(variable, "_Unsigned", "false")).lower() == "true":
      values = raw.view(np.dtype("u" + str(raw.dtype.itemsize)))
   if "scale_factor" in variable.ncattrs() or "add_offset" in variable.ncattrs():
      scale = variable.getncattr("scale_factor") if "scale_factor" in variable.ncattrs() else 1.0
      offset = variable.getncattr("add_offset") if "add_offset" in variable.ncattrs() else 0.0
      values = values * scale + offset
   return np.ma.array(values, mask=mask, copy=False)


def _float_with_nan(values: np.ma.MaskedArray) -> np.ndarray:
   result = np.asarray(np.ma.getdata(values), dtype=np.float64).copy()
   result[np.ma.getmaskarray(values)] = np.nan
   return result


def _numeric_comparison(name: str, reference_variable: Any, candidate_variable: Any, reference_dataset: Dataset,
                        candidate_dataset: Dataset) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
   reference_values = _decoded(reference_variable)
   candidate_values, order = _align_candidate(reference_variable, candidate_variable, _decoded(candidate_variable))
   if reference_values.shape != candidate_values.shape:
      return ({"variable": name, "status": "shape mismatch", "reference_shape": list(reference_values.shape),
               "candidate_shape": list(candidate_values.shape)}, [], {})

   reference_mask = np.ma.getmaskarray(reference_values)
   candidate_mask = np.ma.getmaskarray(candidate_values)
   reference = _float_with_nan(reference_values)
   candidate = _float_with_nan(candidate_values)
   finite_reference = np.isfinite(reference)
   finite_candidate = np.isfinite(candidate)
   valid = ~reference_mask & ~candidate_mask & finite_reference & finite_candidate
   difference = candidate - reference
   absolute = np.abs(difference)
   valid_absolute = absolute[valid]
   differing = valid & (reference != candidate)
   dimensions = tuple(reference_variable.dimensions)
   summary: dict[str, Any] = {
      "variable": name,
      "status": "ok" if np.any(valid) else "no jointly finite unmasked values",
      "reference_dimensions": list(dimensions),
      "candidate_dimensions": list(candidate_variable.dimensions),
      "candidate_transpose_order": order,
      "shape": list(reference.shape),
      "elements": int(reference.size),
      "jointly_valid_elements": int(np.count_nonzero(valid)),
      "differing_elements": int(np.count_nonzero(differing)),
      "exactly_equal": bool(np.array_equal(reference_values, candidate_values)),
      "reference_dtype": str(reference_variable.dtype),
      "candidate_dtype": str(candidate_variable.dtype),
      "reference_masked": int(np.count_nonzero(reference_mask)),
      "candidate_masked": int(np.count_nonzero(candidate_mask)),
      "mask_difference": int(np.count_nonzero(reference_mask != candidate_mask)),
      "reference_nan": int(np.count_nonzero(~reference_mask & np.isnan(reference))),
      "candidate_nan": int(np.count_nonzero(~candidate_mask & np.isnan(candidate))),
      "nan_mask_difference": int(np.count_nonzero(np.isnan(reference) != np.isnan(candidate))),
      "reference_inf": int(np.count_nonzero(~reference_mask & np.isinf(reference))),
      "candidate_inf": int(np.count_nonzero(~candidate_mask & np.isinf(candidate))),
      "packing_reference": _variable_schema(reference_variable)["packing"],
      "packing_candidate": _variable_schema(candidate_variable)["packing"],
   }
   coordinate_rows: list[dict[str, Any]] = []
   diagnostic = {"difference": difference, "absolute": absolute, "valid": valid, "dimensions": dimensions}
   if not np.any(valid):
      return summary, coordinate_rows, diagnostic

   scale = float(np.max(np.abs(reference[valid])))
   near_zero_threshold = max(1.0e-12, scale * 1.0e-6)
   relative_valid = valid & (np.abs(reference) > near_zero_threshold)
   symmetric_valid = valid & ((np.abs(reference) + np.abs(candidate)) > 2.0 * near_zero_threshold)
   near_zero = valid & (np.abs(reference) <= near_zero_threshold)
   max_flat = int(np.nanargmax(np.where(valid, absolute, np.nan)))
   max_index = tuple(int(item) for item in np.unravel_index(max_flat, reference.shape))
   percentiles = np.percentile(valid_absolute, [50.0, 90.0, 95.0, 99.0, 99.9, 100.0])
   summary.update({
      "reference_minimum": float(np.min(reference[valid])),
      "reference_maximum": float(np.max(reference[valid])),
      "candidate_minimum": float(np.min(candidate[valid])),
      "candidate_maximum": float(np.max(candidate[valid])),
      "maximum_absolute_difference": float(np.max(valid_absolute)),
      "mean_absolute_difference": float(np.mean(valid_absolute)),
      "rms_difference": float(np.sqrt(np.mean(difference[valid] ** 2))),
      "absolute_difference_percentiles": {str(p): float(value) for p, value in zip((50, 90, 95, 99, 99.9, 100), percentiles)},
      "near_zero_threshold": near_zero_threshold,
      "near_zero_elements": int(np.count_nonzero(near_zero)),
      "maximum_absolute_difference_near_zero": float(np.max(absolute[near_zero])) if np.any(near_zero) else None,
      "maximum_relative_difference": float(np.max(absolute[relative_valid] / np.abs(reference[relative_valid]))) if np.any(relative_valid) else None,
      "relative_difference_elements": int(np.count_nonzero(relative_valid)),
      "maximum_symmetric_relative_difference": float(np.max(2.0 * absolute[symmetric_valid] /
                                                               (np.abs(reference[symmetric_valid]) + np.abs(candidate[symmetric_valid]))))
                                                if np.any(symmetric_valid) else None,
      "maximum_absolute_index": list(max_index),
      "maximum_absolute_location": _location(reference_dataset, dimensions, reference.shape, max_index),
      "maximum_reference_value": float(reference[max_index]),
      "maximum_candidate_value": float(candidate[max_index]),
   })

   for axis, dimension in enumerate(dimensions):
      coordinates = _coordinate_values(reference_dataset, dimension, reference.shape[axis])
      for coordinate_index in range(reference.shape[axis]):
         selector = [slice(None)] * reference.ndim
         selector[axis] = coordinate_index
         selector_tuple = tuple(selector)
         group_valid = valid[selector_tuple]
         group_absolute = absolute[selector_tuple]
         group_difference = difference[selector_tuple]
         if np.any(group_valid):
            values = group_absolute[group_valid]
            coordinate_rows.append({
               "variable": name,
               "dimension": dimension,
               "coordinate_index": coordinate_index,
               "coordinate_value": _json_value(coordinates[coordinate_index]),
               "elements": int(np.count_nonzero(group_valid)),
               "differing_elements": int(np.count_nonzero(differing[selector_tuple] & group_valid)),
               "maximum_absolute_difference": float(np.max(values)),
               "p99_absolute_difference": float(np.percentile(values, 99.0)),
               "rms_difference": float(np.sqrt(np.mean(group_difference[group_valid] ** 2))),
               "mean_absolute_difference": float(np.mean(values)),
            })
   return summary, coordinate_rows, diagnostic


def _non_numeric_comparison(name: str, reference_variable: Any, candidate_variable: Any) -> dict[str, Any]:
   reference = _decoded(reference_variable)
   try:
      candidate, order = _align_candidate(reference_variable, candidate_variable, _decoded(candidate_variable))
   except ValueError:
      candidate = _decoded(candidate_variable)
      order = None
   equal = reference.shape == candidate.shape and np.array_equal(reference, candidate)
   return {
      "variable": name,
      "status": "ok" if equal else "non-numeric mismatch",
      "reference_dimensions": list(reference_variable.dimensions),
      "candidate_dimensions": list(candidate_variable.dimensions),
      "candidate_transpose_order": order,
      "reference_dtype": str(reference_variable.dtype),
      "candidate_dtype": str(candidate_variable.dtype),
      "elements": int(reference.size),
      "differing_elements": 0 if equal else None,
      "exactly_equal": bool(equal),
   }


def _direct_call_rows(reference_dataset: Dataset, candidate_dataset: Dataset) -> list[dict[str, Any]]:
   if "R_0v" not in reference_dataset.variables or "R_0v" not in candidate_dataset.variables:
      return []
   reference_variable = reference_dataset.variables["R_0v"]
   candidate_variable = candidate_dataset.variables["R_0v"]
   dimensions = tuple(reference_variable.dimensions)
   if not all(name in dimensions for name in CALL_DIMENSIONS):
      return []
   reference = _float_with_nan(_decoded(reference_variable))
   candidate_values, _ = _align_candidate(reference_variable, candidate_variable, _decoded(candidate_variable))
   candidate = _float_with_nan(candidate_values)
   call_axes = [dimensions.index(name) for name in CALL_DIMENSIONS]
   view_axes = [axis for axis in range(reference.ndim) if axis not in call_axes]
   order = call_axes + view_axes
   reference = np.transpose(reference, order)
   candidate = np.transpose(candidate, order)
   difference = candidate - reference
   call_shape = reference.shape[:len(CALL_DIMENSIONS)]
   flattened_shape = call_shape + (-1,)
   reference = reference.reshape(flattened_shape)
   candidate = candidate.reshape(flattened_shape)
   difference = difference.reshape(flattened_shape)
   finite_count = np.sum(np.isfinite(difference), axis=-1)
   safe_difference = np.where(np.isfinite(difference), difference, 0.0)
   maximum = np.max(np.abs(safe_difference), axis=-1)
   rms = np.sqrt(np.sum(safe_difference ** 2, axis=-1) / np.maximum(finite_count, 1))
   p99 = np.zeros(call_shape, dtype=np.float64)
   for index in np.ndindex(call_shape):
      finite = np.abs(difference[index][np.isfinite(difference[index])])
      p99[index] = np.percentile(finite, 99.0) if finite.size else np.nan
      if not finite.size:
         maximum[index] = np.nan
         rms[index] = np.nan
   ranked = np.argsort(np.nan_to_num(maximum, nan=-1.0).ravel())[::-1]
   rows = []
   for rank, flat_index in enumerate(ranked, start=1):
      index = tuple(int(item) for item in np.unravel_index(int(flat_index), call_shape))
      if not finite_count[index]:
         continue
      row: dict[str, Any] = {"rank": len(rows) + 1, "maximum_absolute_R_0v_difference": float(maximum[index]),
                             "p99_absolute_R_0v_difference": float(p99[index]), "rms_R_0v_difference": float(rms[index]),
                             "reference_R_0v_minimum": float(np.nanmin(reference[index])),
                             "reference_R_0v_maximum": float(np.nanmax(reference[index])),
                             "candidate_R_0v_minimum": float(np.nanmin(candidate[index])),
                             "candidate_R_0v_maximum": float(np.nanmax(candidate[index]))}
      for axis, dimension in enumerate(CALL_DIMENSIONS):
         coordinates = _coordinate_values(reference_dataset, dimension, call_shape[axis])
         row[dimension + "_index"] = index[axis]
         row[dimension] = _json_value(coordinates[index[axis]])
      rows.append(row)
   return rows[:1000]


def compare(reference_path: Path, candidate_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
   with Dataset(reference_path, "r") as reference, Dataset(candidate_path, "r") as candidate:
      attribute_differences = _attribute_differences(_attributes(reference), _attributes(candidate), "global")
      reference_names = list(reference.variables)
      candidate_names = list(candidate.variables)
      names = sorted(set(reference_names) | set(candidate_names))
      comparisons = []
      coordinate_rows = []
      diagnostics = {}
      for name in names:
         if name not in reference.variables or name not in candidate.variables:
            comparisons.append({"variable": name, "status": "missing", "in_reference": name in reference.variables,
                                "in_candidate": name in candidate.variables})
            continue
         reference_variable = reference.variables[name]
         candidate_variable = candidate.variables[name]
         attribute_differences.extend(_attribute_differences(_attributes(reference_variable), _attributes(candidate_variable),
                                                              "variable:" + name))
         numeric = np.issubdtype(reference_variable.dtype, np.number) and np.issubdtype(candidate_variable.dtype, np.number)
         if numeric:
            try:
               summary, rows, diagnostic = _numeric_comparison(name, reference_variable, candidate_variable, reference, candidate)
            except ValueError as error:
               summary, rows, diagnostic = ({"variable": name, "status": str(error)}, [], {})
            comparisons.append(summary)
            coordinate_rows.extend(rows)
            diagnostics[name] = diagnostic
         else:
            comparisons.append(_non_numeric_comparison(name, reference_variable, candidate_variable))
      result = {
         "reference": _file_record(reference_path, reference),
         "candidate": _file_record(candidate_path, candidate),
         "dimensions": {"reference": _dimension_record(reference), "candidate": _dimension_record(candidate),
                        "same_names_and_lengths": _dimension_record(reference) == _dimension_record(candidate)},
         "variable_order": {"reference": reference_names, "candidate": candidate_names, "same": reference_names == candidate_names},
         "variables": {"reference": {name: _variable_schema(variable) for name, variable in reference.variables.items()},
                       "candidate": {name: _variable_schema(variable) for name, variable in candidate.variables.items()}},
         "attribute_differences": attribute_differences,
         "comparisons": comparisons,
         "relative_difference_definition": "abs(candidate-reference)/abs(reference), only where abs(reference) exceeds max(1e-12, 1e-6*max(abs(reference)))",
         "symmetric_relative_difference_definition": "2*abs(candidate-reference)/(abs(reference)+abs(candidate)), with the same near-zero exclusion",
      }
      calls = _direct_call_rows(reference, candidate)
      result["direct_beam_calls"] = calls[:100]
   return result, coordinate_rows, calls


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
   path.parent.mkdir(parents=True, exist_ok=True)
   if not rows:
      path.write_text("")
      return
   fieldnames = []
   for row in rows:
      for name in row:
         if name not in fieldnames:
            fieldnames.append(name)
   with path.open("w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=fieldnames)
      writer.writeheader()
      writer.writerows(rows)


def _positive_floor(values: np.ndarray) -> float:
   positive = values[np.isfinite(values) & (values > 0.0)]
   return float(np.min(positive) / 2.0) if positive.size else 1.0e-30


def _make_plots(result: dict[str, Any], coordinate_rows: list[dict[str, Any]], plot_directory: Path, label: str) -> list[str]:
   os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "generated" / "matplotlib-cache"))
   Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
   import matplotlib
   matplotlib.use("Agg")
   import matplotlib.pyplot as plt

   plot_directory.mkdir(parents=True, exist_ok=True)
   written = []
   numerical = [row for row in result["comparisons"] if row.get("maximum_absolute_difference") is not None and
                row.get("differing_elements", 0) > 0]
   numerical.sort(key=lambda row: row["maximum_absolute_difference"])
   if numerical:
      names = [row["variable"] for row in numerical]
      maximum = np.asarray([row["maximum_absolute_difference"] for row in numerical])
      p99 = np.asarray([row["absolute_difference_percentiles"]["99"] for row in numerical])
      rms = np.asarray([row["rms_difference"] for row in numerical])
      floor = _positive_floor(np.concatenate((maximum, p99, rms)))
      y = np.arange(len(names))
      fig, axis = plt.subplots(figsize=(10, max(5, 0.32 * len(names))))
      axis.scatter(np.maximum(rms, floor), y, marker="o", label="RMS")
      axis.scatter(np.maximum(p99, floor), y, marker="s", label="99th percentile |difference|")
      axis.scatter(np.maximum(maximum, floor), y, marker="x", label="maximum |difference|")
      axis.set_xscale("log")
      axis.set_yticks(y, names)
      axis.set_xlabel("absolute decoded-value difference")
      axis.set_title(label + ": all differing numerical variables")
      axis.grid(True, which="both", axis="x", alpha=0.3)
      axis.legend()
      fig.tight_layout()
      path = plot_directory / (label + "_variable_overview.png")
      fig.savefig(path, dpi=180)
      plt.close(fig)
      written.append(str(path))

   for variable in RT_OPERATORS:
      rows = [row for row in coordinate_rows if row["variable"] == variable]
      dimensions = list(dict.fromkeys(row["dimension"] for row in rows))
      if not dimensions:
         continue
      columns = min(3, len(dimensions))
      nrows = (len(dimensions) + columns - 1) // columns
      fig, axes = plt.subplots(nrows, columns, figsize=(5.0 * columns, 3.7 * nrows), squeeze=False)
      for panel, dimension in zip(axes.flat, dimensions):
         group = [row for row in rows if row["dimension"] == dimension]
         maximum = np.asarray([row["maximum_absolute_difference"] for row in group])
         p99 = np.asarray([row["p99_absolute_difference"] for row in group])
         floor = _positive_floor(np.concatenate((maximum, p99)))
         x = np.arange(len(group))
         panel.plot(x, np.maximum(maximum, floor), marker="o", ms=3, label="maximum")
         panel.plot(x, np.maximum(p99, floor), marker=".", ms=3, label="99th percentile")
         panel.set_yscale("log")
         panel.set_title(dimension)
         panel.set_xlabel("coordinate value")
         panel.set_ylabel("absolute difference")
         step = max(1, len(group) // 8)
         ticks = np.arange(0, len(group), step)
         panel.set_xticks(ticks, [str(group[i]["coordinate_value"]) for i in ticks], rotation=35, ha="right")
         panel.grid(True, which="both", alpha=0.25)
      for panel in axes.flat[len(dimensions):]:
         panel.set_visible(False)
      axes.flat[0].legend()
      fig.suptitle(label + ": " + variable + " residual dependence")
      fig.tight_layout()
      path = plot_directory / (label + "_" + variable + "_coordinates.png")
      fig.savefig(path, dpi=180)
      plt.close(fig)
      written.append(str(path))
   return written


def main() -> int:
   parser = argparse.ArgumentParser(description=__doc__)
   parser.add_argument("reference", type=Path, help="reference V21 NetCDF path")
   parser.add_argument("candidate", type=Path, help="candidate V22 NetCDF path")
   parser.add_argument("--output-directory", type=Path, required=True)
   parser.add_argument("--plot-directory", type=Path)
   parser.add_argument("--label", help="output stem; default is the candidate filename without .nc")
   args = parser.parse_args()
   if not args.reference.is_file():
      parser.error("reference file does not exist: " + str(args.reference))
   if not args.candidate.is_file():
      parser.error("candidate file does not exist: " + str(args.candidate))
   label = args.label or args.candidate.stem
   output_directory = args.output_directory
   output_directory.mkdir(parents=True, exist_ok=True)
   result, coordinate_rows, call_rows = compare(args.reference, args.candidate)
   plot_paths = _make_plots(result, coordinate_rows, args.plot_directory, label) if args.plot_directory else []
   result["plots"] = plot_paths
   (output_directory / (label + "_comparison.json")).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
   comparison_rows = [{key: value for key, value in row.items() if not isinstance(value, (dict, list))}
                      for row in result["comparisons"]]
   _write_csv(output_directory / (label + "_variables.csv"), comparison_rows)
   _write_csv(output_directory / (label + "_coordinate_dependence.csv"), coordinate_rows)
   _write_csv(output_directory / (label + "_direct_beam_calls.csv"), call_rows)
   print(json.dumps({"label": label, "reference": str(args.reference.resolve()), "candidate": str(args.candidate.resolve()),
                     "variables": len(result["comparisons"]), "coordinate_rows": len(coordinate_rows),
                     "direct_beam_calls": len(call_rows), "plots": len(plot_paths)}, indent=2))
   return 0


if __name__ == "__main__":
   raise SystemExit(main())
