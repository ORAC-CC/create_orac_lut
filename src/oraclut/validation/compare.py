"""Stage-oriented comparison helpers for compact legacy reference outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from ..io.lut import LutFile, read_lut


def _json_value(value: Any) -> Any:
    """Convert NetCDF/NumPy attribute values to stable JSON values."""

    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, bytes):
        return value.decode("ascii", errors="replace")
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _attribute_differences(
    reference: Mapping[str, Any], candidate: Mapping[str, Any], *, scope: str
) -> list[dict[str, Any]]:
    differences = []
    for name in sorted(set(reference) | set(candidate)):
        reference_value = _json_value(reference.get(name))
        candidate_value = _json_value(candidate.get(name))
        if reference_value != candidate_value:
            differences.append({
                "scope": scope,
                "attribute": name,
                "reference": reference_value,
                "candidate": candidate_value,
            })
    return differences


def _structural_comparison(reference: LutFile, candidate: LutFile) -> dict[str, Any]:
    attribute_differences = _attribute_differences(
        reference.global_attributes, candidate.global_attributes, scope="global"
    )
    for name in sorted(set(reference.variable_attributes) | set(candidate.variable_attributes)):
        attribute_differences.extend(_attribute_differences(
            reference.variable_attributes.get(name, {}),
            candidate.variable_attributes.get(name, {}),
            scope=f"variable:{name}",
        ))
    exact = (
        list(reference.dimensions) == list(candidate.dimensions)
        and reference.dimensions == candidate.dimensions
        and reference.variable_names == candidate.variable_names
        and reference.variable_dimensions == candidate.variable_dimensions
        and reference.variable_dtypes == candidate.variable_dtypes
        and not attribute_differences
    )
    return {
        "status": "exact" if exact else "mismatch",
        "dimensions": {
            "reference": dict(reference.dimensions),
            "candidate": dict(candidate.dimensions),
            "exact": list(reference.dimensions) == list(candidate.dimensions)
            and reference.dimensions == candidate.dimensions,
        },
        "variable_names": {
            "reference": list(reference.variable_names),
            "candidate": list(candidate.variable_names),
            "exact": reference.variable_names == candidate.variable_names,
        },
        "variable_dimensions_exact": reference.variable_dimensions == candidate.variable_dimensions,
        "variable_dtypes_exact": reference.variable_dtypes == candidate.variable_dtypes,
        "attribute_differences": attribute_differences,
    }


def _serialization_comparison(reference: Path, candidate: Path) -> dict[str, Any]:
    reference_bytes = reference.read_bytes()
    candidate_bytes = candidate.read_bytes()
    return {
        "reference_sha256": hashlib.sha256(reference_bytes).hexdigest(),
        "candidate_sha256": hashlib.sha256(candidate_bytes).hexdigest(),
        "byte_equal": reference_bytes == candidate_bytes,
    }


@dataclass(frozen=True)
class ArrayComparison:
    quantity: str
    units: str | None
    reference_shape: tuple[int, ...] | None
    candidate_shape: tuple[int, ...] | None
    reference_dtype: str | None
    candidate_dtype: str | None
    reference_minimum: float | None
    candidate_minimum: float | None
    reference_maximum: float | None
    candidate_maximum: float | None
    maximum_absolute_difference: float | None
    rms_difference: float | None
    maximum_relative_difference: float | None
    maximum_absolute_index: tuple[int, ...] | None
    maximum_relative_index: tuple[int, ...] | None
    nan_inf_difference: dict[str, int]
    status: str
    mean_absolute_difference: float | None = None
    compared_elements: int | None = None
    differing_elements: int | None = None


def compare_arrays(
    quantity: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    units: str | None = None,
) -> ArrayComparison:
    reference = np.asarray(reference)
    candidate = np.asarray(candidate)
    if reference.shape != candidate.shape:
        return ArrayComparison(
            quantity, units, reference.shape, candidate.shape,
            str(reference.dtype), str(candidate.dtype), None, None, None, None,
            None, None, None, None, None,
            {"reference_nan": int(np.isnan(reference).sum()) if np.issubdtype(reference.dtype, np.floating) else 0,
             "candidate_nan": int(np.isnan(candidate).sum()) if np.issubdtype(candidate.dtype, np.floating) else 0,
             "reference_inf": int(np.isinf(reference).sum()) if np.issubdtype(reference.dtype, np.floating) else 0,
             "candidate_inf": int(np.isinf(candidate).sum()) if np.issubdtype(candidate.dtype, np.floating) else 0},
            "shape mismatch",
        )
    if not (np.issubdtype(reference.dtype, np.number) and np.issubdtype(candidate.dtype, np.number)):
        equal = np.array_equal(reference, candidate)
        return ArrayComparison(
            quantity, units, reference.shape, candidate.shape,
            str(reference.dtype), str(candidate.dtype), None, None, None, None,
            0.0 if equal else None, 0.0 if equal else None, 0.0 if equal else None,
            None, None, {"value_difference": 0 if equal else 1},
            "ok" if equal else "non-numeric value mismatch",
        )
    ref = reference.astype(np.float64, copy=False)
    cand = candidate.astype(np.float64, copy=False)
    ref_nan = np.isnan(ref)
    cand_nan = np.isnan(cand)
    ref_inf = np.isinf(ref)
    cand_inf = np.isinf(cand)
    reference_fill = np.isfinite(ref) & (np.abs(ref) >= 0.5 * 9.96921e36)
    candidate_fill = np.isfinite(cand) & (np.abs(cand) >= 0.5 * 9.96921e36)
    finite = np.isfinite(ref) & np.isfinite(cand) & ~reference_fill & ~candidate_fill
    if not np.any(finite):
        return ArrayComparison(
            quantity, units, reference.shape, candidate.shape,
            str(reference.dtype), str(candidate.dtype), None, None, None, None,
            None, None, None, None, None,
            {"reference_nan": int(ref_nan.sum()), "candidate_nan": int(cand_nan.sum()),
             "reference_inf": int(ref_inf.sum()), "candidate_inf": int(cand_inf.sum()),
             "reference_fill": int(reference_fill.sum()),
             "candidate_fill": int(candidate_fill.sum()),
             "fill_mask_difference": int(np.count_nonzero(reference_fill != candidate_fill)),
             "nan_mask_difference": int(np.count_nonzero(ref_nan != cand_nan)),
             "inf_mask_difference": int(np.count_nonzero(ref_inf != cand_inf))},
            "no finite values",
        )
    difference = np.abs(ref[finite] - cand[finite])
    denominator = np.maximum(np.abs(ref[finite]), np.finfo(float).eps)
    relative = difference / denominator
    finite_indices = np.flatnonzero(finite)
    max_position = int(finite_indices[int(np.argmax(difference))])
    index = tuple(int(x) for x in np.unravel_index(max_position, reference.shape))
    relative_position = int(finite_indices[int(np.argmax(relative))])
    relative_index = tuple(int(x) for x in np.unravel_index(relative_position, reference.shape))
    return ArrayComparison(
        quantity, units, reference.shape, candidate.shape,
        str(reference.dtype), str(candidate.dtype),
        float(np.min(ref[finite])), float(np.min(cand[finite])),
        float(np.max(ref[finite])), float(np.max(cand[finite])),
        float(np.max(difference)), float(np.sqrt(np.mean((ref[finite] - cand[finite]) ** 2))),
        float(np.max(relative)), index, relative_index,
        {"reference_nan": int(ref_nan.sum()), "candidate_nan": int(cand_nan.sum()),
         "reference_inf": int(ref_inf.sum()), "candidate_inf": int(cand_inf.sum()),
         "reference_fill": int(reference_fill.sum()), "candidate_fill": int(candidate_fill.sum()),
         "fill_mask_difference": int(np.count_nonzero(reference_fill != candidate_fill)),
         "nan_mask_difference": int(np.count_nonzero(ref_nan != cand_nan)),
         "inf_mask_difference": int(np.count_nonzero(ref_inf != cand_inf))},
        "ok",
        mean_absolute_difference=float(np.mean(difference)),
        compared_elements=int(difference.size),
        differing_elements=int(np.count_nonzero(difference)),
    )


def compare_lut_files(reference_path: str | Path, candidate_path: str | Path) -> dict[str, Any]:
    reference = read_lut(reference_path, validate_v21=False)
    candidate = read_lut(candidate_path, validate_v21=False)
    quantities = sorted(set(reference.variables) | set(candidate.variables))
    comparisons = []
    for quantity in quantities:
        if quantity not in reference.variables or quantity not in candidate.variables:
            comparisons.append({"quantity": quantity, "status": "missing from one file"})
            continue
        units = reference.variable_attributes.get(quantity, {}).get("units")
        comparisons.append(asdict(compare_arrays(quantity, reference.variables[quantity], candidate.variables[quantity], units=units)))
    return {
        "reference": str(reference.path),
        "candidate": str(candidate.path),
        "structural": _structural_comparison(reference, candidate),
        "serialization": _serialization_comparison(reference.path, candidate.path),
        "dimension_order": {
            "reference": reference.variable_dimensions,
            "candidate": candidate.variable_dimensions,
        },
        "comparisons": comparisons,
        "phase_function_normalisation": "not present in V2 NetCDF; inspect captured phs/amom arrays with legacy convention",
        "legendre_convention": "not inferred; captured amom arrays retain current legpexp ordering",
        "optical_depth_reference": "0.55 micron, as recorded by the legacy STG capture manifest",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    result = compare_lut_files(args.reference, args.candidate)
    text = json.dumps(result, indent=2, default=list) + "\n"
    print(text, end="")
    if args.json_out:
        args.json_out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
