"""Compare a Python ORAC V2 product with an independently generated legacy one.

The module reports what the two files actually contain before it classifies
anything: structure first, then per-variable residuals, then formulation
specific checks.  It never modifies either product and never adjusts the
scientific calculation to force agreement.

Two independent results are produced:

* structural: ``PASS`` or ``FAIL``;
* numerical: ``GREEN``, ``AMBER`` or ``RED``.

The numerical thresholds are stated explicitly in :data:`NUMERICAL_BANDS` and
are specific to this comparison; they are calibrated on the already validated
Meteosat-10 SEVIRI liquid-water cloud case, where the legacy IDL and Python
radiative-transfer operators agree to about 1.3e-4 absolute on quantities of
order one.  They are a reporting aid, not a universal scientific tolerance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..io.lut import LutFile, read_lut
from .compare import ArrayComparison, compare_arrays


# Radiative-transfer operator variables of the V2 product.  These carry the
# scientific result and drive the numerical classification.
RT_OPERATORS = ("R_0v", "R_0d", "R_dv", "R_dd", "T_00", "T_0d", "T_dv", "T_dd", "E_md")
# Single-scattering optical properties.  The legacy writer is known to emit
# fill values here for singleton-channel products, so these are reported
# separately rather than driving the classification.
OPTICS_VARIABLES = (
    "extinction_coefficient", "extinction_coefficient_ratio",
    "single_scatter_albedo", "asymmetry_parameter",
    "average_volume_per_particle", "phase_function", "phase_moments",
    "legendre_moments", "phs", "amom",
)
CHANNEL_METADATA = (
    "channel_id", "central_wavelength", "central_wavenumber", "solar_channel_flag",
    "mixed_channel_flag", "thermal_channel_flag", "solar_channel_id", "F0", "oldf0",
    "oldf1", "snr", "oldnefr", "thermal_channel_id", "refbt", "nedt", "B1", "B2", "T1",
    "T2", "oldwvn", "oldb1", "oldb2", "oldt1", "oldt2", "oldnebt", "mixed_channel_id",
)
COORDINATE_VARIABLES = (
    "optical_depth", "effective_radius", "solar_zenith", "satellite_zenith",
    "relative_azimuth", "surface_pressure", "channel_id",
)

NUMERICAL_BANDS = {
    "green_maximum_absolute_difference": 1.0e-3,
    "amber_maximum_absolute_difference": 1.0e-2,
    "systematic_fraction": 0.5,
    "systematic_absolute_difference": 1.0e-3,
    "systematic_relative_difference": 5.0e-3,
    "systematic_relative_floor": 1.0e-2,
    "rationale": (
        "Calibrated on the validated Meteosat-10 SEVIRI liquid-water STG cloud "
        "case, whose legacy-versus-Python residuals peak near 1.3e-4 on "
        "reflectance/transmittance operators of order one. Specific to this "
        "comparison; not a universal scientific tolerance. A run is also held at "
        "AMBER when more than half of an operator's elements with |reference| >= 1e-2 "
        "differ by more than 5e-3 relatively, so a small but systematic bias cannot "
        "hide behind a small absolute residual."
    ),
}

IDL_FILL = 9.96921e36


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_ready(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [_json_ready(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, bytes):
        return value.decode("ascii", errors="replace")
    if isinstance(value, (tuple, list)):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


def _text_value(array: np.ndarray) -> str:
    """Join a NetCDF character array into a string."""

    flat = np.asarray(array).ravel().tolist()
    parts = []
    for item in flat:
        if isinstance(item, bytes):
            parts.append(item.decode("ascii", errors="replace"))
        else:
            parts.append(str(item))
    return "".join(parts).strip()


def structural_comparison(reference: LutFile, candidate: LutFile) -> dict[str, Any]:
    """Compare dimensions, ordering, variables, dtypes and coordinates."""

    mismatches: list[str] = []

    reference_dimensions = dict(reference.dimensions)
    candidate_dimensions = dict(candidate.dimensions)
    if list(reference_dimensions) != list(candidate_dimensions):
        mismatches.append(
            f"dimension order/set differs: reference {list(reference_dimensions)} "
            f"vs candidate {list(candidate_dimensions)}"
        )
    for name in sorted(set(reference_dimensions) | set(candidate_dimensions)):
        if reference_dimensions.get(name) != candidate_dimensions.get(name):
            mismatches.append(
                f"dimension {name}: reference {reference_dimensions.get(name)} "
                f"vs candidate {candidate_dimensions.get(name)}"
            )

    reference_variables = list(reference.variable_names)
    candidate_variables = list(candidate.variable_names)
    only_reference = [name for name in reference_variables if name not in candidate_variables]
    only_candidate = [name for name in candidate_variables if name not in reference_variables]
    if only_reference:
        mismatches.append(f"variables only in reference: {only_reference}")
    if only_candidate:
        mismatches.append(f"variables only in candidate: {only_candidate}")

    dimension_order: dict[str, dict[str, list[str]]] = {}
    dtypes: dict[str, dict[str, str]] = {}
    for name in reference_variables:
        if name not in candidate.variables:
            continue
        reference_order = list(reference.variable_dimensions[name])
        candidate_order = list(candidate.variable_dimensions[name])
        if reference_order != candidate_order:
            dimension_order[name] = {"reference": reference_order, "candidate": candidate_order}
            mismatches.append(
                f"variable {name} dimension order differs: {reference_order} vs {candidate_order}"
            )
        reference_dtype = reference.variable_dtypes[name]
        candidate_dtype = candidate.variable_dtypes[name]
        if reference_dtype != candidate_dtype:
            dtypes[name] = {"reference": reference_dtype, "candidate": candidate_dtype}
            mismatches.append(
                f"variable {name} dtype differs: {reference_dtype} vs {candidate_dtype}"
            )
        if reference.variables[name].shape != candidate.variables[name].shape:
            mismatches.append(
                f"variable {name} shape differs: {reference.variables[name].shape} "
                f"vs {candidate.variables[name].shape}"
            )

    coordinates: dict[str, Any] = {}
    for name in COORDINATE_VARIABLES:
        if name not in reference.variables or name not in candidate.variables:
            continue
        reference_values = np.asarray(reference.variables[name]).ravel()
        candidate_values = np.asarray(candidate.variables[name]).ravel()
        equal = (
            reference_values.shape == candidate_values.shape
            and bool(np.allclose(reference_values, candidate_values, rtol=0.0, atol=0.0, equal_nan=True))
        )
        close = (
            reference_values.shape == candidate_values.shape
            and bool(np.allclose(reference_values, candidate_values, rtol=1e-6, atol=1e-6, equal_nan=True))
        )
        coordinates[name] = {
            "reference": _json_ready(reference_values),
            "candidate": _json_ready(candidate_values),
            "bitwise_equal": equal,
            "close": close,
        }
        if not close:
            mismatches.append(f"coordinate {name} values differ")

    fill_and_finite: dict[str, Any] = {}
    for name in reference_variables:
        if name not in candidate.variables:
            continue
        entry = {}
        for label, source in (("reference", reference), ("candidate", candidate)):
            values = np.asarray(source.variables[name])
            if np.issubdtype(values.dtype, np.floating):
                finite = np.isfinite(values)
                fill = finite & (np.abs(values) >= 0.5 * IDL_FILL)
                entry[label] = {
                    "finite": int(np.count_nonzero(finite & ~fill)),
                    "fill": int(np.count_nonzero(fill)),
                    "nonfinite": int(np.count_nonzero(~finite)),
                    "total": int(values.size),
                }
            else:
                entry[label] = {"finite": int(values.size), "fill": 0, "nonfinite": 0, "total": int(values.size)}
        entry["fill_pattern_differs"] = entry["reference"]["fill"] != entry["candidate"]["fill"]
        fill_and_finite[name] = entry

    attribute_differences = []
    for name in sorted(set(reference.global_attributes) | set(candidate.global_attributes)):
        reference_value = _json_ready(reference.global_attributes.get(name))
        candidate_value = _json_ready(candidate.global_attributes.get(name))
        if reference_value != candidate_value:
            attribute_differences.append(
                {"scope": "global", "attribute": name,
                 "reference": reference_value, "candidate": candidate_value}
            )
    for variable in sorted(set(reference.variable_attributes) | set(candidate.variable_attributes)):
        reference_attributes = reference.variable_attributes.get(variable, {})
        candidate_attributes = candidate.variable_attributes.get(variable, {})
        for name in sorted(set(reference_attributes) | set(candidate_attributes)):
            reference_value = _json_ready(reference_attributes.get(name))
            candidate_value = _json_ready(candidate_attributes.get(name))
            if reference_value != candidate_value:
                attribute_differences.append(
                    {"scope": f"variable:{variable}", "attribute": name,
                     "reference": reference_value, "candidate": candidate_value}
                )

    return {
        "result": "PASS" if not mismatches else "FAIL",
        "mismatches": mismatches,
        "dimensions": {
            "reference": reference_dimensions,
            "candidate": candidate_dimensions,
            "order_reference": list(reference_dimensions),
            "order_candidate": list(candidate_dimensions),
        },
        "variables": {
            "reference": reference_variables,
            "candidate": candidate_variables,
            "only_in_reference": only_reference,
            "only_in_candidate": only_candidate,
        },
        "dimension_order_differences": dimension_order,
        "dtype_differences": dtypes,
        "coordinates": coordinates,
        "fill_and_finite": fill_and_finite,
        "attribute_differences": attribute_differences,
    }


def numerical_comparison(reference: LutFile, candidate: LutFile) -> list[ArrayComparison]:
    """Compare every variable present in both products."""

    comparisons = []
    for name in sorted(set(reference.variables) & set(candidate.variables)):
        units = reference.variable_attributes.get(name, {}).get("units")
        comparisons.append(
            compare_arrays(name, reference.variables[name], candidate.variables[name], units=units)
        )
    return comparisons


def _channel_axis(dims: Sequence[str]) -> tuple[int, str] | None:
    for axis, name in enumerate(dims):
        if name in ("channels", "solar_channels", "thermal_channels", "mixed_channels"):
            return axis, name
    return None


def _channel_ids(product: LutFile, dimension: str) -> list[int]:
    variable = {"channels": "channel_id", "solar_channels": "solar_channel_id",
                "thermal_channels": "thermal_channel_id", "mixed_channels": "mixed_channel_id"}[dimension]
    if variable not in product.variables:
        return []
    return [int(value) for value in np.asarray(product.variables[variable]).ravel()]


def per_channel_comparison(reference: LutFile, candidate: LutFile) -> list[dict[str, Any]]:
    """Residuals of every operator split by instrument channel.

    A two-channel product must not be judged by a single number: a visible
    channel that agrees and an infrared channel that does not would otherwise
    average into a misleading result.
    """

    rows: list[dict[str, Any]] = []
    for name in RT_OPERATORS + OPTICS_VARIABLES:
        if name not in reference.variables or name not in candidate.variables:
            continue
        dims = reference.variable_dimensions[name]
        located = _channel_axis(dims)
        if located is None or reference.variables[name].shape != candidate.variables[name].shape:
            continue
        axis, dimension = located
        ids = _channel_ids(reference, dimension)
        for index in range(reference.variables[name].shape[axis]):
            ref_slice = np.take(np.asarray(reference.variables[name]), index, axis=axis)
            cand_slice = np.take(np.asarray(candidate.variables[name]), index, axis=axis)
            result = compare_arrays(name, ref_slice, cand_slice)
            rows.append({
                "quantity": name,
                "channel_dimension": dimension,
                "channel_index": index,
                "channel_id": ids[index] if index < len(ids) else None,
                "status": result.status,
                "maximum_absolute_difference": result.maximum_absolute_difference,
                "rms_difference": result.rms_difference,
                "mean_absolute_difference": result.mean_absolute_difference,
                "maximum_relative_difference": result.maximum_relative_difference,
                "differing_elements": result.differing_elements,
                "compared_elements": result.compared_elements,
                "maximum_absolute_index": result.maximum_absolute_index,
                "reference_range": [result.reference_minimum, result.reference_maximum],
                "candidate_range": [result.candidate_minimum, result.candidate_maximum],
                "systematic_relative_fraction": systematic_relative_fraction(ref_slice, cand_slice),
            })
    return rows


def systematic_relative_fraction(
    reference: np.ndarray, candidate: np.ndarray, *, bands: Mapping[str, Any] = NUMERICAL_BANDS
) -> float | None:
    """Fraction of well-conditioned elements whose relative difference exceeds the band.

    Only elements with |reference| >= ``systematic_relative_floor`` are
    considered, so near-zero references never inflate the statistic.
    """

    ref = np.asarray(reference, dtype=np.float64)
    cand = np.asarray(candidate, dtype=np.float64)
    if ref.shape != cand.shape or not np.issubdtype(ref.dtype, np.number):
        return None
    usable = np.isfinite(ref) & np.isfinite(cand) & (np.abs(ref) < 0.5 * IDL_FILL) \
        & (np.abs(cand) < 0.5 * IDL_FILL) & (np.abs(ref) >= bands["systematic_relative_floor"])
    if not np.any(usable):
        return None
    relative = np.abs(cand[usable] - ref[usable]) / np.abs(ref[usable])
    return float(np.mean(relative > bands["systematic_relative_difference"]))


def channel_diagnostics(reference: LutFile, candidate: LutFile) -> dict[str, Any]:
    """Channel metadata and the thermal path, reported side by side.

    These are the quantities that distinguish an infrared channel from a
    visible one - solar/thermal flags, central wavelength and wavenumber,
    blackbody fit constants, thermal noise metadata and the emissivity
    operator - and they are checked explicitly rather than assumed to follow
    from the visible-channel result.
    """

    def values(product: LutFile, name: str) -> Any:
        if name not in product.variables:
            return None
        array = np.asarray(product.variables[name])
        if array.dtype.kind in ("S", "U", "O"):
            return _text_value(array)
        return _json_ready(array.ravel())

    metadata: dict[str, Any] = {}
    mismatched: list[str] = []
    for name in CHANNEL_METADATA:
        ref = values(reference, name)
        cand = values(candidate, name)
        if ref is None and cand is None:
            continue
        equal = ref == cand
        close = equal
        if not equal and isinstance(ref, list) and isinstance(cand, list) and len(ref) == len(cand):
            try:
                close = bool(np.allclose(np.asarray(ref, dtype=float), np.asarray(cand, dtype=float),
                                         rtol=1e-5, atol=1e-6))
            except (TypeError, ValueError):
                close = False
        metadata[name] = {"reference": ref, "candidate": cand, "equal": equal, "close": close}
        if not close:
            mismatched.append(name)

    def classify(product: LutFile) -> dict[str, Any]:
        ids = values(product, "channel_id") or []
        solar = values(product, "solar_channel_flag") or []
        thermal = values(product, "thermal_channel_flag") or []
        mixed = values(product, "mixed_channel_flag") or []
        wavelengths = values(product, "central_wavelength") or []
        out = []
        for index, channel in enumerate(ids):
            kind = "mixed (solar + thermal)" if index < len(mixed) and mixed[index] else (
                "solar" if index < len(solar) and solar[index] else (
                    "thermal" if index < len(thermal) and thermal[index] else "unflagged"))
            out.append({"channel_id": channel, "kind": kind,
                        "central_wavelength_micron": wavelengths[index] if index < len(wavelengths) else None})
        return {"channels": out,
                "has_solar_dimension": "solar_channels" in product.dimensions,
                "has_thermal_dimension": "thermal_channels" in product.dimensions,
                "has_E_md": "E_md" in product.variables}

    def filename_codes(product: LutFile) -> dict[str, Any]:
        stem = Path(product.path).name
        parts = stem.split("_")
        code = next((part[1:] for part in parts if part.startswith("a") and part[1:].isdigit()), None)
        mono = next((part for part in parts if part in ("m", "b")), None)
        gas = code is not None and code.startswith("1") and len(code) == 2
        return {"atmospheric_model_code": code,
                "gas_absorption_in_product_name": gas,
                "rayleigh_in_product_name": code not in ("00",) if code else None,
                "srf_treatment_in_product_name": {"m": "monochromatic", "b": "band-integrated"}.get(mono)}

    emissivity = {}
    for label, product in (("reference", reference), ("candidate", candidate)):
        if "E_md" in product.variables:
            array = np.asarray(product.variables["E_md"], dtype=np.float64)
            finite = array[np.isfinite(array) & (np.abs(array) < 0.5 * IDL_FILL)]
            emissivity[label] = {
                "dimensions": list(product.variable_dimensions["E_md"]),
                "shape": list(array.shape),
                "finite": int(finite.size), "total": int(array.size),
                "minimum": float(finite.min()) if finite.size else None,
                "maximum": float(finite.max()) if finite.size else None,
            }
    if "reference" in emissivity and "candidate" in emissivity:
        r_max = emissivity["reference"]["maximum"] or 0.0
        c_max = emissivity["candidate"]["maximum"] or 0.0
        if r_max and c_max:
            emissivity["maximum_ratio_candidate_over_reference"] = c_max / r_max

    return {
        "reference": classify(reference) | {"product_name_codes": filename_codes(reference)},
        "candidate": classify(candidate) | {"product_name_codes": filename_codes(candidate)},
        "metadata": metadata,
        "metadata_mismatches": mismatched,
        "emissivity": emissivity,
        "thermal_path_notes": (
            "The thermal emissivity operator E_md is produced by a DISORT run with the "
            "Planck source on (no beam, no isotropic illumination) over +/-0.5% of the "
            "channel centre wavenumber, normalised by the band-averaged Planck radiance "
            "at 250 K (cloud) or 270 K (aerosol), then SRF-weighted. Gas absorption enters "
            "only through the aerosol formulation's MODTRAN optical-depth profiles."
        ),
    }


def aerosol_checks(
    reference: LutFile,
    candidate: LutFile,
    *,
    expected_surface_pressure: Sequence[float] = (950.0, 1013.0, 1050.0),
    expected_states: int = 96,
    expected_components: int = 2,
    expected_component_order: Sequence[str] = ("waf", "saf"),
    expected_channels: Sequence[int] | None = None,
    expected_atmosphere: int = 2,
    gas: bool = True,
    rayleigh: bool = True,
    srf_quad: int = 1,
) -> dict[str, Any]:
    """Formulation-specific checks for the distributed-profile aerosol product."""

    checks: list[dict[str, Any]] = []

    def record(name: str, ok: bool | None, detail: str) -> None:
        checks.append({"check": name, "ok": ok, "detail": detail})

    for label, product in (("reference", reference), ("candidate", candidate)):
        pressure = product.variables.get("surface_pressure")
        if pressure is None:
            record(f"{label}: surface_pressure present", False, "variable absent")
            continue
        values = [float(value) for value in np.asarray(pressure).ravel()]
        record(
            f"{label}: surface_pressure values", values == list(expected_surface_pressure),
            f"{values} (expected {list(expected_surface_pressure)})",
        )
        operator = product.variable_dimensions.get("R_0v")
        if operator is None:
            record(f"{label}: R_0v present", False, "variable absent")
        else:
            position = operator.index("surface_pressure") if "surface_pressure" in operator else None
            record(
                f"{label}: pressure dimension position in R_0v",
                position is not None,
                f"R_0v dimensions {list(operator)}; surface_pressure at index {position}",
            )
        sizes = product.dimensions
        states = 1
        for name in ("optical_depth", "effective_radius", "solar_zenith",
                     "satellite_zenith", "relative_azimuth", "surface_pressure"):
            states *= int(sizes.get(name, 1))
        record(
            f"{label}: total aerosol grid states", states == expected_states,
            f"{states} (expected {expected_states}) from "
            + " x ".join(f"{name}={sizes.get(name)}" for name in (
                "optical_depth", "effective_radius", "solar_zenith",
                "satellite_zenith", "relative_azimuth", "surface_pressure"))
        )
        channels = product.variables.get("channel_id")
        if channels is not None:
            listed = [int(value) for value in np.asarray(channels).ravel()]
            expected = list(expected_channels) if expected_channels is not None else None
            record(f"{label}: channel selection", listed == expected if expected is not None else None,
                   f"channel_id {listed}" + (f" (expected {expected})" if expected is not None else ""))
        wavelength = product.variables.get("central_wavelength")
        if wavelength is not None:
            record(f"{label}: central wavelength", None,
                   f"{[float(value) for value in np.asarray(wavelength).ravel()]} micron")
        radii = product.variables.get("effective_radius")
        if radii is not None:
            record(f"{label}: effective radius grid", None,
                   f"{[float(value) for value in np.asarray(radii).ravel()]} micron")

    # The component structure is a property of the run configuration rather
    # than of the V2 product, which stores mixed optical properties only.
    record(
        "configuration: aerosol components", None,
        f"aerosol_a79.mm defines {expected_components} components in order "
        f"{list(expected_component_order)}; the V2 product stores the mixed "
        f"optical properties, so component structure is verified from the "
        f"microphysics input, not from either NetCDF file",
    )
    record("configuration: gas absorption", None, f"requested gas={gas}")
    record("configuration: Rayleigh scattering", None, f"requested rayleigh={rayleigh}")
    record("configuration: atmosphere", None, f"MODTRAN code {expected_atmosphere}")
    record("configuration: SRF treatment", None, f"srf_quad={srf_quad}")
    for label, product in (("reference", reference), ("candidate", candidate)):
        name = product.variables.get("instrument_filename")
        if name is not None:
            record(f"{label}: instrument file recorded in product", None, _text_value(name))

    failures = [check for check in checks if check["ok"] is False]
    return {"result": "PASS" if not failures else "FAIL", "checks": checks, "failures": failures}


def classify_numerical(
    comparisons: Sequence[ArrayComparison],
    *,
    bands: Mapping[str, Any] = NUMERICAL_BANDS,
    systematic_relative: Mapping[str, float | None] | None = None,
) -> dict[str, Any]:
    """Classify the radiative-transfer residuals as GREEN, AMBER or RED.

    ``systematic_relative`` maps operator name to the fraction of
    well-conditioned elements whose relative difference exceeds the band (see
    :func:`systematic_relative_fraction`); a fraction above
    ``systematic_fraction`` holds the result at AMBER even when every absolute
    residual is small.
    """

    reasons: list[str] = []
    operators = [item for item in comparisons if item.quantity in RT_OPERATORS]
    if not operators:
        return {"result": "RED", "reasons": ["no radiative-transfer operator variables in common"],
                "bands": dict(bands)}

    worst = 0.0
    for item in operators:
        if item.status == "shape mismatch":
            reasons.append(f"{item.quantity}: shape mismatch {item.reference_shape} vs {item.candidate_shape}")
            continue
        if item.status == "no finite values":
            reasons.append(f"{item.quantity}: no finite values to compare")
            continue
        if item.nan_inf_difference.get("candidate_nan") or item.nan_inf_difference.get("candidate_inf"):
            reasons.append(f"{item.quantity}: candidate contains nonfinite values")
        if item.nan_inf_difference.get("fill_mask_difference"):
            reasons.append(
                f"{item.quantity}: fill-value pattern differs "
                f"({item.nan_inf_difference['fill_mask_difference']} elements)"
            )
        if item.maximum_absolute_difference is not None:
            worst = max(worst, item.maximum_absolute_difference)

    red = any(
        "shape mismatch" in reason or "nonfinite" in reason or "no finite values" in reason
        for reason in reasons
    )
    systematic = []
    for item in operators:
        if item.compared_elements and item.differing_elements is not None and item.maximum_absolute_difference:
            fraction = item.differing_elements / item.compared_elements
            if (
                fraction > bands["systematic_fraction"]
                and item.maximum_absolute_difference > bands["systematic_absolute_difference"]
            ):
                systematic.append(
                    f"{item.quantity}: {fraction:.0%} of elements differ with maximum "
                    f"{item.maximum_absolute_difference:.3e}"
                )

    relative_bias = []
    for name, fraction in (systematic_relative or {}).items():
        if fraction is not None and fraction > bands["systematic_fraction"]:
            relative_bias.append(
                f"{name}: {fraction:.0%} of elements with |reference| >= "
                f"{bands['systematic_relative_floor']:.0e} differ by more than "
                f"{bands['systematic_relative_difference']:.0e} relatively (systematic bias)"
            )

    if red or worst > bands["amber_maximum_absolute_difference"] or systematic:
        result = "RED"
    elif worst > bands["green_maximum_absolute_difference"] or relative_bias:
        result = "AMBER"
    else:
        result = "GREEN"
    if worst:
        reasons.append(
            f"worst radiative-transfer residual {worst:.3e} "
            f"(GREEN <= {bands['green_maximum_absolute_difference']:.0e}, "
            f"AMBER <= {bands['amber_maximum_absolute_difference']:.0e})"
        )
    reasons.extend(systematic)
    reasons.extend(relative_bias)
    return {
        "result": result,
        "worst_operator_absolute_difference": worst,
        "reasons": reasons,
        "bands": dict(bands),
    }


def compare_products(
    reference_path: str | Path,
    candidate_path: str | Path,
    *,
    formulation: str = "aerosol",
    legacy_log: str | Path | None = None,
    legacy_stderr: str | Path | None = None,
    python_log: str | Path | None = None,
    repeatability: Mapping[str, Any] | None = None,
    expected_channels: Sequence[int] | None = None,
) -> dict[str, Any]:
    """Compare a legacy reference product with a Python product."""

    reference_path = Path(reference_path)
    candidate_path = Path(candidate_path)
    reference = read_lut(reference_path, validate_v21=False)
    candidate = read_lut(candidate_path, validate_v21=False)

    structural = structural_comparison(reference, candidate)
    comparisons = numerical_comparison(reference, candidate)
    systematic = {
        name: systematic_relative_fraction(reference.variables[name], candidate.variables[name])
        for name in RT_OPERATORS if name in reference.variables and name in candidate.variables
        and reference.variables[name].shape == candidate.variables[name].shape
    }
    numerical = classify_numerical(comparisons, systematic_relative=systematic)
    per_channel = per_channel_comparison(reference, candidate)
    report: dict[str, Any] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reference": {
            "path": str(reference_path), "sha256": sha256(reference_path),
            "size_bytes": reference_path.stat().st_size, "role": "legacy IDL reference",
        },
        "candidate": {
            "path": str(candidate_path), "sha256": sha256(candidate_path),
            "size_bytes": candidate_path.stat().st_size, "role": "Python generator product",
        },
        "formulation": formulation,
        "structural": structural,
        "numerical": numerical,
        "comparisons": [_json_ready(asdict(item)) for item in comparisons],
        "systematic_relative_fraction": systematic,
        "per_channel": _json_ready(per_channel),
        "channel_diagnostics": channel_diagnostics(reference, candidate),
        "variable_groups": {
            "radiative_transfer_operators": [
                item.quantity for item in comparisons if item.quantity in RT_OPERATORS
            ],
            "optical_properties": [
                item.quantity for item in comparisons if item.quantity in OPTICS_VARIABLES
            ],
        },
    }
    if formulation == "aerosol":
        report["aerosol_checks"] = aerosol_checks(reference, candidate, expected_channels=expected_channels)
    report["warnings"] = collect_warnings(legacy_log, python_log, legacy_stderr=legacy_stderr)
    if repeatability is not None:
        report["repeatability"] = _json_ready(dict(repeatability))
    return report


WARNING_PATTERNS = (
    "UPBEAM--SGECO says matrix near singular",
    "SGECO says matrix near singular",
    "near singular",
    "Floating underflow",
    "Floating overflow",
    "Floating illegal operand",
    "Floating divide by 0",
    "% Program caused arithmetic error",
    "DISORT",
    "Warning",
    "WARNING",
)


def collect_warnings(
    legacy_log: str | Path | None,
    python_log: str | Path | None,
    *,
    legacy_stderr: str | Path | None = None,
) -> dict[str, Any]:
    """Extract notable warning lines from run logs, without judging them."""

    def scan(path: str | Path | None) -> dict[str, Any]:
        if path is None:
            return {"log": None, "lines": [], "counts": {}}
        path = Path(path)
        if not path.is_file():
            return {"log": str(path), "lines": [], "counts": {}, "note": "log not found"}
        lines = []
        counts: dict[str, int] = {}
        for line in path.read_text(errors="replace").splitlines():
            for pattern in WARNING_PATTERNS:
                if pattern in line:
                    counts[pattern] = counts.get(pattern, 0) + 1
                    if len(lines) < 200:
                        lines.append(line.strip())
                    break
        return {"log": str(path), "lines": lines, "counts": counts}

    legacy = scan(legacy_log)
    legacy_err = scan(legacy_stderr)
    python = scan(python_log)
    combined = {}
    for source in (legacy["counts"], legacy_err["counts"]):
        for pattern, count in source.items():
            combined[pattern] = combined.get(pattern, 0) + count
    singular = combined.get("UPBEAM--SGECO says matrix near singular", 0) or combined.get(
        "SGECO says matrix near singular", 0)
    underflow = combined.get("Floating underflow", 0)
    return {
        "legacy": legacy,
        "legacy_stderr": legacy_err,
        "python": python,
        "upbeam_sgeco_near_singular_occurrences": singular,
        "floating_underflow_occurrences": underflow,
        "interpretation": (
            "A DISORT near-singular warning or an IDL floating-underflow report "
            "accompanied by a complete, finite NetCDF product is recorded, not treated "
            "as failure; a missing or nonfinite product is."
        ),
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    """Render the comparison as a readable report."""

    lines: list[str] = []
    add = lines.append
    formulation = report.get("formulation", "unknown")
    add(f"# Legacy versus Python comparison: {formulation} formulation")
    add("")
    add(f"Generated {report['generated_utc']}.")
    add("")
    add("| Role | Path | SHA-256 | Bytes |")
    add("| --- | --- | --- | --- |")
    for key in ("reference", "candidate"):
        entry = report[key]
        add(f"| {entry['role']} | `{entry['path']}` | `{entry['sha256']}` | {entry['size_bytes']:,} |")
    add("")
    add(f"**Structural: {report['structural']['result']}** · "
        f"**Numerical: {report['numerical']['result']}**")
    add("")

    add("## Structure")
    add("")
    structural = report["structural"]
    if structural["mismatches"]:
        add("Mismatches:")
        add("")
        for item in structural["mismatches"]:
            add(f"- {item}")
    else:
        add("Dimensions, dimension order, variable set, shapes, dtypes and coordinate "
            "values are identical.")
    add("")
    add("| Dimension | Reference | Candidate |")
    add("| --- | --- | --- |")
    dimensions = structural["dimensions"]
    for name in dict(dimensions["reference"]) | dict(dimensions["candidate"]):
        add(f"| {name} | {dimensions['reference'].get(name, '-')} | {dimensions['candidate'].get(name, '-')} |")
    add("")
    coordinates = structural.get("coordinates", {})
    if coordinates:
        add("| Coordinate | Values (reference) | Identical |")
        add("| --- | --- | --- |")
        for name, entry in coordinates.items():
            values = entry["reference"]
            shown = values if len(values) <= 8 else values[:4] + ["..."] + values[-2:]
            add(f"| {name} | {shown} | {'yes' if entry['bitwise_equal'] else ('close' if entry['close'] else 'NO')} |")
        add("")

    add("## Numerical residuals")
    add("")
    add("| Variable | Status | max abs | RMS | mean abs | max rel | differing/compared | max-abs index |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- |")

    def number(value: Any) -> str:
        if value is None:
            return "-"
        return f"{value:.3e}"

    groups = report.get("variable_groups", {})
    operators = set(groups.get("radiative_transfer_operators", []))
    optics = set(groups.get("optical_properties", []))
    ordered = (
        [item for item in report["comparisons"] if item["quantity"] in operators]
        + [item for item in report["comparisons"] if item["quantity"] in optics]
        + [item for item in report["comparisons"]
           if item["quantity"] not in operators and item["quantity"] not in optics]
    )
    for item in ordered:
        differing = (
            f"{item.get('differing_elements')}/{item.get('compared_elements')}"
            if item.get("compared_elements") is not None else "-"
        )
        add(
            f"| `{item['quantity']}` | {item['status']} | {number(item.get('maximum_absolute_difference'))} "
            f"| {number(item.get('rms_difference'))} | {number(item.get('mean_absolute_difference'))} "
            f"| {number(item.get('maximum_relative_difference'))} | {differing} "
            f"| {item.get('maximum_absolute_index')} |"
        )
    add("")
    add("Relative differences are computed against the reference magnitude and are "
        "meaningless where the reference is near zero; read the absolute columns first.")
    add("")

    per_channel = report.get("per_channel", [])
    if per_channel:
        add("## Per-channel residuals")
        add("")
        add("Each instrument channel is judged on its own; a two-channel product is never "
            "summarised by one number.")
        add("")
        add("| Variable | Channel | Kind | Status | max abs | RMS | max rel | systematic rel. fraction | reference range |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        kinds = {}
        for side in ("reference",):
            for entry in report.get("channel_diagnostics", {}).get(side, {}).get("channels", []):
                kinds[entry["channel_id"]] = entry["kind"]
        for row in per_channel:
            fraction = row.get("systematic_relative_fraction")
            rng = row.get("reference_range") or [None, None]
            add(
                f"| `{row['quantity']}` | {row.get('channel_id')} | {kinds.get(row.get('channel_id'), '-')} "
                f"| {row['status']} | {number(row.get('maximum_absolute_difference'))} "
                f"| {number(row.get('rms_difference'))} | {number(row.get('maximum_relative_difference'))} "
                f"| {('%.0f%%' % (100 * fraction)) if fraction is not None else '-'} "
                f"| [{number(rng[0])}, {number(rng[1])}] |"
            )
        add("")

    diagnostics = report.get("channel_diagnostics")
    if diagnostics:
        add("## Channel and thermal-path diagnostics")
        add("")
        for side in ("reference", "candidate"):
            entry = diagnostics.get(side, {})
            add(f"**{side}**: " + "; ".join(
                f"channel {c['channel_id']} {c['kind']} at {c['central_wavelength_micron']:.4f} um"
                if c.get("central_wavelength_micron") is not None else f"channel {c['channel_id']} {c['kind']}"
                for c in entry.get("channels", [])
            ) + f"; solar dimension {entry.get('has_solar_dimension')}, thermal dimension "
              f"{entry.get('has_thermal_dimension')}, E_md present {entry.get('has_E_md')}; "
              f"product-name codes {entry.get('product_name_codes')}")
            add("")
        mismatched = diagnostics.get("metadata_mismatches", [])
        add("Channel metadata (flags, wavelengths, wavenumbers, blackbody constants, thermal "
            "noise terms) " + ("**differ** for: " + ", ".join(f"`{m}`" for m in mismatched)
                                if mismatched else "agree between the two products (to 1e-5 relative)."))
        add("")
        if mismatched:
            add("| Variable | Reference | Candidate |")
            add("| --- | --- | --- |")
            for name in mismatched:
                item = diagnostics["metadata"][name]
                add(f"| `{name}` | {item['reference']} | {item['candidate']} |")
            add("")
        emissivity = diagnostics.get("emissivity", {})
        if emissivity:
            add("Emissivity operator `E_md`:")
            add("")
            for side in ("reference", "candidate"):
                if side in emissivity:
                    e = emissivity[side]
                    add(f"- {side}: dims {e['dimensions']}, finite {e['finite']}/{e['total']}, "
                        f"range [{number(e['minimum'])}, {number(e['maximum'])}]")
            if "maximum_ratio_candidate_over_reference" in emissivity:
                add(f"- maximum ratio candidate/reference: {emissivity['maximum_ratio_candidate_over_reference']:.4f}")
            add("")
        add(diagnostics.get("thermal_path_notes", ""))
        add("")

    fill = {
        name: entry for name, entry in structural.get("fill_and_finite", {}).items()
        if entry.get("fill_pattern_differs")
    }
    if fill:
        add("### Fill-value differences")
        add("")
        add("One product carries fill values where the other carries finite values. "
            "This is reported, not classified as a scientific error.")
        add("")
        add("| Variable | Reference finite/fill | Candidate finite/fill |")
        add("| --- | --- | --- |")
        for name, entry in fill.items():
            add(f"| `{name}` | {entry['reference']['finite']}/{entry['reference']['fill']} "
                f"| {entry['candidate']['finite']}/{entry['candidate']['fill']} |")
        add("")

    if "aerosol_checks" in report:
        add("## Aerosol-specific checks")
        add("")
        add(f"Result: **{report['aerosol_checks']['result']}**")
        add("")
        add("| Check | Outcome | Detail |")
        add("| --- | --- | --- |")
        for check in report["aerosol_checks"]["checks"]:
            outcome = {True: "ok", False: "FAIL", None: "info"}[check["ok"]]
            add(f"| {check['check']} | {outcome} | {check['detail']} |")
        add("")

    warnings = report.get("warnings", {})
    add("## Warnings")
    add("")
    occurrences = warnings.get("upbeam_sgeco_near_singular_occurrences", 0)
    add(f"`UPBEAM--SGECO says matrix near singular` occurrences in the legacy logs: {occurrences}.")
    add(f"`Floating underflow` reports in the legacy logs: {warnings.get('floating_underflow_occurrences', 0)}.")
    add("")
    add(warnings.get("interpretation", ""))
    add("")
    for key in ("legacy", "legacy_stderr", "python"):
        entry = warnings.get(key) or {}
        if entry.get("counts"):
            add(f"### {key} log `{entry.get('log')}`")
            add("")
            for pattern, count in entry["counts"].items():
                add(f"- `{pattern}`: {count}")
            add("")

    if "repeatability" in report:
        add("## Legacy repeatability")
        add("")
        repeatability = report["repeatability"]
        for key, value in repeatability.items():
            add(f"- {key}: {value}")
        add("")

    add("## Classification")
    add("")
    add(f"- structural: **{report['structural']['result']}**")
    add(f"- numerical: **{report['numerical']['result']}**")
    add("")
    for reason in report["numerical"]["reasons"]:
        add(f"- {reason}")
    add("")
    add(f"Bands used: {report['numerical']['bands']['rationale']}")
    add("")
    return "\n".join(lines) + "\n"


def write_report(report: Mapping[str, Any], out_dir: str | Path) -> tuple[Path, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "comparison.json"
    markdown_path = out_dir / "comparison.md"
    json_path.write_text(json.dumps(_json_ready(dict(report)), indent=2) + "\n")
    markdown_path.write_text(render_markdown(report))
    return json_path, markdown_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m oraclut.validation.product_comparison",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("reference", type=Path, help="Legacy IDL reference product")
    parser.add_argument("candidate", type=Path, help="Python generator product")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory for comparison.json and comparison.md")
    parser.add_argument("--formulation", choices=("aerosol", "cloud"), default="aerosol")
    parser.add_argument("--legacy-log", type=Path, help="Legacy stdout log to scan for warnings")
    parser.add_argument("--legacy-stderr", type=Path, help="Legacy stderr log to scan for warnings")
    parser.add_argument("--expected-channels", help="Comma-separated channel ids the products should contain")
    parser.add_argument("--python-log", type=Path, help="Python run log to scan for warnings")
    parser.add_argument("--repeatability-json", type=Path, help="JSON file describing a legacy repeat run")
    args = parser.parse_args(argv)

    repeatability = None
    if args.repeatability_json and args.repeatability_json.is_file():
        repeatability = json.loads(args.repeatability_json.read_text())
    expected = None
    if args.expected_channels:
        expected = [int(part) for part in args.expected_channels.split(",") if part.strip()]
    report = compare_products(
        args.reference, args.candidate, formulation=args.formulation,
        legacy_log=args.legacy_log, legacy_stderr=args.legacy_stderr, python_log=args.python_log,
        repeatability=repeatability, expected_channels=expected,
    )
    json_path, markdown_path = write_report(report, args.out_dir)
    print(f"structural: {report['structural']['result']}")
    print(f"numerical:  {report['numerical']['result']}")
    worst = report["numerical"].get("worst_operator_absolute_difference")
    if worst is not None:
        print(f"worst radiative-transfer residual: {worst:.3e}")
    print(f"wrote {json_path}")
    print(f"wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
