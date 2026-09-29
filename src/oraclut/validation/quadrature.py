"""Tools for the LUT sampling (quadrature-convergence) experiments.

Everything here is pure numerics on arrays that a forward-model run has already
produced: candidate-grid construction, off-node point classification,
operational-style interpolation, error metrics and measurement-uncertainty
normalisation.  No radiative transfer happens in this module.

The interpolation "coordinates" mirror what ORAC does with a LUT: piecewise
linear in a transformed coordinate (``linear``, ``log10`` or the experimental
``log1p``).  Candidate grids must be subsets of the truth grid, so candidate
node values are the truth values at those nodes and on-node errors are exactly
zero; only off-node points measure the sampling error.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

C2_M_K = 1.4387769e-2  # second radiation constant, m K


# ---------------------------------------------------------------------------
# Grids
# ---------------------------------------------------------------------------

def linear_nodes(start: float, stop: float, count: int) -> np.ndarray:
    return np.linspace(start, stop, count, dtype=np.float64)


def log_nodes(start: float, stop: float, per_octave: int) -> np.ndarray:
    """Nodes at 2**(k/per_octave) covering [start, stop] (both powers of two included)."""

    k0 = int(round(np.log2(start) * per_octave))
    k1 = int(round(np.log2(stop) * per_octave))
    return 2.0 ** (np.arange(k0, k1 + 1) / per_octave)


def log10_nodes(start: float, stop: float, count: int) -> np.ndarray:
    return 10.0 ** np.linspace(np.log10(start), np.log10(stop), count)


def subset(truth: Sequence[float], every: int, *, keep_last: bool = True) -> np.ndarray:
    """Every ``every``-th truth node, always including the first and (optionally) last."""

    truth = np.asarray(truth, dtype=np.float64)
    picked = list(truth[::every])
    if keep_last and not np.isclose(picked[-1], truth[-1]):
        picked.append(truth[-1])
    return np.asarray(picked)


def assert_subset(candidate: Sequence[float], truth: Sequence[float], rtol: float = 1e-9) -> np.ndarray:
    """Return the truth indices of candidate nodes; raise if any node is not a truth node."""

    truth = np.asarray(truth, dtype=np.float64)
    indices = []
    for value in np.asarray(candidate, dtype=np.float64):
        matches = np.flatnonzero(np.isclose(truth, value, rtol=rtol, atol=0.0))
        if matches.size != 1:
            raise ValueError(f"candidate node {value!r} is not a unique truth node")
        indices.append(int(matches[0]))
    if indices != sorted(indices):
        raise ValueError("candidate nodes must be increasing")
    return np.asarray(indices)


# ---------------------------------------------------------------------------
# Coordinates and interpolation
# ---------------------------------------------------------------------------

def transform(x: np.ndarray, coordinate: str) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if coordinate == "linear":
        return x
    if coordinate == "log10":
        if np.any(x <= 0):
            raise ValueError("log10 coordinate requires positive values")
        return np.log10(x)
    if coordinate == "log1p":
        if np.any(x < 0):
            raise ValueError("log1p coordinate requires non-negative values")
        return np.log1p(x)
    if coordinate == "cos":
        return np.cos(np.deg2rad(x))
    raise ValueError(f"unknown interpolation coordinate {coordinate!r}")


def interpolate(candidate_x: Sequence[float], candidate_y: np.ndarray, truth_x: Sequence[float],
                coordinate: str = "linear", axis: int = 0) -> np.ndarray:
    """Piecewise-linear interpolation in the transformed coordinate along ``axis``.

    ``candidate_y`` holds the LUT values at ``candidate_x`` (any trailing shape);
    the result has the length of ``truth_x`` along ``axis``.  Points outside the
    candidate range are not extrapolated: they take the end value, as a LUT
    lookup would clamp them, and are flagged separately by :func:`classify`.
    """

    cx = transform(candidate_x, coordinate)
    tx = transform(truth_x, coordinate)
    y = np.moveaxis(np.asarray(candidate_y, dtype=np.float64), axis, 0)
    flat = y.reshape(y.shape[0], -1)
    out = np.empty((tx.size, flat.shape[1]))
    for column in range(flat.shape[1]):
        out[:, column] = np.interp(tx, cx, flat[:, column])
    return np.moveaxis(out.reshape((tx.size,) + y.shape[1:]), 0, axis)


# ---------------------------------------------------------------------------
# Point classification
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PointClasses:
    on_node: np.ndarray
    midpoint: np.ndarray
    edge: np.ndarray
    high_curvature: np.ndarray
    off_node: np.ndarray


def classify(truth_x: Sequence[float], candidate_x: Sequence[float], truth_y: np.ndarray | None = None,
             coordinate: str = "linear", *, midpoint_tolerance: float = 0.05,
             curvature_fraction: float = 0.2) -> PointClasses:
    """Boolean masks over the truth points.

    * ``on_node``: coincides with a candidate node;
    * ``midpoint``: within ``midpoint_tolerance`` (fraction of the interval) of the
      centre of a candidate interval, in the interpolation coordinate;
    * ``edge``: lies in the first or last candidate interval (off node);
    * ``high_curvature``: off-node points in the top ``curvature_fraction`` of
      |second difference| of the truth curve (``truth_y`` reduced to 1-D by
      taking the maximum over trailing axes), if ``truth_y`` is given.
    """

    tx = transform(truth_x, coordinate)
    cx = transform(candidate_x, coordinate)
    on_node = np.array([np.any(np.isclose(cx, value, rtol=1e-9, atol=1e-12)) for value in tx])
    segment = np.clip(np.searchsorted(cx, tx, side="right") - 1, 0, cx.size - 2)
    fraction = (tx - cx[segment]) / (cx[segment + 1] - cx[segment])
    midpoint = (~on_node) & (np.abs(fraction - 0.5) <= midpoint_tolerance)
    edge = (~on_node) & ((segment == 0) | (segment == cx.size - 2))
    off_node = ~on_node
    high = np.zeros_like(on_node)
    if truth_y is not None and tx.size >= 3:
        y = np.asarray(truth_y, dtype=np.float64)
        y = y.reshape(y.shape[0], -1)
        # |second difference| per point, normalised by the local spacing, max over observables
        d2 = np.zeros(tx.size)
        for column in range(y.shape[1]):
            g = np.gradient(np.gradient(y[:, column], tx), tx)
            d2 = np.maximum(d2, np.abs(g))
        threshold = np.quantile(d2[off_node], 1.0 - curvature_fraction) if off_node.any() else np.inf
        high = off_node & (d2 >= threshold)
    return PointClasses(on_node, midpoint, edge, high, off_node)


# ---------------------------------------------------------------------------
# Error metrics and measurement-uncertainty normalisation
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ErrorMetrics:
    points: int
    max_abs: float | None
    rms_abs: float | None
    mean_abs: float | None
    max_rel: float | None
    max_abs_location: float | None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def error_metrics(truth: np.ndarray, interpolated: np.ndarray, mask: np.ndarray | None = None,
                  *, truth_x: Sequence[float] | None = None, relative_floor: float = 1e-3,
                  axis: int = 0) -> ErrorMetrics:
    """Errors over the points selected by ``mask`` along ``axis`` (all other axes pooled)."""

    t = np.moveaxis(np.asarray(truth, dtype=np.float64), axis, 0)
    i = np.moveaxis(np.asarray(interpolated, dtype=np.float64), axis, 0)
    if mask is None:
        mask = np.ones(t.shape[0], dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return ErrorMetrics(0, None, None, None, None, None)
    d = np.abs(i[mask] - t[mask])
    finite = np.isfinite(d)
    if not finite.any():
        return ErrorMetrics(int(mask.sum()), None, None, None, None, None)
    flat_index = int(np.nanargmax(np.where(finite, d, -np.inf)))
    point_index = np.unravel_index(flat_index, d.shape)[0]
    location = float(np.asarray(truth_x)[mask][point_index]) if truth_x is not None else None
    meaningful = finite & (np.abs(t[mask]) >= relative_floor)
    max_rel = float(np.max(d[meaningful] / np.abs(t[mask][meaningful]))) if meaningful.any() else None
    return ErrorMetrics(int(mask.sum()), float(np.max(d[finite])), float(np.sqrt(np.mean(d[finite] ** 2))),
                        float(np.mean(d[finite])), max_rel, location)


@dataclass(frozen=True)
class MeasurementUncertainty:
    """Per-channel noise in the units of the operators being compared."""

    channel: int
    solar_reflectance_sigma: float | None      # reflectance units (legacy SAD NeFr, used as ORAC does)
    solar_relative_sigma: float | None         # 1/SNR
    thermal_nedt_k: float | None               # K
    thermal_kelvin_per_unit_operator: float | None  # dT for a unit change of an emissivity/transmission-like operator
    reference_temperature_k: float | None
    source: str

    def thermal_sigma_operator(self) -> float | None:
        """NEdT expressed as an equivalent error in an emissivity-like operator."""

        if self.thermal_nedt_k is None or not self.thermal_kelvin_per_unit_operator:
            return None
        return self.thermal_nedt_k / self.thermal_kelvin_per_unit_operator


def measurement_uncertainty(channel: int, *, wavelength_um: float, solar: bool, thermal: bool,
                            oldnefr: float | None, snr: float | None, nedt: float | None,
                            refbt: float | None) -> MeasurementUncertainty:
    """Build the normalisation from the instrument-file metadata.

    Solar: ``oldnefr`` is the legacy SAD noise-equivalent reflectance the ORAC
    retrieval uses directly, and ``1/snr`` gives a relative uncertainty.
    Thermal: with the Planck function, a change dE of an emissivity-like operator
    at scene temperature T changes brightness temperature by
    dT = dE * lambda * T^2 / c2, so ``nedt`` maps onto an operator error of
    ``nedt / (lambda T^2 / c2)``.
    """

    kelvin_per_unit = None
    if thermal and refbt:
        kelvin_per_unit = (wavelength_um * 1e-6) * refbt ** 2 / C2_M_K
    return MeasurementUncertainty(
        channel,
        oldnefr if (solar and oldnefr) else None,
        (1.0 / snr) if (solar and snr) else None,
        nedt if (thermal and nedt) else None,
        kelvin_per_unit,
        refbt if thermal else None,
        "meteosat-10_seviri_v1.inst: oldnefr (SAD NeFr), snr, nedt, refbt; Planck linearisation dT = dE*lambda*T^2/c2",
    )


def normalised_error(metrics: ErrorMetrics, uncertainty: MeasurementUncertainty, kind: str) -> dict[str, float | None]:
    """Error divided by measurement uncertainty; ``kind`` is 'solar' or 'thermal'."""

    out: dict[str, float | None] = {}
    if metrics.max_abs is None:
        return {"max_over_sigma": None}
    if kind == "solar":
        sigma = uncertainty.solar_reflectance_sigma
        out["max_over_sigma_reflectance"] = (metrics.max_abs / sigma) if sigma else None
        out["max_rel_over_inverse_snr"] = (
            (metrics.max_rel / uncertainty.solar_relative_sigma)
            if (metrics.max_rel is not None and uncertainty.solar_relative_sigma) else None)
    else:
        sigma = uncertainty.thermal_sigma_operator()
        out["max_over_sigma_thermal"] = (metrics.max_abs / sigma) if sigma else None
        out["max_equivalent_kelvin"] = (
            metrics.max_abs * uncertainty.thermal_kelvin_per_unit_operator
            if uncertainty.thermal_kelvin_per_unit_operator else None)
    return out


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------

def to_jsonable(value: Any) -> Any:
    if isinstance(value, ErrorMetrics):
        return value.as_dict()
    if hasattr(value, "__dataclass_fields__"):
        return to_jsonable(asdict(value))
    if isinstance(value, np.ndarray):
        return [to_jsonable(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    return value


def write_json(path: str | Path, payload: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(to_jsonable(payload), indent=1) + "\n")
    return path


# ---------------------------------------------------------------------------
# Refinement-experiment additions: sectioned truth grids, designed candidates,
# 2-D interpolation, acceptance criterion, independent test points, cache keys
# ---------------------------------------------------------------------------

def sectioned_nodes(sections: Sequence[Sequence[float]], decimals: int = 8) -> np.ndarray:
    """Concatenate ``(start, stop, step)`` sections into one increasing grid.

    Each section is inclusive of both ends; shared boundaries appear once.
    """

    nodes: list[float] = []
    for start, stop, step in sections:
        count = int(round((stop - start) / step))
        if count < 1 or not np.isclose(start + count * step, stop, rtol=0, atol=1e-9):
            raise ValueError(f"section {(start, stop, step)} is not an integer number of steps")
        nodes.extend(np.round(start + step * np.arange(count + 1), decimals).tolist())
    unique = np.unique(np.asarray(nodes, dtype=np.float64))
    return unique


def merge_nodes(*groups: Sequence[float], decimals: int = 12) -> np.ndarray:
    """Union of node groups, sorted, duplicates removed."""

    return np.unique(np.round(np.concatenate([np.asarray(g, dtype=np.float64).ravel() for g in groups]), decimals))


def equidistributed_nodes(x: Sequence[float], y: np.ndarray, count: int, coordinate: str = "linear",
                          *, floor: float = 0.0) -> np.ndarray:
    """Design ``count`` nodes so that the piecewise-linear error indicator is equidistributed.

    Linear interpolation error scales with h^2 |f''|; equidistributing
    sqrt(|f''|) (in the interpolation coordinate) over the cells makes the
    expected cell errors equal.  ``y`` may have trailing axes (observables,
    geometries) which are reduced by taking the maximum indicator.  The end
    nodes are always kept and every returned node is a truth node.
    """

    x = np.asarray(x, dtype=np.float64)
    if count < 2 or count > x.size:
        raise ValueError("count must be between 2 and the number of truth nodes")
    t = transform(x, coordinate)
    y2 = np.asarray(y, dtype=np.float64).reshape(x.size, -1)
    indicator = np.zeros(x.size)
    for column in range(y2.shape[1]):
        col = y2[:, column]
        if not np.all(np.isfinite(col)):
            continue
        indicator = np.maximum(indicator, np.abs(np.gradient(np.gradient(col, t), t)))
    density = np.sqrt(indicator) + floor
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(t))])
    targets = np.linspace(0.0, cumulative[-1], count)
    chosen = sorted({int(np.argmin(np.abs(cumulative - target))) for target in targets} | {0, x.size - 1})
    # ensure exactly `count` nodes: fill the largest gaps if rounding merged some
    while len(chosen) < count:
        gaps = np.diff(chosen)
        widest = int(np.argmax(gaps))
        if gaps[widest] < 2:
            break
        chosen.insert(widest + 1, chosen[widest] + gaps[widest] // 2)
    return x[np.asarray(chosen)]


def interpolate_2d(node_x: Sequence[float], node_y: Sequence[float], node_values: np.ndarray,
                   query_x: Sequence[float], query_y: Sequence[float],
                   coordinate_x: str = "linear", coordinate_y: str = "linear") -> np.ndarray:
    """Bilinear interpolation on a tensor grid, in transformed coordinates (LUT semantics).

    ``node_values`` has shape (len(node_x), len(node_y), ...); the result has
    shape (len(query_x), len(query_y), ...) — every query pair is evaluated.
    Queries outside the node range clamp to the edge, as a LUT lookup does.
    """

    first = interpolate(node_x, node_values, query_x, coordinate_x, axis=0)
    return interpolate(node_y, first, query_y, coordinate_y, axis=1)


def cell_midpoints(nodes: Sequence[float], coordinate: str = "linear") -> np.ndarray:
    """Midpoints of every cell in the interpolation coordinate (never a node)."""

    t = transform(nodes, coordinate)
    mid = 0.5 * (t[1:] + t[:-1])
    return _inverse(mid, coordinate)


def cell_fractions(nodes: Sequence[float], fraction: float, coordinate: str = "linear") -> np.ndarray:
    """Points a fixed fraction into every cell (e.g. 0.25 for an independent test set)."""

    if not 0.0 < fraction < 1.0:
        raise ValueError("fraction must be strictly inside (0, 1)")
    t = transform(nodes, coordinate)
    return _inverse(t[:-1] + fraction * np.diff(t), coordinate)


def _inverse(t: np.ndarray, coordinate: str) -> np.ndarray:
    if coordinate == "linear":
        return np.asarray(t, dtype=np.float64)
    if coordinate == "log10":
        return 10.0 ** np.asarray(t, dtype=np.float64)
    if coordinate == "log1p":
        return np.expm1(np.asarray(t, dtype=np.float64))
    if coordinate == "cos":
        return np.rad2deg(np.arccos(np.asarray(t, dtype=np.float64)))
    raise ValueError(coordinate)


@dataclass(frozen=True)
class Acceptance:
    accepted: bool
    limit_fraction: float
    worst_channel: int | None
    worst_ratio: float | None
    ratios: dict[str, float | None]


def acceptance(ratios: Mapping[Any, float | None], *, limit_fraction: float = 0.5) -> Acceptance:
    """Accept only if every channel's (max error / uncertainty) <= limit_fraction.

    A channel with no usable uncertainty (``None``) cannot be judged and is
    excluded from the decision but reported; an empty judgeable set is rejected.
    """

    judgeable = {str(k): v for k, v in ratios.items() if v is not None}
    if not judgeable:
        return Acceptance(False, limit_fraction, None, None, {str(k): v for k, v in ratios.items()})
    worst = max(judgeable, key=judgeable.get)
    return Acceptance(all(v <= limit_fraction for v in judgeable.values()), limit_fraction,
                      int(worst) if worst.isdigit() else None, judgeable[worst], {str(k): v for k, v in ratios.items()})


def cache_key(configuration: Mapping[str, Any], *, channel_set: Sequence[int], effective_radius: float,
              optical_depth: float, solar_zenith: Sequence[float], satellite_zenith: Sequence[float],
              relative_azimuth: Sequence[float]) -> str:
    """Deterministic key for one truth state: physical configuration + state coordinates.

    ``configuration`` must contain everything that changes the numbers
    (instrument, microphysics, atmosphere, srf_quad, streams, phase order,
    Rayleigh/gas, formulation); it is hashed together with the state so that
    a cache can never mix physically different runs.
    """

    import hashlib as _hashlib

    payload = json.dumps({
        "configuration": to_jsonable(dict(sorted(configuration.items()))),
        "channels": [int(c) for c in channel_set],
        "effective_radius": round(float(effective_radius), 10),
        "optical_depth": float(f"{float(optical_depth):.12g}"),
        "solar_zenith": [round(float(v), 6) for v in solar_zenith],
        "satellite_zenith": [round(float(v), 6) for v in satellite_zenith],
        "relative_azimuth": [round(float(v), 6) for v in relative_azimuth],
    }, sort_keys=True, separators=(",", ":"))
    return _hashlib.sha256(payload.encode()).hexdigest()
