#!/usr/bin/env python3
"""Compare two ORAC V2 LUTs with Nakajima–King style diagrams."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from netCDF4 import Dataset


class ComparisonError(ValueError):
    """Raised when two LUTs cannot be compared without interpolation."""


@dataclass(frozen=True)
class LutMetadata:
    path: Path
    platform: str
    instrument: str
    coordinates: dict[str, np.ndarray]
    channel_ids: np.ndarray
    solar_channel_ids: np.ndarray
    channel_wavelengths: dict[int, float]
    r0v_dimensions: tuple[str, ...]


@dataclass(frozen=True)
class Selection:
    x_channel: int
    y_channel: int
    solar_zenith: float
    satellite_zenith: float
    relative_azimuth: float
    surface_pressure: float | None = None


@dataclass(frozen=True)
class NakajimaKingData:
    metadata: LutMetadata
    selection: Selection
    optical_depth: np.ndarray
    effective_radius: np.ndarray
    x_reflectance: np.ndarray
    y_reflectance: np.ndarray


_CORE_COORDINATES = (
    "optical_depth",
    "effective_radius",
    "satellite_zenith",
    "solar_zenith",
    "relative_azimuth",
)


def _decode_text(variable) -> str:
    values = np.asarray(variable[:]).reshape(-1)
    if values.dtype.kind == "S":
        return b"".join(values.tolist()).decode("ascii", errors="replace").strip()
    if values.dtype.kind == "U":
        return "".join(values.tolist()).strip()
    return str(values.tolist()).strip()


def _require_variable(dataset: Dataset, name: str):
    if name not in dataset.variables:
        raise ComparisonError(f"{dataset.filepath()}: missing required variable {name!r}")
    return dataset.variables[name]


def _numeric_vector(dataset: Dataset, name: str) -> np.ndarray:
    variable = _require_variable(dataset, name)
    values = np.ma.asarray(variable[:]).filled(np.nan)
    if values.ndim != 1:
        raise ComparisonError(f"{dataset.filepath()}: variable {name!r} must be one-dimensional")
    return np.asarray(values)


def read_lut_metadata(path: str | Path) -> LutMetadata:
    """Read the small metadata and coordinate portion of one ORAC V2 LUT."""

    lut_path = Path(path)
    if not lut_path.is_file():
        raise FileNotFoundError(f"ORAC LUT does not exist: {lut_path}")
    with Dataset(lut_path, "r") as dataset:
        coordinates = {name: _numeric_vector(dataset, name) for name in _CORE_COORDINATES}
        if "surface_pressure" in dataset.variables:
            coordinates["surface_pressure"] = _numeric_vector(dataset, "surface_pressure")
        channel_ids = _numeric_vector(dataset, "channel_id").astype(int)
        solar_channel_ids = _numeric_vector(dataset, "solar_channel_id").astype(int)
        wavelengths = _numeric_vector(dataset, "central_wavelength").astype(float)
        if channel_ids.size != wavelengths.size:
            raise ComparisonError(
                f"{lut_path}: channel_id and central_wavelength have different lengths"
            )
        if len(set(channel_ids.tolist())) != channel_ids.size:
            raise ComparisonError(f"{lut_path}: channel_id contains duplicates")
        r0v = _require_variable(dataset, "R_0v")
        required_dimensions = set(_CORE_COORDINATES) | {"solar_channels"}
        missing = required_dimensions.difference(r0v.dimensions)
        if missing:
            raise ComparisonError(
                f"{lut_path}: R_0v is missing dimensions {sorted(missing)}; got {r0v.dimensions}"
            )
        r0v_dimensions = tuple(r0v.dimensions)
        platform = _decode_text(_require_variable(dataset, "platform"))
        instrument = _decode_text(_require_variable(dataset, "instrument"))

    return LutMetadata(
        path=lut_path,
        platform=platform,
        instrument=instrument,
        coordinates=coordinates,
        channel_ids=channel_ids,
        solar_channel_ids=solar_channel_ids,
        channel_wavelengths=dict(zip(channel_ids.tolist(), wavelengths.tolist())),
        r0v_dimensions=r0v_dimensions,
    )


def check_compatible(lut_a: LutMetadata, lut_b: LutMetadata) -> None:
    """Reject comparisons that would require interpolation or remapping."""

    problems: list[str] = []
    if (lut_a.platform, lut_a.instrument) != (lut_b.platform, lut_b.instrument):
        problems.append(
            "platform/instrument differ: "
            f"{lut_a.platform}/{lut_a.instrument} versus {lut_b.platform}/{lut_b.instrument}"
        )
    coordinate_names = set(lut_a.coordinates) | set(lut_b.coordinates)
    for name in sorted(coordinate_names):
        if name not in lut_a.coordinates or name not in lut_b.coordinates:
            problems.append(f"coordinate {name!r} is present in only one LUT")
        elif not np.array_equal(lut_a.coordinates[name], lut_b.coordinates[name]):
            problems.append(f"coordinate {name!r} differs")
    if not np.array_equal(lut_a.channel_ids, lut_b.channel_ids):
        problems.append("channel_id arrays differ")
    if not np.array_equal(lut_a.solar_channel_ids, lut_b.solar_channel_ids):
        problems.append("solar_channel_id arrays differ")
    wavelengths_a = np.asarray([lut_a.channel_wavelengths[int(c)] for c in lut_a.channel_ids])
    wavelengths_b = np.asarray([lut_b.channel_wavelengths[int(c)] for c in lut_b.channel_ids])
    # Historical versions can differ by a few 1e-5 microns because the stored
    # centre was recomputed from the same SRF. Channel IDs remain authoritative;
    # reject wavelength changes larger than 0.1 nm.
    if not np.allclose(wavelengths_a, wavelengths_b, rtol=0.0, atol=1.0e-4):
        problems.append("central_wavelength arrays differ by more than 0.0001 microns")
    if lut_a.r0v_dimensions != lut_b.r0v_dimensions:
        problems.append(f"R_0v dimensions differ: {lut_a.r0v_dimensions} versus {lut_b.r0v_dimensions}")
    if problems:
        raise ComparisonError(
            "LUTs are incompatible; no interpolation is performed:\n  - " + "\n  - ".join(problems)
        )


def _coordinate_value(values: np.ndarray, requested: float | None, name: str) -> float:
    if requested is None:
        return float(values[0])
    matches = np.flatnonzero(np.isclose(values.astype(float), requested, rtol=0.0, atol=1.0e-6))
    if matches.size != 1:
        choices = ", ".join(f"{float(value):g}" for value in values)
        raise ComparisonError(f"{name}={requested:g} is unavailable; choose exactly one of: {choices}")
    return float(values[matches[0]])


def resolve_selection(
    metadata: LutMetadata,
    *,
    x_channel: int | None,
    y_channel: int | None,
    solar_zenith: float | None,
    satellite_zenith: float | None,
    relative_azimuth: float | None,
    surface_pressure: float | None,
) -> Selection:
    """Resolve optional CLI choices against coordinates stored in the LUT."""

    solar_channels = metadata.solar_channel_ids.tolist()
    if len(solar_channels) < 2:
        raise ComparisonError(f"{metadata.path}: at least two solar channels are required")
    if (x_channel is None) != (y_channel is None):
        raise ComparisonError("--x-channel and --y-channel must be supplied together")
    if x_channel is None:
        if len(solar_channels) != 2:
            choices = ", ".join(
                f"{channel} ({metadata.channel_wavelengths[channel]:g} µm)"
                for channel in solar_channels
            )
            raise ComparisonError(
                "multiple solar-channel pairs are possible; select a scientifically appropriate pair "
                f"with --x-channel and --y-channel. Available solar channels: {choices}"
            )
        x_channel, y_channel = (int(channel) for channel in solar_channels)
    if x_channel == y_channel:
        raise ComparisonError("x and y channels must be different")
    for channel in (x_channel, y_channel):
        if channel not in solar_channels:
            raise ComparisonError(
                f"channel {channel} is not a solar channel; available solar channels: {solar_channels}"
            )

    pressure_coordinates = metadata.coordinates.get("surface_pressure")
    if pressure_coordinates is None:
        if surface_pressure is not None:
            raise ComparisonError("--surface-pressure was supplied but the LUT has no surface-pressure axis")
        selected_pressure = None
    else:
        selected_pressure = _coordinate_value(pressure_coordinates, surface_pressure, "surface_pressure")

    return Selection(
        x_channel=int(x_channel),
        y_channel=int(y_channel),
        solar_zenith=_coordinate_value(metadata.coordinates["solar_zenith"], solar_zenith, "solar_zenith"),
        satellite_zenith=_coordinate_value(
            metadata.coordinates["satellite_zenith"], satellite_zenith, "satellite_zenith"
        ),
        relative_azimuth=_coordinate_value(
            metadata.coordinates["relative_azimuth"], relative_azimuth, "relative_azimuth"
        ),
        surface_pressure=selected_pressure,
    )


def _index_for_value(values: np.ndarray, value: float) -> int:
    matches = np.flatnonzero(np.isclose(values.astype(float), value, rtol=0.0, atol=1.0e-6))
    if matches.size != 1:  # selection resolution should make this unreachable
        raise ComparisonError(f"could not uniquely locate coordinate value {value:g}")
    return int(matches[0])


def _read_channel_grid(metadata: LutMetadata, selection: Selection, channel: int) -> np.ndarray:
    with Dataset(metadata.path, "r") as dataset:
        variable = dataset.variables["R_0v"]
        indices: list[int | slice] = []
        remaining_dimensions: list[str] = []
        coordinate_values = {
            "relative_azimuth": selection.relative_azimuth,
            "satellite_zenith": selection.satellite_zenith,
            "solar_zenith": selection.solar_zenith,
            "surface_pressure": selection.surface_pressure,
        }
        for dimension in variable.dimensions:
            if dimension in ("optical_depth", "effective_radius"):
                indices.append(slice(None))
                remaining_dimensions.append(dimension)
            elif dimension == "solar_channels":
                matches = np.flatnonzero(metadata.solar_channel_ids == channel)
                if matches.size != 1:
                    raise ComparisonError(f"cannot uniquely map solar channel {channel} in {metadata.path}")
                indices.append(int(matches[0]))
            elif dimension in coordinate_values and coordinate_values[dimension] is not None:
                indices.append(
                    _index_for_value(metadata.coordinates[dimension], float(coordinate_values[dimension]))
                )
            else:
                raise ComparisonError(
                    f"{metadata.path}: unsupported unselected R_0v dimension {dimension!r}; "
                    "this utility does not silently collapse dimensions"
                )
        values = np.ma.asarray(variable[tuple(indices)]).filled(np.nan)
    if set(remaining_dimensions) != {"optical_depth", "effective_radius"}:
        raise ComparisonError(
            f"{metadata.path}: selected R_0v does not reduce to optical_depth/effective_radius"
        )
    order = [remaining_dimensions.index("optical_depth"), remaining_dimensions.index("effective_radius")]
    return np.asarray(values, dtype=float).transpose(order)


def extract_nakajima_king(metadata: LutMetadata, selection: Selection) -> NakajimaKingData:
    """Extract two channel-reflectance grids for one geometry."""

    return NakajimaKingData(
        metadata=metadata,
        selection=selection,
        optical_depth=metadata.coordinates["optical_depth"].astype(float),
        effective_radius=metadata.coordinates["effective_radius"].astype(float),
        x_reflectance=_read_channel_grid(metadata, selection, selection.x_channel),
        y_reflectance=_read_channel_grid(metadata, selection, selection.y_channel),
    )


def calculate_displacement(lut_a: NakajimaKingData, lut_b: NakajimaKingData):
    """Return x/y reflectance displacement and Euclidean magnitude."""

    if lut_a.x_reflectance.shape != lut_b.x_reflectance.shape:
        raise ComparisonError("selected reflectance grids have different shapes")
    dx = lut_b.x_reflectance - lut_a.x_reflectance
    dy = lut_b.y_reflectance - lut_a.y_reflectance
    return dx, dy, np.hypot(dx, dy)


def _channel_label(metadata: LutMetadata, channel: int) -> str:
    wavelength = metadata.channel_wavelengths[channel]
    return f"Channel {channel} reflectance ({wavelength:g} µm)"


def _representative_indices(values: np.ndarray, maximum: int, *, logarithmic: bool) -> tuple[int, ...]:
    """Choose readable labels spanning an actual LUT coordinate grid."""

    values = np.asarray(values, dtype=float)
    if values.size <= maximum:
        return tuple(range(values.size))
    selected: set[int] = set()
    targets_available = maximum
    if logarithmic and values[0] == 0.0:
        selected.add(0)
        targets_available -= 1
    if logarithmic and np.all(values >= 0.0) and np.any(values > 0.0):
        positive_indices = np.flatnonzero(values > 0.0)
        positive_values = values[positive_indices]
        targets = np.geomspace(positive_values.min(), positive_values.max(), targets_available)
        for target in targets:
            offset = int(np.argmin(np.abs(np.log(positive_values) - np.log(target))))
            selected.add(int(positive_indices[offset]))
    else:
        targets = np.linspace(values.min(), values.max(), targets_available)
        for target in targets:
            selected.add(int(np.argmin(np.abs(values - target))))
    return tuple(sorted(selected))


def _label_point(x_values: np.ndarray, y_values: np.ndarray, fraction: float) -> tuple[float, float] | None:
    finite = np.flatnonzero(np.isfinite(x_values) & np.isfinite(y_values))
    if finite.size == 0:
        return None
    location = finite[int(round(fraction * (finite.size - 1)))]
    return float(x_values[location]), float(y_values[location])


def _plot_grid(axis, data: NakajimaKingData, *, linestyle: str, annotate: bool, alpha: float = 1.0):
    tau_color = "#2166ac"
    radius_color = "#b2182b"
    tau_labels = set(_representative_indices(data.optical_depth, 7, logarithmic=True))
    radius_label_indices = _representative_indices(data.effective_radius, 6, logarithmic=False)
    radius_labels = set(radius_label_indices)
    radius_label_fractions = dict(
        zip(radius_label_indices, np.linspace(0.46, 0.9, len(radius_label_indices)))
    )
    label_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.3}
    for index, optical_depth in enumerate(data.optical_depth):
        axis.plot(
            data.x_reflectance[index, :], data.y_reflectance[index, :],
            color=tau_color, linestyle=linestyle, linewidth=1.0, alpha=alpha,
        )
        label_location = _label_point(data.x_reflectance[index, :], data.y_reflectance[index, :], 0.86)
        if annotate and index in tau_labels and label_location is not None:
            axis.annotate(
                f"τ={optical_depth:g}",
                label_location,
                fontsize=6.5,
                color=tau_color,
                xytext=(2, 2),
                textcoords="offset points",
                bbox=label_box,
            )
    for index, radius in enumerate(data.effective_radius):
        axis.plot(
            data.x_reflectance[:, index], data.y_reflectance[:, index],
            color=radius_color, linestyle=linestyle, linewidth=1.0, alpha=alpha,
        )
        label_location = _label_point(
            data.x_reflectance[:, index],
            data.y_reflectance[:, index],
            radius_label_fractions.get(index, 0.68),
        )
        if annotate and index in radius_labels and label_location is not None:
            axis.annotate(
                f"rₑ={radius:g} µm",
                label_location,
                fontsize=6.5,
                color=radius_color,
                xytext=(2, -7),
                textcoords="offset points",
                bbox=label_box,
            )


def _common_limits(*arrays: np.ndarray) -> tuple[float, float]:
    finite = np.concatenate([np.asarray(values)[np.isfinite(values)] for values in arrays])
    if finite.size == 0:
        raise ComparisonError("selected reflectance grids contain no finite values")
    lower, upper = float(np.min(finite)), float(np.max(finite))
    padding = max((upper - lower) * 0.05, 1.0e-4)
    return lower - padding, upper + padding


def create_figure(lut_a: NakajimaKingData, lut_b: NakajimaKingData):
    """Create the four-panel comparison figure."""

    dx, dy, magnitude = calculate_displacement(lut_a, lut_b)
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 10.0), constrained_layout=True)
    panel_a, panel_b, overlay, difference = axes.flat

    _plot_grid(panel_a, lut_a, linestyle="-", annotate=True)
    _plot_grid(panel_b, lut_b, linestyle="-", annotate=True)
    _plot_grid(overlay, lut_a, linestyle="-", annotate=False, alpha=0.9)
    _plot_grid(overlay, lut_b, linestyle="--", annotate=False, alpha=0.9)

    tau_handle = Line2D([0], [0], color="#2166ac", label="constant optical depth")
    radius_handle = Line2D([0], [0], color="#b2182b", label="constant effective radius")
    panel_a.legend(handles=[tau_handle, radius_handle], fontsize=8)
    panel_b.legend(handles=[tau_handle, radius_handle], fontsize=8)
    overlay.legend(
        handles=[
            Line2D([0], [0], color="0.2", linestyle="-", label="LUT A"),
            Line2D([0], [0], color="0.2", linestyle="--", label="LUT B"),
            tau_handle,
            radius_handle,
        ],
        fontsize=8,
    )

    difference.quiver(
        lut_a.x_reflectance,
        lut_a.y_reflectance,
        dx,
        dy,
        angles="xy",
        scale_units="xy",
        scale=1.0,
        color="0.35",
        alpha=0.65,
        width=0.0025,
    )
    points = difference.scatter(
        lut_a.x_reflectance,
        lut_a.y_reflectance,
        c=magnitude,
        s=9,
        cmap="viridis",
        alpha=0.8,
        label="LUT A grid points",
    )
    difference.legend(fontsize=8)
    figure.colorbar(points, ax=difference, label="|ΔR| = √(ΔRₓ² + ΔRᵧ²)")
    finite_magnitude = magnitude[np.isfinite(magnitude)]
    if finite_magnitude.size:
        difference.text(
            0.02,
            0.98,
            f"median |ΔR|={np.median(finite_magnitude):.3g}\nmax |ΔR|={np.max(finite_magnitude):.3g}",
            transform=difference.transAxes,
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.85},
        )

    x_limits = _common_limits(lut_a.x_reflectance, lut_b.x_reflectance)
    y_limits = _common_limits(lut_a.y_reflectance, lut_b.y_reflectance)
    for axis in axes.flat:
        axis.set_xlim(x_limits)
        axis.set_ylim(y_limits)
        axis.set_xlabel(_channel_label(lut_a.metadata, lut_a.selection.x_channel))
        axis.set_ylabel(_channel_label(lut_a.metadata, lut_a.selection.y_channel))
        axis.grid(True, color="0.88", linewidth=0.6)
        axis.set_axisbelow(True)

    panel_a.set_title(f"LUT A: {lut_a.metadata.path.name}", fontsize=10)
    panel_b.set_title(f"LUT B: {lut_b.metadata.path.name}", fontsize=10)
    overlay.set_title("Overlay: LUT A solid, LUT B dashed", fontsize=10)
    difference.set_title("Pointwise displacement A → B", fontsize=10)
    geometry = (
        f"SZA={lut_a.selection.solar_zenith:g}°, VZA={lut_a.selection.satellite_zenith:g}°, "
        f"RAA={lut_a.selection.relative_azimuth:g}°"
    )
    if lut_a.selection.surface_pressure is not None:
        geometry += f", pressure={lut_a.selection.surface_pressure:g} hPa"
    figure.suptitle(
        f"Nakajima–King LUT comparison — {lut_a.metadata.platform} {lut_a.metadata.instrument}\n{geometry}",
        fontsize=12,
    )
    return figure


def compare_luts(
    lut_a_path: str | Path,
    lut_b_path: str | Path,
    *,
    output: str | Path,
    x_channel: int | None = None,
    y_channel: int | None = None,
    solar_zenith: float | None = None,
    satellite_zenith: float | None = None,
    relative_azimuth: float | None = None,
    surface_pressure: float | None = None,
) -> Path:
    """Read, validate, compare, and plot two LUTs."""

    metadata_a = read_lut_metadata(lut_a_path)
    metadata_b = read_lut_metadata(lut_b_path)
    check_compatible(metadata_a, metadata_b)
    selection = resolve_selection(
        metadata_a,
        x_channel=x_channel,
        y_channel=y_channel,
        solar_zenith=solar_zenith,
        satellite_zenith=satellite_zenith,
        relative_azimuth=relative_azimuth,
        surface_pressure=surface_pressure,
    )
    data_a = extract_nakajima_king(metadata_a, selection)
    data_b = extract_nakajima_king(metadata_b, selection)
    figure = create_figure(data_a, data_b)
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return output_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Overlay two compatible ORAC V2 LUTs as Nakajima–King diagrams."
    )
    parser.add_argument("lut_a", type=Path, help="first ORAC V2 NetCDF LUT")
    parser.add_argument("lut_b", type=Path, help="second ORAC V2 NetCDF LUT")
    parser.add_argument("-o", "--output", type=Path, help="output image (default: validation/nakajima_king/results/...png)")
    parser.add_argument(
        "--x-channel", type=int,
        help="solar channel ID for the x axis (required with --y-channel when more than two exist)",
    )
    parser.add_argument(
        "--y-channel", type=int,
        help="solar channel ID for the y axis (required with --x-channel when more than two exist)",
    )
    parser.add_argument("--solar-zenith", type=float, help="exact solar zenith coordinate in degrees")
    parser.add_argument("--satellite-zenith", type=float, help="exact satellite zenith coordinate in degrees")
    parser.add_argument("--relative-azimuth", type=float, help="exact relative azimuth coordinate in degrees")
    parser.add_argument("--surface-pressure", type=float, help="exact pressure coordinate in hPa, for aerosol LUTs")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    arguments = parser.parse_args(argv)
    output = arguments.output
    if output is None:
        output = Path(__file__).resolve().parent / "results" / (
            f"{arguments.lut_a.stem}_vs_{arguments.lut_b.stem}_nakajima_king.png"
        )
    try:
        result = compare_luts(
            arguments.lut_a,
            arguments.lut_b,
            output=output,
            x_channel=arguments.x_channel,
            y_channel=arguments.y_channel,
            solar_zenith=arguments.solar_zenith,
            satellite_zenith=arguments.satellite_zenith,
            relative_azimuth=arguments.relative_azimuth,
            surface_pressure=arguments.surface_pressure,
        )
    except (ComparisonError, FileNotFoundError, OSError, ValueError) as exc:
        parser.exit(2, f"compare_nakajima_king.py: {exc}\n")
    print(f"Wrote {result}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
