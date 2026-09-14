"""Inventory and safety checks for the existing legacy test LUT definitions."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ..config.legacy import read_instrument, read_lut_grid, read_microphysics, reference_configuration


@dataclass(frozen=True)
class LegacyTestCaseDimensions:
    optical_depth: int
    effective_radius: int
    solar_zenith: int
    satellite_zenith: int
    relative_azimuth: int
    selected_channels: int
    source_channels: int
    components: int
    pressure: int | None


@dataclass(frozen=True)
class ExistingTestCase:
    name: str
    lut_file: str
    spacings: tuple[str, ...]
    optical_depth: tuple[float, ...]
    effective_radius: tuple[float, ...]
    solar_zenith: tuple[float, ...]
    satellite_zenith: tuple[float, ...]
    relative_azimuth: tuple[float, ...]
    dimensions: LegacyTestCaseDimensions
    particle_model: str
    instrument_file: str
    atmosphere_code: int
    srf_mode: int
    rayleigh: bool
    gas: bool
    disort_streams: int


def load_existing_test_case(
    repository_root: str | Path,
    lut_file: str | Path | None = None,
    selected_channels: tuple[int, ...] | None = None,
) -> ExistingTestCase:
    root = Path(repository_root).resolve()
    config = reference_configuration(root)
    path = Path(lut_file) if lut_file is not None else (
        root / "create_orac_lut/input_files/lut/liquid-water-cloud_test.lut"
    )
    path = path.resolve()
    grid = read_lut_grid(path)
    instrument = read_instrument(config.instrument_file)
    model = read_microphysics(config.microphysics_file)
    selected = selected_channels if selected_channels is not None else instrument.available_channels
    unknown = sorted(set(selected) - set(instrument.available_channels))
    if unknown:
        raise ValueError(f"Selected channels are not in the instrument definition: {unknown}")
    return ExistingTestCase(
        name=path.stem,
        lut_file=str(path),
        spacings=grid.spacings,
        optical_depth=tuple(float(x) for x in grid.optical_depth),
        effective_radius=tuple(float(x) for x in grid.effective_radius),
        solar_zenith=tuple(float(x) for x in grid.solar_zenith),
        satellite_zenith=tuple(float(x) for x in grid.satellite_zenith),
        relative_azimuth=tuple(float(x) for x in grid.relative_azimuth),
        dimensions=LegacyTestCaseDimensions(
            optical_depth=grid.optical_depth.size,
            effective_radius=grid.effective_radius.size,
            solar_zenith=grid.solar_zenith.size,
            satellite_zenith=grid.satellite_zenith.size,
            relative_azimuth=grid.relative_azimuth.size,
            selected_channels=len(selected),
            source_channels=len(instrument.available_channels),
            components=1,
            pressure=None,
        ),
        particle_model=config.microphysics_file.name,
        instrument_file=str(config.instrument_file),
        atmosphere_code=config.atmosphere_code,
        srf_mode=config.srf_mode,
        rayleigh=config.rayleigh,
        gas=config.gas,
        disort_streams=config.disort_streams,
    )


def enforce_small_case(case: ExistingTestCase) -> None:
    limits = {
        "optical_depth": 2,
        "effective_radius": 3,
        "solar_zenith": 2,
        "satellite_zenith": 2,
        "relative_azimuth": 2,
        "selected_channels": 1,
        "components": 1,
    }
    dimensions = asdict(case.dimensions)
    violations = [
        f"{name}={dimensions[name]} exceeds guard {limit}"
        for name, limit in limits.items()
        if dimensions[name] > limit
    ]
    if case.dimensions.pressure not in (None, 1):
        violations.append(f"pressure={case.dimensions.pressure} is not supported by this guard")
    if violations:
        raise ValueError("Refusing unexpectedly large legacy test case: " + "; ".join(violations))


def case_record(case: ExistingTestCase) -> dict[str, object]:
    result = asdict(case)
    result["dimensions"] = asdict(case.dimensions)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--lut-file", type=Path, required=True)
    parser.add_argument("--channels", type=int, nargs="+", required=True)
    parser.add_argument("--output-manifest", type=Path)
    args = parser.parse_args()
    case = load_existing_test_case(args.repository_root, args.lut_file, tuple(args.channels))
    enforce_small_case(case)
    record = case_record(case)
    print(json.dumps(record, indent=2))
    if args.output_manifest:
        args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
        args.output_manifest.write_text(json.dumps(record, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
