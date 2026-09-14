"""Readers for the existing text inputs used by the STG reference case.

These readers retain the legacy files as the source of truth. They intentionally
do not translate them into a new on-disk format.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np


@dataclass(frozen=True)
class DriverConfig:
    input_root: str
    instrument_file: str
    microphysics_file: str
    lut_file: str
    atmosphere: int
    options: dict[str, str]


@dataclass(frozen=True)
class InstrumentConfig:
    platform: str
    instrument: str
    version: str
    view: int
    available_channels: tuple[int, ...]
    solar_channels: tuple[int, ...]
    thermal_channels: tuple[int, ...]
    maximum_satellite_zenith: float
    srf_files: dict[int, str]
    oldf0: dict[int, float] = field(default_factory=dict)
    oldf1: dict[int, float] = field(default_factory=dict)
    oldnefr: dict[int, float] = field(default_factory=dict)
    oldwvn: dict[int, float] = field(default_factory=dict)
    oldb1: dict[int, float] = field(default_factory=dict)
    oldb2: dict[int, float] = field(default_factory=dict)
    oldt1: dict[int, float] = field(default_factory=dict)
    oldt2: dict[int, float] = field(default_factory=dict)
    oldnebt: dict[int, float] = field(default_factory=dict)
    snr: dict[int, float] = field(default_factory=dict)
    refbt: dict[int, float] = field(default_factory=dict)
    nedt: dict[int, float] = field(default_factory=dict)


@dataclass(frozen=True)
class LutGrid:
    optical_depth: np.ndarray
    effective_radius: np.ndarray
    solar_zenith: np.ndarray
    satellite_zenith: np.ndarray
    relative_azimuth: np.ndarray
    spacings: tuple[str, str, str, str, str]
    surface_pressure: np.ndarray | None = None
    surface_pressure_spacing: str | None = None


@dataclass(frozen=True)
class MicrophysicalComponent:
    """One component from a legacy ``.mm`` definition."""

    component: str
    size_distribution: str
    size_parameters: tuple[float, ...]
    scattering_code: str
    mixing_ratio: float
    refractive_index_file: str


@dataclass(frozen=True)
class MicrophysicalModel:
    substance: str
    shortname: str
    description: str
    profile_height_km: np.ndarray
    profile_relative_amount: np.ndarray
    component: str
    size_distribution: str
    size_parameters: tuple[float, ...]
    scattering_code: str
    mixing_ratio: float
    refractive_index_file: str
    components: tuple[MicrophysicalComponent, ...] = ()


@dataclass(frozen=True)
class AtmosphereProfile:
    height_km: np.ndarray
    pressure_hpa: np.ndarray
    temperature_k: np.ndarray


@dataclass(frozen=True)
class SpectralTable:
    coordinate: np.ndarray
    value: np.ndarray
    coordinate_units: str
    value_units: str


@dataclass(frozen=True)
class ReferenceConfiguration:
    """Effective current configuration for the initial V2 reference case."""

    platform: str
    instrument: str
    particle_model: str
    lut_level: int
    revision: int
    forward_model: str
    driver: Path
    input_root: Path
    instrument_file: Path
    lut_file: Path
    microphysics_file: Path
    atmosphere_file: Path
    solar_spectrum_file: Path
    srf_mode: int
    channels: tuple[int, ...]
    atmosphere_code: int
    gas: bool
    rayleigh: bool
    disort_streams: int


def reference_configuration(repository_root: str | Path) -> ReferenceConfiguration:
    root = Path(repository_root).resolve()
    input_root = root / "create_orac_lut" / "input_files"
    driver = input_root / "driver" / "meteosat-10_seviri_cloud.driver"
    driver_values = read_driver(driver)
    instrument = read_instrument(input_root / "inst" / driver_values.instrument_file)
    return ReferenceConfiguration(
        platform=instrument.platform,
        instrument=instrument.instrument,
        particle_model="liquid-water_stg.mm",
        lut_level=2,
        revision=21,
        forward_model="cloud",
        driver=driver,
        input_root=input_root,
        instrument_file=input_root / "inst" / instrument_path_name(driver_values.instrument_file),
        lut_file=input_root / "lut" / "liquid-water-cloud.lut",
        microphysics_file=input_root / "microphysics" / "liquid-water_stg.mm",
        atmosphere_file=input_root / "atm" / "mls.atm",
        solar_spectrum_file=input_root / "sun" / "Gueymard2018.sssi",
        srf_mode=1,
        channels=instrument.available_channels,
        atmosphere_code=driver_values.atmosphere,
        gas=False,
        rayleigh=True,
        disort_streams=60,
    )


def instrument_path_name(name: str) -> str:
    """Return a driver instrument filename without accepting path traversal."""

    candidate = Path(name)
    if candidate.name != name or candidate.suffix != ".inst":
        raise ValueError(f"Invalid instrument filename: {name!r}")
    return name


def _logical_lines(path: str | Path) -> list[str]:
    lines = []
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            lines.append(line)
    return lines


def read_driver(path: str | Path) -> DriverConfig:
    values: list[str] = []
    options: dict[str, str] = {}
    for line in _logical_lines(path):
        if "=" in line and not line.split("=", 1)[0].strip().replace(".", "").isdigit():
            key, value = line.split("=", 1)
            options[key.strip().lower()] = value.strip()
        else:
            values.append(line)
    if len(values) < 5:
        raise ValueError(f"Driver requires five logical values: {path}")
    try:
        atmosphere = int(values[4].split()[0])
    except ValueError as exc:
        raise ValueError(f"Invalid atmosphere selector in driver: {values[4]!r}") from exc
    return DriverConfig(*values[:4], atmosphere, options)


def _numbers(value: str) -> tuple[int, ...]:
    return tuple(int(x) for x in re.findall(r"\d+", value))


def read_instrument(path: str | Path) -> InstrumentConfig:
    scalar: dict[str, str] = {}
    available: tuple[int, ...] = ()
    solar: tuple[int, ...] = ()
    thermal: tuple[int, ...] = ()
    srfs: dict[int, str] = {}
    per_channel: dict[str, dict[int, float]] = {
        key: {} for key in (
            "oldf0", "oldf1", "oldnefr", "oldwvn", "oldb1", "oldb2",
            "oldt1", "oldt2", "oldnebt", "snr", "refbt", "nedt",
        )
    }
    for line in _logical_lines(path):
        if "=" not in line:
            continue
        key, value = (x.strip() for x in line.split("=", 1))
        key_norm = re.sub(r"\s+", "", key.lower())
        channel_match = re.search(r"\[(\d+)\]", key)
        base_key = key_norm.split("[", 1)[0]
        if base_key == "srf" and channel_match:
            srfs[int(channel_match.group(1))] = value
        elif channel_match and base_key in per_channel:
            per_channel[base_key][int(channel_match.group(1))] = float(value)
        elif key_norm == "availablechannels":
            available = _numbers(value)
        elif key_norm == "solarchannels":
            solar = _numbers(value)
        elif key_norm == "thermalchannels":
            thermal = _numbers(value)
        else:
            scalar[key_norm] = value
    required = ("platform", "instrument", "instrumentversion", "view", "maximumsatellitezenith")
    missing = [key for key in required if key not in scalar]
    if missing:
        raise ValueError(f"Instrument missing required keys {missing}: {path}")
    return InstrumentConfig(
        scalar["platform"], scalar["instrument"], scalar["instrumentversion"],
        int(scalar["view"]), available, solar, thermal,
        float(scalar["maximumsatellitezenith"]), srfs,
        *(per_channel[key] for key in per_channel),
    )


_GRID_SPACINGS = ("linear", "uneven_linear", "logarithmic", "uneven_logarithmic")


def _grid_values(first: str, second: str) -> np.ndarray:
    """Expand one LUT grid exactly as legacy ``load_lutstr.pro`` does.

    ``linear``/``logarithmic`` take two endpoints and generate ``count`` values
    evenly spaced in linear or log10 space; ``uneven_*`` list every value.
    """

    header = first.split()
    count = int(header[0])
    spacing = header[1].lower() if len(header) > 1 else "linear"
    if spacing not in _GRID_SPACINGS:
        raise ValueError(f"Invalid spacing descriptor {spacing!r}; must be one of {_GRID_SPACINGS}")
    values = tuple(float(x) for x in second.split())
    if spacing.startswith("uneven"):
        if len(values) != count:
            raise ValueError(f"Grid contains {len(values)} values, expected {count}")
        return np.asarray(values, dtype=np.float32)
    if len(values) != 2:
        raise ValueError(f"{spacing} grid endpoint line requires two values")
    if spacing == "logarithmic":
        exponents = np.linspace(np.log10(values[0]), np.log10(values[1]), count, dtype=np.float32)
        return (10.0 ** exponents).astype(np.float32)
    return np.linspace(values[0], values[1], count, dtype=np.float32)


def read_lut_grid(path: str | Path) -> LutGrid:
    lines = _logical_lines(path)
    if len(lines) < 10:
        raise ValueError(f"LUT definition is incomplete: {path}")
    grids = []
    spacings = []
    index = 0
    while index + 1 < len(lines) and len(grids) < 5:
        header = lines[index].split()
        spacings.append(header[1].lower() if len(header) > 1 else "unknown")
        grids.append(_grid_values(lines[index], lines[index + 1]))
        index += 2
    if len(grids) != 5:
        raise ValueError(f"LUT definition does not contain five grids: {path}")
    pressure = None
    pressure_spacing = None
    if index + 1 < len(lines):
        header = lines[index].split()
        if len(header) < 2:
            raise ValueError(f"Invalid sixth LUT grid header: {lines[index]!r}")
        pressure = _grid_values(lines[index], lines[index + 1])
        pressure_spacing = header[1].lower()
    return LutGrid(*grids, tuple(spacings), pressure, pressure_spacing)


def read_microphysics(path: str | Path) -> MicrophysicalModel:
    lines = _logical_lines(path)
    scalar: dict[str, str] = {}
    profile_height: list[float] = []
    profile_amount: list[float] = []
    components: list[MicrophysicalComponent] = []
    current: dict[str, object] | None = None
    i = 0
    while i < len(lines):
        line = lines[i]
        lower = line.lower()
        if lower.startswith("substance "):
            scalar["substance"] = line.split(None, 1)[1]
        elif lower.startswith("shortname "):
            scalar["shortname"] = line.split(None, 1)[1]
        elif lower.startswith("description "):
            scalar["description"] = line.split(None, 1)[1]
        elif lower == "profile":
            count = int(lines[i + 1]); i += 2
            for row in lines[i:i + count]:
                h, amount = row.split()[:2]
                profile_height.append(float(h))
                profile_amount.append(float(amount))
            i += count - 1
        elif lower.startswith("component "):
            if current is not None:
                components.append(_component_from_values(current, path))
            fields = line.split(None, 2)
            if len(fields) < 3:
                raise ValueError(f"Component requires a type and name: {line!r}")
            current = {"component": fields[2]}
        elif lower == "*size":
            _require_component(current, path)
            current["size_distribution"] = lines[i + 1].strip().lower()
            current["size_parameters"] = tuple(float(x) for x in lines[i + 2].split())
            i += 2
        elif lower == "*scattering code":
            _require_component(current, path)
            current["scattering_code"] = lines[i + 1].lower()
            i += 1
        elif lower == "*mixing ratio":
            _require_component(current, path)
            current["mixing_ratio"] = float(lines[i + 1])
            i += 1
        elif lower == "*refractive index":
            _require_component(current, path)
            current["refractive_index_file"] = lines[i + 1]
            i += 1
        i += 1
    if current is not None:
        components.append(_component_from_values(current, path))
    missing = [
        key for key, value in {
            "substance": scalar.get("substance"),
            "shortname": scalar.get("shortname"),
        }.items() if value is None
    ]
    if missing or not components:
        raise ValueError(f"Microphysical model missing required fields {missing or ['component']}: {path}")
    first = components[0]
    return MicrophysicalModel(
        scalar["substance"], scalar["shortname"], scalar.get("description", ""),
        np.asarray(profile_height, dtype=np.float32), np.asarray(profile_amount, dtype=np.float32),
        first.component, first.size_distribution, first.size_parameters,
        first.scattering_code, first.mixing_ratio, first.refractive_index_file,
        tuple(components),
    )


def _require_component(current: dict[str, object] | None, path: str | Path) -> None:
    if current is None:
        raise ValueError(f"Microphysical field appears before a component: {path}")


def _component_from_values(values: dict[str, object], path: str | Path) -> MicrophysicalComponent:
    required = (
        "component", "size_distribution", "size_parameters", "scattering_code",
        "mixing_ratio", "refractive_index_file",
    )
    missing = [name for name in required if name not in values]
    if missing:
        raise ValueError(f"Microphysical component missing {missing}: {path}")
    return MicrophysicalComponent(
        str(values["component"]), str(values["size_distribution"]),
        tuple(values["size_parameters"]), str(values["scattering_code"]),
        float(values["mixing_ratio"]), str(values["refractive_index_file"]),
    )


def read_midsatm(path: str | Path) -> AtmosphereProfile:
    rows = []
    for line in _logical_lines(path):
        parts = line.split()
        if len(parts) >= 3:
            rows.append(tuple(float(x) for x in parts[:3]))
    if not rows:
        raise ValueError(f"Atmosphere contains no three-column data: {path}")
    values = np.asarray(rows, dtype=np.float32)
    return AtmosphereProfile(values[:, 0], values[:, 1], values[:, 2])


def read_atmosphere(path: str | Path, atmosphere_code: int) -> AtmosphereProfile:
    """Read the current IDL atmosphere formats for the selected code."""

    if atmosphere_code == 0:
        return read_midsatm(path)
    text = Path(path).read_text()
    sections: dict[str, list[float]] = {}
    current = None
    for line in text.splitlines():
        marker = re.match(r"\s*\*(HGT|PRE|TEM)\b", line, re.IGNORECASE)
        if marker:
            current = marker.group(1).upper()
            sections[current] = []
            continue
        if re.match(r"\s*\*", line):
            current = None
            continue
        if current:
            sections[current].extend(float(token) for token in re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", line))
    if set(sections) != {"HGT", "PRE", "TEM"}:
        raise ValueError(f"MODTRAN/RFM atmosphere lacks HGT, PRE, TEM sections: {path}")
    lengths = {key: len(value) for key, value in sections.items()}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Atmosphere section lengths differ: {lengths}")
    values = np.asarray([sections[key] for key in ("HGT", "PRE", "TEM")], dtype=np.float32)
    values = values[:, values[0] <= 100.0]
    return AtmosphereProfile(values[0, ::-1], values[1, ::-1], values[2, ::-1])


def read_refractive_index(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    for line in _logical_lines(path):
        parts = line.split()
        if len(parts) >= 3:
            rows.append(tuple(float(x) for x in parts[:3]))
    if not rows:
        raise ValueError(f"Refractive-index file contains no data: {path}")
    values = np.asarray(rows, dtype=np.float64)
    return values[:, 0], values[:, 1], values[:, 2]


def read_srf(path: str | Path) -> SpectralTable:
    """Read an ORAC SRF text file as wavenumber and response arrays."""

    rows = []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0].rstrip(",")), float(parts[1].rstrip(","))))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"SRF file contains no numeric pairs: {path}")
    values = np.asarray(rows, dtype=np.float32)
    return SpectralTable(values[:, 0], values[:, 1], "cm^-1", "dimensionless")


def read_solar_spectrum(path: str | Path) -> SpectralTable:
    """Read the two-column Gueymard-style solar spectrum text format."""

    rows = []
    for line in _logical_lines(path):
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"Solar spectrum contains no numeric pairs: {path}")
    values = np.asarray(rows, dtype=np.float32)
    return SpectralTable(values[:, 0], values[:, 1], "microns", "mW cm^-2 micron^-1")
