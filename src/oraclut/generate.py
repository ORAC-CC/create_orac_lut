"""User-facing command line for generating current ORAC V2 LUTs.

This entry point replaces the practical role of the legacy ``makerunfile_v2.pro``
-> ``<platform>_<instrument>_run`` -> ``create_orac_{cloud,aerosol}_lut_wrapper``
chain.  A run is specified either by command-line options or by a readable
driver/configuration file in the legacy driver format (see ``configs/``), and
all inputs are resolved from the normal ``create_orac_lut/input_files``
hierarchy rather than from duplicated Python-only paths.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys

from .config import read_driver, read_instrument, read_lut_grid, read_microphysics, reference_configuration
from .pipeline import generate_aerosol, generate_cloud, read_channels, write_generation
from .provenance import build_record, write_provenance
from .version import LutVersionMismatch, check_lut_version, print_banner


INPUT_ROOT = Path("create_orac_lut") / "input_files"
OUTPUT_ROOT = Path("create_orac_lut") / "luts"
FORWARD_MODELS = ("cloud", "aerosol")

# Legacy wrapper keywords that alter the science and are not implemented here.
UNSUPPORTED_OPTIONS = ("force_n", "force_k", "n_theta", "mie", "srfdat", "tmatrix_path", "opt_prop_luts")
# Legacy workflow-only keywords with no Python counterpart; accepted and ignored.
IGNORED_OPTIONS = ("no_screen", "reuse_scat", "scat_only", "null")
# Every key a configuration file may set after its five positional values.
CONFIGURATION_KEYS = (
    "forward_model", "platform", "instrument", "channelid", "srf_quad", "gas",
    "no_rayleigh", "rayleigh", "streams", "phase_order", "version", "output",
)

# makerunfile_v2.pro derived the formulation and LUT family from the material,
# i.e. the part of the microphysics name before the underscore.
MATERIAL_FORWARD_MODEL = {
    "aerosol": "aerosol", "biomass": "cloud", "liquid-water": "cloud",
    "sulphuric-acid": "cloud", "volcanic-ash": "cloud", "water-ice": "cloud",
}
ATMOSPHERE_NAMES = {
    0: "midlatitude summer, old midsatm.dat profile",
    1: "tropical", 2: "midlatitude summer", 3: "midlatitude winter",
    4: "subarctic summer", 5: "subarctic winter", 6: "1976 US standard",
}
MATERIAL_LUT_FAMILY = {
    "aerosol": "aerosol", "biomass": "biomass-plume", "liquid-water": "liquid-water-cloud",
    "sulphuric-acid": "sulphuric-acid-cloud", "volcanic-ash": "ash-plume", "water-ice": "ice-cloud",
}


@dataclass(frozen=True)
class ResolvedConfiguration:
    root: Path
    forward_model: str
    instrument_path: Path
    lut_path: Path
    microphysics_path: Path
    output_path: Path
    channels: tuple[int, ...]
    atmosphere: int
    srf_quad: int
    rayleigh: bool
    gas: bool
    phase_order: int
    streams: int
    lut_level: int
    revision: int
    instrument: object
    grid: object
    model: object
    config_path: Path | None = None


def _channels(value: str) -> tuple[int, ...]:
    """Parse ``1,2,3``, ``[1, 2, 3]`` (the legacy ``channelid`` form) or ``all``."""

    text = value.strip()
    if text.lower() == "all":
        return ()
    tokens = [token for token in re.split(r"[,\s]+", text.strip("[]() ")) if token]
    bad = [token for token in tokens if not re.fullmatch(r"\d+", token)]
    if bad:
        raise ValueError(f"channels must be integers or 'all'; cannot read {bad} in {value!r}")
    result = tuple(int(token) for token in tokens)
    if not result:
        raise ValueError("channels must contain at least one channel number")
    if len(set(result)) != len(result):
        raise ValueError(f"channels contains repeated entries: {value!r}")
    return result


def _channels_argument(value: str) -> tuple[int, ...]:
    try:
        return _channels(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _under_root(root: Path, value: str | Path, *, label: str) -> Path:
    candidate = Path(value)
    resolved = (candidate if candidate.is_absolute() else root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"{label} must remain inside the repository: {resolved}") from exc
    return resolved


def _input_file(root: Path, value: str, *, subdirectory: str, suffix: str, label: str) -> Path:
    """Resolve a bare legacy input name or an explicit path within the repository.

    ``liquid-water_stg.mm`` resolves to ``input_files/microphysics/liquid-water_stg.mm``
    exactly as the legacy ``mmfile=`` keyword did; an explicit relative or
    absolute path is honoured.  A name that does not exist fails immediately,
    listing nearby candidates, rather than being substituted or deferred to a
    reader further down the pipeline.
    """

    candidate = Path(value)
    if candidate.suffix != suffix:
        raise ValueError(f"{label} must be a {suffix} file name: {value!r}")
    directory = root / INPUT_ROOT / subdirectory
    resolved = (
        _under_root(root, INPUT_ROOT / subdirectory / value, label=label)
        if candidate.name == value
        else _under_root(root, candidate, label=label)
    )
    if not resolved.is_file():
        stem = candidate.stem.split("_", 1)[0].split("-", 1)[0].lower()
        near = sorted(
            path.name for path in directory.glob(f"*{suffix}")
            if not path.name.endswith(".bak") and path.name.lower().startswith(stem[:4])
        )
        hint = f" Similar names in {directory.relative_to(root)}: {', '.join(near[:8])}." if near else ""
        raise FileNotFoundError(f"{label} not found: {resolved}.{hint}")
    return resolved


def resolve_instrument_file(root: Path, platform: str, instrument: str) -> Path:
    """Find ``input_files/inst/<platform>_<instrument>_v*.inst`` as makerunfile_v2 did."""

    directory = root / INPUT_ROOT / "inst"
    pattern = f"{platform.lower()}_{instrument.lower()}_v*.inst"
    matches = sorted(directory.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No instrument definition matches {directory / pattern}")
    if len(matches) > 1:
        names = ", ".join(path.name for path in matches)
        raise ValueError(f"Several instrument definitions match {pattern}; pass --instrument-file ({names})")
    return matches[0]


def material_of(microphysics_name: str) -> str:
    """Return the makerunfile_v2 material key of a microphysics file name."""

    return Path(microphysics_name).name.split("_", 1)[0].lower()


def lut_family(microphysics_name: str) -> str:
    """Return the makerunfile_v2 LUT family stem, e.g. ``liquid-water-cloud``."""

    material = material_of(microphysics_name)
    if material not in MATERIAL_LUT_FAMILY:
        raise ValueError(
            f"No makerunfile_v2 LUT family for material {material!r}; name the LUT grid explicitly"
        )
    return MATERIAL_LUT_FAMILY[material]


def default_forward_model(microphysics_name: str) -> str:
    material = material_of(microphysics_name)
    if material not in MATERIAL_FORWARD_MODEL:
        raise ValueError(
            f"No makerunfile_v2 formulation for material {material!r}; pass --forward-model explicitly"
        )
    return MATERIAL_FORWARD_MODEL[material]


def legacy_lut_filename(
    *, platform: str, instrument: str, srf_quad: int, substance: str, shortname: str,
    atmosphere: int, gas: bool, rayleigh: bool, revision: int,
) -> str:
    """Reproduce the legacy V2 product filename convention.

    ``<platform>_<instrument>_<m|b>_<substance>_a<code>_p<shortname>_v<revision>.nc``
    where the atmospheric-model code is ``1X`` (Rayleigh plus gas for MODTRAN
    atmosphere X), ``00`` (no Rayleigh) or ``01`` (Rayleigh, no gas).
    """

    mono_or_band = "m" if srf_quad == 1 else "b"
    if gas:
        code = f"1{atmosphere}"
    elif not rayleigh:
        code = "00"
    else:
        code = "01"
    return (
        f"{platform.lower()}_{instrument.lower()}_{mono_or_band}_{substance.lower()}"
        f"_a{code}_p{shortname.lower()}_v{revision:02d}.nc"
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m oraclut.generate",
        description=(
            "Generate an ORAC V2 LUT using the legacy-equivalent cloud (discrete particle "
            "layer) or aerosol (particles distributed through the atmospheric column) "
            "forward-model formulation. Inputs are resolved from "
            "create_orac_lut/input_files exactly as the legacy makerunfile_v2 workflow did."
        ),
        epilog=(
            "Settings are combined in the order: built-in defaults < --config file < "
            "explicit command-line options. Bare names such as liquid-water_stg.mm or "
            "aerosol_test.lut are looked up in input_files/microphysics and "
            "input_files/lut; --platform/--instrument select input_files/inst/"
            "<platform>_<instrument>_v*.inst."
        ),
    )
    parser.add_argument("--config", type=Path, help="Driver/configuration file in the legacy driver format")
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--forward-model", choices=FORWARD_MODELS, help="Legacy formulation name (default: derived from the microphysics material as in makerunfile_v2)")
    parser.add_argument("--platform", help="Platform name, e.g. meteosat-10")
    parser.add_argument("--instrument", help="Instrument name, e.g. seviri")
    parser.add_argument("--instrument-file", help="Explicit .inst file; overrides --platform/--instrument lookup")
    parser.add_argument("--lut-file", "--lut", dest="lut_file", help="LUT grid definition (.lut), bare name or path. Required unless --test is given; there is no implicit production grid")
    parser.add_argument("--test", action="store_true", help="Use the compact <family>_test.lut development grid for the material, like Test=1 in makerunfile_v2. Cannot be combined with an explicit LUT grid")
    parser.add_argument("--microphysics-file", "--microphysics", dest="microphysics_file", help="Microphysical model (.mm), bare name or path")
    parser.add_argument("--channels", type=_channels_argument, help="Comma-separated channel ids or 'all' (default: all instrument channels)")
    parser.add_argument("--atmosphere", type=int, choices=tuple(range(7)), help="MODTRAN atmosphere code 0-6 (default: 2)")
    parser.add_argument("--srf-quad", type=int, choices=(1, 2), help="1 = monochromatic at effective centre, 2 = band integration (default: 1)")
    parser.add_argument("--phase-order", type=int, help="Number of phase-function moments (default: 1000)")
    parser.add_argument("--streams", type=int, help="DISORT streams (default: 60)")
    rayleigh = parser.add_mutually_exclusive_group()
    rayleigh.add_argument("--rayleigh", dest="rayleigh", action="store_true", default=None, help="Include Rayleigh scattering (default)")
    rayleigh.add_argument("--no-rayleigh", dest="rayleigh", action="store_false", help="Exclude Rayleigh scattering")
    gas = parser.add_mutually_exclusive_group()
    gas.add_argument("--gas", dest="gas", action="store_true", default=None, help="Include legacy gas optical-depth profiles (aerosol formulation)")
    gas.add_argument("--no-gas", dest="gas", action="store_false", help="Exclude gas absorption (default)")
    parser.add_argument("--lut-level", type=int, default=2, choices=(2,))
    parser.add_argument("--version", "--revision", dest="revision", type=int, help="Output LUT product revision, e.g. 21 (default: 21)")
    parser.add_argument("--output", type=Path, help="Output NetCDF file, or a directory to receive the legacy-named file (default: create_orac_lut/luts/<platform>_<instrument>_<formulation>/)")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Resolve inputs and print the plan without Mie or DISORT")
    return parser


def _truthy(value: str) -> bool:
    text = value.strip().lower()
    if text in ("", "1", "true", "yes", "on"):
        return True
    if text in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"Expected a flag value, got {value!r}")


def read_configuration_file(path: Path) -> dict[str, object]:
    """Read a legacy-format driver file into generator settings.

    The five positional values (input root, instrument file, microphysics file,
    LUT file, atmosphere code) follow the legacy wrapper convention; the
    ``key=value`` lines accept the legacy wrapper keywords (``channelid``,
    ``srf_quad``, ``gas``, ``no_rayleigh``, ``version``) plus the Python-only
    keys ``forward_model``, ``platform``, ``instrument``, ``output``,
    ``streams`` and ``phase_order``.  ``default.mm``/``default.lut``
    placeholders, as used by the legacy per-instrument drivers, leave the value
    to the command line.
    """

    def integer(key: str, value: str) -> int:
        try:
            return int(value.strip())
        except ValueError as exc:
            raise ValueError(f"Configuration key {key!r} expects an integer, got {value.strip()!r}: {path}") from exc

    def flag(key: str, value: str) -> bool:
        try:
            return _truthy(value)
        except ValueError as exc:
            raise ValueError(f"Configuration key {key!r} expects 1/0 (or true/false), got {value.strip()!r}: {path}") from exc

    driver = read_driver(path)
    settings: dict[str, object] = {"atmosphere": driver.atmosphere}
    if driver.input_root.strip() not in ("input_files", str(INPUT_ROOT)):
        raise ValueError(
            f"Configuration input root must be the repository input hierarchy 'input_files', got {driver.input_root!r}: {path}"
        )
    if driver.instrument_file.strip().lower() != "default.inst":
        settings["instrument_file"] = driver.instrument_file.strip()
    if driver.microphysics_file.strip().lower() != "default.mm":
        settings["microphysics_file"] = driver.microphysics_file.strip()
    if driver.lut_file.strip().lower() != "default.lut":
        settings["lut_file"] = driver.lut_file.strip()
    for key, value in driver.options.items():
        if key in ("channelid", "channels"):
            try:
                settings["channels"] = _channels(value)
            except ValueError as exc:
                raise ValueError(f"{exc} in {path}") from exc
        elif key == "srf_quad":
            settings["srf_quad"] = integer(key, value)
        elif key == "gas":
            settings["gas"] = flag(key, value)
        elif key == "no_rayleigh":
            settings["rayleigh"] = not flag(key, value)
        elif key == "rayleigh":
            settings["rayleigh"] = flag(key, value)
        elif key in ("version", "revision"):
            settings["revision"] = integer(key, value)
        elif key == "forward_model":
            if value.strip().lower() not in FORWARD_MODELS:
                raise ValueError(
                    f"Configuration key 'forward_model' must be one of {FORWARD_MODELS}, "
                    f"got {value.strip()!r}: {path}"
                )
            settings["forward_model"] = value.strip().lower()
        elif key in ("platform", "instrument", "output"):
            settings[key] = value.strip()
        elif key in ("streams", "phase_order"):
            settings[key] = integer(key, value)
        elif key in UNSUPPORTED_OPTIONS:
            raise NotImplementedError(f"Legacy option {key!r} is not supported by the Python generator: {path}")
        elif key in IGNORED_OPTIONS:
            print(f"note: legacy workflow option {key!r} has no effect in the Python generator", file=sys.stderr)
        else:
            raise ValueError(
                f"Unrecognised configuration key {key!r}: {path}. Recognised keys: "
                f"{', '.join(CONFIGURATION_KEYS)}."
            )
    if "forward_model" not in settings:
        stem = path.stem.lower()
        for model in FORWARD_MODELS:
            if stem.endswith(f"_{model}") or f"_{model}_" in stem or f"_{model}." in path.name.lower():
                settings["forward_model"] = model
                break
    return settings


def _merge(args: argparse.Namespace, settings: dict[str, object]) -> dict[str, object]:
    """Apply defaults < configuration file < explicit command line."""

    merged: dict[str, object] = {
        "forward_model": None, "channels": (), "atmosphere": 2, "srf_quad": 1,
        "phase_order": 1000, "streams": 60, "rayleigh": True, "gas": False, "revision": 21,
        "platform": None, "instrument": None, "instrument_file": None,
        "lut_file": None, "microphysics_file": None, "output": None,
    }
    merged.update(settings)
    for key in merged:
        value = getattr(args, key, None)
        if value is not None:
            merged[key] = value
    return merged


def resolve_configuration(args: argparse.Namespace) -> ResolvedConfiguration:
    root = args.repository_root.resolve()
    config_path = None
    settings: dict[str, object] = {}
    if args.config is not None:
        config_path = _under_root(root, args.config, label="configuration file")
        settings = read_configuration_file(config_path)
    values = _merge(args, settings)
    if not values["microphysics_file"]:
        raise ValueError(
            "No microphysical model selected. Set the third configuration value, or pass "
            "--microphysics <name>.mm (a file in create_orac_lut/input_files/microphysics)."
        )
    if values["lut_file"] and args.test:
        raise ValueError(
            f"--test cannot be combined with an explicit LUT grid ({values['lut_file']!r}); "
            "the configuration already names the grid to use."
        )
    if not values["lut_file"]:
        # The LUT grid is never chosen implicitly: a missing grid must not turn
        # into a production-sized calculation by default.
        family = lut_family(str(values["microphysics_file"]))
        if not args.test:
            raise ValueError(
                "No LUT grid selected, and there is no default: name the grid explicitly. "
                f"For {values['microphysics_file']} the makerunfile_v2 family is {family!r}, so use "
                f"--lut {family}.lut for the full production grid or --test for the compact "
                f"{family}_test.lut development grid (in a configuration file, set the fourth value)."
            )
        values["lut_file"] = f"{family}_test.lut"
    if not values["forward_model"]:
        values["forward_model"] = default_forward_model(str(values["microphysics_file"]))
    forward_model = str(values["forward_model"])

    if values["instrument_file"]:
        instrument_path = _input_file(root, str(values["instrument_file"]), subdirectory="inst", suffix=".inst", label="instrument file")
    elif values["platform"] and values["instrument"]:
        instrument_path = resolve_instrument_file(root, str(values["platform"]), str(values["instrument"]))
    else:
        raise ValueError("Specify --platform and --instrument (or --instrument-file / a configuration file)")
    lut_path = _input_file(root, str(values["lut_file"]), subdirectory="lut", suffix=".lut", label="LUT file")
    microphysics_path = _input_file(root, str(values["microphysics_file"]), subdirectory="microphysics", suffix=".mm", label="microphysics file")

    instrument = read_instrument(instrument_path)
    if values["platform"] and str(values["platform"]).lower() != instrument.platform.lower():
        raise ValueError(f"platform {values['platform']!r} does not match {instrument.platform!r} in {instrument_path}")
    if values["instrument"] and str(values["instrument"]).lower() != instrument.instrument.lower():
        raise ValueError(f"instrument {values['instrument']!r} does not match {instrument.instrument!r} in {instrument_path}")
    if instrument.view > 0:
        raise NotImplementedError("Dual-view instruments are not yet implemented in the Python V2 generator")
    grid = read_lut_grid(lut_path)
    model = read_microphysics(microphysics_path)
    channels = tuple(values["channels"]) or instrument.available_channels
    unknown = sorted(set(channels) - set(instrument.available_channels))
    if unknown:
        raise ValueError(
            f"Channels {unknown} are not available on {instrument.platform}/{instrument.instrument}; "
            f"{instrument_path.name} defines {list(instrument.available_channels)}"
        )
    atmosphere = int(values["atmosphere"])
    srf_quad = int(values["srf_quad"])
    rayleigh = bool(values["rayleigh"])
    gas = bool(values["gas"])
    phase_order = int(values["phase_order"])
    streams = int(values["streams"])
    revision = int(values["revision"])
    if atmosphere not in range(7):
        raise ValueError(f"atmosphere must be a MODTRAN code 0-6, got {atmosphere}")
    if srf_quad not in (1, 2):
        raise ValueError(f"srf_quad must be 1 (monochromatic) or 2 (band integration), got {srf_quad}")
    if forward_model == "aerosol" and grid.surface_pressure is None:
        raise ValueError(
            f"The aerosol formulation requires a LUT grid with a sixth surface-pressure block, but "
            f"{lut_path.name} defines only five grids. Use an aerosol grid such as aerosol.lut "
            f"or aerosol_test.lut, or set forward_model=cloud."
        )
    if forward_model == "cloud" and grid.surface_pressure is not None:
        raise ValueError(
            f"The cloud formulation expects a five-grid LUT definition, but {lut_path.name} adds a "
            f"surface-pressure grid. Use a cloud grid, or set forward_model=aerosol."
        )
    if forward_model == "cloud" and gas:
        raise ValueError(
            "gas absorption is currently implemented for the aerosol formulation only; "
            "remove gas=1 or use forward_model=aerosol"
        )
    if streams < 2 or streams % 2:
        raise ValueError(f"DISORT streams must be a positive even number, got {streams}")
    if phase_order < 2 * streams:
        raise ValueError(
            f"phase order ({phase_order}) must be at least twice the DISORT stream count ({streams})"
        )

    filename = legacy_lut_filename(
        platform=instrument.platform, instrument=instrument.instrument, srf_quad=srf_quad,
        substance=model.substance, shortname=model.shortname, atmosphere=atmosphere,
        gas=gas, rayleigh=rayleigh, revision=revision,
    )
    if values["output"]:
        output_path = _under_root(root, str(values["output"]), label="output")
        if output_path.is_dir() or output_path.suffix != ".nc":
            output_path = output_path / filename
    else:
        directory = OUTPUT_ROOT / f"{instrument.platform.lower()}_{instrument.instrument.lower()}_{forward_model}"
        output_path = _under_root(root, directory / filename, label="output")

    # This validates SRF files without invoking either scientific kernel.
    read_channels(root, instrument, channels, srf_quad=srf_quad)
    return ResolvedConfiguration(
        root, forward_model, instrument_path, lut_path, microphysics_path, output_path,
        channels, atmosphere, srf_quad, rayleigh, gas, phase_order, streams,
        int(args.lut_level), revision, instrument, grid, model, config_path,
    )


def grid_dimensions(config: ResolvedConfiguration) -> dict[str, int]:
    """Return every LUT dimension of the resolved configuration, in LUT order."""

    grid = config.grid
    dimensions = {
        "channels": len(config.channels),
        "optical_depth": len(grid.optical_depth),
        "effective_radius": len(grid.effective_radius),
        "solar_zenith": len(grid.solar_zenith),
        "satellite_zenith": len(grid.satellite_zenith),
        "relative_azimuth": len(grid.relative_azimuth),
    }
    if grid.surface_pressure is not None:
        dimensions["surface_pressure"] = len(grid.surface_pressure)
    return dimensions


def radiative_transfer_states(config: ResolvedConfiguration) -> int:
    """Estimate the number of DISORT states before angular/direct calls."""

    grid = config.grid
    states = len(config.channels) * len(grid.effective_radius) * len(grid.optical_depth)
    if grid.surface_pressure is not None:
        states *= len(grid.surface_pressure)
    solar = sum(channel in config.instrument.solar_channels for channel in config.channels)
    return states * (1 + solar * len(grid.solar_zenith))


def _print_summary(config: ResolvedConfiguration) -> None:
    grid = config.grid
    instrument = config.instrument
    dimensions = grid_dimensions(config)
    grid_points = 1
    for name, size in dimensions.items():
        if name != "channels":
            grid_points *= size
    states = radiative_transfer_states(config)

    def relative(path: Path) -> Path:
        return path.relative_to(config.root) if path.is_relative_to(config.root) else path

    solar = [channel for channel in config.channels if channel in instrument.solar_channels]
    thermal = [channel for channel in config.channels if channel in instrument.thermal_channels]
    formulation = (
        "discrete particle layer" if config.forward_model == "cloud"
        else "particles distributed through the atmospheric column"
    )
    srf_treatment = (
        "monochromatic at the effective channel centre" if config.srf_quad == 1
        else "integrated over the spectral response function"
    )
    scattering = ", ".join(sorted({component.scattering_code for component in config.model.components}))
    grid_kind = "compact development/test grid" if config.lut_path.stem.endswith("_test") else "full grid"

    print("Resolved ORAC LUT configuration")
    if config.config_path is not None:
        print(f"  configuration file:  {relative(config.config_path)}")
    print(f"  forward model:       {config.forward_model} ({formulation})")
    print(f"  platform/instrument: {instrument.platform} / {instrument.instrument} (instrument version {instrument.version})")
    print(f"  instrument file:     {relative(config.instrument_path)}")
    print(f"  microphysics:        {relative(config.microphysics_path)}")
    print(f"                       {config.model.substance} / {config.model.shortname}, "
          f"{len(config.model.components)} component(s), scattering code {scattering or 'unknown'}")
    print(f"  LUT grid definition: {relative(config.lut_path)}  [{grid_kind}]")
    print(f"  channels:            {len(config.channels)} of {len(instrument.available_channels)} available: "
          f"{', '.join(str(channel) for channel in config.channels)}")
    print(f"                       solar {solar or 'none'}; thermal {thermal or 'none'}")
    print(f"  atmosphere:          {config.atmosphere} ({ATMOSPHERE_NAMES.get(config.atmosphere, 'unknown')})")
    print(f"  Rayleigh scattering: {'on' if config.rayleigh else 'off'}; "
          f"gas absorption: {'on' if config.gas else 'off'}")
    print(f"  SRF treatment:       srf_quad={config.srf_quad} ({srf_treatment})")
    print(f"  DISORT streams:      {config.streams}; phase-function order: {config.phase_order}")
    print(f"  product:             LUT level {config.lut_level}, revision {config.revision}")
    print("  LUT dimensions:")
    for name, size in dimensions.items():
        print(f"      {name:<20s} {size}")
    if "surface_pressure" not in dimensions:
        print(f"      {'surface_pressure':<20s} - (cloud formulation: no pressure dimension)")
    print(f"      {'grid points/channel':<20s} {grid_points:,}")
    print(f"      {'grid points total':<20s} {grid_points * len(config.channels):,}")
    print(f"  estimated workload:  {states:,} DISORT states before angular/direct calls")
    print(f"  output:              {relative(config.output_path)}")
    if config.output_path.exists():
        size = config.output_path.stat().st_size
        print(f"  output status:       ALREADY EXISTS ({size:,} bytes) - a real run needs --overwrite")
    else:
        print("  output status:       does not exist yet")


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        config = resolve_configuration(args)
        # Identify the source (CODE_VERSION + Git) before anything else.  A
        # dry run only reports; a real run must request a LUT version this
        # source release produces (oraclut.version.check_lut_version).
        summary = print_banner(lut_version=config.revision, generator="oraclut.generate (configuration file)")
        _print_summary(config)
        if not args.dry_run:
            check_lut_version(config.revision, summary.code,
                              source=f"configuration {config.config_path or 'command line'}")
        if args.dry_run:
            if config.output_path.exists():
                sys.stdout.flush()
                print(
                    "warning: the output file already exists; a real run would refuse it "
                    "unless --overwrite is given",
                    file=sys.stderr,
                )
            print("dry-run: no Mie or DISORT calculation was executed")
            return 0
        if config.output_path.exists() and not args.overwrite:
            raise FileExistsError(
                f"Refusing to overwrite existing LUT {config.output_path}; pass --overwrite explicitly"
            )
        if config.forward_model == "cloud":
            result = generate_cloud(
                config.root, lut_path=config.lut_path, microphysics_path=config.microphysics_path,
                instrument_path=config.instrument_path, channels=config.channels,
                atmosphere_code=config.atmosphere, srf_quad=config.srf_quad,
                rayleigh=config.rayleigh, phase_order=config.phase_order, nstreams=config.streams,
            )
        else:
            result = generate_aerosol(
                config.root, lut_path=config.lut_path, microphysics_path=config.microphysics_path,
                instrument_path=config.instrument_path, channels=config.channels,
                atmosphere_code=config.atmosphere, srf_quad=config.srf_quad,
                rayleigh=config.rayleigh, gas=config.gas, phase_order=config.phase_order,
                nstreams=config.streams,
            )
        write_generation(
            config.output_path, result, lut_level=config.lut_level, revision=config.revision,
            overwrite=args.overwrite,
        )
        print(f"wrote {config.output_path}")
        print(f"dimensions: {result.dimensions}")
        _write_generation_provenance(config, result)
        return 0
    except (FileExistsError, FileNotFoundError, NotImplementedError, ValueError, LutVersionMismatch) as exc:
        parser.error(str(exc))
    return 2


def _write_generation_provenance(config: ResolvedConfiguration, result) -> None:
    """Provenance sidecar for a product of this configuration-driven path (oraclut.provenance)."""

    grid = config.grid
    axes = {}
    for name in ("optical_depth", "effective_radius", "satellite_zenith", "solar_zenith", "relative_azimuth",
                 "surface_pressure"):
        values = getattr(grid, name, None)
        if values is not None:
            axes[name] = {"n": len(values), "values": [float(v) for v in values]}
    try:
        record = build_record(
            config.output_path, generator="oraclut.generate (development pipeline, oraclut.pipeline)",
            lut_version=config.revision,
            configuration={"forward_model": config.forward_model, "channels": list(config.channels),
                           "atmosphere_code": config.atmosphere, "srf_quad": config.srf_quad,
                           "rayleigh": config.rayleigh, "gas": config.gas, "lut_level": config.lut_level,
                           "instrument": str(config.instrument_path.name), "microphysics": str(config.microphysics_path.name)},
            configuration_file=config.config_path,
            input_files={"instrument": config.instrument_path, "microphysics": config.microphysics_path,
                         "lut_definition": config.lut_path},
            grid={"definition_file": config.lut_path.name, **axes, "dimensions": dict(result.dimensions)},
            numerics={"phase_order": config.phase_order, "srf_quadrature": config.srf_quad,
                      "note": "legacy-reproduction pipeline (oraclut.optics, fixed Legendre order)"},
            radiative_transfer={"solver": "DISORT (oraclut.radiative_transfer.legacy_disort)",
                                "nstreams": config.streams, "rayleigh": config.rayleigh, "gas": config.gas,
                                "atmosphere_code": config.atmosphere})
        print(f"provenance record: {write_provenance(config.output_path, record)}")
    except Exception as exc:   # the product is already written; report, do not hide it
        print(f"WARNING: provenance record for {config.output_path} could not be written: {exc}", file=sys.stderr)


def prepare_reference(repository_root: str | Path) -> None:
    """Validate that the original cloud reference inputs can be read."""

    config = reference_configuration(repository_root)
    read_instrument(config.instrument_file)


if __name__ == "__main__":
    raise SystemExit(main())
