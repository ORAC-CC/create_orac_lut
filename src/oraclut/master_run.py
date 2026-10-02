"""Maker-style master run expansion and read-only preflight.

The master format deliberately stays close to the legacy makerunfile: one
``name = value`` assignment per line, with ``mmfiles`` selecting one or more
microphysical definitions.  Expansion resolves one complete single-LUT run
per model; it never performs Mie, T-matrix, DISORT, or NetCDF work.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_ROOT = ROOT / "create_orac_lut" / "input_files"

REQUIRED = (
    "platform", "instrument", "forward_model", "in_path", "instfile", "mmfiles", "lutfile",
    "atmospheres", "channelid", "srf_quad", "nstreams", "version",
    "out_path",
)
# nmom is deprecated (the IDL's fixed Legendre expansion); it is accepted
# only for Baum / T-matrix classes, see create_orac_luts.py.
OPTIONAL = {
    "gas": 0, "no_rayleigh": 0, "reuse_scat": 0, "scat_only": 0,
    "tmatrix_path": None, "nmom": None,
}
LUT_BY_MATERIAL = {
    "aerosol": "aerosol.lut",
    "biomass": "biomass-cloud.lut",
    "liquid-water": "liquid-water-cloud.lut",
    "sulphuric-acid": "sulphuric-acid-cloud.lut",
    "volcanic-ash": "ash-plume.lut",
    "water-ice": "ice-cloud.lut",
}
SUPPORTED_MATERIALS = set(LUT_BY_MATERIAL)
SRF_QUADRATURE_NAMES = {
    1: "centre-wavelength / monochromatic treatment",
    2: "segmented SRF integration",
}


def read_master(path: str | Path) -> dict:
    """Read a master ``name = literal`` specification."""

    path = Path(path)
    values = {}
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split(";", 1)[0].split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{line_number}: expected 'name = value'")
        name, literal = (part.strip() for part in line.split("=", 1))
        name = name.lower()
        if name not in REQUIRED and name not in OPTIONAL:
            raise ValueError(f"{path}:{line_number}: unknown master setting {name!r}")
        if name in values:
            raise ValueError(f"{path}:{line_number}: {name!r} is set twice")
        try:
            values[name] = ast.literal_eval(literal)
        except (SyntaxError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number}: invalid literal for {name!r}") from exc
    missing = [name for name in REQUIRED if name not in values]
    if missing:
        raise ValueError(f"{path}: required settings missing: {', '.join(missing)}")
    for name, default in OPTIONAL.items():
        values.setdefault(name, default)
    if not isinstance(values["mmfiles"], (list, tuple)) or not values["mmfiles"]:
        raise ValueError("mmfiles must be a non-empty explicit list")
    values["mmfiles"] = [str(item) if str(item).endswith(".mm") else f"{item}.mm" for item in values["mmfiles"]]
    channels = values["channelid"]
    if not isinstance(channels, (list, tuple)) or not channels:
        raise ValueError("channelid must be a non-empty explicit list")
    values["channelid"] = [int(channel) for channel in channels]
    if len(set(values["channelid"])) != len(values["channelid"]):
        raise ValueError("channelid contains duplicate channel IDs")
    if values["lutfile"] != "auto":
        if len(values["mmfiles"]) != 1:
            raise ValueError(
                "an explicit master lutfile is allowed only when exactly one mmfile is requested; "
                "use lutfile = 'auto' for multiple materials"
            )
        if not isinstance(values["lutfile"], str) or not values["lutfile"].endswith(".lut"):
            raise ValueError("explicit master lutfile must be a .lut filename")
    if int(values["version"]) != 22:
        raise ValueError("Python master specifications must use version = 22")
    if int(values["srf_quad"]) not in (1, 2):
        raise ValueError("srf_quad must be 1 or 2")
    return values


def _rooted(path_value: str | Path) -> Path:
    path = Path(path_value)
    return path if path.is_absolute() else ROOT / path


def _mm_metadata(path: Path) -> tuple[str, str, list[str]]:
    text = path.read_text()
    substance = re.search(r"^substance\s+(.+)$", text, re.MULTILINE | re.IGNORECASE)
    shortname = re.search(r"^shortname\s+(.+)$", text, re.MULTILINE | re.IGNORECASE)
    if not substance or not shortname:
        raise ValueError(f"{path}: missing substance or shortname")
    return substance.group(1).strip().lower(), shortname.group(1).strip().lower(), text.splitlines()


def _available_channels(inst_path: Path) -> list[int]:
    match = re.search(r"^available channels\s*=\s*(.+)$", inst_path.read_text(), re.MULTILINE | re.IGNORECASE)
    if not match:
        raise ValueError(f"{inst_path}: missing available channels")
    return [int(item) for item in match.group(1).split()]


def _required_external_data(mm_path: Path, in_path: Path, tmatrix_path: str | None) -> list[Path]:
    required = []
    text = mm_path.read_text()
    if re.search(r"^component\s+baum\s+", text, re.MULTILINE | re.IGNORECASE):
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 3 and parts[0].lower() == "component" and parts[1].lower() == "baum":
                if parts[2].startswith("/"):
                    required.append(Path(parts[2]))
                else:
                    repository_relative = ROOT / parts[2]
                    required.append(repository_relative if repository_relative.exists() else ROOT / "create_orac_lut" / parts[2])
    if re.search(r"^\s*\*scattering code\s*$", text, re.MULTILINE | re.IGNORECASE) and re.search(
        r"^\s*tmatrix\s*$", text, re.MULTILINE | re.IGNORECASE
    ):
        if not tmatrix_path:
            raise ValueError(f"{mm_path.name}: T-matrix scattering requires tmatrix_path")
        required.append(Path(tmatrix_path))
    return required


def _lut_dimensions(lut_path: Path, inst_path: Path, forward_model: str, channels: list[int]) -> dict[str, int]:
    """Read dimensions through the existing instrument/LUT parsers."""

    from .idl_mirror.load_inststr import load_inststr
    from .idl_mirror.load_lutstr import load_lutstr

    instrument = load_inststr(inst_path, requestedchannelid=channels)
    grid = load_lutstr(lut_path, instrument.max_sat_zenith, include_pressure=forward_model == "aerosol")
    dimensions = {
        "optical depth": int(grid.opd_n),
        "effective radius": int(grid.efr_n),
        "solar zenith": int(grid.soz_n),
        "satellite zenith": int(grid.saz_n),
        "relative azimuth": int(grid.raa_n),
    }
    if hasattr(grid, "prs_n"):
        dimensions["pressure"] = int(grid.prs_n)
    return dimensions


def expand_master(path: str | Path) -> list[dict]:
    """Resolve and validate every requested composition without calculating."""

    master = read_master(path)
    in_path = _rooted(master["in_path"])
    inst_path = in_path / "inst" / str(master["instfile"])
    if not inst_path.is_file():
        raise FileNotFoundError(f"instrument file not found: {inst_path}")
    available = _available_channels(inst_path)
    invalid = [channel for channel in master["channelid"] if channel not in available]
    if invalid:
        raise ValueError(f"invalid channel(s) {invalid}; {inst_path.name} provides {available}")
    atmosphere_files = {"0": "midsatm.dat", "1": "tro.atm", "2": "mls.atm", "3": "mlw.atm", "4": "sas.atm", "5": "saw.atm", "6": "std.atm"}
    atmosphere = str(master["atmospheres"])
    if atmosphere not in atmosphere_files or not (in_path / "atm" / atmosphere_files[atmosphere]).is_file():
        raise FileNotFoundError(f"atmosphere input is not available for code {atmosphere}")
    instrument_text = inst_path.read_text()
    for channel in master["channelid"]:
        srf = re.search(rf"^\s*srf\s*\[{channel}\]\s*=\s*(\S+)", instrument_text, re.MULTILINE | re.IGNORECASE)
        if not srf or not (in_path / "srf" / srf.group(1)).is_file():
            raise FileNotFoundError(f"SRF for channel {channel} is not available in {inst_path.name}")
    if not (in_path / "sun" / "Gueymard2018.sssi").is_file():
        raise FileNotFoundError("Gueymard2018.sssi is missing")
    outputs = []
    for mmfile in master["mmfiles"]:
        mm_path = in_path / "microphysics" / mmfile
        if not mm_path.is_file():
            raise FileNotFoundError(f"microphysics file not found: {mm_path}")
        substance, shortname, lines = _mm_metadata(mm_path)
        if substance not in SUPPORTED_MATERIALS:
            raise ValueError(f"{mmfile}: unsupported material {substance!r}")
        expected_forward_model = "aerosol" if substance == "aerosol" else "cloud"
        if master["forward_model"] != expected_forward_model:
            raise ValueError(f"{mmfile}: {substance} requires forward_model = {expected_forward_model!r}")
        lutfile = LUT_BY_MATERIAL[substance] if master["lutfile"] == "auto" else master["lutfile"]
        lut_path = in_path / "lut" / lutfile
        if not lut_path.is_file():
            raise FileNotFoundError(f"LUT definition not found for {mmfile}: {lut_path}")
        grid_dimensions = _lut_dimensions(lut_path, inst_path, master["forward_model"], master["channelid"])
        for required in _required_external_data(mm_path, in_path, master["tmatrix_path"]):
            if not required.exists():
                raise FileNotFoundError(f"external scientific data not found for {mmfile}: {required}")
        mm_lines = iter(enumerate(lines))
        for line_number, line in mm_lines:
            if line.strip().lower() == "*refractive index":
                _, ri_name = next(mm_lines, (line_number, ""))
                ri_path = in_path / "ri" / ri_name.strip()
                if not ri_path.is_file():
                    raise FileNotFoundError(f"refractive-index file not found for {mmfile}: {ri_path}")
        if master["gas"]:
            for channel in master["channelid"]:
                gas_path = in_path / "gas" / f"ModtranGasOpd_A{atmosphere}_{master['platform']}_{master['instrument']}_ch{channel:02d}.gas"
                if not gas_path.is_file():
                    raise FileNotFoundError(f"gas file not found for {mmfile}: {gas_path}")
        if any("opac" in line.lower() for line in lines):
            raise ValueError(f"{mmfile}: OPAC optical-property branch is not ported")
        if any("baran" in line.lower() for line in lines):
            raise ValueError(f"{mmfile}: Baran optical-property branch is not ported")
        atmospheric_code = "1" + str(master["atmospheres"])[-1] if master["gas"] else ("00" if master["no_rayleigh"] else "01")
        basename = (
            f"{master['platform'].lower()}_{master['instrument'].lower()}_m_{substance}"
            f"_a{atmospheric_code}_p{shortname}_v{int(master['version']):02d}.nc"
        )
        output = _rooted(master["out_path"]) / basename
        values = dict(master)
        values.update({
            "mmfile": mmfile,
            "lutfile": lutfile,
            "output": output,
            "inst_path": inst_path,
            "mm_path": mm_path,
            "grid_dimensions": grid_dimensions,
        })
        outputs.append(values)
    names = [item["output"] for item in outputs]
    if len(set(names)) != len(names):
        raise ValueError("expanded output filenames are not unique")
    return outputs


def concrete_text(settings: dict) -> str:
    """Render one transient concrete run file for create_orac_luts.py."""

    ordered = (
        "platform", "instrument", "forward_model", "in_path", "instfile", "mmfile", "lutfile",
        "atmospheres", "channelid", "gas", "no_rayleigh", "srf_quad", "reuse_scat", "scat_only",
        "nstreams", "nmom", "version", "out_path", "tmatrix_path",
    )
    lines = ["; Generated transiently from a maker-style master specification; do not edit."]
    for key in ordered:
        value = settings.get(key)
        if value is not None:
            lines.append(f"{key} = {value!r}")
    return "\n".join(lines) + "\n"


def write_manifest(master_path: str | Path, manifest_path: str | Path) -> list[dict]:
    """Write the validated concrete task specifications for one array."""

    expanded = expand_master(master_path)
    manifest = {
        "master": str(Path(master_path)),
        "tasks": [
            {
                "index": index,
                "mmfile": item["mmfile"],
                "output": str(item["output"]),
                "run": concrete_text(item),
            }
            for index, item in enumerate(expanded)
        ],
    }
    destination = Path(manifest_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    return expanded


def write_manifest_task(manifest_path: str | Path, index: int, output: str | Path) -> None:
    """Materialise one previously validated manifest task."""

    manifest = json.loads(Path(manifest_path).read_text())
    tasks = manifest.get("tasks", [])
    if index < 0 or index >= len(tasks):
        raise IndexError(f"manifest task index {index} is outside 0..{len(tasks) - 1}")
    task = tasks[index]
    if task.get("index") != index or not isinstance(task.get("run"), str):
        raise ValueError(f"manifest task {index} is malformed")
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(task["run"])


def print_preflight(path: str | Path) -> None:
    master = read_master(path)
    expanded = expand_master(path)
    print(f"Master expansion preflight: {Path(path)}")
    print()
    print(f"Platform:             {master['platform']}")
    print(f"Instrument:           {master['instrument']}")
    print(f"Forward model:        {master['forward_model']}")
    print(f"Python LUT version:   {master['version']}")
    print()
    print(f"Instrument file:      {master['instfile']}")
    print(f"Microphysics:         {', '.join(item['mmfile'] for item in expanded)}")
    print(f"LUT definition:       {', '.join(sorted({item['lutfile'] for item in expanded}))}")
    print()
    print(f"Channels:             {master['channelid']}")
    print(f"Atmosphere:           {master['atmospheres']}")
    print(f"Gas absorption:       {'on' if master['gas'] else 'off'}")
    print(f"Rayleigh scattering:  {'off' if master['no_rayleigh'] else 'on'}")
    print(f"Scattering only:      {'yes' if master['scat_only'] else 'no'}")
    print(f"SRF quadrature:       {master['srf_quad']} ({SRF_QUADRATURE_NAMES[master['srf_quad']]})")
    print(f"DISORT streams:       {master['nstreams']}")
    if master["nmom"] is None:
        print("Legendre moments:     adaptive (King's criterion on each averaged Mie phase function)")
    else:
        print(f"Legendre moments:     {master['nmom']} (fixed; Baum / T-matrix tabulated phase functions)")
    print()
    print("LUT grid:")
    for name, count in expanded[0]["grid_dimensions"].items():
        print(f"  {name.title():<21} {count}")
    print()
    print(f"Output directory:     {master['out_path']}")
    print("Calculations:")
    for index, item in enumerate(expanded):
        existing = " [EXISTS; single-LUT generator will refuse overwrite]" if item["output"].exists() else ""
        print(f"  task {index}:")
        print(f"    Microphysics:     {item['mmfile']}")
        print(f"    LUT definition:   {item['lutfile']}")
        print(f"    Output:            {item['output'].relative_to(ROOT)}{existing}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Expand or preflight a maker-style ORAC master run")
    parser.add_argument("master", type=Path)
    parser.add_argument("--preflight", action="store_true", help="validate and list without calculation")
    parser.add_argument("--count", action="store_true", help="print the number of expanded calculations")
    parser.add_argument("--index", type=int, help="write one transient concrete run to --output")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--write-manifest", type=Path, help="write the validated array task manifest")
    parser.add_argument("--manifest-task", type=Path, help="materialise one task from a validated manifest")
    args = parser.parse_args(argv)
    if args.write_manifest is not None:
        write_manifest(args.master, args.write_manifest)
    elif args.manifest_task is not None:
        if args.index is None or args.output is None:
            parser.error("--manifest-task requires --index and --output")
        write_manifest_task(args.manifest_task, args.index, args.output)
    else:
        expanded = expand_master(args.master)
        if args.count:
            print(len(expanded))
        elif args.index is not None:
            if args.output is None:
                parser.error("--index requires --output")
            args.output.write_text(concrete_text(expanded[args.index]))
        else:
            print_preflight(args.master)
    return 0


if __name__ == "__main__":
    sys.exit(main())
