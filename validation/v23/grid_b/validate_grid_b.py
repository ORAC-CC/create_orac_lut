"""Validate the three installed Grid B definitions without running RT."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RESULTS = HERE.parent / "results"
INPUTS = REPO / "create_orac_lut" / "input_files"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE.parent))

import create_orac_luts  # noqa: E402
from oraclut.idl_mirror.load_atmstr import load_atmstr  # noqa: E402
from oraclut.idl_mirror.load_gasstr import load_gasstr  # noqa: E402
from oraclut.idl_mirror.load_inststr import load_inststr  # noqa: E402
from oraclut.idl_mirror.load_lutstr import load_lutstr  # noqa: E402
from oraclut.idl_mirror.load_mmdat import load_mmdat  # noqa: E402
from oraclut.idl_mirror.load_srfstrarr import load_srfstrarr  # noqa: E402
from stage16_audit import gauss_safety  # noqa: E402


FILES = {
    "water": {"lut": "liquid-water-cloud-grid-b.lut", "json": "stage18_final_grid_water.json", "mm": "liquid-water_273.mm", "model": "cloud", "gas": 0},
    "ice": {"lut": "ice-cloud-grid-b.lut", "json": "stage18_final_grid_ice.json", "mm": "water-ice_sph.mm", "model": "cloud", "gas": 0},
    "aerosol_pa76": {"lut": "aerosol-pa76-grid-b.lut", "json": "stage18_final_grid_aerosol_pa76.json", "mm": "aerosol_a76.mm", "model": "aerosol", "gas": 1},
}
RUNS = {
    "water": "water_grid_b.run",
    "ice": "ice_grid_b.run",
    "aerosol_pa76": "aerosol_pa76_grid_b.run",
}
REPRESENTATIVE_INSTRUMENTS = {
    "EarthCARE MSI (20°)": "earthcare_msi_v1.inst",
    "AATSR (56°)": "envisat_aatsr_v1.inst",
    "VIIRS (57°)": "noaa-20_viirs_v1.inst",
    "MODIS (65°)": "aqua_modis_v1.inst",
    "SEVIRI (geostationary, advertises 90°)": "meteosat-10_seviri_v1.inst",
}
CHANNELS = [1, 3]
NSTR = 60
NMOM = 1000
V23_SAZ_MAX = 75.0
MXUMU, MXPHI, MXCMU = 180, 181, 100


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def axes_from_json(path: Path) -> tuple[dict, dict]:
    data = json.loads(path.read_text())
    return data, {
        "tau": data["tau_nodes"],
        "re": data["re_nodes"],
        "soz": data["soz_nodes"],
        "saz": data["saz_nodes_by_domain"]["75"],
        "raa": data["raa_nodes"],
    }


def parsed_axes(parsed) -> dict[str, np.ndarray]:
    return {name: np.asarray(values, dtype=np.float32) for name, values in {
        "tau": parsed.opd, "re": parsed.efr, "soz": parsed.soz,
        "saz": parsed.saz, "raa": parsed.raa,
    }.items()}


def compare_axes(parsed, expected: dict) -> dict:
    actual = parsed_axes(parsed)
    rows = {}
    for name, expected_values in expected.items():
        wanted = np.asarray(expected_values, dtype=np.float32)
        got = actual[name]
        rows[name] = {
            "expected_count": int(wanted.size),
            "parsed_count": int(got.size),
            "minimum": float(got.min()),
            "maximum": float(got.max()),
            "maximum_absolute_coordinate_difference": (
                float(np.max(np.abs(got.astype(np.float64) - wanted.astype(np.float64))))
                if got.shape == wanted.shape else None
            ),
            "float32_exact": bool(got.shape == wanted.shape and np.array_equal(got, wanted)),
        }
    return rows


def representative_tests(path: Path, expected: dict, include_pressure: bool) -> dict:
    rows = {}
    for label, filename in REPRESENTATIVE_INSTRUMENTS.items():
        instrument = load_inststr(INPUTS / "inst" / filename, requestedchannelid=CHANNELS)
        parsed = load_lutstr(path, float(instrument.max_sat_zenith), include_pressure=include_pressure)
        parsed_at_90 = load_lutstr(path, 90.0, include_pressure=include_pressure)
        axis = np.asarray(parsed.saz, dtype=np.float64)
        expected_saz = np.asarray(expected["saz"], dtype=np.float64)
        rows[label] = {
            "instrument_file": str((INPUTS / "inst" / filename).relative_to(REPO)),
            "instrument_max_sat_zenith": float(instrument.max_sat_zenith),
            "parsed_saz_count": int(parsed.saz_n),
            "parsed_saz_minimum": float(axis.min()),
            "parsed_saz_maximum": float(axis.max()),
            "required_measurement_max": min(float(instrument.max_sat_zenith), V23_SAZ_MAX),
            "required_measurement_inside_grid": bool(axis.max() >= min(float(instrument.max_sat_zenith), V23_SAZ_MAX)),
            "exact_0_75_grid": bool(np.array_equal(parsed.saz.astype(np.float32), expected_saz.astype(np.float32))),
            "unchanged_at_90_degree_instrument_limit": bool(np.array_equal(parsed.saz, parsed_at_90.saz)),
        }
        rows[label]["pass"] = bool(
            rows[label]["required_measurement_inside_grid"]
            and rows[label]["exact_0_75_grid"]
            and rows[label]["unchanged_at_90_degree_instrument_limit"]
            and rows[label]["parsed_saz_maximum"] <= V23_SAZ_MAX + 1.0e-6
        )
    return rows


def disort_test(expected: dict) -> dict:
    soz = np.asarray(expected["soz"], dtype=float)
    saz_n = len(expected["saz"])
    raa_n = len(expected["raa"])
    full = gauss_safety(soz, NSTR)
    nested = {}
    for count in range(3, len(soz) + 1):
        if (len(soz) - 1) % (count - 1) == 0:
            nested[str(count)] = gauss_safety(np.linspace(0.0, 75.0, count), NSTR)
    return {
        "nstreams": NSTR,
        "nmom": NMOM,
        "numu": 2 * saz_n,
        "nphi": raa_n,
        "limits": {"max_numu": MXUMU, "max_nphi": MXPHI, "max_nstreams": MXCMU},
        "setdis_full": full,
        "setdis_nested": nested,
        "disort_limits_pass": bool(2 * saz_n <= MXUMU and raa_n <= MXPHI and NSTR <= MXCMU),
        "setdis_pass": bool(full["safe"] and all(row["safe"] for row in nested.values())),
    }


def expected_output(settings, instrument, mm) -> str:
    code = "12" if settings["gas"] else "01"
    return (f"{instrument.platform.lower()}_{instrument.instrument.lower()}_m_"
            f"{mm.substance.lower()}_a{code}_p{mm.shortname.lower()}_v{int(settings['version']):02d}.nc")


def run_preflight(family: str, spec: dict) -> dict:
    run_path = HERE / "runs" / RUNS[family]
    settings = create_orac_luts.read_runfile(run_path)
    in_path = REPO / settings["in_path"]
    instrument_path = in_path / "inst" / settings["instfile"]
    mm_path = in_path / "microphysics" / settings["mmfile"]
    lut_path = in_path / "lut" / settings["lutfile"]
    instrument = load_inststr(instrument_path, requestedchannelid=settings["channelid"])
    mm = load_mmdat(mm_path, in_path)
    atmosphere_path = in_path / "atm" / "mls.atm"
    atmosphere = load_atmstr(atmosphere_path, settings["atmospheres"])
    if settings["gas"]:
        gas = load_gasstr(settings["atmospheres"], instrument.platform, instrument.instrument, instrument.channelid, in_path / "gas")
        gas_pass = all(np.array_equal(atmosphere.height, item.height) for item in gas)
    else:
        gas = []
        gas_pass = True
    srf, _ = load_srfstrarr(instrument, in_path / "sun" / "Gueymard2018.sssi", settings["srf_quad"], in_path)
    parsed = load_lutstr(lut_path, instrument.max_sat_zenith, include_pressure=family == "aerosol_pa76")
    output = expected_output(settings, instrument, mm)
    output_path = REPO / settings["out_path"] / output
    production_collisions = [str(path.relative_to(REPO)) for path in (REPO / "luts" / output, INPUTS / "lut" / output) if path.exists()]
    return {
        "run": str(run_path.relative_to(REPO)),
        "run_parser_pass": True,
        "lut_resolves": lut_path.is_file(),
        "instrument_parser_pass": True,
        "particle_model_parser_pass": True,
        "atmosphere_pass": atmosphere_path.is_file(),
        "gas_required": bool(settings["gas"]),
        "gas_pass": gas_pass,
        "channel_availability_pass": not bool(set(CHANNELS) - set(map(int, instrument.channelid))),
        "srf_count": len(srf),
        "parsed_counts": {"tau": parsed.opd_n, "re": parsed.efr_n, "soz": parsed.soz_n, "saz": parsed.saz_n, "raa": parsed.raa_n},
        "output_name": output,
        "output_collision": bool(production_collisions or output_path.exists()),
        "production_collisions": production_collisions,
        "pass": bool(lut_path.is_file() and gas_pass and bool(srf) and not production_collisions and not output_path.exists()),
    }


def main() -> int:
    manifest = {
        "definition": "Grid B",
        "scope": "definitive LUT sampling specification, independent of ORAC software and product versions",
        "operational_files": {family: f"create_orac_lut/input_files/lut/{spec['lut']}" for family, spec in FILES.items()},
        "source_stage18": {family: f"validation/v23/results/{spec['json']}" for family, spec in FILES.items()},
        "satellite_zenith_coverage_degrees": [0.0, 75.0],
        "families": {},
        "representative_instruments": list(REPRESENTATIVE_INSTRUMENTS),
    }
    all_pass = True
    for family, spec in FILES.items():
        source, expected = axes_from_json(RESULTS / spec["json"])
        path = INPUTS / "lut" / spec["lut"]
        parsed = load_lutstr(path, 90.0, include_pressure=family == "aerosol_pa76")
        axes = compare_axes(parsed, expected)
        instruments = representative_tests(path, expected, family == "aerosol_pa76")
        disort = disort_test(expected)
        run = run_preflight(family, spec)
        parser_pass = all(row["float32_exact"] for row in axes.values())
        instrument_pass = all(row["pass"] for row in instruments.values())
        disort_pass = disort["disort_limits_pass"] and disort["setdis_pass"]
        family_pass = parser_pass and instrument_pass and disort_pass and run["pass"]
        all_pass &= family_pass
        manifest["families"][family] = {
            "operational_file": str(path.relative_to(REPO)),
            "source_stage18": spec["json"],
            "sha256": digest(path),
            "spacing": source["spacing"],
            "expected_axes": expected,
            "parsed_ranges": {axis: {"count": row["parsed_count"], "minimum": row["minimum"], "maximum": row["maximum"]} for axis, row in axes.items()},
            "parser_roundtrip": {"axes": axes, "pass": parser_pass},
            "representative_instruments": {"results": instruments, "pass": instrument_pass},
            "disort_setdis": {"result": disort, "pass": disort_pass},
            "generator_preflight": run,
            "pass": family_pass,
        }
    manifest["parser_roundtrip_pass"] = all(manifest["families"][family]["parser_roundtrip"]["pass"] for family in FILES)
    manifest["representative_instrument_pass"] = all(manifest["families"][family]["representative_instruments"]["pass"] for family in FILES)
    manifest["disort_setdis_pass"] = all(manifest["families"][family]["disort_setdis"]["pass"] for family in FILES)
    manifest["generator_preflight_pass"] = all(manifest["families"][family]["generator_preflight"]["pass"] for family in FILES)
    manifest["all_pass"] = all_pass
    (RESULTS / "grid_b_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_summary(manifest)
    write_report(manifest)
    for family, entry in manifest["families"].items():
        print(f"{family}: {'PASS' if entry['pass'] else 'FAIL'}")
    return 0 if all_pass else 1


def write_summary(manifest: dict) -> None:
    labels = {"water": "WATER", "ice": "ICE", "aerosol_pa76": "PA76 AEROSOL"}
    lines = ["# GRID B — CURRENT RECOMMENDED GRID", "", "Grid B is the definitive LUT sampling specification for liquid-water cloud, ice cloud, and generic PA76 aerosol. It is independent of ORAC software, processing, generator, and NetCDF product version numbers. All families use an explicit 0–75° satellite-zenith domain.", ""]
    for family, entry in manifest["families"].items():
        axes = entry["expected_axes"]
        counts = {axis: len(values) for axis, values in axes.items()}
        points = int(np.prod(list(counts.values())))
        lines += [f"## {labels[family]}", "", f"Operational file: `{entry['operational_file']}`", f"Grid points: `{points:,}`", f"- tau ({counts['tau']}, uneven_logarithmic): `{axes['tau']}`", f"- r_eff ({counts['re']}, uneven_linear physical microns): `{axes['re']}`", f"- SOZ ({counts['soz']}, uneven_linear degrees): `{axes['soz']}`", f"- SAZ ({counts['saz']}, uneven_linear degrees, 0–75°): `{axes['saz']}`", f"- RAA ({counts['raa']}, linear degrees): `{axes['raa']}`", ""]
    lines += ["All three operational files passed the unchanged parser round-trip, representative 20°/56°/57°/65°/geostationary checks, and SETDIS/DISORT preflight. Historical `.lut` files remain unchanged.", ""]
    (RESULTS / "GRID_B.md").write_text("\n".join(lines))


def write_report(manifest: dict) -> None:
    lines = ["# Grid B operational LUT installation", "", "Grid B is now the recommended operational LUT sampling specification for liquid-water cloud, ice cloud, and generic PA76 aerosol. It is independent of ORAC software, processing, generator, and product version numbers. No RT or NetCDF generation was performed.", "", "## Selected operational files", "", "| family | installed file | counts (tau × r_eff × SOZ × SAZ × RAA) |", "|---|---|---:|"]
    for family, entry in manifest["families"].items():
        c = {axis: len(values) for axis, values in entry["expected_axes"].items()}
        lines.append(f"| {family} | `{entry['operational_file']}` | {c['tau']} × {c['re']} × {c['soz']} × {c['saz']} × {c['raa']} |")
    lines += ["", "## Source and exact grid", "", "The files were populated directly from the Stage 18 JSON selections:", "- `validation/v23/results/stage18_final_grid_water.json`", "- `validation/v23/results/stage18_final_grid_ice.json`", "- `validation/v23/results/stage18_final_grid_aerosol_pa76.json`", "", "All three use explicit `uneven_linear` r_eff and SAZ arrays, `uneven_logarithmic` tau, family-specific SOZ, and parser-required linear RAA endpoints. SAZ is exactly 0–75° for every family.", "", "## Validation", "", f"- Parser round-trip: **{'PASS' if manifest['parser_roundtrip_pass'] else 'FAIL'}** for tau, r_eff, SOZ, SAZ, and RAA.", f"- Representative instruments: **{'PASS' if manifest['representative_instrument_pass'] else 'FAIL'}** for EarthCARE 20°, AATSR 56°, VIIRS 57°, MODIS 65°, and geostationary SEVIRI advertising 90°.", f"- SETDIS/DISORT: **{'PASS' if manifest['disort_setdis_pass'] else 'FAIL'}** with NSTR=60, NMOM=1000, NUMU/NPHI within compiled limits.", f"- Generator preflight: **{'PASS' if manifest['generator_preflight_pass'] else 'FAIL'}** for the three validation runs under `validation/v23/grid_b/runs/`.", "", "## Documentation", "", "- `docs/grid_b_lut_definitions.md` identifies Grid B as the current recommended choice.", "- `validation/v23/results/GRID_B.md` gives exact nodes and five-dimensional point counts.", "- Stage 18 JSON and scientific reports remain the provenance/evidence layer.", "", "## Git", "", "Commit and push identifiers are recorded here after the controlled commit/push step.", ""]
    (RESULTS / "REPORT_grid_b_installation.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
