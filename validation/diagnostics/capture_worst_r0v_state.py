"""Capture one exact Python DISORT state for the worst full-LUT R_0v point."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from oraclut.config import (  # noqa: E402
    read_atmosphere,
    read_instrument,
    read_lut_grid,
    read_microphysics,
    read_refractive_index,
)
from oraclut.pipeline import _layer_inputs, _view_cosines, read_channels  # noqa: E402
from oraclut.radiative_transfer import call_disort  # noqa: E402
from oraclut.optics import legacy_stg_optics  # noqa: E402


def main() -> None:
    print("capture: start", flush=True)
    output = ROOT / "validation" / "diagnostics"
    grid = read_lut_grid(ROOT / "create_orac_lut/input_files/lut/liquid-water-cloud.lut")
    instrument = read_instrument(ROOT / "create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst")
    model = read_microphysics(ROOT / "create_orac_lut/input_files/microphysics/liquid-water_stg.mm")
    atmosphere = read_atmosphere(ROOT / "create_orac_lut/input_files/atm/mls.atm", 2)
    channels = read_channels(ROOT, instrument, (1,), srf_quad=1)
    wavelength, real, imaginary = read_refractive_index(
        ROOT / "create_orac_lut/input_files/ri" / model.refractive_index_file
    )
    wavelengths = np.asarray([0.55, float(channels[0].wavelength_microns)], dtype=np.float64)
    refractive_indices = (
        np.interp(wavelengths, wavelength, real).astype(np.float32)
        - 1j * np.interp(wavelengths, wavelength, imaginary).astype(np.float32)
    ).astype(np.complex64)
    target_radius = np.asarray([grid.effective_radius[11]], dtype=np.float64)
    optics = legacy_stg_optics(target_radius, wavelengths, refractive_indices, phase_order=1000)
    print("capture: optics", flush=True)
    extinction_ratio = (
        optics.extinction_coefficient[:, 1:] / optics.reference_extinction_coefficient[:, None]
    ).astype(np.float32)
    ssa_value = optics.single_scatter_albedo[0, 1].astype(np.float32)
    asymmetry_value = optics.asymmetry_parameter[0, 1].astype(np.float32)
    moments = optics.phase_moments[:, 0, 1].astype(np.float32)
    relative_tau, tau_rayleigh, _ = _layer_inputs(
        atmosphere,
        model.profile_height_km,
        model.profile_relative_amount,
        no_rayleigh=False,
        wavelength_microns=float(channels[0].wavelength_microns),
    )
    optical_depth = np.float32(grid.optical_depth[9])
    tauscat = (optical_depth * relative_tau * extinction_ratio[0, 0]).astype(np.float32)
    dtau = (tau_rayleigh + tauscat).astype(np.float32)
    total_tau = np.sum(dtau, dtype=np.float32)
    ssalb = np.divide(
        tau_rayleigh + ssa_value * tauscat,
        dtau,
        out=np.zeros_like(dtau),
        where=dtau != 0,
    ).astype(np.float32)
    ssalb[ssalb > 1.0] = 1.0
    pmom = np.empty((1000, dtau.size), dtype=np.float32, order="F")
    molecular = __import__("oraclut.radiative_transfer", fromlist=["getmom"]).getmom(2, 0.0, 999)
    for layer in range(dtau.size):
        pmom[:, layer] = (
            molecular * tau_rayleigh[layer] + moments * ssa_value * tauscat[layer]
        ) / (tau_rayleigh[layer] + ssa_value * tauscat[layer])
    pmom[pmom > 1.0] = 1.0
    utau = np.asarray([0.0, total_tau], dtype=np.float32)
    umu = _view_cosines(grid.satellite_zenith)
    phi = np.asarray(grid.relative_azimuth, dtype=np.float32)
    umu0 = np.float32(np.cos(np.deg2rad(np.float32(50.0))))
    # DISORT's delta-M setup can modify its array arguments in place. Snapshot
    # the actual pre-call state for reproducible cross-build comparisons.
    input_arrays = {
        "dtau": dtau.copy(), "ssalb": ssalb.copy(), "pmom": pmom.copy(order="F"),
        "utau": utau.copy(), "umu": umu.copy(), "phi": phi.copy(),
    }
    result = call_disort(
        dtau, ssalb, pmom, utau, umu, phi, 100.0, umu0, 0.0, nstreams=60
    )
    print("capture: disort", flush=True)
    print("capture: target_r0v", float(result["uu"][11, 0, 10] * np.pi * 0.01), flush=True)
    repeat = call_disort(
        input_arrays["dtau"].copy(), input_arrays["ssalb"].copy(),
        input_arrays["pmom"].copy(order="F"), input_arrays["utau"].copy(),
        input_arrays["umu"].copy(), input_arrays["phi"].copy(),
        100.0, umu0, 0.0, nstreams=60,
    )
    print("capture: repeat_r0v", float(repeat["uu"][11, 0, 10] * np.pi * 0.01), flush=True)
    tag = os.environ.get("ORACLUT_DIAGNOSTIC_TAG", "default")
    prefix = output / f"worst_r0v_disort_state_{tag}"
    arrays = {
        **input_arrays,
        "rfldir": result["rfldir"],
        "rfldn": result["rfldn"],
        "flup": result["flup"],
        "uu": result["uu"],
    }
    np.savez(prefix.with_suffix(".npz"), **arrays)
    # The raw files are column-major where that is relevant to IDL READU.
    for name, value in arrays.items():
        path = prefix.parent / f"{prefix.name}_{name}.bin"
        order = "F" if name in {"pmom", "uu"} else "C"
        path.write_bytes(np.asarray(value, dtype=np.float32).ravel(order=order).tobytes())
    summary = {
        "point": {
            "channel_id": 1,
            "solar_channel_index": 0,
            "wavelength_microns": float(channels[0].wavelength_microns),
            "optical_depth_index": 9,
            "optical_depth": float(optical_depth),
            "effective_radius_index": 11,
            "effective_radius_microns": float(grid.effective_radius[11]),
            "solar_zenith_index": 5,
            "solar_zenith": float(grid.solar_zenith[5]),
            "satellite_zenith_index": 8,
            "satellite_zenith": float(grid.satellite_zenith[8]),
            "relative_azimuth_index": 0,
            "relative_azimuth": float(grid.relative_azimuth[0]),
        },
        "disort": {
            "nlyr": int(dtau.size), "nstr": 60, "nmom": 999,
            "ntau": 2, "numu": int(umu.size), "nphi": int(phi.size),
            "fbeam": 100.0, "umu0": float(umu0), "fisot": 0.0,
            "phi0": 0.0, "usrtau": 1, "usrang": 1,
            "plank": 0, "onlyfl": 0, "accur": 1.0e-8,
        },
        "compiler": os.environ.get("ORACLUT_FORTRAN", "gfortran"),
        "arrays": {
            name: {"shape": list(value.shape), "dtype": str(value.dtype),
                   "minimum": float(np.min(value)), "maximum": float(np.max(value))}
            for name, value in arrays.items()
        },
        "optics": {
            "ssa": float(ssa_value), "asymmetry": float(asymmetry_value),
            "extinction_ratio": float(extinction_ratio[0, 0]),
            "phase_moment_max": float(np.max(moments)),
        },
    }
    (prefix.with_suffix(".json")).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
