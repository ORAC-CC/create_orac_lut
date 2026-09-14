"""Capture the per-layer DISORT inputs the Python aerosol pipeline assembles for one state.

Diagnostic only. It calls the same production helpers the pipeline uses
(``read_*``, ``_interpolate_profile``, ``_gas_levels``, ``legacy_mixed_mie_optics``,
``getmom``) and repeats the layer arithmetic of ``calculate_aerosol_operators``
verbatim for a single (channel, pressure, radius, optical depth) state, so the
arrays can be written out without instrumenting the scientific code.  Nothing
in ``src/`` is modified or monkey-patched.

Usage:
    PYTHONPATH=src python validation/diagnostics/aerosol_layer_compare/capture_python_state.py \
        --channel 1 --effective-radius 0.01 --optical-depth 1.0 --surface-pressure 950 \
        --out validation/diagnostics/aerosol_layer_compare/python_state.npz
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from oraclut.config import (
    read_atmosphere, read_instrument, read_lut_grid, read_microphysics,
)
from oraclut.optics import legacy_mixed_mie_optics
from oraclut.pipeline import (
    _atmosphere_path, _component_refractive_indices, _gas_levels, _interpolate_profile,
    read_channels,
)
from oraclut.radiative_transfer import getmom


def capture(root: Path, *, channel: int, effective_radius: float, optical_depth: float,
            surface_pressure: float, atmosphere_code: int = 2, phase_order: int = 1000) -> dict:
    inputs = root / "create_orac_lut" / "input_files"
    instrument = read_instrument(inputs / "inst" / "meteosat-10_seviri_v1.inst")
    grid = read_lut_grid(inputs / "lut" / "aerosol_test.lut")
    model = read_microphysics(inputs / "microphysics" / "aerosol_a79.mm")
    atmosphere = read_atmosphere(_atmosphere_path(inputs, atmosphere_code), atmosphere_code)
    channels = read_channels(root, instrument, (channel,), srf_quad=1)
    wavelengths = np.asarray([0.55] + [float(c.wavelength_microns) for c in channels], dtype=np.float64)
    optics = legacy_mixed_mie_optics(
        grid.effective_radius, wavelengths, model.components,
        _component_refractive_indices(root, model, wavelengths), phase_order=phase_order,
    )
    radius_index = int(np.argmin(np.abs(grid.effective_radius - effective_radius)))
    extinction_ratio = (optics.extinction_coefficient[:, 1:] / optics.reference_extinction_coefficient[:, None]).astype(np.float32)
    ssa_particle = optics.single_scatter_albedo[:, 1:].astype(np.float32)
    phase_moments = optics.phase_moments[:, :, 1:].astype(np.float32)
    ch = channels[0]

    # --- identical arithmetic to calculate_aerosol_operators ---------------
    n_moments = phase_moments.shape[0]
    molecular_moments = getmom(2, 0.0, n_moments - 1)
    n_layers = atmosphere.height_km.size - 1
    layer_height = (atmosphere.height_km[:n_layers] + atmosphere.height_km[1:]) / np.float32(2.0)
    relative_tau_raw = _interpolate_profile(model.profile_height_km, model.profile_relative_amount, layer_height)
    relative_tau = relative_tau_raw / np.sum(relative_tau_raw, dtype=np.float32)
    gas = np.asarray(_gas_levels(root, instrument, atmosphere_code, channel), dtype=np.float32)
    tau_gas = gas[1:] - gas[:-1]
    column_tau = np.float32(atmosphere.pressure_hpa[-1] / np.float32(surface_pressure)) / (
        np.float32(117.03) * ch.wavelength_microns ** np.float32(4.0)
        - np.float32(1.316) * ch.wavelength_microns ** np.float32(2.0)
    )
    rayleigh_level = column_tau * np.exp(
        np.float32(-0.1188) * atmosphere.height_km - np.float32(0.00116) * atmosphere.height_km ** np.float32(2.0)
    )
    tau_rayleigh = (rayleigh_level[1:] - rayleigh_level[:-1]).astype(np.float32)
    w = ssa_particle[radius_index, 0]
    tauscat = np.asarray(np.float32(optical_depth) * relative_tau * extinction_ratio[radius_index, 0], dtype=np.float32)
    dtau = (tau_gas + tau_rayleigh + tauscat).astype(np.float32)
    ssa = np.divide(tau_rayleigh + w * tauscat, dtau, out=np.zeros_like(dtau), where=dtau != 0).astype(np.float32)
    ssa[ssa > 1.0] = 1.0
    pmom = np.empty((n_moments, dtau.size), dtype=np.float32, order="F")
    for layer in range(dtau.size):
        if tauscat[layer] <= 0.0:
            pmom[:, layer] = molecular_moments
            continue
        scattering_tau = tau_rayleigh[layer] + w * tauscat[layer]
        pmom[:, layer] = molecular_moments if scattering_tau == 0.0 else (
            molecular_moments * tau_rayleigh[layer] + phase_moments[:, radius_index, 0] * w * tauscat[layer]
        ) / scattering_tau
    pmom[pmom > 1.0] = 1.0
    # -----------------------------------------------------------------------

    return {
        "state": json.dumps({
            "channel": channel, "effective_radius_um": float(grid.effective_radius[radius_index]),
            "optical_depth": optical_depth, "surface_pressure_hpa": surface_pressure,
            "atmosphere_code": atmosphere_code, "wavelength_um": float(ch.wavelength_microns),
            "particle_ssa": float(w), "extinction_ratio": float(extinction_ratio[radius_index, 0]),
            "particle_asymmetry": float(optics.asymmetry_parameter[radius_index, 1]),
        }),
        "level_height_km": atmosphere.height_km.astype(np.float32),
        "level_pressure_hpa": atmosphere.pressure_hpa.astype(np.float32),
        "layer_height_km": layer_height.astype(np.float32),
        "profile_height_km": model.profile_height_km.astype(np.float32),
        "profile_relative_amount": model.profile_relative_amount.astype(np.float32),
        "relative_tau_raw": relative_tau_raw.astype(np.float32),
        "relative_tau": relative_tau.astype(np.float32),
        "gas_level": gas, "tau_gas": tau_gas.astype(np.float32),
        "column_tau_rayleigh": np.float32(column_tau), "rayleigh_level": rayleigh_level.astype(np.float32),
        "tau_rayleigh": tau_rayleigh, "tauscat": tauscat,
        "tau_aerosol_scattering": (w * tauscat).astype(np.float32),
        "tau_aerosol_absorption": ((np.float32(1.0) - w) * tauscat).astype(np.float32),
        "dtau": dtau, "ssa": ssa, "pmom": np.ascontiguousarray(pmom),
        "molecular_moments": molecular_moments.astype(np.float32),
        "particle_moments": phase_moments[:, radius_index, 0].astype(np.float32),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--channel", type=int, default=1)
    parser.add_argument("--effective-radius", type=float, default=0.01)
    parser.add_argument("--optical-depth", type=float, default=1.0)
    parser.add_argument("--surface-pressure", type=float, default=950.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    arrays = capture(root, channel=args.channel, effective_radius=args.effective_radius,
                     optical_depth=args.optical_depth, surface_pressure=args.surface_pressure)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, **arrays)
    print(f"wrote {args.out}: {arrays['dtau'].size} layers, {arrays['pmom'].shape[0]} moments; state {arrays['state']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
