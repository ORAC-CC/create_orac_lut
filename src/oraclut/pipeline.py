"""Executable legacy-equivalent ORAC cloud LUT pipeline.

This module keeps the historical cloud calculation order explicit.  The Mie
and DISORT scientific kernels are the preserved repository sources; Python
replaces the IDL orchestration and NetCDF assembly around them.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from .config import (
    AtmosphereProfile,
    InstrumentConfig,
    LutGrid,
    ReferenceConfiguration,
    read_atmosphere,
    read_instrument,
    read_lut_grid,
    read_microphysics,
    read_refractive_index,
    read_solar_spectrum,
    read_srf,
)
from .io.v2 import write_v2_lut
from .optics import legacy_mixed_mie_optics, legacy_stg_optics
from .radiative_transfer import call_disort, getmom, plkavg


@dataclass(frozen=True)
class Channel:
    channel_id: int
    srf_filename: str
    wavelength_microns: np.float32
    wavenumber_cm_inverse: np.float32
    solar: bool
    mixed: bool
    thermal: bool
    f0: np.float32
    b1: np.float32
    b2: np.float32
    t1: np.float32
    t2: np.float32


@dataclass(frozen=True)
class GenerationResult:
    variables: dict[str, np.ndarray]
    variable_dimensions: dict[str, tuple[str, ...]]
    variable_attributes: dict[str, dict[str, object]]
    dimensions: dict[str, int]
    channels: tuple[Channel, ...]


def _legacy_integral(x: np.ndarray, y: np.ndarray) -> np.float32:
    x = np.asarray(x, dtype=np.float32)
    y = np.asarray(y, dtype=np.float32)
    return np.sum(
        (x[1:] - x[:-1]) * (y[1:] + y[:-1]) * np.float32(0.5),
        dtype=np.float32,
    )


def _legacy_simpson_bug(x: np.ndarray, y: np.ndarray) -> np.float32:
    """Match the current ``simpson_integral.pro`` implementation exactly."""

    return np.sum(np.asarray(x) * np.asarray(y), dtype=np.float32)


def _solar_f0(
    srf_wavenumber: np.ndarray,
    srf_response: np.ndarray,
    centre_wavelength: np.float32,
    solar_wavelength: np.ndarray,
    solar_value: np.ndarray,
) -> np.float32:
    wvl = np.float32(1.0e4) / np.asarray(srf_wavenumber, dtype=np.float32)
    response = np.asarray(srf_response, dtype=np.float32)
    solar_wavelength = np.asarray(solar_wavelength, dtype=np.float32)
    solar_value = np.asarray(solar_value, dtype=np.float32)
    if centre_wavelength > 3.0:
        solar_wavenumber = np.float64(1.0e4) / solar_wavelength.astype(np.float64)
        irradiance = np.float64(1.0e8) * solar_value.astype(np.float64) / solar_wavenumber**2
        order = np.argsort(solar_wavenumber)
        interpolated = np.interp(
            srf_wavenumber.astype(np.float64), solar_wavenumber[order], irradiance[order]
        ).astype(np.float32)
        numerator = _legacy_simpson_bug(srf_wavenumber, interpolated * response)
        denominator = _legacy_simpson_bug(srf_wavenumber, response)
    else:
        interpolated = np.interp(
            wvl.astype(np.float64), solar_wavelength.astype(np.float64), solar_value.astype(np.float64)
        ).astype(np.float32)
        numerator = _legacy_simpson_bug(wvl, interpolated * np.float32(10.0) * response)
        denominator = _legacy_simpson_bug(wvl, response)
    return np.float32(numerator / denominator / np.float32(np.pi))


def _blackbody_constants(root: Path) -> dict[str, tuple[float, float, float, float]]:
    path = root / "create_orac_lut" / "input_files" / "bbconstants.pro"
    text = path.read_text()
    result: dict[str, tuple[float, float, float, float]] = {}
    pattern = re.compile(
        r"'([^']+)':\s*begin\s*\n\s*b1=\s*([^&]+)&\s*b2=\s*([^&]+)&\s*t1=\s*([^&]+)&\s*t2=\s*([^\n]+)",
        re.IGNORECASE,
    )
    for match in pattern.finditer(text):
        result[match.group(1)] = tuple(float(match.group(index)) for index in range(2, 6))
    return result


def read_channels(
    repository_root: str | Path,
    instrument: InstrumentConfig,
    channel_ids: tuple[int, ...],
    *,
    srf_quad: int,
) -> tuple[Channel, ...]:
    """Read current SRFs and reproduce ``load_srfstrarr`` channel metadata."""

    root = Path(repository_root).resolve()
    solar = read_solar_spectrum(root / "create_orac_lut" / "input_files" / "sun" / "Gueymard2018.sssi")
    constants = _blackbody_constants(root)
    result = []
    for channel_id in channel_ids:
        filename = instrument.srf_files[channel_id]
        table = read_srf(root / "create_orac_lut" / "input_files" / "srf" / filename)
        wavenumber = np.asarray(table.coordinate, dtype=np.float32)
        response = np.asarray(table.value, dtype=np.float32)
        wavelength = np.float32(1.0e4) / wavenumber
        centre_wavelength = _legacy_integral(wavelength, response * wavelength) / _legacy_integral(wavelength, response)
        centre_wavenumber = _legacy_integral(wavenumber, response * wavenumber) / _legacy_integral(wavenumber, response)
        if srf_quad not in (0, 1, 2):
            raise ValueError("srf_quad must be 0, 1, or 2")
        if srf_quad != 1:
            raise NotImplementedError("Python reproduction currently requires srf_quad=1")
        b1, b2, t1, t2 = constants.get(filename, (0.0, 0.0, 0.0, 0.0))
        result.append(Channel(
            channel_id=channel_id,
            srf_filename=filename,
            wavelength_microns=centre_wavelength,
            wavenumber_cm_inverse=centre_wavenumber,
            solar=channel_id in instrument.solar_channels,
            mixed=channel_id in instrument.solar_channels and channel_id in instrument.thermal_channels,
            thermal=channel_id in instrument.thermal_channels,
            f0=_solar_f0(wavenumber, response, centre_wavelength, solar.coordinate, solar.value),
            b1=np.float32(b1), b2=np.float32(b2), t1=np.float32(t1), t2=np.float32(t2),
        ))
    return tuple(result)


def _interpolate_profile(profile_x: np.ndarray, profile_y: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Reproduce IDL ``INTERPOL(V, X, XOUT)`` as used for the particle profile.

    The legacy generators interpolate the ``.mm`` relative-amount profile onto
    layer mid-heights with plain ``INTERPOL``: piecewise linear inside the
    tabulated range and, because the segment index is clamped to the end
    segments, **linear extrapolation from the two nearest nodes** below the
    first node and above the last one.  Ascending or descending abscissae are
    accepted.  Arithmetic is single precision, as in the legacy run.
    """

    x = np.asarray(profile_x, dtype=np.float32)
    y = np.asarray(profile_y, dtype=np.float32)
    q = np.asarray(query, dtype=np.float32)
    if x.ndim != 1 or x.size < 2 or y.shape != x.shape:
        raise ValueError("Profile interpolation requires at least two (height, amount) nodes")
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError("Profile heights must be strictly monotonic")
    segment = np.clip(np.searchsorted(x, q, side="right") - 1, 0, x.size - 2)
    x0, x1 = x[segment], x[segment + 1]
    y0, y1 = y[segment], y[segment + 1]
    return (y0 + (q - x0) * (y1 - y0) / (x1 - x0)).astype(np.float32)


def _layer_inputs(
    atmosphere: AtmosphereProfile,
    model_height: np.ndarray,
    model_relative_amount: np.ndarray,
    *,
    no_rayleigh: bool,
    wavelength_microns: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_layers = atmosphere.height_km.size - 1
    layer_height = (atmosphere.height_km[:n_layers] + atmosphere.height_km[1:]) / np.float32(2.0)
    relative_tau = _interpolate_profile(model_height, model_relative_amount, layer_height)
    relative_tau /= np.sum(relative_tau, dtype=np.float32)
    if no_rayleigh:
        column_tau = np.float32(1.0e-6)
    else:
        column_tau = np.float32(atmosphere.pressure_hpa[-1] / np.float32(1013.0)) / (
            np.float32(117.03) * np.float32(wavelength_microns) ** np.float32(4.0)
            - np.float32(1.316) * np.float32(wavelength_microns) ** np.float32(2.0)
        )
    rayleigh_level = column_tau * np.exp(
        np.float32(-0.1188) * atmosphere.height_km
        - np.float32(0.00116) * atmosphere.height_km ** np.float32(2.0)
    )
    rayleigh_layer = rayleigh_level[1:] - rayleigh_level[:-1]
    return relative_tau, rayleigh_layer.astype(np.float32), layer_height.astype(np.float32)


def _view_cosines(satellite_zenith: np.ndarray) -> np.ndarray:
    values = np.asarray(satellite_zenith, dtype=np.float32).copy()
    values[values == np.float32(90.0)] = np.float32(89.99)
    down = -np.cos(np.deg2rad(values)).astype(np.float32)
    down[np.abs(down) == np.float32(1.0)] *= np.float32(0.99999)
    return np.concatenate((down, -down[::-1])).astype(np.float32)


def _operator_arrays(grid: LutGrid, channels: tuple[Channel, ...]) -> dict[str, np.ndarray]:
    nc, nr, no, nz, ns, na = len(channels), grid.effective_radius.size, grid.optical_depth.size, grid.solar_zenith.size, grid.satellite_zenith.size, grid.relative_azimuth.size
    return {
        "rfd": np.zeros((nc, nr, no), dtype=np.float32),
        "tfd": np.zeros((nc, nr, no), dtype=np.float32),
        "rd": np.zeros((nc, nr, no, ns), dtype=np.float32),
        "td": np.zeros((nc, nr, no, ns), dtype=np.float32),
        "tb": np.zeros((nc, nr, no, nz), dtype=np.float32),
        "rfbd": np.zeros((nc, nr, no, nz), dtype=np.float32),
        "tfbd": np.zeros((nc, nr, no, nz), dtype=np.float32),
        "rbd": np.zeros((nc, nr, no, nz, ns, na), dtype=np.float32),
        "em": np.zeros((nc, nr, no, ns), dtype=np.float32),
    }


def _pressure_operator_arrays(grid: LutGrid, channels: tuple[Channel, ...]) -> dict[str, np.ndarray]:
    if grid.surface_pressure is None:
        raise ValueError("The aerosol formulation requires a surface-pressure LUT grid")
    nc, npres, nr, no, nz, ns, na = (
        len(channels), grid.surface_pressure.size, grid.effective_radius.size,
        grid.optical_depth.size, grid.solar_zenith.size,
        grid.satellite_zenith.size, grid.relative_azimuth.size,
    )
    return {
        "rfd": np.zeros((nc, npres, nr, no), dtype=np.float32),
        "tfd": np.zeros((nc, npres, nr, no), dtype=np.float32),
        "rd": np.zeros((nc, npres, nr, no, ns), dtype=np.float32),
        "td": np.zeros((nc, npres, nr, no, ns), dtype=np.float32),
        "tb": np.zeros((nc, npres, nr, no, nz), dtype=np.float32),
        "rfbd": np.zeros((nc, npres, nr, no, nz), dtype=np.float32),
        "tfbd": np.zeros((nc, npres, nr, no, nz), dtype=np.float32),
        "rbd": np.zeros((nc, npres, nr, no, nz, ns, na), dtype=np.float32),
        "em": np.zeros((nc, npres, nr, no, ns), dtype=np.float32),
    }


def _gas_levels(root: Path, instrument: InstrumentConfig, atmosphere_code: int, channel_id: int) -> np.ndarray:
    filename = (
        f"ModtranGasOpd_A{atmosphere_code}_{instrument.platform}_"
        f"{instrument.instrument}_ch{channel_id:02d}.gas"
    )
    path = root / "create_orac_lut" / "input_files" / "gas" / filename
    if not path.is_file():
        raise FileNotFoundError(f"Gas profile for channel {channel_id} does not exist: {path}")
    lines = path.read_text().splitlines()
    nlevels = None
    data_start = None
    for index, line in enumerate(lines):
        if line.strip().lower() == "*nlevels":
            nlevels = int(lines[index + 1].strip())
            data_start = index + 2
            break
    if nlevels is None or data_start is None:
        raise ValueError(f"Gas profile lacks an nlevels section: {path}")
    rows = []
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
        if len(rows) == nlevels:
            break
    if len(rows) != nlevels:
        raise ValueError(f"Gas profile contains {len(rows)} levels, expected {nlevels}: {path}")
    return np.asarray(rows, dtype=np.float32)[:, 1]


def calculate_cloud_operators(
    grid: LutGrid,
    atmosphere: AtmosphereProfile,
    model_height: np.ndarray,
    model_relative_amount: np.ndarray,
    channels: tuple[Channel, ...],
    extinction_ratio: np.ndarray,
    single_scatter_albedo: np.ndarray,
    asymmetry: np.ndarray,
    phase_moments: np.ndarray,
    *,
    rayleigh: bool = True,
    nstreams: int = 60,
) -> dict[str, np.ndarray]:
    """Port the current cloud ``setup_disort``/``call_disort`` loops."""

    arrays = _operator_arrays(grid, channels)
    n_moments = phase_moments.shape[0]
    molecular_moments = getmom(2, 0.0, n_moments - 1)
    umu = _view_cosines(grid.satellite_zenith)
    phi = np.asarray(grid.relative_azimuth, dtype=np.float32)
    diffuse_umu0 = np.float32(np.cos(np.deg2rad(np.float32(50.0))))
    for channel_index, channel in enumerate(channels):
        relative_tau, tau_rayleigh, _ = _layer_inputs(
            atmosphere, model_height, model_relative_amount,
            no_rayleigh=not rayleigh, wavelength_microns=float(channel.wavelength_microns),
        )
        for radius_index in range(grid.effective_radius.size):
            for optical_depth_index, optical_depth in enumerate(grid.optical_depth):
                tauscat = np.asarray(
                    np.float32(optical_depth) * relative_tau * extinction_ratio[radius_index, channel_index],
                    dtype=np.float32,
                )
                dtau = (tau_rayleigh + tauscat).astype(np.float32)
                total_tau = np.sum(dtau, dtype=np.float32)
                denominator = dtau.copy()
                ssa = np.divide(
                    tau_rayleigh + single_scatter_albedo[radius_index, channel_index] * tauscat,
                    denominator,
                    out=np.zeros_like(denominator), where=denominator != 0,
                ).astype(np.float32)
                ssa[ssa > 1.0] = 1.0
                aerosol_present = tauscat > 0.0
                pmom = np.empty((n_moments, dtau.size), dtype=np.float32, order="F")
                for layer in range(dtau.size):
                    if not aerosol_present[layer]:
                        pmom[:, layer] = molecular_moments
                    else:
                        pmom[:, layer] = (
                            molecular_moments * tau_rayleigh[layer]
                            + phase_moments[:, radius_index, channel_index]
                            * single_scatter_albedo[radius_index, channel_index]
                            * tauscat[layer]
                        ) / (
                            tau_rayleigh[layer]
                            + single_scatter_albedo[radius_index, channel_index] * tauscat[layer]
                        )
                pmom[pmom > 1.0] = 1.0
                utau = np.asarray([0.0, total_tau], dtype=np.float32)
                diffuse = call_disort(
                    dtau, ssa, pmom, utau, umu, phi, 0.0, diffuse_umu0, 100.0,
                    nstreams=nstreams,
                )
                arrays["rfd"][channel_index, radius_index, optical_depth_index] = diffuse["flup"][0] / np.float32(np.pi)
                arrays["tfd"][channel_index, radius_index, optical_depth_index] = diffuse["rfldn"][1] / np.float32(np.pi)
                up_indices = slice(2 * grid.satellite_zenith.size - 1, grid.satellite_zenith.size - 1, -1)
                arrays["rd"][channel_index, radius_index, optical_depth_index, :] = diffuse["uu"][up_indices, 0, 0]
                down = diffuse["uu"][:grid.satellite_zenith.size, 1, 0]
                if channel.solar:
                    down = down - np.float32(100.0) * np.exp(total_tau / umu[:grid.satellite_zenith.size])
                arrays["td"][channel_index, radius_index, optical_depth_index, :] = down
                if channel.thermal:
                    layers = np.flatnonzero(tauscat > 0.0)
                    if layers.size:
                        em_tau = dtau[layers]
                        em_ssa = ssa[layers]
                        em_pmom = np.asfortranarray(pmom[:, layers])
                        em_utau = np.asarray([0.0, np.sum(em_tau, dtype=np.float32)], dtype=np.float32)
                        wavenumber = np.float32(1.0e4) / channel.wavelength_microns
                        wnlo = np.float32(0.995) * wavenumber
                        wnhi = np.float32(1.005) * wavenumber
                        emission = call_disort(
                            em_tau, em_ssa, em_pmom, em_utau, umu, phi,
                            0.0, diffuse_umu0, 0.0, nstreams=nstreams, plank=True,
                            wavenumber_low=float(wnlo), wavenumber_high=float(wnhi),
                            temperature=np.full(layers.size + 1, 250.0, dtype=np.float32),
                        )
                        bbe = plkavg(float(wnlo), float(wnhi), 250.0)
                        arrays["em"][channel_index, radius_index, optical_depth_index, :] = emission["uu"][up_indices, 0, 0] / np.float32(bbe)
                if channel.solar:
                    for solar_zenith_index, solar_zenith in enumerate(grid.solar_zenith):
                        solar_angle = np.float32(89.99 if solar_zenith == 90.0 else solar_zenith)
                        umu0 = np.float32(np.cos(np.deg2rad(solar_angle)))
                        direct = call_disort(
                            dtau, ssa, pmom, utau, umu, phi, 100.0, umu0, 0.0,
                            nstreams=nstreams,
                        )
                        arrays["tb"][channel_index, radius_index, optical_depth_index, solar_zenith_index] = np.float32(100.0) * direct["rfldir"][1] / direct["rfldir"][0]
                        arrays["rfbd"][channel_index, radius_index, optical_depth_index, solar_zenith_index] = np.float32(100.0) * direct["flup"][0] / direct["rfldir"][0]
                        arrays["tfbd"][channel_index, radius_index, optical_depth_index, solar_zenith_index] = np.float32(100.0) * direct["rfldn"][1] / direct["rfldir"][0]
                        up = direct["uu"][up_indices, 0, :] * np.float32(np.pi)
                        for azimuth_index in range(grid.relative_azimuth.size):
                            arrays["rbd"][channel_index, radius_index, optical_depth_index, solar_zenith_index, :, grid.relative_azimuth.size - azimuth_index - 1] = up[:, azimuth_index]
        # srf_quad=1 has a unit weight, but retain the normalization operation.
        weight = np.float32(1.0)
        for name in arrays:
            if arrays[name].shape[0] == len(channels):
                arrays[name][channel_index] /= weight
    return arrays


def calculate_aerosol_operators(
    grid: LutGrid,
    atmosphere: AtmosphereProfile,
    model_height: np.ndarray,
    model_relative_amount: np.ndarray,
    channels: tuple[Channel, ...],
    extinction_ratio: np.ndarray,
    single_scatter_albedo: np.ndarray,
    asymmetry: np.ndarray,
    phase_moments: np.ndarray,
    *,
    rayleigh: bool = True,
    gas_levels: dict[int, np.ndarray] | None = None,
    nstreams: int = 60,
) -> dict[str, np.ndarray]:
    """Port the pressure-aware distributed-profile aerosol formulation.

    This follows ``create_orac_aerosol_lut.pro``: pressure is an output
    dimension, gas and Rayleigh optical depths are added per atmospheric layer,
    and particle optical depth is distributed using the model profile.
    """

    arrays = _pressure_operator_arrays(grid, channels)
    n_moments = phase_moments.shape[0]
    molecular_moments = getmom(2, 0.0, n_moments - 1)
    umu = _view_cosines(grid.satellite_zenith)
    phi = np.asarray(grid.relative_azimuth, dtype=np.float32)
    diffuse_umu0 = np.float32(np.cos(np.deg2rad(np.float32(50.0))))
    n_layers = atmosphere.height_km.size - 1
    layer_height = (atmosphere.height_km[:n_layers] + atmosphere.height_km[1:]) / np.float32(2.0)
    relative_tau = _interpolate_profile(model_height, model_relative_amount, layer_height)
    relative_tau /= np.sum(relative_tau, dtype=np.float32)
    up_indices = slice(2 * grid.satellite_zenith.size - 1, grid.satellite_zenith.size - 1, -1)
    for channel_index, channel in enumerate(channels):
        if gas_levels is None:
            gas = np.zeros(atmosphere.height_km.size, dtype=np.float32)
        else:
            gas = np.asarray(gas_levels.get(channel.channel_id, np.zeros(atmosphere.height_km.size)), dtype=np.float32)
            if gas.size != atmosphere.height_km.size:
                raise ValueError(f"Gas profile for channel {channel.channel_id} has the wrong number of levels")
        tau_gas = gas[1:] - gas[:-1]
        for pressure_index, pressure in enumerate(grid.surface_pressure):
            if rayleigh:
                column_tau = np.float32(atmosphere.pressure_hpa[-1] / np.float32(pressure)) / (
                    np.float32(117.03) * channel.wavelength_microns ** np.float32(4.0)
                    - np.float32(1.316) * channel.wavelength_microns ** np.float32(2.0)
                )
            else:
                column_tau = np.float32(1.0e-6)
            rayleigh_level = column_tau * np.exp(
                np.float32(-0.1188) * atmosphere.height_km
                - np.float32(0.00116) * atmosphere.height_km ** np.float32(2.0)
            )
            tau_rayleigh = (rayleigh_level[1:] - rayleigh_level[:-1]).astype(np.float32)
            for radius_index in range(grid.effective_radius.size):
                for optical_depth_index, optical_depth in enumerate(grid.optical_depth):
                    tauscat = np.asarray(
                        np.float32(optical_depth) * relative_tau * extinction_ratio[radius_index, channel_index],
                        dtype=np.float32,
                    )
                    dtau = (tau_gas + tau_rayleigh + tauscat).astype(np.float32)
                    total_tau = np.sum(dtau, dtype=np.float32)
                    denominator = dtau.copy()
                    ssa = np.divide(
                        tau_rayleigh + single_scatter_albedo[radius_index, channel_index] * tauscat,
                        denominator, out=np.zeros_like(denominator), where=denominator != 0,
                    ).astype(np.float32)
                    ssa[ssa > 1.0] = 1.0
                    pmom = np.empty((n_moments, dtau.size), dtype=np.float32, order="F")
                    for layer in range(dtau.size):
                        if tauscat[layer] <= 0.0:
                            pmom[:, layer] = molecular_moments
                            continue
                        scattering_tau = tau_rayleigh[layer] + single_scatter_albedo[radius_index, channel_index] * tauscat[layer]
                        if scattering_tau == 0.0:
                            pmom[:, layer] = molecular_moments
                        else:
                            pmom[:, layer] = (
                                molecular_moments * tau_rayleigh[layer]
                                + phase_moments[:, radius_index, channel_index]
                                * single_scatter_albedo[radius_index, channel_index] * tauscat[layer]
                            ) / scattering_tau
                    pmom[pmom > 1.0] = 1.0
                    utau = np.asarray([0.0, total_tau], dtype=np.float32)
                    diffuse = call_disort(
                        dtau, ssa, pmom, utau, umu, phi, 0.0, diffuse_umu0, 100.0,
                        nstreams=nstreams,
                    )
                    arrays["rfd"][channel_index, pressure_index, radius_index, optical_depth_index] = diffuse["flup"][0] / np.float32(np.pi)
                    arrays["tfd"][channel_index, pressure_index, radius_index, optical_depth_index] = diffuse["rfldn"][1] / np.float32(np.pi)
                    arrays["rd"][channel_index, pressure_index, radius_index, optical_depth_index, :] = diffuse["uu"][up_indices, 0, 0]
                    down = diffuse["uu"][:grid.satellite_zenith.size, 1, 0]
                    if channel.solar:
                        down = down - np.float32(100.0) * np.exp(total_tau / umu[:grid.satellite_zenith.size])
                    arrays["td"][channel_index, pressure_index, radius_index, optical_depth_index, :] = down
                    if channel.thermal:
                        layers = np.flatnonzero(tauscat > 0.0)
                        if layers.size:
                            em_tau = dtau[layers]
                            em_ssa = ssa[layers]
                            em_pmom = np.asfortranarray(pmom[:, layers])
                            em_utau = np.asarray([0.0, np.sum(em_tau, dtype=np.float32)], dtype=np.float32)
                            wavenumber = np.float32(1.0e4) / channel.wavelength_microns
                            wnlo = np.float32(0.995) * wavenumber
                            wnhi = np.float32(1.005) * wavenumber
                            emission = call_disort(
                                em_tau, em_ssa, em_pmom, em_utau, umu, phi,
                                0.0, diffuse_umu0, 0.0, nstreams=nstreams, plank=True,
                                wavenumber_low=float(wnlo), wavenumber_high=float(wnhi),
                                temperature=np.full(layers.size + 1, 270.0, dtype=np.float32),
                            )
                            bbe = plkavg(float(wnlo), float(wnhi), 270.0)
                            arrays["em"][channel_index, pressure_index, radius_index, optical_depth_index, :] = emission["uu"][up_indices, 0, 0] / np.float32(bbe)
                    if channel.solar:
                        for solar_zenith_index, solar_zenith in enumerate(grid.solar_zenith):
                            solar_angle = np.float32(89.99 if solar_zenith == 90.0 else solar_zenith)
                            umu0 = np.float32(np.cos(np.deg2rad(solar_angle)))
                            direct = call_disort(
                                dtau, ssa, pmom, utau, umu, phi, 100.0, umu0, 0.0,
                                nstreams=nstreams,
                            )
                            arrays["tb"][channel_index, pressure_index, radius_index, optical_depth_index, solar_zenith_index] = np.float32(100.0) * direct["rfldir"][1] / direct["rfldir"][0]
                            arrays["rfbd"][channel_index, pressure_index, radius_index, optical_depth_index, solar_zenith_index] = np.float32(100.0) * direct["flup"][0] / direct["rfldir"][0]
                            arrays["tfbd"][channel_index, pressure_index, radius_index, optical_depth_index, solar_zenith_index] = np.float32(100.0) * direct["rfldn"][1] / direct["rfldir"][0]
                            up = direct["uu"][up_indices, 0, :] * np.float32(np.pi)
                            for azimuth_index in range(grid.relative_azimuth.size):
                                arrays["rbd"][channel_index, pressure_index, radius_index, optical_depth_index, solar_zenith_index, :, grid.relative_azimuth.size - azimuth_index - 1] = up[:, azimuth_index]
        weight = np.float32(1.0)
        for name in arrays:
            arrays[name][channel_index] /= weight
    return arrays


def _chars(text: str, length: int) -> np.ndarray:
    encoded = text.encode("ascii", errors="replace")[:length]
    values = np.full(length, b" ", dtype="S1")
    values[:len(encoded)] = np.frombuffer(encoded, dtype="S1")
    return values


def _metadata_variables(
    root: Path,
    instrument: InstrumentConfig,
    channels: tuple[Channel, ...],
    grid: LutGrid,
    optics,
    operators: dict[str, np.ndarray],
    *,
    instrument_filename: str = "meteosat-10_seviri_v1.inst",
    surface_pressure: np.ndarray | None = None,
) -> GenerationResult:
    channel_ids = [channel.channel_id for channel in channels]
    solar_channels = [channel for channel in channels if channel.solar]
    thermal_channels = [channel for channel in channels if channel.thermal]
    mixed_channels = [channel for channel in channels if channel.mixed]
    nr, no, ns, nz, na, nc = (
        grid.effective_radius.size, grid.optical_depth.size, grid.satellite_zenith.size,
        grid.solar_zenith.size, grid.relative_azimuth.size, len(channels),
    )
    dimensions = {
        "optical_depth": no, "effective_radius": nr, "satellite_zenith": ns,
        "solar_zenith": nz, "relative_azimuth": na, "channels": nc,
        "length": 32, "st1": 26,
        "st2": len(instrument.platform), "st3": len(instrument.instrument), "st4": len(instrument.version),
    }
    if surface_pressure is not None:
        dimensions["surface_pressure"] = int(np.asarray(surface_pressure).size)
    if solar_channels:
        dimensions["solar_channels"] = len(solar_channels)
    if thermal_channels:
        dimensions["thermal_channels"] = len(thermal_channels)
    if mixed_channels:
        dimensions["mixed_channels"] = len(mixed_channels)
    variables: dict[str, np.ndarray] = {
        "instrument_filename": _chars(instrument_filename, dimensions["st1"]),
        "platform": _chars(instrument.platform, dimensions["st2"]),
        "instrument": _chars(instrument.instrument, dimensions["st3"]),
        "instrument_version": _chars(instrument.version, dimensions["st4"]),
        "max_sat_zenith": np.float32(instrument.maximum_satellite_zenith),
        "number_of_channels": np.int32(nc),
        "SRF_file": np.stack([_chars(channel.srf_filename, dimensions["length"]) for channel in channels]),
        "channel_id": np.asarray(channel_ids, dtype=np.int16),
        "central_wavelength": np.asarray([channel.wavelength_microns for channel in channels], dtype=np.float32),
        "central_wavenumber": np.asarray([channel.wavenumber_cm_inverse for channel in channels], dtype=np.float32),
        "solar_channel_flag": np.asarray([channel.solar for channel in channels], dtype=np.int16),
        "mixed_channel_flag": np.asarray([channel.mixed for channel in channels], dtype=np.int16),
        "thermal_channel_flag": np.asarray([channel.thermal for channel in channels], dtype=np.int16),
        "average_volume_per_particle": np.asarray(optics.average_volume_per_particle, dtype=np.float32),
        "extinction_coefficient": np.asarray(optics.extinction_coefficient[:, 1:], dtype=np.float32),
        "extinction_coefficient_ratio": np.asarray(operators["extinction_ratio"], dtype=np.float32),
        "single_scatter_albedo": np.asarray(optics.single_scatter_albedo[:, 1:], dtype=np.float32),
        "asymmetry_parameter": np.asarray(optics.asymmetry_parameter[:, 1:], dtype=np.float32),
        "optical_depth": np.asarray(grid.optical_depth, dtype=np.float32),
        "effective_radius": np.asarray(grid.effective_radius, dtype=np.float32),
        "satellite_zenith": np.asarray(grid.satellite_zenith, dtype=np.float32),
        "solar_zenith": np.asarray(grid.solar_zenith, dtype=np.float32),
        "relative_azimuth": np.asarray(grid.relative_azimuth, dtype=np.float32),
    }
    if surface_pressure is not None:
        variables["surface_pressure"] = np.asarray(surface_pressure, dtype=np.float32)
    if solar_channels:
        oldnefr = np.full(nc, np.float32(9.96921e36), dtype=np.float32)
        oldnefr[:len(solar_channels)] = np.asarray(
            [instrument.oldnefr.get(c.channel_id, 0.0) for c in solar_channels],
            dtype=np.float32,
        )
        variables.update({
            "solar_channel_id": np.asarray([c.channel_id for c in solar_channels], dtype=np.int16),
            "oldf0": np.asarray([instrument.oldf0.get(c.channel_id, 0.0) for c in solar_channels], dtype=np.float32),
            "oldf1": np.asarray([instrument.oldf1.get(c.channel_id, 0.0) for c in solar_channels], dtype=np.float32),
            "F0": np.asarray([c.f0 for c in solar_channels], dtype=np.float32),
            "snr": np.asarray([instrument.snr.get(c.channel_id, 0.0) for c in solar_channels], dtype=np.float32),
            "oldnefr": oldnefr,
        })
    if thermal_channels:
        def legacy_thermal_vector(values: dict[int, float]) -> np.ndarray:
            result = np.full(nc, np.float32(9.96921e36), dtype=np.float32)
            result[:len(thermal_channels)] = np.asarray(
                [values.get(c.channel_id, 0.0) for c in thermal_channels],
                dtype=np.float32,
            )
            return result

        variables.update({
            "thermal_channel_id": np.asarray([c.channel_id for c in thermal_channels], dtype=np.int16),
            "refbt": np.asarray([instrument.refbt.get(c.channel_id, 0.0) for c in thermal_channels], dtype=np.float32),
            "nedt": np.asarray([instrument.nedt.get(c.channel_id, 0.0) for c in thermal_channels], dtype=np.float32),
            "B1": np.asarray([c.b1 for c in thermal_channels], dtype=np.float32),
            "B2": np.asarray([c.b2 for c in thermal_channels], dtype=np.float32),
            "T1": np.asarray([c.t1 for c in thermal_channels], dtype=np.float32),
            "T2": np.asarray([c.t2 for c in thermal_channels], dtype=np.float32),
            "oldwvn": legacy_thermal_vector(instrument.oldwvn),
            "oldb1": legacy_thermal_vector(instrument.oldb1),
            "oldb2": legacy_thermal_vector(instrument.oldb2),
            "oldt1": legacy_thermal_vector(instrument.oldt1),
            "oldt2": legacy_thermal_vector(instrument.oldt2),
            "oldnebt": legacy_thermal_vector(instrument.oldnebt),
        })
    if mixed_channels:
        variables["mixed_channel_id"] = np.asarray([c.channel_id for c in mixed_channels], dtype=np.int16)
    solar_dim = "solar_channels"
    thermal_dim = "thermal_channels"
    variable_dimensions: dict[str, tuple[str, ...]] = {
        "instrument_filename": ("st1",), "platform": ("st2",), "instrument": ("st3",),
        "instrument_version": ("st4",), "max_sat_zenith": (), "number_of_channels": (),
        "SRF_file": ("channels", "length"), "channel_id": ("channels",),
        "central_wavelength": ("channels",), "central_wavenumber": ("channels",),
        "solar_channel_flag": ("channels",), "mixed_channel_flag": ("channels",),
        "thermal_channel_flag": ("channels",), "average_volume_per_particle": ("effective_radius",),
        "extinction_coefficient": ("effective_radius", "channels"),
        "extinction_coefficient_ratio": ("effective_radius", "channels"),
        "single_scatter_albedo": ("effective_radius", "channels"),
        "asymmetry_parameter": ("effective_radius", "channels"),
        "optical_depth": ("optical_depth",), "effective_radius": ("effective_radius",),
        "satellite_zenith": ("satellite_zenith",), "solar_zenith": ("solar_zenith",),
        "relative_azimuth": ("relative_azimuth",),
    }
    if surface_pressure is not None:
        variable_dimensions["surface_pressure"] = ("surface_pressure",)
    if solar_channels:
        for name in ("solar_channel_id", "oldf0", "oldf1", "F0", "snr", "oldnefr"):
            variable_dimensions[name] = (solar_dim,) if name != "oldnefr" else ("channels",)
    if thermal_channels:
        for name in ("thermal_channel_id", "refbt", "nedt", "B1", "B2", "T1", "T2"):
            variable_dimensions[name] = (thermal_dim,)
        for name in ("oldwvn", "oldb1", "oldb2", "oldt1", "oldt2", "oldnebt"):
            variable_dimensions[name] = ("channels",)
    if mixed_channels:
        variable_dimensions["mixed_channel_id"] = ("mixed_channels",)
    pressure_dim = ("surface_pressure",) if surface_pressure is not None else ()
    variable_dimensions.update({
        "T_dv": ("satellite_zenith", "optical_depth", "effective_radius", *pressure_dim, "channels"),
        "T_dd": ("optical_depth", "effective_radius", *pressure_dim, "channels"),
        "R_dv": ("satellite_zenith", "optical_depth", "effective_radius", *pressure_dim, "channels"),
        "R_dd": ("optical_depth", "effective_radius", *pressure_dim, "channels"),
    })
    if solar_channels:
        variable_dimensions.update({
            "R_0v": ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth", "effective_radius", *pressure_dim, solar_dim),
            "R_0d": ("solar_zenith", "optical_depth", "effective_radius", *pressure_dim, solar_dim),
            "T_0d": ("solar_zenith", "optical_depth", "effective_radius", *pressure_dim, solar_dim),
            "T_00": ("solar_zenith", "optical_depth", "effective_radius", *pressure_dim, solar_dim),
        })
    if thermal_channels:
        variable_dimensions["E_md"] = ("satellite_zenith", "optical_depth", "effective_radius", *pressure_dim, thermal_dim)
    scale = np.float32(0.01)
    solar_indices = [index for index, channel in enumerate(channels) if channel.solar]
    thermal_indices = [index for index, channel in enumerate(channels) if channel.thermal]
    if surface_pressure is None:
        td_values = np.transpose(operators["td"], (3, 2, 1, 0))
        tfd_values = np.transpose(operators["tfd"], (2, 1, 0))
        rd_values = np.transpose(operators["rd"], (3, 2, 1, 0))
        rfd_values = np.transpose(operators["rfd"], (2, 1, 0))
    else:
        td_values = np.transpose(operators["td"], (4, 3, 2, 1, 0))
        tfd_values = np.transpose(operators["tfd"], (3, 2, 1, 0))
        rd_values = np.transpose(operators["rd"], (4, 3, 2, 1, 0))
        rfd_values = np.transpose(operators["rfd"], (3, 2, 1, 0))
    variables.update({
        "T_dv": scale * td_values,
        "T_dd": scale * tfd_values,
        "R_dv": scale * rd_values,
        "R_dd": scale * rfd_values,
    })
    if solar_channels:
        rbd = operators["rbd"][solar_indices]
        rfbd = operators["rfbd"][solar_indices]
        tfbd = operators["tfbd"][solar_indices]
        tb = operators["tb"][solar_indices]
        if surface_pressure is None:
            rbd_values = np.transpose(rbd, (5, 4, 3, 2, 1, 0))
            rfbd_values = np.transpose(rfbd, (3, 2, 1, 0))
            tfbd_values = np.transpose(tfbd, (3, 2, 1, 0))
            tb_values = np.transpose(tb, (3, 2, 1, 0))
        else:
            rbd_values = np.transpose(rbd, (6, 5, 4, 3, 2, 1, 0))
            rfbd_values = np.transpose(rfbd, (4, 3, 2, 1, 0))
            tfbd_values = np.transpose(tfbd, (4, 3, 2, 1, 0))
            tb_values = np.transpose(tb, (4, 3, 2, 1, 0))
        variables.update({
            "R_0v": scale * rbd_values,
            "R_0d": scale * rfbd_values,
            "T_0d": scale * tfbd_values,
            "T_00": scale * tb_values,
        })
    if thermal_channels:
        # Both operator calculators store UU/BBE for emission, i.e. the
        # emissivity as a fraction in 0-1, which is the convention of the legacy
        # V2 product (the legacy writer divides its percentage Em by 100).
        # The 0.01 scale applied to the other operators must not be applied here.
        em = operators["em"][thermal_indices]
        variables["E_md"] = np.transpose(em, (4, 3, 2, 1, 0) if surface_pressure is not None else (3, 2, 1, 0))
    reference_variable_order = (
        "instrument_filename", "platform", "instrument", "instrument_version",
        "max_sat_zenith", "number_of_channels", "SRF_file", "channel_id",
        "central_wavelength", "central_wavenumber", "solar_channel_flag",
        "mixed_channel_flag", "thermal_channel_flag", "solar_channel_id",
        "oldf0", "oldf1", "F0", "snr", "oldnefr", "mixed_channel_id",
        "thermal_channel_id", "refbt", "nedt", "B1", "B2", "T1", "T2",
        "oldwvn", "oldb1", "oldb2", "oldt1", "oldt2", "oldnebt",
        "average_volume_per_particle", "extinction_coefficient",
        "extinction_coefficient_ratio", "single_scatter_albedo",
        "asymmetry_parameter", "optical_depth", "effective_radius",
        "satellite_zenith", "solar_zenith", "relative_azimuth", "surface_pressure", "T_dv",
        "T_dd", "R_dv", "R_dd", "R_0v", "R_0d", "T_0d", "T_00", "E_md",
    )
    variables = {name: variables[name] for name in reference_variable_order if name in variables}
    positive_float = np.asarray([0.0, np.finfo(np.float32).max], dtype=np.float32)
    channel_range = np.asarray([0, nc], dtype=np.int32)
    attrs = {
        "max_sat_zenith": {
            "units": "degrees", "valid_range": np.asarray([0.0, 90.0], dtype=np.float32),
        },
        "number_of_channels": {"units": "dimensionless", "valid_range": channel_range},
        "SRF_file": {"long_name": "file containing the spectral response for the channel"},
        "channel_id": {
            "long_name": "Instrument channel identifier", "units": "dimensionless",
            "valid_range": channel_range,
        },
        "central_wavelength": {
            "long_name": "effective central wavelength for the channel", "units": "microns",
            "valid_range": positive_float,
        },
        "central_wavenumber": {
            "long_name": "effective central wavenumber for the channel", "units": "cm^{-1}",
            "valid_range": positive_float,
        },
        "solar_channel_flag": {
            "long_name": "Flag set to 1 if channel measures reflected solar radiation otherwise 0",
            "units": "dimensionless", "valid_range": np.asarray([0, 1], dtype=np.int16),
        },
        "mixed_channel_flag": {
            "long_name": "Flag set to 1 if channel measures reflected solar and emitted infrared radiation otherwise 0",
            "units": "dimensionless", "valid_range": np.asarray([0, 1], dtype=np.int16),
        },
        "thermal_channel_flag": {
            "long_name": "Flag set to 1 if channel measures emitted infrared radiation otherwise 0",
            "units": "dimensionless", "valid_range": np.asarray([0, 1], dtype=np.int16),
        },
        "solar_channel_id": {
            "long_name": "Instrument channel identifier for solar channels",
            "units": "dimensionless", "valid_range": channel_range,
        },
        "oldf0": {"long_name": "SAD file F0", "units": "W/(m^2 um)", "valid_range": positive_float},
        "oldf1": {"long_name": "SAD file F1", "units": "W/(m^2 um)", "valid_range": positive_float},
        "F0": {"long_name": "in-band solar radiance", "units": "W/(m^2 um)", "valid_range": positive_float},
        "snr": {"long_name": "signal-to-noise ratio", "units": "dimensionless", "valid_range": positive_float},
        "oldnefr": {"long_name": "SAD file NeFr (measurement uncertainty)", "units": "W/(m^2 sr um)", "valid_range": positive_float},
        "mixed_channel_id": {
            "long_name": "Instrument channel identifier for mixed channels",
            "units": "dimensionless", "valid_range": channel_range,
        },
        "thermal_channel_id": {
            "long_name": "Instrument channel identifier for thermal channels",
            "units": "dimensionless", "valid_range": channel_range,
        },
        "refbt": {
            "long_name": "reference brightness temperate at neDT has been calculated",
            "units": "K", "valid_range": np.asarray([100, 500], dtype=np.int16),
        },
        "nedt": {"long_name": "noise equivalent delta temperature", "units": "K", "valid_range": positive_float},
        "oldwvn": {"long_name": "SAD file Wvn", "units": "1/cm", "valid_range": positive_float},
        "oldb1": {"long_name": "SAD file B1"},
        "oldb2": {"long_name": "SAD file B2"},
        "oldt1": {"long_name": "SAD file T1"},
        "oldt2": {"long_name": "SAD file T2"},
        "oldnebt": {"long_name": "SAD file NeBT", "units": "K", "valid_range": positive_float},
        "average_volume_per_particle": {
            "long_name": "average volume per particle", "units": "to be investigated",
            "valid_range": positive_float,
        },
        "extinction_coefficient": {
            "long_name": "volume extinction coefficient", "units": "to be investigated",
            "valid_range": positive_float,
        },
        "extinction_coefficient_ratio": {
            "long_name": "ratio of volume extinction coefficient to the volume extinction coefficient at 550 nm",
            "units": "dimensionless", "valid_range": positive_float,
        },
        "single_scatter_albedo": {
            "long_name": "single scatter albedo", "units": "dimensionless",
            "valid_range": np.asarray([0.0, 1.0], dtype=np.float32),
        },
        "asymmetry_parameter": {
            "long_name": "asymmetry parameter", "units": "dimensionless",
            "valid_range": np.asarray([-1.0, 1.0], dtype=np.float32),
        },
        "optical_depth": {
            "long_name": "optical depth", "spacing": grid.spacings[0],
            "units": "dimensionless", "valid_range": positive_float,
        },
        "effective_radius": {
            "long_name": "particle effective radius", "spacing": grid.spacings[1],
            "units": "microns", "valid_range": positive_float,
        },
        "satellite_zenith": {
            "long_name": "satellite zenith angle", "spacing": grid.spacings[3],
            "units": "degrees", "valid_range": np.asarray([0.0, 180.0], dtype=np.float32),
        },
        "solar_zenith": {
            "long_name": "solar zenith angle", "spacing": grid.spacings[2],
            "units": "degrees", "valid_range": np.asarray([0.0, 180.0], dtype=np.float32),
        },
        "relative_azimuth": {
            "long_name": "satellite azimuth relative to the Sun", "spacing": grid.spacings[4],
            "units": "degrees", "valid_range": np.asarray([0.0, 180.0], dtype=np.float32),
        },
    }
    operator_attributes = {
        "T_dv": "diffuse transmission of direct light", "T_dd": "diffuse transmission",
        "R_dv": "direct reflection of diffuse light", "R_dd": "diffuse reflection of diffuse light",
        "R_0v": "bi-directional reflectance", "R_0d": "diffuse reflectance of direct beam",
        "T_0d": "diffuse transmission of diffuse light", "T_00": "direct transmission",
        "E_md": "diffuse emissivity",
    }
    for name, long_name in operator_attributes.items():
        attrs[name] = {
            "long_name": long_name, "units": "dimensionless",
            "valid_range": np.asarray([0.0, 1.0], dtype=np.float32),
        }
    return GenerationResult(variables, variable_dimensions, attrs, dimensions, channels)


def _atmosphere_path(input_root: Path, atmosphere_code: int) -> Path:
    names = {
        0: "midsatm.dat", 1: "tro.atm", 2: "mls.atm", 3: "mlw.atm",
        4: "sas.atm", 5: "saw.atm", 6: "std.atm",
    }
    try:
        return input_root / "atm" / names[atmosphere_code]
    except KeyError as exc:
        raise ValueError(f"Unsupported atmosphere selector: {atmosphere_code}") from exc


def _component_refractive_indices(root: Path, model, wavelengths: np.ndarray) -> np.ndarray:
    input_root = root / "create_orac_lut" / "input_files"
    values = []
    for component in model.components:
        ri_wavelength, ri_real, ri_imaginary = read_refractive_index(
            input_root / "ri" / component.refractive_index_file
        )
        real = np.interp(wavelengths, ri_wavelength, ri_real)
        imaginary = np.interp(wavelengths, ri_wavelength, ri_imaginary)
        values.append(real.astype(np.float32) - 1j * imaginary.astype(np.float32))
    return np.asarray(values, dtype=np.complex64)


def generate_cloud(
    repository_root: str | Path,
    *,
    lut_path: str | Path,
    microphysics_path: str | Path,
    channels: tuple[int, ...],
    atmosphere_code: int = 2,
    srf_quad: int = 1,
    rayleigh: bool = True,
    phase_order: int = 1000,
    nstreams: int = 60,
    instrument_path: str | Path | None = None,
) -> GenerationResult:
    root = Path(repository_root).resolve()
    input_root = root / "create_orac_lut" / "input_files"
    instrument_file = Path(instrument_path) if instrument_path is not None else input_root / "inst" / "meteosat-10_seviri_v1.inst"
    instrument = read_instrument(instrument_file)
    grid = read_lut_grid(lut_path)
    model = read_microphysics(microphysics_path)
    atmosphere = read_atmosphere(_atmosphere_path(input_root, atmosphere_code), atmosphere_code)
    channel_defs = read_channels(root, instrument, channels, srf_quad=srf_quad)
    wavelengths = np.asarray([0.55] + [float(channel.wavelength_microns) for channel in channel_defs], dtype=np.float64)
    if len(model.components) == 1 and model.components[0].size_distribution == "modified_gamma":
        ri_wavelength, ri_real, ri_imaginary = read_refractive_index(input_root / "ri" / model.refractive_index_file)
        real = np.interp(wavelengths, ri_wavelength, ri_real).astype(np.float32)
        imaginary = np.interp(wavelengths, ri_wavelength, ri_imaginary).astype(np.float32)
        refractive_indices = real.astype(np.complex64) - 1j * imaginary.astype(np.complex64)
        optics = legacy_stg_optics(
            grid.effective_radius, wavelengths, refractive_indices, phase_order=phase_order,
        )
    else:
        optics = legacy_mixed_mie_optics(
            grid.effective_radius, wavelengths, model.components,
            _component_refractive_indices(root, model, wavelengths), phase_order=phase_order,
        )
    extinction_ratio = (optics.extinction_coefficient[:, 1:] / optics.reference_extinction_coefficient[:, None]).astype(np.float32)
    operators = calculate_cloud_operators(
        grid, atmosphere, model.profile_height_km, model.profile_relative_amount,
        channel_defs, extinction_ratio, optics.single_scatter_albedo[:, 1:].astype(np.float32),
        optics.asymmetry_parameter[:, 1:].astype(np.float32), optics.phase_moments[:, :, 1:].astype(np.float32),
        rayleigh=rayleigh, nstreams=nstreams,
    )
    operators["extinction_ratio"] = extinction_ratio
    return _metadata_variables(
        root, instrument, channel_defs, grid, optics, operators,
        instrument_filename=instrument_file.name,
    )


def generate_aerosol(
    repository_root: str | Path,
    *,
    lut_path: str | Path,
    microphysics_path: str | Path,
    instrument_path: str | Path,
    channels: tuple[int, ...],
    atmosphere_code: int = 2,
    srf_quad: int = 1,
    rayleigh: bool = True,
    gas: bool = False,
    phase_order: int = 1000,
    nstreams: int = 60,
) -> GenerationResult:
    """Generate the current pressure-aware distributed-profile aerosol LUT."""

    root = Path(repository_root).resolve()
    input_root = root / "create_orac_lut" / "input_files"
    instrument_file = Path(instrument_path)
    instrument = read_instrument(instrument_file)
    grid = read_lut_grid(lut_path)
    model = read_microphysics(microphysics_path)
    if grid.surface_pressure is None:
        raise ValueError("Aerosol LUT definitions must contain a surface-pressure grid")
    if not model.profile_height_km.size or not model.profile_relative_amount.size:
        raise ValueError("Aerosol microphysics must define a non-empty vertical profile")
    channel_defs = read_channels(root, instrument, channels, srf_quad=srf_quad)
    wavelengths = np.asarray([0.55] + [float(channel.wavelength_microns) for channel in channel_defs], dtype=np.float64)
    refractive_indices = _component_refractive_indices(root, model, wavelengths)
    optics = legacy_mixed_mie_optics(
        grid.effective_radius, wavelengths, model.components, refractive_indices,
        phase_order=phase_order,
    )
    extinction_ratio = (
        optics.extinction_coefficient[:, 1:] / optics.reference_extinction_coefficient[:, None]
    ).astype(np.float32)
    gas_profiles = None
    if gas:
        gas_profiles = {
            channel.channel_id: _gas_levels(root, instrument, atmosphere_code, channel.channel_id)
            for channel in channel_defs
        }
    atmosphere = read_atmosphere(_atmosphere_path(input_root, atmosphere_code), atmosphere_code)
    operators = calculate_aerosol_operators(
        grid, atmosphere, model.profile_height_km, model.profile_relative_amount,
        channel_defs, extinction_ratio, optics.single_scatter_albedo[:, 1:].astype(np.float32),
        optics.asymmetry_parameter[:, 1:].astype(np.float32), optics.phase_moments[:, :, 1:].astype(np.float32),
        rayleigh=rayleigh, gas_levels=gas_profiles, nstreams=nstreams,
    )
    operators["extinction_ratio"] = extinction_ratio
    return _metadata_variables(
        root, instrument, channel_defs, grid, optics, operators,
        instrument_filename=instrument_file.name, surface_pressure=grid.surface_pressure,
    )


def write_generation(
    path: str | Path,
    result: GenerationResult,
    *,
    lut_level: int = 2,
    revision: int = 21,
    overwrite: bool = False,
) -> None:
    output = Path(path)
    if output.exists() and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite existing LUT {output}; pass overwrite=True explicitly"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    write_v2_lut(
        output, lut_level=lut_level, revision=revision,
        dimensions=result.dimensions, variables=result.variables,
        variable_dimensions=result.variable_dimensions,
        variable_attributes=result.variable_attributes,
    )
