"""Legacy-equivalent particle optical properties.

The current IDL STG path calls the repository's Mie/DLM machinery through
``create_bwgp``. This module defines the downstream scientific contract without
substituting a different Mie implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from statistics import NormalDist

import numpy as np

from .legacy_mie import mie_single_batch


@dataclass(frozen=True)
class OpticalProperties:
    """Particle properties required by the current ORAC RT stage."""

    wavelength_microns: np.ndarray
    effective_radius_microns: np.ndarray
    extinction_coefficient: np.ndarray
    scattering_coefficient: np.ndarray | None
    single_scatter_albedo: np.ndarray
    asymmetry_parameter: np.ndarray
    phase_moments: np.ndarray
    average_volume_per_particle: np.ndarray | None = None
    absorption_coefficient: np.ndarray | None = None
    reference_wavelength_microns: float = 0.55
    reference_extinction_coefficient: np.ndarray | None = None


def _mie_distribution(
    effective_radius_microns: float,
    effective_variance: float,
    wavelength_microns: float,
    refractive_index: complex,
    scattering_cosines: np.ndarray,
    *,
    xres: float = 0.4,
    distribution: str = "modified_gamma",
    size_parameters: tuple[float, ...] | None = None,
) -> tuple[float, float, float, np.ndarray, float]:
    """Port the Mie size integration in ``mie_size_dist_new.pro``."""

    if distribution == "log_normal":
        if size_parameters is None or len(size_parameters) < 2:
            raise ValueError("log_normal requires mode radius and spread")
        # IDL GAUSS_CVF(0.999) is the negative lower-tail deviate used by
        # the legacy routine; retaining that sign makes Rl < Ru.
        tq = -NormalDist().inv_cdf(0.999)
        radius_lower = math.exp(math.log(size_parameters[0]) + tq * math.log(size_parameters[1]))
        radius_upper = math.exp(math.log(size_parameters[0]) - tq * math.log(size_parameters[1]) + math.log(4.0))
    elif distribution == "modified_gamma":
        radius_lower = 0.001
        radius_upper = 100.0
    else:
        raise NotImplementedError(f"Unsupported legacy size distribution: {distribution}")
    wavenumber = 1.0 / float(wavelength_microns)
    npts = max(int(2.0 * np.pi * (radius_upper - radius_lower) * wavenumber / xres), 200)
    abscissa = np.linspace(-1.0, 1.0, npts, dtype=np.float64)
    weights = np.full(npts, 2.0 / (npts - 1), dtype=np.float64)
    weights[[0, -1]] = 1.0 / (npts - 1)
    radius = ((radius_lower + radius_upper) + (radius_upper - radius_lower) * abscissa) / 2.0
    radius_weights = (radius_upper - radius_lower) * weights / 2.0

    if distribution == "modified_gamma":
        alpha = (1.0 - 3.0 * effective_variance) / effective_variance
        beta = 1.0 / (effective_radius_microns * effective_variance)
        normalization = beta ** (-alpha - 1.0) * math.gamma(alpha + 1.0)
        distribution_weight = radius_weights / normalization * radius ** alpha * np.exp(-beta * radius)
    else:
        mode_radius, spread = size_parameters[:2]
        log_spread = math.log(spread)
        distribution_weight = (
            radius_weights * np.exp(-0.5 * (np.log(radius / mode_radius) / log_spread) ** 2)
            / (math.sqrt(2.0) * math.sqrt(math.pi) * radius * log_spread)
        )
    size_parameters = 2.0 * np.pi * radius * wavenumber
    particle = mie_single_batch(size_parameters, refractive_index, scattering_cosines)

    area_weight = distribution_weight * np.pi * radius ** 2
    volume_weight = area_weight * (4.0 / 3.0) * radius
    extinction = float(np.sum(area_weight * particle["qext"]))
    scattering = float(np.sum(area_weight * particle["qsca"]))
    asymmetry = float(np.sum(area_weight * particle["g"] * particle["qsca"]) / scattering)
    phase = np.sum(
        area_weight[:, None] * particle["f11"] * particle["qsca"][:, None], axis=0
    ) / scattering
    return extinction, scattering / extinction, asymmetry, phase, float(np.sum(volume_weight))


def legendre_moments(
    scattering_cosines: np.ndarray,
    quadrature_weights: np.ndarray,
    phase: np.ndarray,
    *,
    minimum_absolute_coefficient: float = 1.0e-9,
) -> np.ndarray:
    """Match the current ``legpexp.pro`` coefficient and stopping convention."""

    qv = np.asarray(scattering_cosines, dtype=np.float64)
    qw = np.asarray(quadrature_weights, dtype=np.float64)
    values = np.asarray(phase, dtype=np.float64)
    if not (qv.shape == qw.shape == values.shape):
        raise ValueError("qv, qw and phase must have equal shapes")
    coefficients = np.zeros(qv.size, dtype=np.float64)
    coefficients[0] = np.sum(values * qw) / 2.0
    if qv.size == 1:
        return coefficients
    coefficients[1] = 3.0 * np.sum(values * qv * qw) / 2.0
    previous = np.ones_like(qv)
    current = qv.copy()
    for order in range(2, qv.size):
        polynomial = ((2.0 * order - 1.0) / order) * qv * current - ((order - 1.0) / order) * previous
        coefficients[order] = (2.0 * order + 1.0) * np.sum(values * polynomial * qw) / 2.0
        if abs(coefficients[order]) < minimum_absolute_coefficient:
            break
        previous, current = current, polynomial
    return coefficients


def legacy_stg_optics(
    effective_radii_microns: np.ndarray,
    wavelengths_microns: np.ndarray,
    refractive_indices: np.ndarray,
    *,
    phase_order: int = 1000,
    size_parameter_resolution: float = 0.4,
) -> OpticalProperties:
    """Calculate the current single-component STG optical-property arrays."""

    radii = np.asarray(effective_radii_microns, dtype=np.float64).ravel()
    wavelengths = np.asarray(wavelengths_microns, dtype=np.float64).ravel()
    indices = np.asarray(refractive_indices, dtype=np.complex128).ravel()
    if wavelengths.size != indices.size:
        raise ValueError("wavelengths_microns and refractive_indices must have equal size")
    if phase_order < 2:
        raise ValueError("phase_order must be at least 2")
    # IDL quadrature('g', NMom) returns Gaussian nodes in ascending order and
    # create_bwgp reverses their cosine convention with qv=-abscissa.
    abscissa, phase_weights = np.polynomial.legendre.leggauss(phase_order)
    scattering_cosines = -abscissa
    extinction = np.empty((radii.size, wavelengths.size), dtype=np.float64)
    ssa = np.empty_like(extinction)
    asymmetry = np.empty_like(extinction)
    moments = np.empty((phase_order, radii.size, wavelengths.size), dtype=np.float64)
    volume = np.empty(radii.size, dtype=np.float64)
    for radius_index, radius in enumerate(radii):
        for wavelength_index, (wavelength, refractive_index) in enumerate(zip(wavelengths, indices)):
            result = _mie_distribution(
                radius, 0.1111111, wavelength, refractive_index,
                scattering_cosines, xres=size_parameter_resolution,
            )
            extinction[radius_index, wavelength_index] = result[0]
            ssa[radius_index, wavelength_index] = result[1]
            asymmetry[radius_index, wavelength_index] = result[2]
            moments[:, radius_index, wavelength_index] = legendre_moments(
                scattering_cosines, phase_weights, result[3]
            ) / (2.0 * np.arange(phase_order) + 1.0)
            if wavelength_index == 0:
                volume[radius_index] = result[4]
    reference_extinction = np.empty(radii.size, dtype=np.float64)
    reference_wavelength = 0.55
    # The caller normally includes the reference wavelength as the first point;
    # if it does not, calculate it explicitly with the interpolated RI supplied
    # by the caller.
    if wavelengths.size and np.isclose(wavelengths[0], reference_wavelength):
        reference_extinction[:] = extinction[:, 0]
    else:
        raise ValueError("wavelength grid must begin with the 0.55 micron reference point")
    return OpticalProperties(
        wavelength_microns=wavelengths,
        effective_radius_microns=radii,
        extinction_coefficient=extinction,
        scattering_coefficient=extinction * ssa,
        single_scatter_albedo=ssa,
        asymmetry_parameter=asymmetry,
        phase_moments=moments,
        average_volume_per_particle=volume,
        absorption_coefficient=extinction * (1.0 - ssa),
        reference_wavelength_microns=reference_wavelength,
        reference_extinction_coefficient=reference_extinction,
    )


def _component_optics(
    mode_radii_microns: np.ndarray,
    wavelengths_microns: np.ndarray,
    refractive_indices: np.ndarray,
    distribution: str,
    size_parameters: tuple[float, ...],
    *,
    phase_order: int,
    size_parameter_resolution: float,
) -> OpticalProperties:
    """Calculate one legacy Mie component over the requested modes."""

    radii = np.asarray(mode_radii_microns, dtype=np.float64).ravel()
    wavelengths = np.asarray(wavelengths_microns, dtype=np.float64).ravel()
    indices = np.asarray(refractive_indices, dtype=np.complex128).ravel()
    abscissa, phase_weights = np.polynomial.legendre.leggauss(phase_order)
    scattering_cosines = -abscissa
    extinction = np.empty((radii.size, wavelengths.size), dtype=np.float64)
    ssa = np.empty_like(extinction)
    asymmetry = np.empty_like(extinction)
    moments = np.empty((phase_order, radii.size, wavelengths.size), dtype=np.float64)
    volume = np.empty(radii.size, dtype=np.float64)
    for radius_index, radius in enumerate(radii):
        for wavelength_index, (wavelength, refractive_index) in enumerate(zip(wavelengths, indices)):
            if distribution == "modified_gamma":
                if len(size_parameters) < 2:
                    raise ValueError("modified_gamma requires effective radius and variance")
                result = _mie_distribution(
                    radius, size_parameters[1], wavelength, refractive_index,
                    scattering_cosines, xres=size_parameter_resolution,
                )
            else:
                result = _mie_distribution(
                    radius, 0.0, wavelength, refractive_index, scattering_cosines,
                    xres=size_parameter_resolution, distribution=distribution,
                    size_parameters=(radius, size_parameters[1]),
                )
            extinction[radius_index, wavelength_index] = result[0]
            ssa[radius_index, wavelength_index] = result[1]
            asymmetry[radius_index, wavelength_index] = result[2]
            moments[:, radius_index, wavelength_index] = legendre_moments(
                scattering_cosines, phase_weights, result[3]
            ) / (2.0 * np.arange(phase_order) + 1.0)
            if wavelength_index == 0:
                volume[radius_index] = result[4]
    return OpticalProperties(
        wavelength_microns=wavelengths.astype(np.float64),
        effective_radius_microns=radii,
        extinction_coefficient=extinction,
        scattering_coefficient=extinction * ssa,
        single_scatter_albedo=ssa,
        asymmetry_parameter=asymmetry,
        phase_moments=moments,
        average_volume_per_particle=volume,
        absorption_coefficient=extinction * (1.0 - ssa),
        reference_wavelength_microns=0.55,
        reference_extinction_coefficient=extinction[:, 0],
    )


def legacy_mixed_mie_optics(
    effective_radii_microns: np.ndarray,
    wavelengths_microns: np.ndarray,
    components,
    refractive_indices: np.ndarray,
    *,
    phase_order: int = 1000,
    size_parameter_resolution: float = 0.4,
) -> OpticalProperties:
    """Combine legacy Mie components using ``generate_scattering_properties`` weights.

    The current implementation covers the production aerosol family in which
    all components share one log-normal size mode (for example ``aerosol_a79``)
    and the existing single-component modified-gamma cloud path.
    """

    radii = np.asarray(effective_radii_microns, dtype=np.float64).ravel()
    wavelengths = np.asarray(wavelengths_microns, dtype=np.float64).ravel()
    components = tuple(components)
    if not components:
        raise ValueError("At least one microphysical component is required")
    refractive_indices = np.asarray(refractive_indices, dtype=np.complex128)
    if refractive_indices.shape != (len(components), wavelengths.size):
        raise ValueError("refractive_indices must have shape (component, wavelength)")
    distributions = {component.size_distribution for component in components}
    if distributions == {"modified_gamma"} and len(components) == 1:
        return _component_optics(
            radii, wavelengths, refractive_indices[0], "modified_gamma",
            components[0].size_parameters, phase_order=phase_order,
            size_parameter_resolution=size_parameter_resolution,
        )
    if distributions != {"log_normal"}:
        raise NotImplementedError(
            "Python aerosol optics currently supports shared-mode log_normal Mie components"
        )
    effective_component_radii = np.asarray(
        [p[0] * math.exp(2.5 * math.log(p[1]) ** 2) for p in (c.size_parameters for c in components)],
        dtype=np.float64,
    )
    if not np.allclose(effective_component_radii, effective_component_radii[0], rtol=1e-7, atol=1e-12):
        raise NotImplementedError("Multiple distinct aerosol size modes are not yet supported")
    mode_radii = radii[:, None] * np.asarray(
        [components[0].size_parameters[0] / effective_component_radii[0]], dtype=np.float64
    )
    mode_radii = mode_radii[:, 0]
    component_optics = tuple(
        _component_optics(
            mode_radii, wavelengths, refractive_indices[index], component.size_distribution,
            component.size_parameters, phase_order=phase_order,
            size_parameter_resolution=size_parameter_resolution,
        )
        for index, component in enumerate(components)
    )
    ratios = np.asarray([component.mixing_ratio for component in components], dtype=np.float64)
    if np.any(ratios < 0) or not np.any(ratios > 0):
        raise ValueError("Microphysical mixing ratios must contain a positive value")
    extinction = np.zeros((radii.size, wavelengths.size), dtype=np.float64)
    ssa = np.zeros_like(extinction)
    asymmetry = np.zeros_like(extinction)
    moments = np.zeros((phase_order, radii.size, wavelengths.size), dtype=np.float64)
    for radius_index in range(radii.size):
        for wavelength_index in range(wavelengths.size):
            bext = np.asarray([
                ratio * optics.extinction_coefficient[radius_index, wavelength_index]
                for ratio, optics in zip(ratios, component_optics)
            ])
            bextw = np.asarray([
                ratio * optics.extinction_coefficient[radius_index, wavelength_index]
                * optics.single_scatter_albedo[radius_index, wavelength_index]
                for ratio, optics in zip(ratios, component_optics)
            ])
            total_bext = np.sum(bext)
            total_bextw = np.sum(bextw)
            extinction[radius_index, wavelength_index] = total_bext / np.sum(ratios)
            ssa[radius_index, wavelength_index] = total_bextw / total_bext
            asymmetry[radius_index, wavelength_index] = np.sum(np.asarray([
                value * optics.asymmetry_parameter[radius_index, wavelength_index]
                for value, optics in zip(bextw, component_optics)
            ])) / total_bextw
            moments[:, radius_index, wavelength_index] = np.sum(np.asarray([
                value * optics.phase_moments[:, radius_index, wavelength_index]
                for value, optics in zip(bextw, component_optics)
            ]), axis=0) / total_bextw
    return OpticalProperties(
        wavelength_microns=wavelengths,
        effective_radius_microns=radii,
        extinction_coefficient=extinction,
        scattering_coefficient=extinction * ssa,
        single_scatter_albedo=ssa,
        asymmetry_parameter=asymmetry,
        phase_moments=moments,
        average_volume_per_particle=component_optics[0].average_volume_per_particle,
        absorption_coefficient=extinction * (1.0 - ssa),
        reference_wavelength_microns=0.55,
        reference_extinction_coefficient=extinction[:, 0],
    )
