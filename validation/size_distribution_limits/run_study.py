#!/usr/bin/env python3
"""Convergence study for ORAC particle-size integration domains.

This is validation-only code.  It reproduces the radius quadrature and bulk
averaging in ``mie_size_dist_new`` while making the radius bounds explicit.
The preserved production Mie kernel is called without modification.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oraclut.idl_mirror.create_bwgp import gauss_cvf, quadrature  # noqa: E402
from oraclut.idl_mirror.load_mmdat import read_ri  # noqa: E402
from oraclut.optics.legacy_mie import mie_single_batch  # noqa: E402


WAVELENGTHS = np.asarray(
    [0.55, 0.64583063, 0.86778653, 1.62807465, 3.78025293, 11.02639866],
    dtype=np.float64,
)
UPPER_LIMITS = (5.0, 10.0, 20.0, 40.0, 60.0, 100.0)
LOWER_LIMITS = (0.0001, 0.001, 0.005, 0.01, 0.05)
REFERENCE_BOUNDS = (0.001, 100.0)
EXTENDED_UPPER_LIMITS = (60.0, 80.0, 100.0, 120.0, 150.0, 200.0)
SCALED_UPPER_FACTORS = (2.0, 2.5, 3.0, 3.5, 4.0)


@dataclass(frozen=True)
class DistributionCase:
    name: str
    family: str
    distribution: str
    radius_parameter: float
    spread: float
    refractive_index_file: str
    provenance: str

    @property
    def effective_radius(self) -> float:
        if self.distribution == "modified_gamma":
            return self.radius_parameter
        return self.radius_parameter * math.exp(2.5 * math.log(self.spread) ** 2)

    @property
    def legacy_bounds(self) -> tuple[float, float]:
        if self.distribution == "modified_gamma":
            return REFERENCE_BOUNDS
        tq = gauss_cvf(0.999)
        lower = math.exp(math.log(self.radius_parameter) + tq * math.log(self.spread))
        upper = math.exp(
            math.log(self.radius_parameter) - tq * math.log(self.spread) + math.log(4.0)
        )
        return lower, upper


CASES = tuple(
    DistributionCase(
        f"cloud_re{radius:g}", "liquid_cloud", "modified_gamma", radius,
        0.1111111, "H2O_Segelstein_1981.ri",
        "liquid-water_stg.mm; effective variance 0.1111111",
    )
    for radius in (5.0, 10.0, 15.0, 20.0, 30.0, 40.0)
) + (
    DistributionCase(
        "aerosol_small", "aerosol", "log_normal", 0.070, 1.700,
        "waf.ri", "aerosol_a79.mm small water-soluble mode",
    ),
    DistributionCase(
        "aerosol_accumulation", "aerosol", "log_normal", 0.217, 1.770,
        "chaiten_ash_Deguine_2020.ri", "volcanic-ash_ctn.mm accumulation mode",
    ),
    DistributionCase(
        "aerosol_coarse", "aerosol", "log_normal", 0.788, 1.822,
        "ssc.ri", "aerosol_a75.mm coarse sea-salt mode",
    ),
)


@dataclass
class KernelResult:
    radius: np.ndarray
    radius_weight: np.ndarray
    qext: np.ndarray
    qsca: np.ndarray
    g: np.ndarray
    f11: np.ndarray
    seconds: float


@dataclass
class BulkResult:
    extinction: float
    scattering: float
    ssa: float
    asymmetry: float
    phase: np.ndarray
    moments: np.ndarray
    number_fraction: float
    radius: np.ndarray
    number_terms: np.ndarray
    extinction_terms: np.ndarray
    scattering_terms: np.ndarray
    moment_terms: dict[int, np.ndarray]


def refractive_indices(filename: str, wavelengths: np.ndarray) -> np.ndarray:
    """Interpolate an existing ORAC RI table with the legacy sign convention."""

    table = read_ri(ROOT / "create_orac_lut" / "input_files" / "ri" / filename)
    order = np.argsort(table.wavl, kind="stable")
    real = np.interp(wavelengths, table.wavl[order], table.n[order])
    imaginary = np.interp(wavelengths, table.wavl[order], table.k[order])
    return real - 1j * imaginary


def radius_quadrature(
    lower: float, upper: float, wavelength: float, *, xres: float = 0.4
) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce the legacy linear-radius trapezoidal quadrature exactly."""

    if not 0.0 < lower < upper:
        raise ValueError(f"invalid radius bounds [{lower}, {upper}]")
    npoints = max(int(2.0 * np.pi * (upper - lower) / wavelength / xres), 200)
    radius = np.linspace(lower, upper, npoints, dtype=np.float64)
    radius_weight = np.full(npoints, (upper - lower) / (npoints - 1), dtype=np.float64)
    radius_weight[[0, -1]] *= 0.5
    return radius, radius_weight


def evaluate_kernel(
    lower: float,
    upper: float,
    wavelength: float,
    refractive_index: complex,
    scattering_cosines: np.ndarray,
    *,
    xres: float = 0.4,
    minimum_points: int = 200,
) -> KernelResult:
    radius, radius_weight = radius_quadrature(lower, upper, wavelength, xres=xres)
    if radius.size < minimum_points:
        radius = np.linspace(lower, upper, minimum_points, dtype=np.float64)
        radius_weight = np.full(
            minimum_points, (upper - lower) / (minimum_points - 1), dtype=np.float64
        )
        radius_weight[[0, -1]] *= 0.5
    size_parameter = 2.0 * np.pi * radius / wavelength
    started = time.perf_counter()
    particle = mie_single_batch(size_parameter, refractive_index, scattering_cosines)
    seconds = time.perf_counter() - started
    return KernelResult(
        radius, radius_weight, particle["qext"], particle["qsca"],
        particle["g"], particle["f11"], seconds,
    )


def restrict_kernel(kernel: KernelResult, lower: float, upper: float) -> KernelResult:
    """Select one domain without moving the common reference-grid nodes."""

    selected = (kernel.radius >= lower) & (kernel.radius <= upper)
    if np.count_nonzero(selected) < 2:
        raise ValueError(f"domain [{lower}, {upper}] has fewer than two reference-grid nodes")
    return KernelResult(
        kernel.radius[selected], kernel.radius_weight[selected], kernel.qext[selected],
        kernel.qsca[selected], kernel.g[selected], kernel.f11[selected], 0.0,
    )


def distribution_density(case: DistributionCase, radius: np.ndarray) -> np.ndarray:
    if case.distribution == "modified_gamma":
        variance = case.spread
        alpha = (1.0 - 3.0 * variance) / variance
        beta = 1.0 / (case.radius_parameter * variance)
        normalization = beta ** (-alpha - 1.0) * math.gamma(alpha + 1.0)
        return radius**alpha * np.exp(-beta * radius) / normalization
    log_spread = math.log(case.spread)
    return np.exp(-0.5 * (np.log(radius / case.radius_parameter) / log_spread) ** 2) / (
        math.sqrt(2.0 * math.pi) * radius * log_spread
    )


def legendre_polynomials(mu: np.ndarray, orders: tuple[int, ...]) -> dict[int, np.ndarray]:
    wanted = set(orders)
    output: dict[int, np.ndarray] = {}
    previous = np.ones_like(mu)
    if 0 in wanted:
        output[0] = previous.copy()
    if max(wanted, default=0) == 0:
        return output
    current = mu.copy()
    if 1 in wanted:
        output[1] = current.copy()
    for order in range(2, max(wanted) + 1):
        following = (
            (2.0 * order - 1.0) / order * mu * current
            - (order - 1.0) / order * previous
        )
        if order in wanted:
            output[order] = following.copy()
        previous, current = current, following
    return output


def phase_moments(
    phase: np.ndarray, mu: np.ndarray, angle_weights: np.ndarray
) -> np.ndarray:
    """Return ORAC-stored beta_l = 0.5 integral(P phase_l) dmu."""

    moments = np.empty(mu.size, dtype=np.float64)
    previous = np.ones_like(mu)
    moments[0] = 0.5 * np.sum(phase * angle_weights)
    if mu.size == 1:
        return moments
    current = mu.copy()
    moments[1] = 0.5 * np.sum(phase * current * angle_weights)
    for order in range(2, mu.size):
        following = (
            (2.0 * order - 1.0) / order * mu * current
            - (order - 1.0) / order * previous
        )
        moments[order] = 0.5 * np.sum(phase * following * angle_weights)
        previous, current = current, following
    return moments


def average_distribution(
    case: DistributionCase,
    kernel: KernelResult,
    mu: np.ndarray,
    angle_weights: np.ndarray,
    selected_moment_orders: tuple[int, ...],
) -> BulkResult:
    number_terms = kernel.radius_weight * distribution_density(case, kernel.radius)
    area_terms = number_terms * np.pi * kernel.radius**2
    extinction_terms = area_terms * kernel.qext
    scattering_terms = area_terms * kernel.qsca
    extinction = float(np.sum(extinction_terms))
    scattering = float(np.sum(scattering_terms))
    phase = np.sum(
        scattering_terms[:, None] * kernel.f11, axis=0
    ) / scattering
    moments = phase_moments(phase, mu, angle_weights)
    polynomials = legendre_polynomials(mu, selected_moment_orders)
    moment_terms = {
        order: scattering_terms
        * (0.5 * np.sum(kernel.f11 * polynomial[None, :] * angle_weights[None, :], axis=1))
        / scattering
        for order, polynomial in polynomials.items()
    }
    return BulkResult(
        extinction=extinction,
        scattering=scattering,
        ssa=scattering / extinction,
        asymmetry=float(np.sum(scattering_terms * kernel.g) / scattering),
        phase=phase,
        moments=moments,
        number_fraction=float(np.sum(number_terms)),
        radius=kernel.radius,
        number_terms=number_terms,
        extinction_terms=extinction_terms,
        scattering_terms=scattering_terms,
        moment_terms=moment_terms,
    )


def comparison_metrics(
    value: BulkResult, reference: BulkResult, angle_weights: np.ndarray
) -> dict[str, float]:
    return {
        "extinction_relative_error": (value.extinction - reference.extinction) / reference.extinction,
        "scattering_relative_error": (value.scattering - reference.scattering) / reference.scattering,
        "ssa_absolute_error": value.ssa - reference.ssa,
        "asymmetry_absolute_error": value.asymmetry - reference.asymmetry,
        "phase_l1_error": 0.5 * float(np.sum(np.abs(value.phase - reference.phase) * angle_weights)),
        "phase_max_absolute_error": float(np.max(np.abs(value.phase - reference.phase))),
        "moment_max_absolute_error": float(np.max(np.abs(value.moments - reference.moments))),
        "moment_1_absolute_error": float(value.moments[1] - reference.moments[1]),
        "moment_10_absolute_error": float(value.moments[10] - reference.moments[10]),
        "moment_50_absolute_error": float(value.moments[50] - reference.moments[50]),
        "moment_100_absolute_error": float(value.moments[100] - reference.moments[100]),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def crossing_radius(radius: np.ndarray, terms: np.ndarray, fraction: float) -> float:
    cumulative = np.cumsum(terms)
    target = fraction * cumulative[-1]
    return float(radius[int(np.searchsorted(cumulative, target, side="left"))])


def plot_cumulative(
    case: DistributionCase,
    wavelength: float,
    result: BulkResult,
    output: Path,
) -> None:
    import matplotlib.pyplot as plt

    radius = result.radius
    density = distribution_density(case, radius)
    figure, axes = plt.subplots(1, 3, figsize=(13.0, 4.1), constrained_layout=True)
    axes[0].loglog(radius, density)
    axes[0].set(xlabel="particle radius (µm)", ylabel="number density (normalised)")
    for terms, label in (
        (result.number_terms, "number"),
        (result.extinction_terms, "extinction"),
        (result.scattering_terms, "scattering"),
    ):
        axes[1].semilogx(radius, np.cumsum(terms) / np.sum(terms), label=label)
    axes[1].set(xlabel="particle radius (µm)", ylabel="cumulative fraction", ylim=(-0.02, 1.02))
    axes[1].legend()
    for order, terms in result.moment_terms.items():
        axes[2].semilogx(radius, np.cumsum(terms), label=rf"$\beta_{{{order}}}$")
    axes[2].set(xlabel="particle radius (µm)", ylabel="cumulative ORAC Legendre moment")
    axes[2].legend()
    figure.suptitle(
        f"{case.name}: {case.provenance}\n{wavelength:g} µm, "
        f"domain {radius[0]:g}–{radius[-1]:g} µm"
    )
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_envelopes(rows: list[dict[str, object]], output: Path, study: str) -> None:
    import matplotlib.pyplot as plt

    limit_key = "upper_microns" if study == "upper" else "lower_microns"
    filtered = [row for row in rows if row["study"] == study]
    limits = sorted({float(row[limit_key]) for row in filtered})
    metrics = (
        ("extinction_relative_error", "max |relative extinction error|"),
        ("phase_l1_error", "max phase-function L1 error"),
        ("moment_max_absolute_error", "max |Legendre-moment error|"),
        ("speedup_vs_reference", "median measured Mie speed-up"),
    )
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), constrained_layout=True)
    for axis, (metric, label) in zip(axes.ravel(), metrics):
        for family in ("liquid_cloud", "aerosol"):
            values = []
            for limit in limits:
                subset = [
                    abs(float(row[metric]))
                    for row in filtered
                    if row["family"] == family and float(row[limit_key]) == limit
                ]
                values.append(np.median(subset) if metric.startswith("speedup") else max(subset))
            axis.plot(limits, values, marker="o", label=family.replace("_", " "))
        if metric != "speedup_vs_reference":
            axis.set_yscale("log")
        if study == "lower":
            axis.set_xscale("log")
        axis.set(xlabel=f"{study} radius limit (µm)", ylabel=label)
        axis.grid(True, alpha=0.25)
        axis.legend()
    figure.suptitle(f"Worst-case size-domain convergence: {study}-limit study")
    figure.savefig(output, dpi=180)
    plt.close(figure)


def run(results: Path, phase_order: int) -> None:
    if phase_order < 101:
        raise ValueError("phase_order must be at least 101 to assess moment 100")
    results.mkdir(parents=True, exist_ok=True)
    abscissa, angle_weights = quadrature("g", phase_order)
    mu = -abscissa
    selected_moments = (1, 10, 50, 100)
    convergence_rows: list[dict[str, object]] = []
    timing_rows: list[dict[str, object]] = []
    cumulative_rows: list[dict[str, object]] = []
    fixed_grid_rows: list[dict[str, object]] = []
    extended_cloud_rows: list[dict[str, object]] = []
    scaled_cloud_rows: list[dict[str, object]] = []

    groups: dict[str, list[DistributionCase]] = {}
    for case in CASES:
        groups.setdefault(case.refractive_index_file, []).append(case)

    for ri_file, cases in groups.items():
        indices = refractive_indices(ri_file, WAVELENGTHS)
        case_domains = {case.name: case.legacy_bounds for case in cases}
        common_domains = {
            REFERENCE_BOUNDS,
            *((case.legacy_bounds[0], upper) for case in cases for upper in UPPER_LIMITS),
            *((lower, case.legacy_bounds[1]) for case in cases for lower in LOWER_LIMITS),
            *((LOWER_LIMITS[0], case.legacy_bounds[1]) for case in cases),
            *case_domains.values(),
        }
        for wavelength, refractive_index in zip(WAVELENGTHS, indices):
            kernels: dict[tuple[float, float], KernelResult] = {}
            physical_kernels: dict[tuple[float, float], KernelResult] = {}
            for lower, upper in sorted(common_domains):
                kernel = evaluate_kernel(lower, upper, wavelength, refractive_index, mu)
                kernels[(lower, upper)] = kernel
                timing_rows.append({
                    "refractive_index_file": ri_file,
                    "wavelength_microns": wavelength,
                    "lower_microns": lower,
                    "upper_microns": upper,
                    "radius_points": kernel.radius.size,
                    "mie_seconds": kernel.seconds,
                })

            for case in cases:
                averaged = {
                    bounds: average_distribution(
                        case, kernel, mu, angle_weights, selected_moments
                    )
                    for bounds, kernel in kernels.items()
                }
                # The fixed domain is the production reference only for the
                # modified-gamma cloud path.  The production aerosol path
                # ignores those two parameters and derives adaptive log-normal
                # bounds from mode radius and spread.
                reference_bounds = case.legacy_bounds
                reference = averaged[reference_bounds]
                reference_seconds = kernels[reference_bounds].seconds
                specifications = [
                    ("upper", reference_bounds[0], upper)
                    for upper in UPPER_LIMITS
                ] + [
                    ("lower", lower, reference_bounds[1])
                    for lower in LOWER_LIMITS
                ] + [
                    ("production_reference", *reference_bounds),
                    ("common_fixed_domain", *REFERENCE_BOUNDS),
                ]
                for study, lower, upper in specifications:
                    result = averaged[(lower, upper)]
                    row: dict[str, object] = {
                        "case": case.name,
                        "family": case.family,
                        "distribution": case.distribution,
                        "provenance": case.provenance,
                        "effective_radius_microns": case.effective_radius,
                        "mode_radius_microns": case.radius_parameter,
                        "spread_or_effective_variance": case.spread,
                        "refractive_index_file": case.refractive_index_file,
                        "wavelength_microns": wavelength,
                        "study": study,
                        "reference_lower_microns": reference_bounds[0],
                        "reference_upper_microns": reference_bounds[1],
                        "lower_microns": lower,
                        "upper_microns": upper,
                        "radius_points": kernels[(lower, upper)].radius.size,
                        "mie_seconds": kernels[(lower, upper)].seconds,
                        "speedup_vs_reference": reference_seconds / kernels[(lower, upper)].seconds,
                        "number_fraction": result.number_fraction,
                        "extinction": result.extinction,
                        "scattering": result.scattering,
                        "single_scattering_albedo": result.ssa,
                        "asymmetry_parameter": result.asymmetry,
                    }
                    row.update(comparison_metrics(result, reference, angle_weights))
                    convergence_rows.append(row)

                # Isolate physical tail omission from changes in quadrature
                # sampling.  All candidates below use identical radius nodes
                # from a denser xres=0.1 master grid; unlike the operational
                # rows above, moving a bound cannot shift every Mie sample.
                physical_bounds = (LOWER_LIMITS[0], reference_bounds[1])
                if physical_bounds not in physical_kernels:
                    physical_kernels[physical_bounds] = evaluate_kernel(
                        physical_bounds[0], physical_bounds[1], wavelength,
                        refractive_index, mu, xres=0.1, minimum_points=800,
                    )
                physical_kernel = physical_kernels[physical_bounds]
                physical_reference = average_distribution(
                    case,
                    restrict_kernel(physical_kernel, *reference_bounds),
                    mu,
                    angle_weights,
                    selected_moments,
                )
                physical_specifications = [
                    ("upper", reference_bounds[0], min(upper, reference_bounds[1]))
                    for upper in UPPER_LIMITS
                ] + [
                    ("lower", lower, reference_bounds[1]) for lower in LOWER_LIMITS
                ]
                for study, lower, upper in physical_specifications:
                    physical_result = average_distribution(
                        case,
                        restrict_kernel(physical_kernel, lower, upper),
                        mu,
                        angle_weights,
                        selected_moments,
                    )
                    physical_row: dict[str, object] = {
                        "case": case.name,
                        "family": case.family,
                        "distribution": case.distribution,
                        "effective_radius_microns": case.effective_radius,
                        "wavelength_microns": wavelength,
                        "study": study,
                        "lower_microns": lower,
                        "upper_microns": upper,
                        "reference_lower_microns": reference_bounds[0],
                        "reference_upper_microns": reference_bounds[1],
                        "master_grid_xres": 0.1,
                        "master_grid_points": physical_kernel.radius.size,
                        "number_fraction": physical_result.number_fraction,
                        "extinction": physical_result.extinction,
                        "scattering": physical_result.scattering,
                        "single_scattering_albedo": physical_result.ssa,
                        "asymmetry_parameter": physical_result.asymmetry,
                    }
                    physical_row.update(
                        comparison_metrics(physical_result, physical_reference, angle_weights)
                    )
                    fixed_grid_rows.append(physical_row)

                if np.isclose(wavelength, WAVELENGTHS[0]):
                    cumulative_bounds = (LOWER_LIMITS[0], reference_bounds[1])
                    full = averaged[cumulative_bounds]
                    cumulative_rows.append({
                        "case": case.name,
                        "family": case.family,
                        "distribution": case.distribution,
                        "effective_radius_microns": case.effective_radius,
                        "mode_radius_microns": case.radius_parameter,
                        "spread_or_effective_variance": case.spread,
                        "wavelength_microns": wavelength,
                        **{
                            f"number_r{100 * fraction:g}_microns": crossing_radius(
                                full.radius, full.number_terms, fraction
                            )
                            for fraction in (0.9, 0.99, 0.999, 0.9999)
                        },
                        **{
                            f"extinction_r{100 * fraction:g}_microns": crossing_radius(
                                full.radius, full.extinction_terms, fraction
                            )
                            for fraction in (0.9, 0.99, 0.999, 0.9999)
                        },
                        **{
                            f"scattering_r{100 * fraction:g}_microns": crossing_radius(
                                full.radius, full.scattering_terms, fraction
                            )
                            for fraction in (0.9, 0.99, 0.999, 0.9999)
                        },
                        **{f"beta_{order}": full.moments[order] for order in selected_moments},
                    })
                    plot_cumulative(
                        case, wavelength, full, results / f"cumulative_{case.name}.png"
                    )

            cloud_cases = [case for case in cases if case.family == "liquid_cloud"]
            if cloud_cases:
                extended_kernel = evaluate_kernel(
                    LOWER_LIMITS[0], EXTENDED_UPPER_LIMITS[-1], wavelength,
                    refractive_index, mu, xres=0.1, minimum_points=800,
                )
                for case in cloud_cases:
                    extended_reference = average_distribution(
                        case,
                        restrict_kernel(
                            extended_kernel, REFERENCE_BOUNDS[0], EXTENDED_UPPER_LIMITS[-1]
                        ),
                        mu,
                        angle_weights,
                        selected_moments,
                    )
                    for upper in EXTENDED_UPPER_LIMITS:
                        value = average_distribution(
                            case,
                            restrict_kernel(extended_kernel, REFERENCE_BOUNDS[0], upper),
                            mu,
                            angle_weights,
                            selected_moments,
                        )
                        row: dict[str, object] = {
                            "case": case.name,
                            "effective_radius_microns": case.effective_radius,
                            "wavelength_microns": wavelength,
                            "lower_microns": REFERENCE_BOUNDS[0],
                            "upper_microns": upper,
                            "reference_upper_microns": EXTENDED_UPPER_LIMITS[-1],
                        }
                        row.update(comparison_metrics(value, extended_reference, angle_weights))
                        extended_cloud_rows.append(row)
                    for factor in SCALED_UPPER_FACTORS:
                        upper = min(factor * case.effective_radius, EXTENDED_UPPER_LIMITS[-1])
                        value = average_distribution(
                            case,
                            restrict_kernel(extended_kernel, REFERENCE_BOUNDS[0], upper),
                            mu,
                            angle_weights,
                            selected_moments,
                        )
                        row = {
                            "case": case.name,
                            "effective_radius_microns": case.effective_radius,
                            "wavelength_microns": wavelength,
                            "upper_factor_times_effective_radius": factor,
                            "upper_microns": upper,
                            "reference_upper_microns": EXTENDED_UPPER_LIMITS[-1],
                        }
                        row.update(comparison_metrics(value, extended_reference, angle_weights))
                        scaled_cloud_rows.append(row)

    write_csv(results / "convergence.csv", convergence_rows)
    write_csv(results / "kernel_timing.csv", timing_rows)
    write_csv(results / "cumulative_quantiles.csv", cumulative_rows)
    write_csv(results / "fixed_grid_tail_convergence.csv", fixed_grid_rows)
    write_csv(results / "cloud_extended_upper_convergence.csv", extended_cloud_rows)
    write_csv(results / "cloud_scaled_upper_convergence.csv", scaled_cloud_rows)
    plot_envelopes(convergence_rows, results / "upper_limit_convergence.png", "upper")
    plot_envelopes(convergence_rows, results / "lower_limit_convergence.png", "lower")

    metadata = {
        "phase_order": phase_order,
        "xres": 0.4,
        "reference_bounds_microns": REFERENCE_BOUNDS,
        "reference_note": (
            "0.001–100 is the production modified-gamma reference. Log-normal "
            "aerosols use their legacy mode-radius/spread-dependent bounds."
        ),
        "wavelengths_microns": WAVELENGTHS.tolist(),
        "wavelength_provenance": {
            "0.55": "ORAC reference wavelength and near SLSTR S1 (0.5541 µm)",
            "0.64583063": "Aqua MODIS channel 1 centre from repository SRF",
            "0.86778653": "Sentinel-3A SLSTR channel 3 centre from repository SRF",
            "1.62807465": "Aqua MODIS channel 6 centre from repository SRF",
            "3.78025293": "Aqua MODIS channel 20 centre from repository SRF",
            "11.02639866": "Aqua MODIS channel 31 centre from repository SRF",
        },
        "cases": [case.__dict__ | {"legacy_bounds": case.legacy_bounds} for case in CASES],
    }
    (results / "study_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", type=Path, default=HERE / "results",
        help="output directory (default: validation/size_distribution_limits/results)",
    )
    parser.add_argument(
        "--phase-order", type=int, default=128,
        help="Gauss-Legendre phase angles / moments (default: 128)",
    )
    args = parser.parse_args()
    output = args.results.resolve()
    try:
        output.relative_to(ROOT)
    except ValueError as exc:
        raise SystemExit("results must remain inside the repository") from exc
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "validation" / "tmp" / "matplotlib-cache"))
    run(output, args.phase_order)
    print(f"Results written to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
