#!/usr/bin/env python3
"""Targeted 1000-moment tail check for the size-domain study."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))

from run_study import (  # noqa: E402
    CASES,
    average_distribution,
    evaluate_kernel,
    refractive_indices,
    restrict_kernel,
)
from oraclut.idl_mirror.create_bwgp import quadrature  # noqa: E402


def moment_band_errors(value: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    difference = np.abs(value - reference)
    return {
        "moment_max_0_127": float(np.max(difference[:128])),
        "moment_max_128_511": float(np.max(difference[128:512])),
        "moment_max_512_999": float(np.max(difference[512:])),
        "moment_250_absolute_error": float(value[250] - reference[250]),
        "moment_500_absolute_error": float(value[500] - reference[500]),
        "moment_999_absolute_error": float(value[999] - reference[999]),
    }


def main() -> int:
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    phase_order = 1000
    wavelength = 0.55
    abscissa, angle_weights = quadrature("g", phase_order)
    mu = -abscissa
    selected = (1, 10, 50, 100)
    rows: list[dict[str, object]] = []

    specifications = (
        (next(case for case in CASES if case.name == "cloud_re40"), 0.001, 200.0,
         (60.0, 80.0, 100.0, 120.0, 150.0, 200.0)),
        (next(case for case in CASES if case.name == "aerosol_coarse"), 0.0001,
         next(case for case in CASES if case.name == "aerosol_coarse").legacy_bounds[1],
         (10.0, 15.0, 20.0)),
    )
    for case, master_lower, master_upper, upper_limits in specifications:
        refractive_index = refractive_indices(
            case.refractive_index_file, np.asarray([wavelength])
        )[0]
        kernel = evaluate_kernel(
            master_lower, master_upper, wavelength, refractive_index, mu, xres=0.4
        )
        reference_lower = 0.001 if case.family == "liquid_cloud" else case.legacy_bounds[0]
        reference = average_distribution(
            case,
            restrict_kernel(kernel, reference_lower, master_upper),
            mu,
            angle_weights,
            selected,
        )
        for upper in upper_limits:
            value = average_distribution(
                case,
                restrict_kernel(kernel, reference_lower, upper),
                mu,
                angle_weights,
                selected,
            )
            row: dict[str, object] = {
                "case": case.name,
                "wavelength_microns": wavelength,
                "phase_order": phase_order,
                "lower_microns": reference_lower,
                "upper_microns": upper,
                "reference_upper_microns": master_upper,
                "extinction_relative_error": (
                    value.extinction - reference.extinction
                ) / reference.extinction,
                "scattering_relative_error": (
                    value.scattering - reference.scattering
                ) / reference.scattering,
                "phase_l1_error": 0.5 * float(
                    np.sum(np.abs(value.phase - reference.phase) * angle_weights)
                ),
            }
            row.update(moment_band_errors(value.moments, reference.moments))
            rows.append(row)

    path = results / "high_order_1000_moment_check.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Results written to {path}")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "validation" / "tmp" / "matplotlib-cache"))
    raise SystemExit(main())
