# VALIDATION.md

## Purpose

Validation is part of the implementation, not an end-stage activity.

The first objective is to demonstrate that Python reproduces a selected legacy ORAC LUT for the same inputs and scientific assumptions.

## Validation hierarchy

Compare at the earliest available stage.

If a final LUT differs, do not begin by tuning the final output. Compare successively:

1. input configuration;
2. coordinate grids;
3. atmospheric quantities;
4. particle optical properties;
5. optical-depth scaling;
6. phase-function representation;
7. radiative-transfer layer inputs;
8. monochromatic radiative-transfer outputs;
9. spectral-response-integrated quantities;
10. LUT assembly and dimension ordering;
11. stored output.

Stop at the first scientifically meaningful divergence and diagnose it there.

For the first observed V21 product, validate the complete Meteosat-10 SEVIRI
liquid-water `stg` case without pressure, then validate the Meteosat-10 aerosol
`a79` case with its pressure dimension and command-line gas setting. Preserve
the exact wrapper, driver, model/LUT overrides, and source revision in each
case record; the output version field alone does not identify the legacy
calculation.

## Required checks

For arrays, record:

- shape;
- dtype/precision;
- dimension names/order;
- units;
- minimum and maximum;
- NaN/Inf counts;
- absolute difference;
- relative difference where meaningful;
- RMS difference;
- maximum absolute and relative differences;
- location/index of extrema.

For coordinate grids, compare values directly rather than only dimensions.

## Optical-property checks

Verify explicitly:

- wavelength units;
- size/effective-radius definition;
- extinction/scattering convention;
- single-scattering albedo;
- reference wavelength for optical depth;
- phase-function angular convention;
- phase-function normalisation;
- Legendre-moment convention if used;
- interpolation/extrapolation behaviour;
- mixture weighting for multi-component particles.

## Geometry checks

Verify explicitly:

- solar zenith angle;
- viewing/satellite zenith angle;
- relative azimuth definition and range;
- degrees versus radians;
- cosine conventions;
- ordering of geometry dimensions.

## Spectral checks

For every channel establish whether the legacy code uses:

- centre wavelength only;
- monochromatic sampling;
- spectral response convolution;
- solar-weighted response;
- thermal/Planck weighting;
- another channel-specific rule.

Visible and infrared paths must not be assumed identical.

## Radiative-transfer checks

Record the exact legacy DISORT setup needed for reproduction, including where relevant:

- number of streams;
- delta-M or other corrections;
- layer ordering;
- boundary conditions;
- thermal source treatment;
- solar source treatment;
- phase-function moments;
- Rayleigh treatment;
- cloud/aerosol layer placement;
- numerical precision.

Do not change solver settings during equivalence testing.

## LUT comparison artefacts

Each validation case should produce a compact machine-readable summary and, where useful, diagnostic plots.

Suggested outputs:

```text
validation/comparisons/<case>/
   summary.txt
   summary.json
   coordinates.txt
   residual_<quantity>.npy
   plot_<quantity>_slice.png
```

Large generated LUTs should normally remain outside version control.

The executable Python reproduction writes compact and full candidates under
`validation/generated/`. Its comparison records dimensions, ordering, dtypes,
finite minimum/maximum values, fill masks, NaN/Inf counts, absolute/RMS/
relative residuals, and the indices of both the maximum absolute and maximum
relative residuals. It also records structural equivalence, typed metadata
differences, and byte-level SHA-256 identity separately. The full comparison against the read-only Meteosat-10
SEVIRI liquid-water reference is
`validation/generated/python_full_comparison.json`.

The current forensic follow-up is in
`validation/diagnostics/full_r0v_residual_summary.json`. The largest `R_0v`
residual is at channel 1, wavelength 0.6381768 microns, optical depth 2.0,
effective radius 23 microns, solar zenith 50 degrees, satellite zenith 80
degrees, relative azimuth 0 degrees: reference 0.25314441 versus Python
0.15585500. Forty-eight of the top 50 residuals share that channel and
optical/radius/solar-angle state. For the explicit 110-point suspected
sensitive cluster the RMS is 0.041996; outside it the maximum is 0.013135 and
RMS is `5.1170e-5`.

The targeted Python DISORT state is captured before the Fortran call in
`validation/diagnostics/worst_r0v_disort_state_ifort.json` (45 layers, 60
streams, 1000 float32 phase coefficients, two output levels, 20 views, 11
azimuths). The ifort and gfortran input arrays are bitwise equal. Their
outputs differ substantially at this warning state: target `R_0v`
0.15590107 versus 0.05717356, with maximum `UU` difference 4.6931 and RMS
1.1147. A fresh-process serialized Python state gives 0.28200659. The user’s
licensed IDL 8.9 production `DISORTTOIDL` call gives `R_0v = 0.288875`,
`RFldir = 64.2788, 2.61441`, `RFldn = -0.000579834, 48.5968`, and
`Flup = 13.0656, 0.000335587`, while emitting the same
`UPBEAM--SGECO says matrix near singular` warning. The archived `0.25314441`
value is consequently documented numerical evidence, not a bit-for-bit
acceptance target. Repeated calls after initialization are repeatable, and
the preserved DISORT source contains `SAVE` state. The runnable diagnostic is
`validation/diagnostics/run_worst_r0v_idl.pro`.

The tiny fill difference is not treated as a Python mask rule. The legacy
scattering generator populates all four optical arrays, and the full
reference is finite; the 12 fills are classified as a singleton-channel
legacy V2 NetCDF writer artifact in the direct `reform`/`ncdf_varput` path.
The Python writer retains finite optics, with the comparator and regression
test recording the observed 12-element mask difference.

The tiny captured reference agrees with the Python cloud operators at the
expected single-precision/compiler scale. Its captured microphysical arrays
are fill-valued while the Python legacy-Mie arrays are finite, and the
comparison records that fill-mask discrepancy explicitly. The full case must be
validated with the temporary `intel-compilers/2022` build: gfortran can take a
different branch in the preserved DISORT `UPBEAM--SGECO` near-singular solve.
This is recorded as a kernel/compiler reproducibility issue, not hidden by a
loose tolerance. See `docs/python_legacy_reproduction.md` for the generator
and comparison commands.

## Aerosol/distributed-profile implementation

The compact Python aerosol case uses the existing `aerosol_test.lut` grid,
`aerosol_a79.mm`, Meteosat-10 SEVIRI channel 1, atmosphere 2, Rayleigh, gas,
60 DISORT streams, and phase order 1000. It writes the pressure-aware V2
dimensions and finite operators; a gas-enabled run is retained as
`validation/generated/python_aerosol_test_v21_gas.nc`. The selected `a79`
model has two log-normal Mie components, `waf`/`saf`, with number-mixing
ratios `0.625`/`0.375`; it is selected because the complete current input set
and an archived production V21 product are present.

The compact legacy aerosol reference has not yet been captured in this shell:
`command -v idl` is empty and the IDL licence is unavailable. The safe runner
and driver are `scripts/generate_aerosol_test_reference.sh` and
`validation/reference_cases/meteosat10_seviri_aerosol_a79_test.driver`; the
manual licensed command is documented in
`docs/generating_luts_with_python.md`. The Python optics nevertheless agree
with the archived full `a79` product at shared radius `0.01` microns, channel
1: extinction absolute difference `1.42e-13`, extinction-ratio difference
`1.19e-7`, SSA difference `1.86e-9`, and asymmetry difference `0`. This is an
optics check, not a substitute for the pending pressure-aware RT comparison;
the latter should be run when the compact legacy capture is available.

## Acceptance

There is no single universal tolerance.

Use exact comparison for integers, identifiers, dimensions, and coordinate definitions where exact equivalence is expected.

For floating-point calculations:

- establish the precision of the legacy calculation;
- determine whether algorithms and operation ordering are identical;
- quantify residual structure;
- choose tolerances only after understanding expected numerical differences.

A result is not considered validated merely because it passes a loose relative tolerance.

All remaining systematic differences must be documented.

## Regression tests

Once a component is matched, add a regression test before moving on.

Prefer compact reference slices or intermediate quantities rather than committing very large complete LUT files.

Tests should make accidental changes to scientific conventions fail loudly.

## POM-stage validation

POM validation is distinct from legacy reproduction.

When POM is introduced:

- first verify interface consistency;
- then compare POM with the legacy optical model for cases where equivalent physics is expected;
- finally quantify intentional scientific differences for new particle shapes or microphysical models.

Never describe an intentional physics change as a porting discrepancy.
