# Visible + infrared validation matrix

Validating the Python ORAC LUT generator against SEVIRI channel 1 alone is not
sufficient. Channel 1 is a solar channel: it never exercises the thermal
emission route (Planck source, band-averaged blackbody normalisation, the
`E_md` emissivity operator, the thermal-channel metadata) or the gas optical
depths of an infrared window. A formulation is therefore only described as
broadly validated after **both** a visible and an infrared case agree with the
legacy IDL implementation, and the two are always reported separately — never
averaged into one number.

## The infrared channel: SEVIRI channel 9 (IR10.8)

Chosen from `create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst`, not
guessed:

| Property | Channel 1 | Channel 9 |
| --- | --- | --- |
| instrument-file classification | `solar channels` | `thermal channels` only (not in `solar channels`) |
| SRF file / points | `rtcoef_msg_3_seviri_srf_ch01.txt`, 51 | `rtcoef_msg_3_seviri_srf_ch09.txt` ("Channel IR10.8"), 101 |
| central wavelength / wavenumber | 0.6382 µm / 15 731.5 cm⁻¹ | 10.7963 µm / 928.72 cm⁻¹ |
| blackbody fit constants (`bbconstants.pro`) | none (0, 0, 0, 0) | B1 9577.68, B2 1337.94, T1 0.618241, T2 0.998317 |
| thermal noise metadata | none | refbt 300 K, nedt 0.06 K |
| solar terms (F0, `R_0v`, `R_0d`, `T_0d`, `T_00`) | present | absent |
| emissivity operator `E_md` | absent | present |
| gas optical depths available | `ModtranGasOpd_A0-6_meteosat-10_seviri_ch01.gas` | `…_ch09.gas` |

Channel 9 is the canonical ORAC 10.8 µm window channel and is thermal-only,
so it isolates the emission route. Channel 4 (3.9 µm) is flagged both solar
and thermal ("mixed") and would not isolate it; it remains a further test for
later. Gas absorption follows the legacy products: off for the cloud
formulation (the `a01` product family never used it; ORAC applies gas
corrections separately), on for the aerosol formulation (`a12`).

Channel selection lives in the driver/configuration, not in the `.lut` grid
file, so the existing compact grids (`liquid-water-cloud_test.lut`,
`aerosol_test.lut`) are reused unchanged and only the channel list differs.
Every validation product sets `output=` explicitly because the legacy file name
does not encode the channel set.

## Cases

| Case | Formulation | Channels | Python configuration | Legacy driver | Outputs |
| --- | --- | --- | --- | --- | --- |
| cloud visible | cloud | 1 | `configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver` | (existing capture) | `validation/generated/meteosat10_seviri_liquid_water_stg_cloud_test/` |
| cloud IR | cloud | 9 | `configs/validation_…_cloud_…_ch09.driver` | `…_cloud_test_ch09.driver` | `validation/generated/cloud_ir_comparison/` |
| cloud visible+IR | cloud | 1, 9 | `…_cloud_…_ch01_ch09.driver` | `…_cloud_test_ch01_ch09.driver` | `validation/generated/cloud_visible_ir_comparison/` |
| aerosol visible | aerosol | 1 | `configs/meteosat-10_seviri_aerosol_a79_test.driver` | `…_aerosol_a79_test.driver` | `validation/generated/aerosol_legacy_reference/`, `aerosol_comparison/` |
| aerosol IR | aerosol | 9 | `…_aerosol_a79_test_ch09.driver` | `…_aerosol_a79_test_ch09.driver` | `validation/generated/aerosol_ir_comparison/` |
| aerosol visible+IR | aerosol | 1, 9 | `…_aerosol_a79_test_ch01_ch09.driver` | `…_aerosol_a79_test_ch01_ch09.driver` | `validation/generated/aerosol_visible_ir_comparison/` |

Each comparison directory holds `python/` (Python NetCDF and run log),
`legacy/` (legacy NetCDF, `capture_manifest.json`, `provenance.log`,
`stdout.log`, `stderr.log`, `repeatability.json`, `scatfile.sav`),
`comparison.json` and `comparison.md`. All legacy references were generated
with `scripts/generate_legacy_reference.sh` from the coherent extracted working
tree; the provenance gate confirmed 26 routines resolved from it and none from
the active V2.1 tree in every run.

## Results

Residuals are maximum absolute differences on the radiative-transfer operators
(quantities of order one, `E_md` in the range 0–1). Bands: GREEN ≤ 1e-3,
AMBER ≤ 1e-2, RED beyond — plus AMBER for a systematic relative bias — as
calibrated on the original cloud channel-1 case; case-specific, not universal.

### First pass (2026-09-14, before the E_md correction)

| Case | Structural | Numerical | Worst operator residual | Notes |
| --- | --- | --- | --- | --- |
| cloud visible (ch 1) | PASS | **GREEN** | 1.26e-4 (`R_dd`) | benchmark |
| cloud IR (ch 9) | PASS | **GREEN** | 6.6e-5 (`R_dd`); `E_md` 1.5e-6 | thermal metadata identical |
| cloud visible+IR (ch 1, 9) | PASS | **GREEN** | 1.26e-4 (ch 1), 6.6e-5 (ch 9) | one product, both regimes |
| aerosol visible (ch 1) | FAIL* | **GREEN** (band) | 9.4e-4 (`T_dv`), 7.8e-4 (`R_0d`) | localised 0.3–0.9 % low bias at r_eff = 0.01 µm, see Finding 2 |
| aerosol IR (ch 9) | FAIL* | **RED** | `E_md` 99.0; all other operators ≤ 3e-5 | Python `E_md` exactly 100 × legacy — Finding 1 |
| aerosol visible+IR (ch 1, 9) | FAIL* | **RED** | inherits both aerosol findings | |

The pre-correction aerosol products and comparison reports are kept under
`validation/generated/aerosol_ir_comparison/before_emd_fix/` and
`…/aerosol_visible_ir_comparison/before_emd_fix/` as the record that the
matrix found the bug.

### Finding 1 — aerosol thermal emissivity scaled by 100 (found, corrected)

In both aerosol IR cases the Python `E_md` was 99.999997 × the legacy value at
every grid point (legacy range 0.284–1.000, Python 28.4–100.0), while `T_dd`,
`T_dv` and `R_dv` at 10.8 µm agreed to 1e-6–1e-7 and the cloud `E_md` agreed
to 1.5e-6. Tracing the convention: the legacy code accumulates
`Em = 100·UU/BBE` in both formulations and `write_v2_lut` divides by 100 on
output, so the legacy product stores the emissivity as a fraction; the Python
`E_md` assembly applies no scale to `arrays["em"]`, so the cloud block
(`UU/BBE`) matched and the aerosol block (`100·UU/BBE`, `pipeline.py` line 490)
was a factor 100 high with no compensation anywhere. This is invisible in any
visible-channel validation.

Correction (2026-09-14, the only source change): the aerosol emission
assignment in `calculate_aerosol_operators` was changed from
`emission["uu"][…] * np.float32(100.0) / np.float32(bbe)` to
`emission["uu"][…] / np.float32(bbe)`, making it identical to the cloud block.
No other operator, the cloud path, Mie, DISORT, gas, profile physics, SRF
handling or the writer was touched. A regression test compares the aerosol
channel-9 `E_md` directly with the captured legacy reference and a static
guard requires both emission blocks to end in `/ np.float32(bbe)` with no
percentage factor; both fail on the pre-correction implementation.

### After the correction

Python aerosol channel-9 and channel-1+9 products regenerated with the
existing validation configurations; legacy references untouched.

| Case | Structural | Numerical | `E_md` residual (max / RMS / mean, ratio) | Worst other operator |
| --- | --- | --- | --- | --- |
| aerosol IR (ch 9) | FAIL* | **GREEN** | 1.16e-5 / 4.55e-6 / 2.23e-6, ratio 1.000000 | 2.97e-5 (`R_dd`, noise-level values) |
| aerosol visible+IR (ch 1, 9) | FAIL* | **GREEN** (band) | as above for ch 9 | 9.4e-4 (`T_dv` ch 1 — Finding 2, unchanged) |

Cloud products were not regenerated and are unchanged (cloud IR re-compared:
PASS / GREEN, 6.6e-5).

\* The aerosol structural FAIL is NetCDF layout only: the `surface_pressure`
dimension is declared in a different position in the file, and the Python
`surface_pressure` coordinate lacks the legacy attributes `long_name`,
`spacing`, `units` and `valid_range`. Every variable's own dimension order,
shape, dtype and coordinate values match. Deliberately left for a separate
writer clean-up.

**Overall status: cloud GREEN (visible and IR); aerosol GREEN (visible, IR and
combined) after the two corrections this matrix found.** Both formulations now
reproduce the legacy radiative transfer on the compact grids for a solar and a
thermal channel at the kernel-residual level; the aerosol structural FAIL is
the NetCDF layout item only.

### Finding 2 — aerosol channel-1 bias at the smallest effective radius (found, corrected)

First pass: Python was 0.66–0.88 % low in `R_dd`, 0.3–0.4 % low in `T_dv` and
~1 % off in `R_0d`/`R_dv`/`T_0d` for r_eff = 0.01 µm at all three pressures,
within 5e-4 for r_eff = 0.1 µm, with `T_dd` agreeing to 1e-5 — a
scattering-partition effect, not a column-total effect. The single-scattering
optics were excluded as the cause (legacy `scatfile.sav` b_ext, ratio, SSA, g
and all 1000 Legendre moments agree with Python to 1e-7–1e-8).

**Per-layer diagnostic** (`validation/diagnostics/aerosol_layer_compare/`): for
the maximum-discrepancy state (ch 1, r_eff 0.01 µm, τ 1.0, 950 hPa) the actual
`dtau`/`ssa`/`pmom` arrays assembled for DISORT were captured from both codes.
Level heights, pressures, gas terms, Rayleigh terms, molecular and particle
moments were bitwise (or float32-rounding) identical. The **first material
difference was the vertical-profile interpolation** of the `.mm`
relative-amount profile onto layer mid-heights: at the 0.5 km layer legacy had
1085.23 and Python 778.80, giving normalised lowest-layer fractions of
46.4 / 33.3 / 20.2 % (legacy) versus 38.4 / 38.4 / 23.3 % (Python). A controlled
run of the production DISORT kernel from Python gave `R_dd` 0.0377878 on the
Python inputs, 0.0381495 on the legacy inputs and 0.0381610 with only the
profile weights swapped — the interpolation convention accounted for the whole
discrepancy. The earlier "different Rayleigh/absorber mixing" hypothesis was
not supported; the mixing arithmetic is identical.

**IDL `INTERPOL` semantics** (IDL 8.9 `lib/interpol.pro`, resolved by the
legacy run; verified by direct evaluation): piecewise linear inside the
tabulated range; the segment index is clamped to the end segments
(`s = VALUE_LOCATE(x, xout) > 0L < (m-2)`), so it **extrapolates linearly from
the two nearest nodes both below the first node and above the last**
(`INTERPOL([1,2,4,8],[1,2,3,4],[0,5])` → 0.0, 12.0 — the two slopes differ);
ascending and descending abscissae give identical results; two-point input
works; FLOAT input gives FLOAT output. Python's `_interpolate_profile` used
`np.interp`, which **clamps** to the end values outside the range.

**Correction** (2026-09-14): `_interpolate_profile` in `src/oraclut/pipeline.py`
now reproduces `INTERPOL` — sorted nodes, end-segment linear extrapolation on
both sides, single-precision arithmetic, and a `ValueError` for fewer than two
nodes or non-monotonic heights. It is the only production change. The routine
has exactly two callers, the cloud `_layer_inputs` and the aerosol operator
loop, matching the identical `INTERPOL(mmstr.rext, mmstr.height, hlayers)`
call in both legacy generators; the other `np.interp` uses (SRF/solar
integration, refractive-index interpolation) are different quantities and are
unchanged. For the liquid-water STG profile both end-node pairs are zero, so
extrapolation returns exactly the clamped value and cloud products are
unchanged by construction (asserted by a test).

**Post-correction layer arrays** (same state): `relative_tau_raw` and
`relative_tau` bitwise identical; `tauscat`/`dtau` max 3.0e-8; `ssa` 3.3e-7;
`pmom` 8.4e-9; Python kernel on Python inputs now 0.0381610 vs 0.0381495 on
legacy inputs (legacy LUT 0.038124). Pre-correction arrays and report kept in
`aerosol_layer_compare/before_profile_fix/`.

**Post-correction aerosol matrix** (Python candidates regenerated to explicit
validation paths; legacy references untouched; pre-correction products and
reports kept under each case's `before_profile_fix/`):

| Case | Structural | Numerical | Worst operator | Key operators |
| --- | --- | --- | --- | --- |
| aerosol visible (ch 1) | FAIL* | **GREEN** | 9.45e-5 (`R_0d`) | `R_dd` 8.9e-5 (was 3.4e-4), `T_dv` 2.5e-6 (was 9.4e-4), `R_0d` 9.5e-5 (was 7.8e-4); systematic-bias fractions all 0 |
| aerosol IR (ch 9) | FAIL* | **GREEN** | 2.25e-5 (`R_dd`, noise-level values) | `E_md` 6.0e-8 (was 1.2e-5), `T_dv` 1.2e-7 |
| aerosol visible+IR (ch 1, 9) | FAIL* | **GREEN** | 9.45e-5 | per-channel residuals identical to the single-channel cases |

The visible systematic bias has disappeared: every channel-1 operator now sits
at or below the cloud benchmark level (1.3e-4), with no localised structure.

### Warnings

Every legacy run reported one `% Program caused arithmetic error: Floating
underflow` (IDL, at exit) and produced a complete finite product; no
`UPBEAM--SGECO says matrix near singular` occurred, so no repeatability run
was triggered. Both codes return noise-level, partly negative `R_dd` at
10.8 µm for the optically thin aerosol (|R_dd| < 3e-5); that is common
behaviour, not a discrepancy.

## Reproducing

```shell
# Python products (tiny; interactive is fine)
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python -m oraclut.generate \
    --config configs/validation_meteosat-10_seviri_cloud_liquid-water_stg_test_ch01_ch09.driver
# legacy reference + comparison (licensed IDL shell; tiny; interactive)
scripts/generate_legacy_reference.sh --formulation cloud \
    --driver validation/reference_cases/meteosat10_seviri_liquid_water_stg_cloud_test_ch01_ch09.driver \
    --case-id cloud_visible_ir \
    --python-product validation/generated/cloud_visible_ir_comparison/python \
    --reference-dir validation/generated/cloud_visible_ir_comparison/legacy \
    --comparison-dir validation/generated/cloud_visible_ir_comparison
```

Anything larger than these grids goes through `scripts/submit_oraclut_lut.sh`
from `atmlxint7`.
