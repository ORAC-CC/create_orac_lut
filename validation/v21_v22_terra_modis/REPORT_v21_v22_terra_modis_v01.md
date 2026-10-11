# Terra MODIS ORAC LUT V21/V22 numerical comparison

## Executive conclusion

Follow-up: the active working implementation now uses an IDL-compatible
sequential `TOTAL` mirror at all corresponding call sites.  See
[the TOTAL compatibility report](../idl_total/REPORT_idl_total.md).  This makes
all 36 Terra MODIS centre wavelengths and wavenumbers bitwise identical to
V21; the sealed historical V22 `_v01` source and existing LUTs remain
unchanged.

The three newly generated V22 files are unambiguous counterparts of the V21
IDL products: they use the same instrument definition, microphysics file, LUT
grid, atmosphere, channel set, spectral-quadrature setting, Rayleigh setting,
DISORT stream count and phase-moment limit.  Their NetCDF dimensions, variable
schemas, data types, packing attributes, global and variable attributes, masks,
NaNs, infinities and decoded coordinate values agree exactly.

Their numerical fields do **not** agree exactly.  Most differences are small
single-precision effects, but sparse direct-beam radiative-transfer calls have
much larger errors.  The clearest failure is liquid-water 240 K, MODIS channel
6, optical depth 256, effective radius 29 um and solar zenith 80 degrees: V22
`R_0v` ranges from -0.473901 to 0.742851 where V21 ranges from 0.011446 to
0.454744, with a maximum absolute error of 0.541448.  The negative reflectance
is unphysical.  In the representative operational propagation used here, the
largest black-surface solar error is 847.9 times the production ORAC
measurement uncertainty.  V22 therefore does not faithfully reproduce V21,
despite close bulk statistics away from a sparse set of unstable calls.

There are two distinct findings:

1. Small, widespread changes begin before scattering or DISORT.  The V22
   centre wavelength and wavenumber are produced by NumPy's reduction tree;
   evaluating the same float32 terms sequentially, as IDL `TOTAL` does,
   reproduces every V21 centre value bit for bit.  This is a confirmed cause.
   The `F0` differences are also confined to float32 rounding scale and arise
   in a path where V21 uses double solar-spectrum arrays and IDL interpolation
   while V22 deliberately loads the spectrum as float32 and uses NumPy
   interpolation/reduction.  Tested combinations reduce but do not exactly
   reproduce V21, so the individual contribution of interpolation versus
   precision/reduction remains unresolved.
2. Small input perturbations expose the known isolated instabilities in the
   production single-precision DISORT 2.0 calculation at 60 streams.  An
   earlier targeted investigation independently reproduced four other failed
   calls: double precision at 60 streams, or single precision at 56, 68 or 100
   streams, removes those outliers and agrees, while only the production
   single-precision 60-stream result is anomalous.  The exact V22 calls found
   here were not recomputed under those alternatives, so assigning them to the
   same mechanism is strongly supported rather than independently confirmed.

No LUT, IDL reference source, or production generator was modified.  No LUT
was regenerated or repaired.

## 1. Files and provenance

All six LUTs were read in place under
`/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS`; that directory
was not modified.  Timestamps are UTC and hashes are SHA-256.

| Case | Version | Exact basename | Bytes | Timestamp | SHA-256 |
|---|---:|---|---:|---|---|
| liquid-water old (`pold`) | 21 | `terra_modis_m_liquid-water_a01_pold_v21.nc` | 41,318,050 | 2026-10-09 22:39:20.303 | `d65a05113d7ca9675154639648e0e9f011f75f21b53d6ce71abd087f94feeb45` |
| | 22 | `terra_modis_m_liquid-water_a01_pold_v22.nc` | 41,330,203 | 2026-10-10 04:33:37.745 | `7426ad01f68a65e69c7b552e1e1aeaea376d27e44a124baed7c7ff404ea84e76` |
| liquid-water 240 K (`p240`) | 21 | `terra_modis_m_liquid-water_a01_p240_v21.nc` | 41,318,050 | 2026-10-10 01:37:59.137 | `ef39386bc5c4b2f280f6419c3f548d9f6d5b8b5f9fa55c64425a783a5ebb0521` |
| | 22 | `terra_modis_m_liquid-water_a01_p240_v22.nc` | 41,330,203 | 2026-10-10 04:32:00.391 | `8055e7e1e9980f40a7c0f045d96366c148e1620e759c84fe7bf52eeaf106a69c` |
| water-ice aggregate (`pagg`) | 21 | `terra_modis_m_water-ice_a01_pagg_v21.nc` | 49,570,690 | 2026-10-10 04:41:18.371 | `5b423406e1f5f4f4ecbabde965d2d6dfc612f3bbb3781a70a130a97d6fe50630` |
| | 22 | `terra_modis_m_water-ice_a01_pagg_v22.nc` | 49,582,811 | 2026-10-10 04:55:23.750 | `a07735f88bc5c95a5c1e859537f323dfb0ac512ca5907148c305f3d05f304fdd` |

The V22 submission log identifies jobs 490632--490634.  The products' sidecar
provenance records source commit `aa62864`; submission was prepared from
`4fee1b4`.  The sealed [V22 release manifest](../../releases/lut_v22.release)
identifies release `lut-code-v22.1`, every versioned Python source, the shared
Fortran/C kernels, and the exact compiled Mie and DISORT libraries by hash.

## 2. Configuration and NetCDF structure

The corresponding inputs are:

| Setting | `pold` | `p240` | `pagg` |
|---|---|---|---|
| instrument | `terra_modis_v1.inst` | same | same |
| microphysics | `liquid-water_old.mm` | `liquid-water_240.mm` | `water-ice_agg.mm` |
| LUT grid | `liquid-water-cloud.lut` | same | `ice-cloud.lut` |
| atmosphere | 2 | 2 | 2 |
| channels | 1--36 | 1--36 | 1--36 |
| gas / Rayleigh | off / on | off / on | off / on |
| `srf_quad` | 1 | 1 | 1 |
| DISORT streams / moments | 60 / 1000 | 60 / 1000 | 60 / 1000 |

Input hashes recorded during discovery are:

- instrument: `8b1f1adcac24e8b066a5735edbffedd2f942157196ad6dc09c5ef9c5e4693d05`;
- liquid grid: `5b233f8f7709bd7bf21bcb4800cf33812c4f5fa15ba9c86597f9156dec54b612`;
- ice grid: `3142e32fcd77279c5189b109385c239214b8b0d72c1d8f1d82fb30c30e3989bd`;
- `pold` microphysics: `0102aec5918ba60ecd44ccb06f07ea7cace849c420ec62d5397509a10f46ffaf`;
- `p240` microphysics: `8145e8c5be0cb565a727eecb827d5def620205dcf9f39b15e722f1fe0269cbf2`;
- `pagg` microphysics: `9cbbcb55743636f89d672198bb84d13b8e0998699db70bba73009e9eb78e5f24`;
- legacy cloud driver: `877e893fe142b0620234cc8c53c7e11a33a75ed0c40dc20fc84ec157e4835e12`.

Full V22 run-file and source hashes are retained in the product provenance
sidecars and release manifest.

The V21 and V22 files have identical dimension names and lengths.  Both liquid
cases have 17 optical depths, 20 radii, 10 view zeniths, 10 solar zeniths, 11
relative azimuths and 36 channels; the ice case has 24 radii.  The spectral
subsets contain 26 solar, 16 thermal and 6 mixed channels.  The decoded grids
are bitwise identical:

- optical depth: `1e-10, 0.0078125, ..., 128, 256` (17 points);
- liquid radius: 1--39 um in steps of 2; ice radius: 1--93 um in steps of 4;
- view zenith: 0--65 degrees (10 points); solar zenith: 0--90 degrees in
  10-degree steps; relative azimuth: 0--180 degrees in 18-degree steps.

There are no global-attribute or variable-attribute differences, no dimension
or variable-schema differences, no dtype differences, and no packing,
fill-mask, NaN or infinity differences.  No numerical variable uses a differing
scale factor or offset.  Therefore the larger V22 file sizes are HDF5/NetCDF
representation differences, not evidence of different physical coordinates or
schemas.

## 3. Numerical comparison

The utility first tests exact equality and then computes absolute, RMS,
reference-relative and symmetric-relative residuals over jointly finite,
unmasked decoded values.  Its relative statistic excludes reference values
within `max(1e-12, 1e-6 * max(abs(reference)))`; near-zero values are counted
and summarized separately.  This prevents division by nearly zero from being
mistaken for physical significance without hiding those residuals.

All numerical variables not listed below are exactly equal.  Full results for
every variable, including element count, difference count, maximum absolute
and relative difference, RMS, near-zero statistics and the exact coordinate of
the maximum, are in the linked CSV/JSON files.

### `pold`

| Variable | Elements different | Maximum absolute | Maximum relative | RMS |
|---|---:|---:|---:|---:|
| `central_wavelength` | 34/36 | 7.62939e-6 um | 1.12985e-6 | 2.22470e-6 |
| `central_wavenumber` | 28/36 | 5.85938e-3 cm-1 | 1.28733e-6 | 2.50006e-3 |
| `F0` | 22/26 | 3.05176e-4 | 6.86295e-7 | 9.25422e-5 |
| `average_volume_per_particle` | 3/20 | 1.95312e-3 | 9.99159e-8 | 4.88284e-4 |
| `extinction_coefficient` | 541/720 | 3.12500e-2 | 1.27104e-5 | 4.09891e-3 |
| `extinction_coefficient_ratio` | 671/720 | 1.23978e-5 | 1.24310e-5 | 1.43720e-6 |
| `single_scatter_albedo` | 409/720 | 1.66893e-6 | 2.11243e-6 | 2.29788e-7 |
| `asymmetry_parameter` | 537/720 | 1.10269e-5 | 1.26713e-5 | 1.21805e-6 |
| `R_0v` | 9,713,200/9,724,000 | 3.72261e-2 | 11.7900 | 4.81426e-5 |
| `R_0d`, `T_0d`, `T_00` | 88,339; 85,619; 55,420 | 3.41159e-3; 1.31256e-3; 3.93391e-6 | see CSV | 6.94799e-5; 4.36177e-5; 2.04734e-7 |
| diffuse/thermal operators (`R_dd`, `R_dv`, `T_dd`, `T_dv`, `E_md`) | 12,236--122,294 | 6.56396e-6--8.68559e-4 | see CSV | 2.30368e-7--6.95898e-5 |

Detailed: [variables CSV](results/terra_modis_pold_v21_v22_variables.csv),
[full JSON](results/terra_modis_pold_v21_v22_comparison.json),
[coordinate profiles](results/terra_modis_pold_v21_v22_coordinate_dependence.csv),
[ranked direct calls](results/terra_modis_pold_v21_v22_direct_beam_calls.csv).

### `p240`

| Variable | Elements different | Maximum absolute | Maximum relative | RMS |
|---|---:|---:|---:|---:|
| spectral centres / `F0` | same counts as `pold` | 7.62939e-6; 5.85938e-3; 3.05176e-4 | <=1.28734e-6 | 2.22470e-6; 2.50006e-3; 9.25422e-5 |
| `average_volume_per_particle` | 3/20 | 1.95312e-3 | 9.99159e-8 | 4.88284e-4 |
| optical properties | 458--607/720 | 2.68221e-6--2.05078e-2 | <=7.03335e-6 | 3.29725e-7--2.63731e-3 |
| `R_0v` | 9,712,123/9,724,000 | **0.541448** | 236.834 | 1.05252e-3 |
| `R_0d` | 88,337/88,400 | **1.20160** | 58.4054 | 4.47667e-3 |
| `T_0d` | 85,604/88,400 | 8.33879e-3 | 54.8824 | 5.75495e-5 |
| other seven RT operators | 47,993--122,291 | 2.44379e-6--1.08105e-3 | see CSV | 1.31984e-7--6.95781e-5 |

Detailed: [variables CSV](results/terra_modis_p240_v21_v22_variables.csv),
[full JSON](results/terra_modis_p240_v21_v22_comparison.json),
[coordinate profiles](results/terra_modis_p240_v21_v22_coordinate_dependence.csv),
[ranked direct calls](results/terra_modis_p240_v21_v22_direct_beam_calls.csv).

### `pagg`

| Variable | Elements different | Maximum absolute | Maximum relative | RMS |
|---|---:|---:|---:|---:|
| spectral centres / `F0` | same counts as `pold` | 7.62939e-6; 5.85938e-3; 3.05176e-4 | <=1.28734e-6 | 2.22470e-6; 2.50006e-3; 9.25422e-5 |
| optical properties | 401--484/864 | 5.36442e-7--1.66893e-6 | <=1.24715e-6 | 1.18916e-7--3.07557e-7 |
| `R_0v` | 11,649,796/11,668,800 | 5.18614e-3 | 28.7170 | 1.57452e-5 |
| `R_0d` | 105,942/106,080 | 1.28968e-2 | 132.988 | 7.08552e-5 |
| `T_0d` | 100,286/106,080 | **5.28509e-2** | 91.3333 | 1.69259e-4 |
| other six RT operators | 13,868--146,697 | 4.47035e-7--1.13983e-3 | see CSV | 5.67089e-8--6.65919e-5 |
| `E_md` | 56,352/65,280 | 4.93228e-6 | 2.80978e-4 | 1.62139e-7 |

Detailed: [variables CSV](results/terra_modis_pagg_v21_v22_variables.csv),
[full JSON](results/terra_modis_pagg_v21_v22_comparison.json),
[coordinate profiles](results/terra_modis_pagg_v21_v22_coordinate_dependence.csv),
[ranked direct calls](results/terra_modis_pagg_v21_v22_direct_beam_calls.csv).

The large quoted relative maxima occur where the V21 value is small; they are
not used alone to judge importance.  Absolute residuals, complete coordinate
locations, symmetric relative residuals and explicit near-zero subsets remain
available in the machine-readable results.

## 4. Coordinate structure and plots

Residual profiles were computed along every dimension of every differing
variable.  For `R_0v`, calls were also ranked after collapsing only the two
view-angle dimensions; this distinguishes a whole unstable DISORT call from an
isolated array element.

The leading direct-call discrepancies are:

| Case | Channel | tau | radius (um) | SZA | max abs `R_0v` | call RMS |
|---|---:|---:|---:|---:|---:|---:|
| `pold` | 13 | 8 | 5 | 40 | 0.0372261 | 0.0117921 |
| `pold` | 16 | 2 | 11 | 60 | 0.0135524 | 0.00380141 |
| `p240` | 6 | 256 | 29 | 80 | **0.541448** | **0.297720** |
| `p240` | 6 | 128 | 29 | 80 | 0.130887 | 0.0749833 |
| `pagg` | 9 | 16 | 93 | 50 | 0.00518614 | 0.00174077 |
| `pagg` | 11 | 0.015625 | 33 | 70 | 0.00304554 | 0.00158190 |

The 33 generated figures are directly under the project [plots](../../plots/)
directory: one variable overview, nine RT-operator coordinate profiles and one
radiance-impact plot for each case.  Useful entry points are:

- `pold`: [overview](../../plots/terra_modis_pold_v21_v22_variable_overview.png),
  [`R_0v` coordinates](../../plots/terra_modis_pold_v21_v22_R_0v_coordinates.png),
  [radiance impact](../../plots/terra_modis_pold_v21_v22_radiance_impact.png);
- `p240`: [overview](../../plots/terra_modis_p240_v21_v22_variable_overview.png),
  [`R_0v` coordinates](../../plots/terra_modis_p240_v21_v22_R_0v_coordinates.png),
  [radiance impact](../../plots/terra_modis_p240_v21_v22_radiance_impact.png);
- `pagg`: [overview](../../plots/terra_modis_pagg_v21_v22_variable_overview.png),
  [`R_0v` coordinates](../../plots/terra_modis_pagg_v21_v22_R_0v_coordinates.png),
  [radiance impact](../../plots/terra_modis_pagg_v21_v22_radiance_impact.png).

## 5. Source tracing and earliest divergence

The V21 cloud path is
`create_orac_cloud_lut_wrapper.pro -> create_orac_cloud_lut.pro ->
generate_scattering_properties.pro -> create_bwgp.pro` (Mie) or the Baum
reader, followed by `setup_disort.pro -> call_disort.pro ->` the Fortran
DISORT library and `write_v2_lut.pro`.  The sealed V22 path is
`create_orac_luts_v01.py` plus the corresponding
`src/oraclut/idl_mirror/*_v01.py` modules.  Its manifest pins the same legacy
Mie and DISORT kernel sources and exact shared libraries.

The first difference is earlier than either Fortran kernel.  Running the
sealed `load_srfstrarr_v01.py` logic from the common SRFs reproduces V22 centre
wavelengths, centre wavenumbers and `F0`.  For a centre, IDL evaluates the
float32 trapezoid terms and performs a sequential `TOTAL`.  V22 evaluates the
same float32 terms but `np.sum(..., dtype=np.float32)` uses a different
reduction order.  A targeted sequential float32 sum reproduces all V21 centre
wavelengths and wavenumbers exactly.  Thus the wavelength passed into the
single-point (`srf_quad=1`) scattering calculation already differs.

For `F0`, the implementations additionally differ in solar-spectrum handling:
IDL's `loadxy/loadxyz` path yields double arrays, whereas V22 explicitly reads
the table into float32; IDL `INTERPOL` and NumPy `interp` do not necessarily
round each arithmetic stage identically.  Sequential/dtype diagnostic
combinations reduced the maximum error but none reproduced all 26 V21 values.
It is therefore confirmed that the divergence is within this input,
interpolation and reduction path, but it is not yet justified to assign the
full residual to one operation.

The liquid Mie path subsequently receives the changed centre wavelength and
performs bulk size-distribution reductions in Python/NumPy, so differences in
volume and single-scattering quantities are expected at approximately 1e-5
relative.  The aggregate-ice optical properties are tabulated/interpolated and
show smaller residuals (approximately 1e-6 relative).  These statements locate
the propagation path; they do not claim that every last-bit optics residual has
been decomposed operation by operation.

The atmospheric selection, gas flag, Rayleigh flag, surface-free LUT
operators, 60-stream setting, 1000-moment limit and Fortran library identity
were checked and do not explain the discrepancy.  NetCDF writing is also
downstream of the first divergence, and identical dtypes/packing plus the
structured coordinate dependence rule out output conversion as the cause.

As an independent reproducibility check, same-named Terra MODIS V11 IDL
products were compared with V21.  Their grids and optical/RT variables are
bitwise equal for `pagg`; `pold` has only three `R_0v` last-bit changes
(maximum 1.819e-12), and `p240` only one (maximum 9.095e-13).  `F0` is not
equal because the historical solar metadata changed by as much as 637.324, so
the files as wholes are not identical.  Nevertheless, the independent
optical/RT agreement shows that the sparse V22 RT failures are not ordinary
V21 run-to-run variability.  These supplemental results are beside the primary
files as `*_v11_v21_*`; V21 remains the authoritative baseline.

## 6. ORAC radiance/measurement impact

The propagation uses the validated ORAC equation-form-3 mirror and the
McGarragh et al. surface multiple-scattering coupling.  It evaluates:

- black-surface solar reflectance over every LUT geometry;
- Lambertian albedos 0.05, 0.2 and 0.8 at nadir, where `T_vv` is an exact LUT
  node;
- thermal and mixed-channel brightness temperatures using the repository
  mid-latitude-summer atmosphere and MODTRAN channel gas profiles.

The production uncertainty includes MODIS instrument noise and the ORAC 1%
homogeneity and 2% coregistration terms.  Solar ratios exclude reference
reflectance below 0.01 to avoid meaningless ratios at zero signal.  Results
are representative calculations, not a frequency-weighted global scene
distribution.

| Case | Measurement giving largest impact | Channel | Maximum absolute impact | Maximum production sigma | 99th-percentile sigma | Fraction >1 sigma |
|---|---|---:|---:|---:|---:|---:|
| `pold` | black-surface solar reflectance | 16 | 0.0135524 | 5.3046 | 0.1208 | 1.31e-4 |
| `pold` | mixed brightness temperature | 20 | 0.2758 K | 0.3286 | lower than maximum | 0 |
| `p240` | black-surface solar reflectance | 6 | 0.541448 | **847.915** | 0.1491 | 3.27e-3 |
| `p240` | mixed brightness temperature | 20 | 15.4296 K | **22.348** | 0.0371 | 2.27e-3 |
| `pagg` | black-surface solar reflectance | 11 | 0.00304554 | 3.061 | 0.0333 | 2.15e-4 |
| `pagg` | mixed brightness temperature | 20 | 2.7663 K | 2.122 | 0.0162 | sparse |

Pure thermal-channel impacts are negligible in these tests: maxima are about
0.00165--0.00190 sigma and below 0.001 K.  The problem is therefore primarily
solar/direct-beam, with mixed thermal channels affected through their solar
contribution.  Nadir Lambertian maxima are 0.268 sigma (`pold`), 356.2 sigma
(`p240`, albedo 0.05) and 2.040 sigma (`pagg`, albedo 0.8).  Nadir is not where
the worst `pold R_0v` view occurs, explaining its smaller Lambertian result.

All channels, uncertainty variants and maximum indices are in:
[pold CSV](results/terra_modis_pold_v21_v22_radiance_impact.csv),
[p240 CSV](results/terra_modis_p240_v21_v22_radiance_impact.csv), and
[pagg CSV](results/terra_modis_pagg_v21_v22_radiance_impact.csv).

## 7. Confirmed causes, hypotheses and unresolved questions

Confirmed:

- V21/V22 scientific configuration and coordinate grids match.
- Spectral centres first diverge in Python versus IDL float32 reduction order;
  sequential float32 summation exactly restores the V21 centres.
- These changed centres feed `srf_quad=1` scattering, before Mie/Baum and
  DISORT.
- V22 contains an unphysical `p240` direct-beam call with operationally large
  radiance impact.
- Prior focused DISORT experiments show that other failures with this
  whole-call signature are single-precision, 60-stream numerical outliers, not
  robust physical solutions.  The unphysical V22 call itself was not rerun.

Strongly supported but not individually reconfirmed here:

- The sparse, call-shaped V22 direct-beam failures, including the worst `p240`
  call, are manifestations of the same DISORT numerical sensitivity.  Their
  structure matches independently reproduced failures, but these exact calls
  were not rerun across precision and stream count during this read-only
  study.

Unresolved:

- The exact stage-by-stage decomposition of `F0` last-bit differences among
  table input precision, IDL versus NumPy interpolation arithmetic, product
  rounding and summation order.
- The contribution of each NumPy reduction in the Mie bulk average to the
  small optical-property residuals after the wavelength difference is removed.
- Whether every sparse RT outlier disappears under the same precision/stream
  alternatives as the confirmed `p240` failure.

## 8. Recommendations

1. Do not accept these V22 files as numerically faithful replacements for V21,
   and specifically do not use the affected calls for operational retrievals.
2. In a new scientific release, not by editing sealed V22, reproduce IDL
   sequential float32 reductions where bitwise V21 reproduction is required.
   Add regression tests for all 36 spectral centres and all 26 `F0` values.
3. Isolate the solar-spectrum path with saved IDL intermediates immediately
   after `loadxy`, `INTERPOL`, elementwise multiplication and each `TOTAL`.
   Those intermediates are not present in the LUTs and are the only essential
   missing evidence for exact `F0` attribution.
4. Add per-DISORT-call physical guards (finite results, plausible
   reflectance/transmittance bounds and energy checks) that fail loudly rather
   than writing an anomalous LUT.  Validate any correction against V21 and
   stable precision/stream variants; do not merely clip values.
5. Evaluate the already demonstrated double-precision solver or another
   validated stable solver configuration for the next release, with targeted
   regressions at every ranked outlier coordinate and representative ordinary
   points.  Changing solver precision/settings is a scientific change and
   must not be backported silently into V22.
6. Use the comparison and radiance utilities as release gates for other
   instruments and particle models before production LUTs are promoted.

## 9. Reproduction

The generic [comparison utility](compare_luts.py) accepts explicit reference
and candidate paths, aligns variables by dimension name, decodes packing,
records masks and schemas, writes JSON/CSV diagnostics and produces coordinate
plots.  The [radiance utility](radiance_impact.py) separately performs the
ORAC measurement propagation.  Example:

```shell
python validation/v21_v22_terra_modis/compare_luts.py \
  /path/to/reference_v21.nc /path/to/candidate_v22.nc \
  --output-directory validation/v21_v22_terra_modis/results \
  --plot-directory plots --label instrument_case_v21_v22

python validation/v21_v22_terra_modis/radiance_impact.py \
  /path/to/reference_v21.nc /path/to/candidate_v22.nc \
  --cloud liquid \
  --output-directory validation/v21_v22_terra_modis/results \
  --plot-directory plots --label instrument_case_v21_v22
```

The synthetic [regression tests](../../tests/test_v21_v22_lut_compare.py)
exercise named-dimension transposition, fill-mask handling, coordinate
locations and mask differences.  Comparisons are against decoded physical
values; file hashes are recorded only for provenance, never used as a proxy
for numerical equality.
