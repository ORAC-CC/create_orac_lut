# Python legacy-reproduction pipeline

The selected reproduction target is the current Meteosat-10 SEVIRI liquid-water
STG cloud LUT: ORAC LUT level 2, product revision v21. The Python generator
keeps the legacy orchestration and output ordering, but does not require IDL at
runtime.

## Pipeline

```text
current ORAC text inputs
  -> src/oraclut/config
  -> legacy Mie kernel (mie/dlm-code/mieint.f)
  -> STG distribution and LEGPEXP-compatible moments
  -> production DISORT kernel (create_orac_lut/disort2/src)
  -> current cloud setup_disort/call_disort semantics
  -> V2 NetCDF writer
```

`src/oraclut/optics/legacy_mie.py` and
`src/oraclut/radiative_transfer/legacy_disort.py` build repository-local shared
libraries from the preserved Fortran kernels and small C ABI adapters. The
default compiler is GNU Fortran 12.2.0 with `-O3 -fPIC
-ffixed-line-length-0` for Mie and `-O2 -fPIC -std=legacy
-ffixed-line-length-0` for DISORT. The temporary Intel variant is selected
with `ORACLUT_FORTRAN=ifort` after loading `intel-compilers/2022` (ifort
2021.6.0, using `-O3 -fPIC -extend_source -D__GFORTRAN__`); neither choice
changes the preserved scientific sources.

The optical-property boundary is intentionally independent of the solver and
is the future insertion point for a separately validated POM backend. POM is
not used by this reproduction.

## Commands

The safe default is the compact validation grid and channel 1:

```shell
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python -m oraclut.generate --repository-root . \
  --output validation/generated/python_liquid_water_cloud_test_v21.nc
```

The full current Meteosat-10 reference configuration is explicit:

```shell
nice -n 19 env ORACLUT_FORTRAN=ifort PYTHONPATH=src \
  /home/g/grainger/miniforge3/envs/science/bin/python -m oraclut.generate \
  --repository-root . \
  --lut-file create_orac_lut/input_files/lut/liquid-water-cloud.lut \
  --microphysics-file create_orac_lut/input_files/microphysics/liquid-water_stg.mm \
  --channels all --atmosphere 2 --srf-quad 1 \
  --output validation/generated/python_meteosat10_seviri_liquid_water_stg_cloud_v21_ifort.nc
```

This is a low-priority interactive development command. Future production
runs should follow the AOPP policy in `AGENTS.md`; the measured full-case
runtime determines whether SLURM submission from `atmlxint7` is preferable.

## Validation and caveats

The captured tiny reference is in
`validation/generated/meteosat10_seviri_liquid_water_stg_cloud_test/`.
`/home/g/grainger/miniforge3/envs/science/bin/python -m oraclut.validation.compare`
reports dimensions, ordering, dtypes,
finite ranges, NaN/Inf and fill masks, absolute/RMS/relative residuals, and
the indices of the largest absolute and relative residuals. It separately
reports structural equivalence, typed metadata differences, and SHA-256/byte
identity; NetCDF byte identity is not required.

The tiny case uses the current `srf_quad=1` centre-wavelength path and matches
the captured legacy operators at single-precision/compiler scale: the largest
operator absolute residual is `1.2622e-4`. The captured tiny file contains
fill values for all four microphysical optical-property arrays, whereas the
Python legacy-Mie path computes finite values for them; this is recorded as a
12-element fill-mask discrepancy rather than replacing valid optics with
placeholders. Dimensions, variable order, dtypes, and attributes are exact.
The default gfortran build can take a different numerical branch in difficult
full-case DISORT `UPBEAM--SGECO` near-singular solves. The temporary Intel 2022
build is therefore retained as one full-reference validation build. The
complete production run was performed on allowed node `atmlxint5` with
`nice -n 19` and took longer than 30 minutes. Its 17x20x10x10x11, 11-channel
output is structurally exact, but the preserved DISORT branch leaves a largest
full-case operator residual of `9.7289e-2` absolute (`R_0v`, RMS `3.6373e-4`);
this is not hidden behind a broad tolerance. The `UPBEAM--SGECO` warning and
complete per-variable residuals are retained in
`validation/generated/python_full_comparison.json`.

The writer preserves the historical percentage conversion for radiative
operators, the unscaled `E_md = UU/PLKAVG` result after the legacy writer’s
conversion, and the historical packed calibration arrays. LUT level/schema
and product revision remain separate: this is a V2 LUT at revision v21, not a
“V21” format.

## Forensic discrepancy analysis

The compact reference has one selected channel and three radii. Its four
microphysical variables `extinction_coefficient`,
`extinction_coefficient_ratio`, `single_scatter_albedo`, and
`asymmetry_parameter` contain the NetCDF fill value at all 12 elements,
while `average_volume_per_particle` is finite. The full production reference
contains finite values for all five quantities. `generate_scattering_properties.pro`
allocates and fills `Bext`, `BextRat`, `w`, and `g` for every channel/radius,
and `create_orac_cloud_lut.pro` passes their centre spectral slices to
`write_v2_lut.pro`; there is no source branch intentionally omitting these
values. The writer declares the arrays with a channel dimension and passes
direct `reform` results. The evidence therefore classifies the 12 fills as a
one-channel legacy NetCDF writer artifact (singleton channel/dimension
handling), not meaningful missing optics. Python continues to write the
finite calculated values; it does not mask them. A live IDL restore/write
reproduction could not be completed in this shell because the IDL licence was
unavailable.

The full `R_0v` forensic summary is stored in
`validation/diagnostics/full_r0v_residual_summary.json`. The worst point is
channel 1 (solar-channel index 0, 0.6381768 microns), optical depth 2.0
(index 9), effective radius 23 microns (index 11), solar zenith 50 degrees
(index 5), satellite zenith 80 degrees (index 8), and relative azimuth 0
degrees (index 0). The reference/Python values are 0.25314441/0.15585500,
an absolute residual of 0.09728941 (relative residual 0.3843). Forty-eight
of the top 50 points share this channel, optical-depth, radius, and solar
zenith state. The 110-point suspected cluster has RMS residual 0.041996;
outside it, the maximum is 0.013135 and RMS is `5.1170e-5`.

The exact pre-call Python state is machine-readable in
`validation/diagnostics/worst_r0v_disort_state_ifort.json` and its NPZ/raw
companions: 45 layers, 60 streams, highest moment 999 (1000 coefficients),
two `UTAU` levels, 20 `UMU` values, 11 `PHI` values, float32 arrays, solar
beam 100, `UMU0=cos(50 degrees)`, `FISOT=0`, Lambertian albedo 0, and the
non-Planck wrapper flags from `call_disort.pro`. The gfortran and ifort
pre-call arrays are bitwise identical. The targeted ifort call emits the
known `UPBEAM--SGECO says matrix near singular` warning; the Python target
value is 0.15590107 when the state is captured after the production optics
sequence. Repeated calls after the state is established are bitwise
repeatable. The preserved DISORT source contains multiple `SAVE` variables,
so this behaviour is recorded rather than “fixed”.

For identical captured inputs, the ifort and gfortran targets are 0.15590107
and 0.05717356 respectively; the largest `UU` difference is 4.6931 and the
RMS is 1.1147. A fresh-process serialized Python state gives 0.28200659. The
user subsequently ran the same targeted state through the licensed IDL 8.9
production `DISORTTOIDL` DLM and obtained `R_0v = 0.288875`, with
`RFldir = 64.2788, 2.61441`, `RFldn = -0.000579834, 48.5968`, and
`Flup = 13.0656, 0.000335587`. That actual IDL call emitted
`UPBEAM--SGECO says matrix near singular`, confirming compiler,
process/workspace, and IDL-DLM sensitivity at this ill-conditioned state. The
archived value `0.25314441` is therefore not an appropriate bit-for-bit
acceptance target. The warning is retained as historical numerical sensitivity;
DISORT streams, tolerances, and solver behaviour were not changed.
