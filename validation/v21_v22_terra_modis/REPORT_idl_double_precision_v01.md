# Double-precision DISORT in the V21 IDL LUT generator

## Summary

A double-precision DISORT DLM (`DISORTtoIDL_DP`) has been built from the
unmodified production DISORT sources and validated through the actual IDL
LUT generator on four representative Terra MODIS `p240` calculations.  All
four pass every criterion set for the decision gate:

1. **It eliminates the single-precision failures.**  The historical V21
   kernel is in error by 36.5, 70.2 and 17.5 times the production MODIS
   uncertainty at the three known failure nodes (A, B, C), including
   reflected flux 132 % of incident, negative intensities down to -26.8 and
   negative LUT reflectances.  The double-precision results have none of
   these defects.
2. **It is physically consistent.**  No negative intensity or flux in any of
   the 90 double-precision calls at any stream count; top-of-atmosphere
   downward diffuse flux at most 1.6e-7 (exact value 0; single precision up
   to 0.67); absorbed and reflected fractions in [0, 1].
3. **It converges with stream count.**  Across NSTR 48, 60, 80 and 100 the
   LUT coefficients agree to at most 8.5e-4 (absolute, `R_0v`), at most
   0.08 sigma at the V21 stream count of 60 and 0.11 sigma at 48.
4. **It is stable and repeatable.**  Repeat calls are bitwise identical; the
   IDL DLM is bitwise identical to the same kernel called directly.
5. **It agrees with independent double- and quadruple-precision builds.**
   gfortran double precision and a full REAL*16 gfortran build agree with
   the IDL result to at most 3.1e-8 in every LUT coefficient of every call.
6. **It retains the V21 physical formulation.**  The compact single-precision
   runs reproduce the production V21 LUT bitwise in every RT coefficient; only
   the solver precision changes.

The double-precision generator backend is prepared and has produced small
end-to-end LUTs (version label 90) that match the validated solver output to
float32 rounding.  It costs about 1.3 times the historical solver at NSTR 60;
a full Terra MODIS product is expected to take about 3.6 h (liquid) or 4 h
(ice) against about 3 h historically.

**Ready for full LUT testing**, subject to your approval to run
`idl_double/run_dp_generator_full.sh`.  Nothing has been promoted: V21 sources
and LUTs are untouched, and no output carries a V21 label.

## 1. Design

### 1.1 The historical implementation

| Component | Precision and role |
|---|---|
| `create_orac_lut/disort2/src/DISORTfunctions.f` | DISORT 2.1 (RCS 2000-04-04) with ORAC additions (`GETMOM`, `PLKAVG`, integer error returns) and enlarged limits `MXCLY = 46`, `MXCMU = 100`, `MXUMU = 180`, `MXPHI = 181`.  Default `REAL` (4 bytes); `ASYMTX` eigen-analysis already in `DOUBLE PRECISION`, the linear solves (`SGBCO`/`SGECO`) single |
| `disort/LINPAK.f`, `ErrPack.f`, `RDI1MACH.f`, `PRTFIN.f`, `BDREF.f` | LINPACK, messages, machine constants (from Fortran intrinsics, so they follow the promoted kind), test printing, BRDF |
| `DISORTwrapper.c`, `DISORTtoIDL.c` | IDL interface; requires FLOAT arrays, reads scalars as float/long, forces `MAXCLY/MAXULV/MAXUMU/MAXPHI/MAXMOM` to `NLYR/NTAU/NUMU/NPHI/NMOM` |
| `DISORTtoIDL.so` | historical binary imported 2026-09-13, built on Ubuntu 18.04 with an unrecorded ifort version (`-O3 -fpic -D__GFORTRAN__`) |

Observations relevant to a promoted-precision build:

- `BDREF` calls the C routine `ambralsfortran_` with `float *` arguments,
  which would be mis-passed once `REAL` is promoted.  ORAC always uses
  `LAMBER = .TRUE.`, so the BRDF path is never reached.
- A missing comma in `EXTERNAL PLKAVG, R1MACH, RATIO, CHEKIN`
  (`DISORTfunctions.f:607`) creates a phantom external `CHEKINALBTRN`; it is
  never called and is tolerated by IDL's lazy binding in both DLMs.
- DISORT overwrites `SSALB = 1` with `1 - DITHER` in the caller's array, and
  its first-call self-test temporarily overwrites scalar arguments.
- `DISORTfunctions.h` declares `disort_` as returning `short`; the Fortran is
  `INTEGER`.
- The historical wrapper reads `PRNT` as 32-bit integers, while V21 passes an
  IDL INT (16-bit) array.

### 1.2 The double-precision DLM

`validation/v21_v22_terra_modis/idl_double/`:

| File | Purpose |
|---|---|
| `src/DISORTtoIDL_DP.c` | new IDL interface, procedure `DISORT_DP` (same 44 positional arguments as `DISORT`) and function `DISORT_DP_INFO()` |
| `src/dp_support.f` | `DPREALBYTES` (kernel REAL size; the DLM refuses to load unless 8) and a fatal `BDREF` stub replacing `BDREF.f` |
| `src/DISORTtoIDL_DP.dlm`, `src/exports.map` | DLM descriptor; linker map exporting only `IDL_Load` |
| `build_dlm.sh` | reproducible build, writes `build/` and `build/BUILD_INFO.txt` |

Interface rules (`DISORTtoIDL_DP.c`):

- every REAL scalar is read with `IDL_DoubleScalar`; every REAL array of any
  numeric IDL type is converted element by element to double (float to
  double is exact, so V21 float32 inputs reach the kernel unchanged);
- integers via `IDL_LongScalar` into 32-bit `int32_t`; logical arguments
  normalised to 0/1; `PRNT` converted from any numeric type;
- array lengths checked against `NLYR`, `NMOM+1`, `NTAU`, `NUMU`, `NPHI`; the
  `PMOM` leading dimension must equal `NMOM+1`; `SSALB` must lie in [0, 1];
  `TEMPER` copied into a zero-padded `NLYR+1` buffer (required when `PLANK`);
- `LAMBER`, `USRTAU` and `USRANG` must be true (as in every ORAC call);
- `MAXCLY` etc. set to the actual sizes exactly as the historical wrapper;
- all eight outputs returned as IDL DOUBLE arrays; a non-zero kernel status
  is an IDL error.

Fortran ABI: default `REAL` 8 bytes (`-real-size 64`, which also promotes
real literal constants and intrinsics), `DOUBLE PRECISION` 8 bytes
(`-double-size 64`), `INTEGER` and default `LOGICAL` 4 bytes
(`-integer-size 32`).  The only COMMON block in the kernel sources
(`DOCHEK`, self-test printing) is private to the DLM.  Because the DLM
exports only `IDL_Load` and is linked `-Bsymbolic`, its kernel symbols cannot
collide with the historical `DISORTtoIDL.so`, which remains loaded in the
same session for `GETMOM` (Rayleigh moments, exact in single precision).

### 1.3 Build provenance

| Item | Value |
|---|---|
| Compiler | ifort (IFORT) 2021.6.0 20220226, module `intel-compilers/2022` |
| Fortran flags | `-O2 -fp-model precise -fpic -real-size 64 -double-size 64 -integer-size 32` |
| C compiler / flags | gcc 12.2.0 (Debian 12.2.0-14+deb12u1), `-O2 -fPIC -std=gnu99 -Wall -Wextra -Wno-cast-function-type`, IDL 8.9 `idl_export.h` |
| Kernel sources (SHA-256) | `DISORTfunctions.f` 3f45489e…5f72, `PRTFIN.f` d3846c63…0770, `ErrPack.f` d1b90b20…baaf, `LINPAK.f` 8f0ec48b…96a0, `RDI1MACH.f` c59dadbf…3738 (all unmodified production files) |
| `DISORTtoIDL_DP.so` | 3742858aa0dcd6e2a03dbc9f008931a5e5d6c0261da16d27007905c63b37d2f9 |
| Exported symbols | `IDL_Load` only |

`-fp-model precise` was chosen deliberately for a reference implementation,
in place of the historical `-O3` default floating-point model.

### 1.4 Independent references

`build_reference.sh` builds the same unmodified sources with gfortran 12.2.0
in double (`-fdefault-real-8 -fdefault-double-8`) and full quadruple precision
(`-freal-4-real-16 -freal-8-real-16`, REAL and DOUBLE PRECISION both 16
bytes), each linked to `src/disort_ref_driver.f90`, which reads every input
as an exact integer mantissa and binary exponent.  No independently
*implemented* RT solver (libRadtran, cdisort, VLIDORT) is installed; none was
installed.  The quad build is therefore the reference for arithmetic
precision, and the gfortran double build the cross-compiler check.

## 2. Execution through IDL

`run_idl_double_validation.sh` (run by you on 2026-10-10, atmlxint5, IDL 8.9;
results `idl_double/results/run_20261010T214143/`, checksums in
`outputs.sha256`, provenance in `provenance.txt`):

1. **Interface test** (`idl/dp_interface_test.pro`): all 11 checks passed —
   8-byte kernel; historical DISORT reproduces the V21 capture **bitwise**
   with both DLMs loaded; DOUBLE outputs of the right dimensions; inputs
   unchanged; float32 and float64 inputs (and INT/LONG scalars) give
   bitwise-identical results; repeat call bitwise identical; `LAMBER = 0`,
   wrong `PMOM` layout, `SSALB > 1` and short `UMU` raise IDL errors.
2. **Validation runs**: the unmodified V21 generator with
   `idl/validation/call_disort.pro`, whose argument checks and layer
   compression are V21 lines 11-92 verbatim.  The V21 single-precision result
   is used for the LUT; every call is also solved by `DISORT_DP` at NSTR
   48/60/80/100 (+ a repeat at 60) from pre-call copies of the inputs.
3. **Generator runs**: the same generator with `idl/generator/call_disort.pro`
   (`DISORT_DP` at NSTR 60) writing small LUTs labelled version 90.

## 3. Test cases

All Terra MODIS liquid-water 240 K (`liquid-water_240.mm`), atmosphere 2,
`srf_quad = 1`, NMOM 999, full 10 x 11 view grid, with the production
optical-property and atmosphere processing unchanged.  Each case includes
the adjacent optical-depth and effective-radius nodes for Jacobians.

| Case | Channel | Centre node | Grid computed | Why |
|---|---:|---|---|---|
| A | 6 | τ 256, r_e 29, SZA 80 | τ 64/128/256, r_e 27/29/31, SZA 70/80 | original severe failure; SZA 80 supports interpolation to SZA < 75 |
| B | 20 | τ 1, r_e 11, SZA 30 | τ 0.5/1/2, r_e 9/11/13 | largest mixed-channel BT error in the operational domain (VZA 36.1, RAA 126) |
| C | 6 | τ 0.125, r_e 21, SZA 30 | τ 0.0625/0.125/0.25, r_e 19/21/23 | V21 R_0v = -0.00578 at VZA 43.3, RAA 0 where V22 is +0.00567 |
| D | 1 | τ 8, r_e 11, SZA 40 | τ 4/8/16, r_e 9/11/13 | ordinary node; V21 and V22 agree to 1.8e-5 |

90 DISORT calls in total (diffuse, direct-beam and, for channel 20,
emission).  For every case the compact single-precision LUT equals the
production V21 LUT **bitwise** at the corresponding nodes for all nine RT
coefficients.

## 4. Results

### 4.1 Precision: historical single precision versus DISORT_DP

Maximum |SP - DP(60)| over each case's LUT coefficients, by call type:

| Case | Direct beam | Diffuse | Emission |
|---|---:|---:|---:|
| A | 0.687 | 1.4e-5 | — |
| B | 0.049 | 9.1e-5 | 9.0e-7 |
| C | 0.011 | 1.1e-4 | — |
| D | 3.9e-5 | 6.3e-5 | — |

Centre node, `R_0v` against the NSTR 100 double-precision result, in units
of the production solar uncertainty sigma = R sqrt(noise^2 + 0.01^2 + 0.02^2)
(mixed channel 20 has no SNR in the LUT, so only the calibration terms):

| Case | R_0v at the reported geometry: SP / DP60 / DP100 | max sigma, SP | max sigma, DP60 | max sigma, DP48 |
|---|---|---:|---:|---:|
| A (VZA 0, RAA 180) | 0.06760 / 0.03939 / 0.03940 | **36.5** | 0.046 | 0.050 |
| B (VZA 36.1, RAA 126) | **-0.01419** / 0.03253 / 0.03253 | **70.2** | 0.015 | 0.041 |
| C (VZA 43.3, RAA 0) | **-0.00578** / 0.00531 / 0.00530 | **17.5** | 0.037 | 0.109 |
| D (VZA 28.9, RAA 90) | 0.28694 / 0.28693 / 0.28692 | 0.10 | 0.077 | 0.057 |

At the ordinary node D the single-precision result is already accurate, and
its difference from double precision is the same size as the stream-count
spread.

### 4.2 Stream-count convergence

Maximum |DP(N) - DP(100)| and successive differences over all coefficients:

| Case | Call type | 48 | 60 | 80 | 48-60 | 60-80 | 80-100 |
|---|---|---:|---:|---:|---:|---:|---:|
| A | direct | 1.9e-4 | 5.9e-4 | 3.2e-4 | 4.2e-4 | 3.8e-4 | 3.2e-4 |
| A | diffuse | 1.1e-6 | 5.8e-7 | 8.6e-7 | 1.3e-6 | 1.3e-6 | 8.6e-7 |
| B | direct | 1.2e-4 | 7.3e-5 | 1.2e-5 | 1.7e-4 | 6.8e-5 | 1.2e-5 |
| B | diffuse / emission | 9.4e-7 / 7.5e-7 | 2.4e-7 / 1.6e-7 | 1.5e-8 / 9.3e-9 | | | |
| C | direct | 1.0e-4 | 5.3e-5 | 8.5e-5 | 1.0e-4 | 6.7e-5 | 8.5e-5 |
| D | direct | 5.8e-4 | 7.3e-4 | 7.4e-4 | 7.9e-4 | 8.5e-4 | 7.4e-4 |

Diffuse and emission coefficients converge monotonically.  The direct-beam
`R_0v` spread is bounded (<= 8.5e-4) but not monotonic in the large-particle
cases A and D: it is the residual of the delta-M truncation and single-
scattering intensity correction at sharp angular features of the 1000-moment
phase function, so NSTR 100 is not automatically the most accurate.  It is the
same discretisation error the V21 configuration has always carried, and at
NSTR 60 it is below 0.08 sigma at every node tested.  **NSTR 60, the V21
value, is retained**; changing it would change the numerical formulation and
cost (x2.1 at 80, x4.2 at 100) for no demonstrated operational benefit.

### 4.3 Arithmetic precision and independent builds

| Comparison | Result |
|---|---|
| DISORT_DP in IDL vs the same ifort kernel called directly | bitwise identical, all 90 calls |
| Repeat call in the same IDL process | bitwise identical, all 90 calls |
| DISORT_DP in IDL vs gfortran double | <= 3.1e-8 in any coefficient (0 for case A) |
| DISORT_DP vs gfortran REAL*16, NSTR 60, all 90 calls | <= 3.1e-8 (A 4.8e-9, B 1.8e-8, C 3.1e-8, D 3.3e-9) |
| DISORT_DP vs REAL*16 at the centre nodes, NSTR 48-100 | <= 1.7e-8 at every stream count |

The arithmetic-precision error of double precision (<= 3e-8) is four orders
of magnitude below the stream-count spread and six below the measurement
uncertainty; the single-precision error at the failure nodes is up to 0.69.

### 4.4 Physical consistency

| Check (direct-beam calls) | Single precision, worst | Double precision, worst (all NSTR) |
|---|---|---|
| Downward diffuse flux at TOA (exact 0) | 0.672 (A) | 1.6e-7 |
| Minimum intensity | -26.8 (C), -2.26 (B) | 0 |
| Reflected fraction of incident flux | 1.323 (A) | 0.647 |
| Absorbed fraction | -0.323 (A) | > 0 (min 3.4e-5, D) |

Diffuse calls: TOA downward flux equals pi FISOT to <= 1.3e-9 relative in
double precision (up to 5.5e-5 in single); no negative intensities; absorbed
fraction positive in every call (>= 3.4e-5, at the nearly conservative
channel-1 node D).  Emission calls:
zero downward intensity at TOA, no negative values.

### 4.5 Jacobians

Finite differences on the computed neighbours, as ORAC's bicubic scheme uses
(d/dlog10 τ and d/dr_e of `R_0v`, maximum over the view grid at the centre):

| Case | |SP - DP100| dR/dlog10τ | |DP60 - DP100| dR/dlog10τ | |SP - DP100| dR/dr_e | |DP60 - DP100| dR/dr_e | DP100 scale (τ / r_e) |
|---|---:|---:|---:|---:|---|
| A | **0.68** | 7.0e-8 | 9.7e-5 | 5.2e-6 | 2.3e-6 / 1.4e-3 |
| B | **0.014** | 5.0e-5 | 1.5e-5 | 1.2e-5 | 0.18 / 9.3e-3 |
| C | 9.5e-5 | 8.2e-5 | 7.4e-6 | 2.3e-6 | 0.039 / 3.7e-4 |
| D | 2.1e-5 | 1.2e-5 | 1.1e-4 | 9.8e-5 | 0.48 / 7.8e-3 |

At the saturated case A the true optical-depth derivative is ~2e-6 while the
single-precision value is in error by 0.68, five orders of magnitude; the
double-precision Jacobians are stable with stream count to well below their
own magnitude except where the stream-count spread of section 4.2 applies
(D, r_e derivative, 1.3 %).

### 4.6 Differences from historical V21

The double-precision generator changes the solver only: the compact version-90
LUTs differ from the compact V21 LUTs by the SP-DP differences of 4.1 (up to
0.69 in `R_0v` at A, ~1e-4 or less in diffuse/thermal coefficients and at the
ordinary node D), and the version-90 `R_0v` equals the captured DP(60) result
to float32 rounding (<= 6e-8).

## 5. Output precision

ORAC holds LUT coefficients in default `REAL` (`src/sad_lut.F90`, `real,
pointer :: Rbd(...)` etc.; `common/common_constants.F90:34` `sreal = 4`; no
real promotion in the build configuration) at commit `eb9a1233`.  Storing
float64 coefficients would therefore be rounded on read and gain nothing; the
float32 storage error (~6e-8 relative) is far below both the stream-count
spread and the measurement uncertainty.  The generator keeps double precision
inside the solver and the existing NC_FLOAT variables, names, dimensions and
`/100` normalisation; `create_orac_cloud_lut` rounds each DOUBLE result once
when adding it to its FLTARR arrays.

## 6. Cost

Per DISORT call, measured inside IDL over the 90 calls (mean, wall clock):

| Solver | Time per call | Relative to V21 |
|---|---:|---:|
| historical single precision (NSTR 60) | 0.061 s | 1.00 |
| DISORT_DP NSTR 48 | 0.053 s | 0.87 |
| DISORT_DP NSTR 60 | 0.081 s | 1.33 |
| DISORT_DP NSTR 80 | 0.170 s | 2.80 |
| DISORT_DP NSTR 100 | 0.337 s | 5.56 |

A full Terra MODIS product needs 312 DISORT calls per (τ, r_e) node — about
106,000 calls for each liquid LUT and 127,000 for the ice LUT.  The V21
production timestamps imply about 3 h per product.  With the solver taking
~1.8 h of that, double precision at NSTR 60 adds ~0.6 h (liquid) to ~0.7 h
(ice): **about 3.6 h for pold or p240 and 4 h for pagg**, ~11.5 h for all
three against ~9 h historically.  These are estimates from a 90-call sample;
the per-call cost varies little with optical depth because layer compression
leaves three layers.

## 7. Limitations

- Four nodes of one LUT (`p240`); `pold` and `pagg` and other instruments are
  covered only by the shared code path.  The full LUTs will provide the
  population test.
- No independently implemented solver was available; the references are
  independent builds (compiler, precision) of the same DISORT algorithm, so
  algorithmic errors common to DISORT 2.1 are not excluded.
- Stream-count convergence of direct-beam `R_0v` at sharp angular features is
  bounded but not monotonic (section 4.2).
- Measurement significance for mixed channel 20 uses the calibration terms of
  the solar uncertainty only; brightness-temperature propagation of the
  double-precision LUTs is left to the full-LUT comparison.
- The remaining SGECO warnings in the run log come only from the historical
  single-precision calls in the validation runs; the double-precision
  generator runs printed none.

## 8. Readiness and next steps

The implementation is ready for full LUT testing.  Prepared but **not run**:

```shell
bash validation/v21_v22_terra_modis/idl_double/run_dp_generator_full.sh          # pold p240 pagg
bash validation/v21_v22_terra_modis/idl_double/run_dp_generator_full.sh p240     # one product
```

It uses the V21 production `terra_modis_run` arguments with version label 90
(it refuses 21) and writes only to `idl_double/results/full_<timestamp>/`.
After it runs: compare the version-90 LUTs with V21 and V22 over the whole
grid and the operational domain, repeat the operational radiance-impact
analysis, and only then decide on promoting double precision to the
reference.  V21 sources and LUTs stay archived and unchanged.

## 9. Reproduction

```shell
bash validation/v21_v22_terra_modis/idl_double/build_dlm.sh
bash validation/v21_v22_terra_modis/idl_double/build_reference.sh
bash validation/v21_v22_terra_modis/idl_double/run_idl_double_validation.sh   # needs IDL
python validation/v21_v22_terra_modis/idl_double/analyse_idl_double.py \
   validation/v21_v22_terra_modis/idl_double/results/run_20261010T214143 --quad-all
```

Analysis output: `results/run_20261010T214143/analysis/summary.json`.
