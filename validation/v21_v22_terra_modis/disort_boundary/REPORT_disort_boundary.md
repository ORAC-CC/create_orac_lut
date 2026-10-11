# V21/V22 DISORT boundary: Terra MODIS p240, channel 6

> **Historical (single-precision V21/V22).**  This report analyses the V21 IDL
> and V22 Python products made with the historical single-precision DISORT.
> Since 2026-10-10 both generators use double-precision DISORT by default;
> see `DOUBLE_PRECISION_BASELINE.md` (repository root) for the current
> baseline.  The findings below remain the record of the single-precision
> products, which stay archived.

Selected call: liquid-water 240 K (`liquid-water_240.mm`), Terra MODIS
channel 6 (1.629 um, `srf_quad=1`), optical depth 256, effective radius 29 um,
solar zenith 80 deg, direct beam (FBEAM = 100), with the full 10 satellite
zenith x 11 relative azimuth grid of that call.  This is the rank-1 call in
[`terra_modis_p240_v21_v22_direct_beam_calls.csv`](../results/terra_modis_p240_v21_v22_direct_beam_calls.csv)
(maximum `R_0v` difference 0.541).

## Conclusion

The LUT values differ because this DISORT solve is numerically ill-conditioned
in single precision.  Every single-precision build prints
`UPBEAM--SGECO says matrix near singular` for this call, and its result is then
decided by rounding: 1-ulp changes in the inputs, or a different compiler and
optimisation level for the same Fortran source, change the top-of-atmosphere
radiances by O(1).  V21 and V22 differ in **both** respects:

1. **Inputs.** The corrected V22 Python inputs differ from V21 by 1-2 ulp in
   layer optical depths, in 837 of the 1000 cloud phase-function moments and in
   4 of the 20 view cosines.  With the exact V21 solver binary, these input
   differences alone change top-of-atmosphere upwelling `UU` by up to 1.42
   (the cloud moments dominate).
2. **Solver build.** V21 calls the preserved historical IDL DLM
   (`DISORTtoIDL.so`, Intel Fortran `-O3`, built on Ubuntu 18.04 with an
   unrecorded ifort version).  V22 calls `liboraclut_disort_gfortran.so`
   (gfortran 12.2 `-O2`) built from the same source.  With identical V21
   inputs the gfortran build differs from V21 by up to 2.75.

Combined, V22 differs from V21 by up to 6.54 in `UU` (`R_0v` 0.54).
**Neither result is correct.**  The double-precision build of the same
DISORT source is well conditioned here: it gives the same answer (to 1.7e-5)
from either input set, and V21 itself is up to 82 % in error against it (V22
up to 400 %, including negative radiances).

So bitwise V21 reproduction requires *both* the identical solver binary *and*
bitwise-identical inputs, and would reproduce a numerically wrong V21 value at
this call.  Rebuilding the DLM source with today's ifort and the same flags is
not sufficient (see below).

## 1. Captures

| Side | Capture | How |
|---|---|---|
| V21 | `results/v21_capture_20261010T175154/v21_disort_boundary_{input,output}.sav` | [`run_v21_capture.sh`](run_v21_capture.sh) with the unmodified V21 IDL sources and a `call_disort` shadow ([`capture_call_disort.pro`](capture_call_disort.pro), numerically identical statements) saving the arguments before and the results after the DLM call; IDL 8.9, intel-compilers/2022, atmlxint5.  Provenance and SHA-256 in `provenance.txt` there. |
| V22 | `results/corrected_v22_disort_boundary.npz` | [`capture_corrected.py`](capture_corrected.py) on the corrected V22 production path (IDL-compatible `TOTAL`). |

Both captures select the same call: μ0 = 0.17364822 is bit-identical, and the
compact LUT gives exactly the production grid values (τ = 256, r_e = 29.0,
SZA = 90·8/9 = 80.0 in float32).

## 2. Inputs: argument-by-argument

All 44 arguments compared bitwise (V21 arrays truncated to NLYR = 3 after
`call_disort` layer compression):

| Argument | Result |
|---|---|
| NLYR, NMOM (999), NSTR (60), NTAU, NUMU (20), NPHI (11), USRTAU, USRANG, IBCND, LAMBER, PLANK, ONLYFL, PRNT, MAXULV, MAXUMU, MAXPHI, MAXMOM | identical |
| FBEAM, UMU0, PHI0, FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, WVNMLO, WVNMHI, ACCUR, TEMPER, PHI | bitwise identical |
| MAXCLY | 45 (IDL) vs 3; **no effect**: the DLM C wrapper overrides it with NLYR, and direct calls with 3 or 45 give identical results in every build |
| SSALB | Clear layers 1.0 (V21) vs 0.9999988 (V22): **capture artefact**.  DISORT overwrites SSALB = 1 with 1 − DITHER in place (`DISORTfunctions.f:681`), and the Python capture recorded the array after the call.  The true inputs are identical; the cloud layer (0.98322916) is bitwise identical. |
| DTAUC | all 3 layers differ by 1-2 ulp (cloud 261.9816 vs 261.98163) |
| UTAU | total optical depth differs by 1 ulp |
| PMOM | clear (Rayleigh) layers identical; cloud layer differs in 837/1000 moments, max 6e-8 |
| UMU | 4 of 20 view cosines differ by 1 ulp |

Upstream origin (first differing intermediates, from the two `scatfile`s):
at the channel wavelength `bext`, `w` and `g` are bitwise identical, but the
phase function `phs` differs by up to 39 ulp and its Legendre moments `amom` in
819/1000 values.  At 550 nm (the optical-depth reference) `bext550` differs by
2 ulp, so `bextrat` differs by 1 ulp, hence the cloud DTAUC.  The clear-layer
DTAUC (Rayleigh) and UMU (cosine of the satellite zeniths) differences arise
separately; their individual arithmetic origins were not traced further.

## 3. Outputs and solver experiments

[`diagnose_solver.py`](diagnose_solver.py) passes both input sets directly to
the Fortran `disort_` entry point of each library, with every argument under
explicit control, one library per process.  Results:
`results/diagnose_solver_{dlm_v21,gfortran_v22,ifort2022}.json`.  `UU_TOA↑`
is `UU` at UTAU = 0 for the 10 upward directions (the quantity behind `RBD`);
the double-precision reference ranges 1.25-16.16.

| Solver | Inputs | FLUP top | RFLDN top | max abs ΔUU_TOA↑ vs double | outputs bitwise = V21 capture |
|---|---|---|---|---|---|
| V21 IDL capture | V21 | 7.4107 | +0.0498 | 1.68 | — |
| V21 DLM binary (`DISORTtoIDL.so`) | V21 | 7.4107 | +0.0498 | 1.68 | 472/490 (rest ≤ 1.4e-5) |
| V21 DLM binary | V22 | 8.8203 | +0.0828 | 0.64 | 152/490 |
| ifort 2022 `-O3` rebuild (DLM flags) | V21 | 7.4107 | +0.0498 | 2.31 | 182/490 (UU up to 0.87 off) |
| ifort 2022 `-O3` rebuild | V22 | 8.8203 | +0.0828 | 1.37 | 152/490 |
| V22 gfortran `-O2` (production) | V21 | 7.0657 | +1.3161 | 2.75 | 152/490 |
| V22 gfortran `-O2` (production) | V22 | 2.5134 | −0.9980 | 8.23 | = V22 capture 490/490 |
| **double precision** (gfortran `-fdefault-real-8`, MXCMU 100) | V21 or V22 | **11.0488** | **≈ 0** | — | V21-vs-V22 inputs agree to 1.7e-5 |

TOA downward diffuse flux (RFLDN) should be ~0 with FISOT = 0; only double
precision gives that.  The gfortran call reproduces the V22 capture bitwise,
confirming the capture.

Attribution with the V21 DLM binary (swap one V21 input group for V22):
DTAUC+UTAU changes `UU_TOA↑` by < 1e-4, UMU by 0.25 (only in the 4 changed
directions), PMOM by 1.42, all together by 1.50.  With the gfortran build the
same swaps give 2.7-4.1 individually and 6.5 together; the effects are not
additive, as expected for an ill-conditioned solve.

## 4. Matching the Python build to the IDL build

- The DLM's Fortran objects carry no compiler identification; the C parts were
  built with GCC 7.5 on Ubuntu 18.04, and the files arrived with the
  2026-09-13 history import.  The original ifort version is unknown.
- A fresh build of the same sources with intel-compilers/2022 and the DLM
  Makefile flags (`-O3 -fpic -D__GFORTRAN__`) does **not** reproduce the DLM
  (182/490 bitwise, `UU` up to 0.87 different, V21 inputs).
- Calling the preserved DLM binary's `disort_` from Python reproduces the IDL
  capture for 472/490 values; the remaining 18 are `UU` at UTAU = 0 in the
  three most nadir upward directions, ≤ 1.4e-5 (≈ 1e-5 relative).  They are
  deterministic and unaffected by preceding the call with the IDL diffuse
  call, MAXCLY padding, array alignment, stack contents, FTZ/DAZ, x87
  precision control, or glibc versus Intel `libimf` symbol binding.  Their
  cause is unresolved; it does not affect the conclusions.

## 5. Limitations

- One call only.  The other ranked V22 outliers were not rerun here, though
  the mechanism matches the earlier census in `validation/v24`.
- The double-precision library is the existing validation build
  `validation/tmp/disort_options/build/double100/libdisort_double100.so`
  (production MXCMU = 100, production DELTAM/CORINT/ACCUR).  It is a
  reference, not a production calculation.
- The ifort 2022 test library is kept at
  `validation/tmp/v21_v22_disort_boundary/ifort2022_O3/` (SHA-256
  `d76b7029…cfdb2f`), built with `legacy_disort._build_library` and
  `ORACLUT_FORTRAN=ifort`.

Follow-up: the double-precision DISORT DLM validated through IDL in
[REPORT_idl_double_precision.md](../REPORT_idl_double_precision.md) reproduces
the double-precision reference for this call (FLUP 11.0488) and agrees with a
REAL*16 build to 5e-9.

No production code, IDL source, DISORT source or LUT was modified.

## 6. Reproduction

```shell
# V21 capture (needs IDL; run by the user)
bash validation/v21_v22_terra_modis/disort_boundary/run_v21_capture.sh

# V22 capture
python validation/v21_v22_terra_modis/disort_boundary/capture_corrected.py

# Solver experiments (intel-compilers/2022 loaded for the DLM's runtime)
ORACLUT_DISORT_DOUBLE=validation/tmp/disort_options/build/double100/libdisort_double100.so \
python validation/v21_v22_terra_modis/disort_boundary/diagnose_solver.py \
  validation/v21_v22_terra_modis/disort_boundary/results/diagnose_solver_dlm_v21.json \
  create_orac_lut/disort2/src/DISORTtoIDL.so
```
