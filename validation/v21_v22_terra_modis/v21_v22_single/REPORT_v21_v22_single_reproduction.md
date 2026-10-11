# Single-precision V21 (IDL) / V22 (Python) DISORT reproduction, stages 1-5

Date 2026-10-11.  Scope: reproduce the historical single-precision V21 IDL
results with V22 Python, including the behaviour of the original
single-precision DISORT.  The double-precision work is suspended, not
removed: its sources, libraries, manifests and products are untouched.  Every
experiment here selects single precision explicitly
(`ORAC_LUT_DISORT_PRECISION=single` for the Python generator; the original
`DISORTtoIDL.so` for IDL); the production default (double) is not relied on.

## Summary

1. **The same compiled kernel gives bitwise-identical results from IDL and
   Python, provided the memory layout matches.**  Called from Python, the
   original Intel DLM reproduces all 90 IDL-captured single-precision DISORT
   calls (four Terra MODIS p240 cases: diffuse, direct-beam and emission)
   **bitwise**, when the argument arrays are placed as IDL places them,
   8 bytes past a 16-byte boundary.  With 64-byte-aligned (NumPy-default)
   arrays only 45/90 match: every direct-beam call differs in `UU` (up to
   1.5 at the failure node).  The controlling array is `UMU`: the `-O3`
   Intel code vectorises loops over user polar angles and peels a different
   number of elements depending on its alignment, which changes the
   summation order.  This explains the 18 "unexplained" `UU` values of
   `disort_boundary/REPORT_disort_boundary.md`.
2. **The two historical kernels are not interchangeable.**  The V22 Python
   kernel (gfortran -O2 build of the same source) differs from the Intel
   DLM for identical inputs in every call: up to 0.06 in diffuse, 1e-9 in
   emission and 29 in direct-beam outputs.  Bitwise V21 reproduction
   therefore requires the Intel DLM kernel, not the gfortran build.
3. **Input differences.**  Integer, logical and configuration arguments
   (NLYR, NMOM = 999, NSTR = 60, NTAU, NUMU, NPHI, USRTAU, USRANG, IBCND,
   LAMBER, PLANK, ONLYFL, ACCUR, ALBEDO, BTEMP, TTEMP, TEMIS, PHI0) are equal
   in every call; `PHI`, `TEMPER`, `FBEAM`, `FISOT`, `WVNMLO/HI` and the
   direct-beam `UMU0` are bitwise equal.  Floating-point inputs differ by
   1-3 ulp: `DTAUC` (90/90 calls, <= 3 ulp), `UTAU` (78, <= 3 ulp),
   `PMOM` (90, <= 1.2e-7 absolute), `UMU` (90, 1 ulp), `SSALB` (14 cloud-layer
   values, 1 ulp) and the nominal `UMU0` = cos 50 deg of diffuse/emission calls
   (1 ulp; unused there because FBEAM = 0).
4. **At the failing p240 call (channel 6, tau 256, r_e 29, SZA 80) only
   `PMOM` and `UMU` matter.**  With the Intel kernel, V21 and V22 inputs
   differ by 1.50 in top-of-atmosphere upwelling `UU`.  Substituting the V21
   phase moments into the V22 inputs leaves 0.25; adding the V21 view
   cosines gives **bitwise V21**.  `DTAUC`/`UTAU`, `SSALB` and the scalars
   change nothing (<= 1e-6) in either direction.
5. **`UMU` can be made bitwise.**  IDL evaluates `cos(tmpSat*!DTOR)` as a
   float32 product with float32 `!DTOR` = 0.017453292 followed by a
   single-precision cosine; the Python `np.cos(np.deg2rad(tmpsat))` differs in
   3 of the 10 view cosines.  Reproducing the IDL arithmetic (float32 product,
   `cosf` or a double cosine rounded to float32) matches all nine non-nadir
   values (the nadir value is fixed by the |UMU| = 1 rule).  The `PMOM`
   differences originate upstream, in the phase function computed by the Mie
   and Legendre-expansion arithmetic (`disort_boundary` report, section 2).

So V22 can reproduce V21 bitwise at the worst node if (a) it computes `UMU`
and `PMOM` with IDL's arithmetic and (b) it calls the original Intel DLM
kernel with IDL's memory layout.  No change was made to any generator,
kernel, LUT or numerical setting in this work.

## 1. Historical reference configurations

| | V21 (IDL) | V22 (Python) |
|---|---|---|
| Source set | `releases/idl_lut_v21_single.release` (161 entries; `call_disort_v01.pro`) | historical: `releases/lut_v22.release` (`*_v01.py`); inputs here from the active TOTAL-corrected revisions (`create_orac_luts_v03.py`, `idl_mirror/*_v02.py`, `legacy_disort_v02.py` single path = `legacy_disort_v01.py` body); hashes in `results/python_single_sources.json` |
| DISORT source | `create_orac_lut/disort2/src/DISORTfunctions.f` (sha256 3f45489e...) | same file |
| Compiled kernel | `create_orac_lut/disort2/src/DISORTtoIDL.so` (sha256 59307e42..., preserved as `DISORTtoIDL_v01.so`); ifort, version unrecorded, Makefile target `ifort`: `-O3 -fpic -D__GFORTRAN__`; C parts GCC 7.5 (Ubuntu 18.04) | `build/oraclut/liboraclut_disort_gfortran.so` (sha256 2e4740e3...); gfortran 12.2.0 `-O2 -fPIC -std=legacy -ffixed-line-length-0`, C adapter `src/oraclut/radiative_transfer/legacy_disort.c` |
| Run-time libraries | Intel 2022 runtime (`libifcore`, `libimf`, `libsvml`, `libirc`, `libintlc`; module `intel-compilers/2022`); math functions (`cosf`, `expf`, `sincosf`, ...) from glibc `libm` in IDL and in Python alike | libgfortran, glibc libm |
| Entry point | `INTEGER FUNCTION DISORT` (`disort_`), 44 arguments by reference; called by `DISORTwrapper.c` | same Fortran entry through `oraclut_disort` (fixed flags in C) |
| Argument precision | REAL = float32, INTEGER = int32, LOGICAL = int32 0/1; MAXCLY/ULV/UMU/PHI/MOM forced to the actual sizes by the wrapper | identical |
| NSTR / NMOM | 60 / 999 (1000 moments) | 60 / 999 |

Both use single precision, 60 streams, the fixed 1000-moment Legendre
expansion, the original coarse grids and the same physical inputs.  One
wrapper quirk with no numerical effect: IDL passes `PRNT` as a 16-bit INT
array that `DISORTwrapper.c` reads as 32-bit, so the last print flags are
read from memory beyond the array; they control printing only.

The original Intel kernel can be called from Python: `ctypes.CDLL(...,
RTLD_LAZY)` resolves `disort_` (the `IDL_*` imports are needed only by the C
wrapper), arguments passed by reference with the types above.  The harness is
`sp_kernels.py`; it places every array at a chosen offset modulo 64.

## 2. Captures used

* V21: the 90 single-precision calls captured inside IDL (pre-call copies of
  all arguments, outputs of the original DLM) by the validation shadow in
  `idl_double/results/run_20261010T214143` — cases A (p240 ch6, tau 64/128/256,
  r_e 27/29/31, SZA 70/80), B (ch20), C (ch6, tau 0.0625-0.25), D (ch1).  They
  are the historical V21 configuration: the compact LUTs of that run equal the
  production V21 products bitwise in every RT coefficient.
* V22: `capture_python_single.py` reruns the active Python generator on the
  same four cases with `ORAC_LUT_DISORT_PRECISION=single`, copying every
  argument before the call (the older `disort_boundary` V22 capture recorded
  `SSALB` after DISORT had dithered it in place).  The same 90 calls, in the
  same order, kinds and sizes.

## 3. Four-way experiment (90 calls)

| Comparison | Bitwise calls | Max abs difference |
|---|---|---|
| A (V21 in, Intel) vs IDL capture | **90/90** | 0 |
| D (V22 in, gfortran) vs Python capture | **90/90** | 0 |
| A vs C (V21 inputs, Intel vs gfortran) | 0/90 | diffuse 0.061, direct 29.2, emission 1.4e-9 |
| B vs D (V22 inputs, Intel vs gfortran) | 0/90 | diffuse 0.063, direct 51.0, emission 1.6e-9 |
| A vs B (Intel, V21 vs V22 inputs) | 0/90 | diffuse 0.039, direct 21.8, emission 1.7e-9 |

Layout: with all arrays 64-byte aligned instead, A matches IDL in 45/90 calls
(all diffuse and emission calls; no direct-beam call).  Shifting only `UMU`
to offset 8 restores the direct-beam calls; `PHI`, `UTAU` and the other
arrays are layout-insensitive.  Ruled out as causes (experiments on the
residual): FTZ/DAZ and x87 precision control, preceding calls and saved state,
stack contents, glibc versus Intel math symbols, library dispatch (IDL loads
the same Intel runtime and exports no math or Intel symbols).

## 4. First differing inputs

From `results/four_way.json` (per argument: calls differing, max ulp, max
absolute difference, first differing element):

| Argument | Calls differing | Max ulp | Max abs | First difference (case, call, element: IDL / Python) |
|---|---|---|---|---|
| DTAUC | 90 | 3 | 3.05e-5 | A, 1 (diffuse), [0]: 7.4304215e-4 / 7.4304221e-4 |
| UTAU | 78 | 3 | 3.05e-5 | A, 1, [1]: 65.571526 / 65.571533 |
| PMOM | 90 | — (sign changes near zero) | 1.19e-7 | A, 1, [1,1]: 0.87199736 / 0.87199730 |
| UMU | 90 | 1 | 5.96e-8 | A, 1, [8]: -0.53320444 / -0.53320438 |
| SSALB | 14 | 1 | 5.96e-8 | B, 4, [1]: 0.92424041 / 0.92424035 |
| UMU0 (diffuse/emission nominal) | 45 | 1 | 5.96e-8 | A, 1: 0.64278764 / 0.64278758 |
| PHI, TEMPER, FBEAM, FISOT, WVNMLO, WVNMHI, direct UMU0 | 0 | 0 | 0 | — |

Integer and logical arguments and all fixed solver settings are identical.

## 5. Substitution at the failing call (case A capture 24)

Intel kernel, IDL layout; `UU_TOA` = upwelling intensity at the top over the
10 x 11 view grid.  V21 vs V22 inputs: max |dUU_TOA| = 1.501.

| Group substituted | V22 base + V21 group: max abs vs V21 | V21 base + V22 group: max abs vs V21 |
|---|---|---|
| DTAUC + UTAU | 1.501 (no effect) | 9.5e-7 |
| SSALB | 1.501 (no effect) | 0 |
| PMOM | **0.250** | **1.422** |
| geometry (UMU, PHI, UMU0) | **1.422** | **0.250** |
| scalars (FBEAM, FISOT, WVNMLO/HI) | 1.501 (no effect) | 0 |

Cumulative from V22: DTAUC, SSALB 1.501; + PMOM 0.250; + geometry **0
(bitwise V21)**.  The two causes interact (the ill-conditioned call amplifies
each), and both must be reproduced exactly.

## Files

`sp_kernels.py` (layout-controlled calls to both historical kernels),
`capture_python_single.py`, `four_way.py`, `results/four_way.json`,
`results/python_single_case_*.npz`, `results/python_single_sources.json`.
