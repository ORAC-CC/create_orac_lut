# DISORT options in ORAC LUT generation: effect on the ORAC TOA measurement, its Jacobians and cost

Written for the owner of the ORAC LUT generator and for reviewers of future
production changes. This was a read-only investigation: no production code,
configuration or V24 calculation was changed.

## 0. Summary and recommendation

**Decision basis.** The decision rests on the complete ORAC simulated
measurement, not on any single LUT operator:
- the LUT operators of each DISORT configuration, propagated together through
  the production ORAC fast forward model (§9–§10);
- its retrieval Jacobians;
- comparison with the measurement uncertainty ORAC itself uses (§8);
- the computational cost (§14).

Operator-level results are kept as diagnostics (§13).

**Numerical reference.** Double-precision DISORT 2.0 with NSTR = 236, all other
settings as in production. This is a convergence reference within the same
radiative-transfer formulation, not truth (§11).

**Principal finding: a production defect.** The production single-precision
DISORT 2.0 sporadically returns a completely wrong direct-beam solution for an
individual (channel, τ, r_e, solar zenith) call (§12).
- **Prevalence:** whole-call failures were found in 17 of 40 V23 and 14 of 35
  V24 products checked: 44 and 32 cells, about 5×10⁻⁶ of direct-beam cells,
  up to 12 per product. They include ORAC retrieval channels (MODIS 0.86 µm;
  SLSTR 0.66, 1.6, 2.25 and 3.7 µm).
- **Partial failures:** some calls are wrong in only part of their views; five
  such calls push R_0v outside the physical range. The census is therefore a
  lower bound (§12.2).
- **Effect on the measurement:** every view of the affected call is wrong.
  R_0v ranges from −4.8 to +5.7 where the correct range is 0.38–0.68. In the
  ORAC measurement this is 13–244σ (median over views) and up to 856σ at the
  LUT node. ORAC's bicubic and trilinear interpolation spreads the error into
  the neighbouring cells.
- **Cause:** confirmed by reproduction. All four reproduced cases are restored
  by double precision at the same NSTR = 60, and equally by single precision at
  NSTR = 56, 68 or 100 (all agree). The production value is the only outlier.
  Only R_0v is affected; R_0d, T_0d and T_00 from the same call agree with
  double precision to ≤ 10⁻⁴.

This is the only DISORT-related change found that is worth considering for a
future version.

**Answers to the questions:**

1. **Is NSTR = 60 sufficiently converged for the complete ORAC measurement?**
   Yes, against ORAC's production uncertainty, everywhere except in a narrow
   region around exact backscatter.
   - Away from backscatter (scattering angle Θ < 170°, over 5,100 points per
     channel and surface), NSTR = 60 differs from the reference by at most
     0.38σ in every solar channel, case and surface. No point exceeds 0.5σ.
     Here σ is the production ORAC measurement uncertainty.
   - Thermal channels differ by at most 0.002σ.
   - The 3.7 µm mixed channel is at most 0.26σ (liquid) and 0.61σ (Baum ice).
     For large ice spheres it reaches 1.1σ, but that is a non-converging stream
     oscillation, not an NSTR = 60 deficiency (§11.3).
   - **Against instrument noise alone** (without ORAC's homogeneity and
     co-registration terms) the same differences reach, away from
     backscatter:
     - 1.5–2.3σ in solar channels (worst at MODIS 0.86 µm, where the SNR is
       highest);
     - 4–11σ in the 3.7 µm channel.

     NSTR = 60 is converged relative to the uncertainty ORAC actually assigns,
     not relative to instrument noise (§10.1).
2. **Where does NSTR = 60 differ by a significant fraction of σ?**
   - Only at and near exact backscatter (solar zenith = view zenith, relative
     azimuth 0; the glory of droplets and the backscatter peak of ice).
   - There, in visible to SWIR and 3.7 µm channels, 10–70% of cloud states
     exceed 0.5σ. The maximum is 5σ (liquid), 6.5σ (ice spheres) and 14σ
     (Baum ice, at small signal).
   - Angular width: for droplets and ice spheres the error exceeds 0.5σ only at Θ = 180° itself, and is ≤ 0.31σ at 177.5° (≤ 0.7σ for Baum, ≤ 1σ for the 3.7 µm channel). Through ORAC's angular interpolation, the exact-backscatter LUT nodes affect measurements within one Grid B cell (5° in RAA, 3.75° in SZA and VZA) (§10.4).
   - Surface conditions do not change this: the effect sits in R_0v and passes
     unchanged into the measurement.
3. **Does increasing NSTR improve the Jacobians where the radiance is already
   close?** No.
   - Away from backscatter the Jacobian errors of NSTR = 60, expressed as the
     measurement error over half a LUT cell, are ≤ 0.28σ (dlog₁₀τ) and
     ≤ 0.18σ (dr_e).
   - At backscatter they follow the measurement error, up to 3.3–3.6σ (liquid,
     ice spheres), and higher NSTR reduces both together (§10.3).
4. **Do operator errors cancel or amplify in the forward model?**
   - Neither, for the solar and mixed channels. Wherever the effect matters,
     R_0v supplies essentially all of it, and the total equals the sum of the
     single-operator effects (median ratio 1.00).
   - In thermal channels the E_md and T_dv errors partly cancel, but the totals
     are ≤ 0.002 K anyway (§10.5).
5. **Do CORINT, DELTAM, ACCUR, moment treatment or precision give a
   meaningful improvement?**
   - **CORINT, DELTAM, moments:** the production settings are necessary.
     Disabling CORINT gives errors of up to 6σ away from backscatter and up to
     170σ elsewhere. Disabling it together with DELTAM, or truncating the phase
     moments to NSTR + 1, gives errors of thousands of σ.
   - **ACCUR** has no effect at any value from 0 to 10⁻², bitwise and in run
     time.
   - **Precision:** a meaningful improvement in robustness, not in accuracy.
     - Systematically, double precision at NSTR = 60 changes the measurement
       by ≤ 0.28σ away from backscatter and does not reduce the backscatter
       error.
     - But it removes the catastrophic single-call failures of §12.
     - It also removes the sporadic partial call failures in single
       precision (up to 17σ) seen at several other NSTR (§10.2).
6. **Is any improvement worth its cost?**
   - **Backscatter:** bringing it below 1σ needs NSTR ≳ 200 in double
     precision, about 60–90× the NSTR = 60 radiative-transfer cost (§14). That
     is not justified as a global change for a defect confined to
     the exact-backscatter LUT nodes and the Grid B cells adjacent to them.
   - **Failed calls:** removing them is worth a modest cost. Double
     precision at NSTR = 60 costs at most 3.8× (an upper bound, §14). A
     detect-and-recompute guard would cost almost nothing but catches only
     whole-call and extreme partial failures (§15).
7. **Is any current setting unnecessarily expensive?**
   - No. ACCUR never terminates the azimuthal series early, so it costs
     nothing.
   - Lower NSTR saves time (NSTR = 48 or 56: 0.46–0.79×) but degrades the
     backscatter region further.
8. **Should any DISORT setting be considered for a future version?** Yes, one:
   the floating-point precision, or an equivalent guard against failed calls
   (§15).
   - No other setting should change on accuracy-per-cost grounds.
   - Two items outside the DISORT option set deserve follow-up:
     - the exact-backscatter treatment;
     - the negative backscatter values produced by the 1000-moment Baum phase
       functions (§12.4).
9. **Which settings should remain?** NSTR = 60, DELTAM and CORINT on,
   ACCUR = 10⁻⁸, and the full adaptive (Mie) or 1000 (Baum) moment set.
   - Away from backscatter, the NSTR = 60 error in the complete measurement is
     ≤ 0.38σ (≤ 0.61σ for Baum 3.7 µm), and in the Jacobians ≤ 0.28σ per half
     cell.
   - Precision is the only open item.
10. **What was implemented?** Nothing. The candidate change, and the evidence
    it would need, are in §15.

## 1. Objective and decision criterion

**Question.** Does any DISORT option available to the ORAC LUT generator
change the final ORAC top-of-atmosphere simulated measurement, or its retrieval
Jacobians, by enough relative to the measurement uncertainty to justify its
cost?

**Not the decision basis.** An individual LUT operator (R_0v, R_dd, T_dd, T_dv,
E_md, …) is not the TOA measurement. The ORAC fast forward model combines
these operators with:
- surface reflectance operators;
- above- and below-cloud clear-sky transmittances;
- atmospheric emission, for thermal channels.

Operator differences are reported only as diagnostics (§13). An earlier draft
of this report treated R_0v as the TOA radiance; that is corrected here.

**Metric.** For every channel, cloud state, geometry and surface case:

  r = |y(candidate) − y(reference)| / σ_y,

where:
- y is the ORAC simulated measurement: reflectance for solar channels,
  brightness temperature for thermal and mixed channels, as in production;
- σ_y is the production ORAC measurement uncertainty (§8).

**Classification:**

| r | Class |
|---|---|
| < 0.1 | negligible |
| 0.1–0.5 | potentially material, but below the measurement uncertainty |
| 0.5–2 | comparable with the measurement uncertainty |
| ≥ 2 | clearly larger |

**Jacobians.** For dy/dlog₁₀τ and dy/dr_e, the error is reported as
|ΔJ| × (half the LUT cell width) / σ_y. That is the measurement error the
Jacobian error implies across half a cell, on the same scale as r.

## 2. Repository state and baseline

- **Start:** 2026-10-04 19:52 BST, host atmlxint5.
- **Commit:** branch `main` at `51717b82dece213d0f93ff1476ff1806983bbbe4`
  (the V24 production commit, equal to `origin/main`).
- **Uncommitted at the start, all untouched:** ` M AGENTS.md`,
  ` M tests/test_baum.py`, `?? documents/`,
  `?? scripts/submit_v23_slstr_cloud_luts.sh`, `?? tests/test_nakajima_king.py`,
  and 40 `?? runs/*_v23.run`.
- **V24 production** was running throughout. Those jobs execute this working
  tree, so no production path was modified.

**Baseline.** The production configuration at `51717b8`:
- V24 numerics: 3.5 r_e limit, adaptive Legendre expansion, refined radius
  grid;
- the production DISORT library;
- `nstreams = 60`.

The V24 expansion changes are part of every experiment, and only the DISORT
settings of §6 vary.

**Baseline reproduction.** The harness baseline, on a compact grid of Grid B
nodes, is **bitwise identical** to the finished V24 product
`terra_modis_m_water-ice_a01_pagg_v24.nc` at the matching nodes in every
radiative-transfer variable (R_0v, R_0d, R_dd, R_dv, T_dd, T_dv, T_0d, T_00,
E_md) and every optical property. The reproduction cases of §12 are bitwise
identical to their products as well; the 1.24 µm GHM case agrees to 1.5×10⁻⁶ in
R_0v (`check_baseline.py`).

This holds once numpy's AVX-512 kernels are disabled
(`NPY_DISABLE_CPU_FEATURES`):
- the production compute node uses the AVX2 path, whose float32 `exp`/`cos`
  differ in the last bit;
- single-precision DISORT turns that into differences of up to 3.6×10⁻⁴
  relative in R_dd (thin cloud, 0.65 µm);
- that is a host-dependent reproducibility floor of the production
  calculation, far below σ.

All experiments were run with AVX-512 disabled.

## 3. The production DISORT implementation (trace)

**Path, for each channel, τ and r_e:**

1. `create_orac_luts.create_orac_cloud_lut` builds the layer inputs.
   - Atmosphere 2 (`mls.atm`): 46 levels, 45 layers.
   - Cloud: one layer at 3.5 km (the `.mm` `profile`, `scatreltau`).
   - Rayleigh: column τ from the centre wavelength, phase moments from DISORT
     `GETMOM`.
   - V24 cloud LUTs have no gas absorption.
   - Moments: NMOM + 1 = max(L, NSTR + 1); L is adaptive for Mie classes and
     1000 for Baum.
2. `oraclut.idl_mirror.setup_disort` (mirror of `setup_disort.pro`) fills the
   `disort_vars` namespace. Only NSTR, NLYR, NUMU, NPHI and NMOM are used from
   it.
3. `oraclut.idl_mirror.call_disort` (mirror of `call_disort.pro`) checks
   dimensions, merges adjacent layers with identical SSALB and PMOM (exact),
   and calls `oraclut.radiative_transfer.legacy_disort.disort`.
4. That module compiles and loads `liboraclut_disort_gfortran.so` from:
   - `create_orac_lut/disort2/src/DISORTfunctions.f`, with BDREF, PRTFIN,
     ErrPack, LINPAK and RDI1MACH;
   - the C adapter `src/oraclut/radiative_transfer/legacy_disort.c`.

   Flags: gfortran `-O2 -fPIC -std=legacy -ffixed-line-length-0`.
5. **DISORT version.** DISORT 2.0 beta (Stamnes et al.; `disort/README.DISORT2.0`),
   localised for ORAC:
   - G. Thomas 2013: error returns instead of aborts; the UTAU check made
     relative;
   - G. McGarragh 2016: gfortran fixes.

   The numerics are those of DISORT 2.0: single precision (REAL, with some
   double-precision internals in the eigensolver) and plane-parallel.
6. **Calls per (channel, τ, r_e).** Every call returns TOA and BOA outputs
   (UTAU = 0 and the total τ) at all LUT view zeniths (±μ) and azimuths.

| Call | Settings | LUT outputs |
|---|---|---|
| Diffuse | FISOT = 100, FBEAM = 0 | R_dd, T_dd, R_dv, T_dv |
| Direct (per solar zenith, solar channels) | FBEAM = 100, μ₀ = cos SZA | R_0v, R_0d, T_0d, T_00 |
| Emission (thermal channels) | PLANK, in-cloud layers only at a uniform 250 K, ±0.5% wavenumber window | E_md = UU/B |

## 4. Inventory of DISORT options

The categories are:
- **A** — used in the production calculation;
- **B** — available in DISORT but irrelevant to this calculation;
- **C** — present in source or configuration but inactive;
- **D** — changeable only by changing production code (here investigated
  through validation copies, §6).

| Option | Production value | Where set | Category |
|---|---|---|---|
| NSTR (streams) | 60 | run file `nstreams` (all 40 V24 run files: 60); generator default 60; IDL hard-wired 60 | A, configurable; never varied in production |
| NMOM / PMOM (moments) | max(L, NSTR + 1) − 1; L adaptive for Mie, 1000 for Baum | `create_orac_luts` | A (the V24 Legendre change; not a DISORT switch) |
| DELTAM (δ-M scaling) | .TRUE. | hard-wired inside `DISORTfunctions.f` (line 622) | A; D to change |
| CORINT (Nakajima–Tanaka TMS/IMS intensity correction) | .TRUE. | hard-wired inside `DISORTfunctions.f` (line 623) | A; D to change |
| ACCUR (azimuthal-series convergence) | 1e-8 | C adapter (`legacy_disort.c`); `setup_disort` records it for reference only | A; D to change |
| Azimuthal terms | up to m = NSTR − 1, stopped early by ACCUR; one term (m = 0) for diffuse, emission and overhead-sun calls | DISORT | A |
| USRTAU, NTAU, UTAU | 1, 2, [0, total τ] | adapter / `create_orac_luts` | A |
| USRANG, NUMU, UMU | 1, 2 × number of VZA, ±cos VZA (μ = ±1 → ±0.99999) | `create_orac_luts` | A |
| NPHI, PHI | Grid B relative azimuths (DISORT convention, reversed to ORAC's) | `create_orac_luts` | A |
| IBCND | 0 (general boundary problem) | adapter | A; IBCND = 1 (albedo/transmissivity for many beams) C |
| FBEAM, UMU0, PHI0, FISOT | 100 or 0; cos SZA (90° → 89.99°); 0; 100 or 0 | `create_orac_luts`, adapter | A |
| LAMBER, ALBEDO | 1, 0 (black Lambertian; the surface is ORAC's) | adapter; `setup_disort(alb≠0)` raises | A; BRDF (LAMBER = 0, `BDREF`) B; `alb` argument C |
| PLANK, TEMPER, WVNMLO/HI | emission calls only; 250 K (cloud), 270 K (aerosol); ±0.5% | `create_orac_luts` | A |
| BTEMP, TTEMP, TEMIS | 0, 0, 0 (no surface or top emission) | adapter | A; non-zero C |
| ONLYFL | 0 (intensities computed) | adapter | A |
| PRNT | all off | adapter | B |
| Array limits | MXCMU = 100 (so NSTR ≤ 100), MXCLY = 46, MXUMU = 180, MXPHI = 181 | `DISORTfunctions.f` PARAMETER | A; raising needs a rebuild (D) |
| Floating-point precision | single (REAL) | DISORT 2.0 source | A; double needs a rebuild (D) |
| SSALB dither | 1 → 1 − 10 ε | DISORT | A |
| Pseudo-spherical geometry | not available (plane-parallel only) | – | not in DISORT 2.0 |
| δ-M+ and the improved intensity corrections of DISORT 3/4 | not available | – | not in DISORT 2.0 |
| Layering | 45 layers; identical Rayleigh layers merged (exact), so effectively Rayleigh / cloud / Rayleigh; DISORT is exact in τ within a layer | `create_orac_luts`, `call_disort` | A; not a numerical accuracy lever |
| Output angles | user angles are evaluated by DISORT's source-function integration, not interpolated; fluxes are bitwise independent of them (verified in this work) | – | A |

**Constraint on NSTR (found in this work).** DISORT stops with a fatal error
when the beam cosine coincides with a computational node:
`SETDIS: |μ₀ − μ_i|/μ₀ < 10⁻⁴`, followed by a Fortran STOP.
- Every NSTR ≡ 2 (mod 4) puts a node at μ = 0.5 (SZA 60°, a Grid B node).
- Several other values clash with other Grid B zeniths (64, 128 and 160 among
  them), and NSTR ≥ 240 clashes with SZA 0°.
- Usable on the full Grid B: 4, 8, …, 48, 56, 60, 68, 72, 76, 80, 88, 96, 100,
  104, 108, 112, 132, 136, 144, 148, 152, 156, 168, 172, 192, 196, 200, 204,
  216, 220, 224 and 236.
- The production value 60 is usable. Any future change must come from this list.

## 5. Production values

As in the table of §4:

| Setting | Production value |
|---|---|
| NSTR | 60 |
| DELTAM, CORINT | on |
| ACCUR | 1e-8 |
| Moments | adaptive L (Mie) or 1000 (Baum), at least NSTR + 1 |
| Surface | black Lambertian |
| Geometry | plane-parallel |
| Precision | single |

## 6. Experiments

All experiments use the production generator
(`create_orac_luts.create_orac_cloud_lut`) on compact grids of Grid B nodes.
The scattering properties are computed once by the production code and reused,
so only the radiative transfer differs.

**Varying the settings:**
- NSTR through the generator's own `nstreams` argument.
- The other settings through experimental copies of the production DISORT
  sources (`build_variants.py`, in `validation/tmp`). They differ only in:
  - MXCMU = 256 (100 for timing);
  - DELTAM/CORINT taken from a COMMON block;
  - a settable ACCUR;
  - optionally `-fdefault-real-8` for double precision.
- The MXCMU = 256 copy reproduces production bitwise at NSTR = 60
  (`kernel_check`).

**Configurations:**

| Group | Configurations |
|---|---|
| Streams | NSTR = 8, 16, 24, 32, 40, 48, 56, **60**, 68, 72, 80, 88, 96, 100, 112, 136, 152, 168, 200, 236 (all usable on Grid B, §4) |
| ACCUR at NSTR 60 | 10⁻², 10⁻³, 10⁻⁴, 10⁻⁶, 0 (production 10⁻⁸) |
| Intensity correction / scaling | CORINT off; DELTAM and CORINT off (DISORT 2.0 applies the correction only with δ-M) |
| Phase moments | truncated to NSTR + 1 (production: the full adaptive L, or 1000 for Baum) |
| Precision | double at NSTR 60, 136, 236 |
| Interactions | CORINT off at NSTR 136, 236; truncated moments at NSTR 136; ACCUR 10⁻³ at NSTR 32, 40 |
| Near backscatter (§10.4) | NSTR 60 (production), 100, 136, 200, 236, double 136, double 236 on SZA, VZA = 30, 33.75°, RAA 0–180° every 5° |
| Failed-call reproductions (§12) | production, double, NSTR 56, 68, 100 on small Grid B neighbourhoods of four failed calls; NSTR 100, 136, double 136 and CORINT off for the GHM 4.5 µm backscatter case |
| Timing (§14) | 3 repeats each of the production-library NSTR ladder up to 100, NSTR 136, the options and double precision |

**Reproducibility:**
- Every run writes `record.json` (settings, call counts, CPU times) next to its
  LUT under `validation/tmp/disort_options/runs/<case>/<experiment>/`.
- Command: `validation/v24/disort_options/launch.sh N < jobs.txt`, which runs
  `run_variant.py CASE EXPERIMENT [REPEAT]` in the environment of
  `validation/tmp/disort_options/env.sh` (minimal PATH, one thread, numpy
  AVX-512 off).

## 7. Test cases

Four cases cover both instrument families and all three V24 particle paths:

| Case | Instrument | Microphysics | Channels | Grid (all Grid B nodes) |
|---|---|---|---|---|
| modis_liquid | Terra MODIS | liquid stg (Mie, adaptive L) | 8, 3, 1, 2, 6, 7, 20 (mixed), 31, 32 | τ 0.25, 1, 4, 16, 64, 256; r_e 3, 7, 15, 25, 35 µm |
| modis_ice_agg | Terra MODIS | Baum aggregate (tabulated, nmom 1000) | as above | τ as above; r_e 5, 15, 33, 61, 93 µm |
| slstr_liquid | Sentinel-3A SLSTR | liquid stg | 1, 3, 5, 6, 7 (mixed), 8, 9 | as liquid |
| slstr_ice_sph | Sentinel-3A SLSTR | ice spheres (Mie) | as above | as ice |

All cases share the geometry SZA = VZA ∈ {0, 30, 60, 75}° and RAA ∈ {0, 30, …,
180}°.

- **Wavelength coverage.** The channels span weakly absorbing visible and NIR
  (0.41–0.87 µm), absorbing SWIR (1.6, 2.1–2.25 µm), the mixed 3.7 µm and the
  thermal window (11, 12 µm).
- **Phase functions.** Large droplets and ice particles (r_e up to 35 and 93 µm)
  give the most strongly forward-peaked phase functions in the V24 products.
- **Validation case.** This is the established compact validation framework of
  the V24 numerics work, extended to the V24 geometry limits (75°).
- **Mixed channels.** Treated as mixed: thermal radiance plus f₀ × solar
  reflectance, then converted to brightness temperature.

## 8. Measurement uncertainty: ORAC's own definition

Traced in ORAC `src/get_measurements.F90` lines 170–404 and
`src/read_sad_lut.F90` lines 1195–1238. The ORAC defaults for cloud retrievals
are `SySelm = SelmAux` and `Homog = Coreg = .true.`
(`src/read_driver.F90` lines 521–523).

**Solar channels.** σ_L = L/SNR (MODIS: the LUT's `snr`), or
σ_L = √(rua²L² + rub²L + ruc²) (SLSTR/AATSR: `rua`, `rub`, `ruc`, with
rub = ruc = 0 here). Both are relative. ORAC adds 1% (homogeneity) and 2%
(co-registration) of the signal in quadrature. For a non-positive measurement
it sets Sy = 10⁶ (lines 270–271), which is mirrored here.

**Thermal channels.** σ_BT = NEBT × (dB/dT at T₀)/(dB/dT at the measured BT),
with NEBT = the LUT's `nedt` and T₀ = `refbt`. ORAC adds 0.5 K and 0.15 K in
quadrature.

**Mixed channels (day).** As thermal, plus 1% and 2% of the total radiance
converted to BT (`USE_OLD_MIXED_UNCERTAINTY = .false.`, line 106).

These are the values ORAC actually uses. Results are also given with
instrument noise alone (no homogeneity or co-registration terms). For solar
channels this is stricter by a factor of 5–10 at high SNR (§10.1).

The project's earlier normalisation, `oraclut.validation.quadrature` (SAD NeFR
0.005, 1/SNR), is consistent with the instrument-noise part.

## 9. The ORAC fast forward model

### 9.1 McGarragh et al. (2018) formulation

Source: AMT 11, 3397–3431; Sect. 3.4–3.7. The paper was downloaded read-only
from the publisher (open access, SHA-256
`de81fafa0825e74a83f213be010fc152e5f4e0b650d982884b79dcd73ca46229`). It was not
stored in this repository.

**Cloud operators.** DISORT computes seven cloud operators per channel:
- reflection: R_bb, R_db, R_dd;
- transmission: T↓_bb, T↑_bb (the same table), T↓_bd, T↑_db;
- emissivity: ε.

Four surface operators (ρ_bb, ρ_bd, ρ_db, ρ_dd) come from the BRDF.

**Solar.**
- Eq. (36) is the series of cloud–surface reflections.
- Eq. (38) is its closed form, with 1/(1 − ρ_dd·R_dd).
- Eq. (39) adds below-cloud gas transmission. Diffuse paths use T_d =
  2∫T(θ) sinθ cosθ dθ (Eq. 46).
- Eq. (40) applies above-cloud transmission: R_TOA = T_ac(θ₀)·T_ac(θ_v)·R_TOC.

**Thermal.** Eq. (47):

L_TOA = L↑_ac + [L↓_ac·R_db + B(T_c)·ε + L↑_bc·T↑_db]·T_ac.

Mixed channels convert the solar reflectance to radiance and add it before the
conversion to brightness temperature.

### 9.2 The production implementation

Source: ORAC-CC/orac `eb9a1233c1136ef95fa0c37f650fe519e9134ed4` (master,
2026-08-03). It was cloned read-only to the session scratch area outside this
repository on 2026-10-04 21:24 BST. The tree was clean, and `git ls-remote`
confirmed upstream `master` at the same commit.

**Production choices for cloud retrievals with V2 NetCDF LUTs:**
- **One layer, overcast:** cloud fraction X(IFr) = 1 (`src/fm.F90` FM; lines
  287–288 call `FM_Solar` with `i_layer` = 0).
- **Solar, `i_equation_form` = 3** (default for cloud classes,
  `read_driver.F90` line 620) with **full BRDF** (`use_full_brdf = .true.`,
  line 481). `src/fm_solar.F90` lines 883–945:

      a = 1 − ρ_dd R_dd Tbc_dd
      b = T_00 ρ_0d Tbc_0 + T_0d ρ_dd Tbc_d
      c = T_dv Tbc_d + R_dd ρ_dv T_vv Tbc_dd Tbc_v
      e = R_0v + (T_00 ρ_0v T_vv Tbc_0v + T_0d ρ_dv T_vv Tbc_dv + b c / a) / sec θ₀
      Ref = Tac_0v e

- **Transmittances:** Tac_0 = Tac^secθ₀ and Tac_v = Tac^secθ_v, from the RTTOV
  nadir transmittance at the cloud pressure. Diffuse paths use the nadir value
  itself (`dif_trans_fac` = 1, line 643), not Eq. (46). Likewise Tbc.
- **Reciprocity-obeying form.** Equation form 3 is algebraically Eq. (39)–(40)
  for a Lambertian surface, once normalisation is accounted for:
  - ORAC's Ref, like the LUT's R_0v, includes cos θ₀ (πL/F₀);
  - the paper's reflectance is μ₀-normalised;
  - hence the surface terms are multiplied by cos θ₀.

  Verified numerically to 10⁻¹⁰ (`check_fm.py`).
- **Thermal:** `src/fm_thermal.F90` lines 253–260, then `R2T`
  (`src/r2t.F90`):

      R = (Rbc_up T_dv + B(T_c) E_md + Rac_dwn R_dv) Tac + Rac_up,  BT = R2T(R)

  The production form uses the diffuse reflectance R_dv (the paper's R_db) and
  transmission T_dv, as in Eq. (47).
- **Mixed (3.7 µm):** `src/fm.F90` lines 348–355:
  BT = R2T(R_thermal + f₀·Ref), with f₀ = the LUT's `F0` (= F/π, in the band's
  Planck units).
- **Measured quantities:** reflectance for solar channels, brightness
  temperature for thermal and mixed channels.
- **LUT interpolation:**
  - (θ₀, θ_v, Δφ) trilinear at the 4 × 4 (τ, r_e) stencil
    (`src/int_lut_tausatsolazire.F90` lines 130–176);
  - then bicubic in (log₁₀τ, r_e) (`src/int_lut_routines.F90` lines 163–213).
    This is the default `LUTIntSelm` (`read_driver.F90` line 624) in the conda
    build, which defines `INCLUDE_NR` (`config/lib.conda.inc` lines 37–39);
  - corner derivatives are centred finite differences on the LUT grid,
    one-sided at its edges (`src/set_gzero.F90` lines 95–174);
  - τ is log₁₀ because our LUTs carry `spacing = uneven_logarithmic`
    (`read_sad_lut.F90` lines 870–878).
- **Jacobians:** analytic, from the interpolated operator gradients through the
  chain rule (`fm_solar.F90` lines 269–286 for form 3; `fm_thermal.F90` lines
  268–275; `fm.F90` lines 356–360 for mixed channels). Retrieval state: X(ITau)
  = log₁₀τ and X(IRe) = r_e.
- **Surface operators:** from the pre-processor, `ross_thick_li_sparse_r.F90`
  (land, Ross-thick/Li-sparse-R, 4 × 4 Gauss quadrature, lines 524–804) or
  `cox_munk.F90` (sea). ORAC divides them by a solar factor only for
  equation forms 2 and 4 (`src/get_surface.F90` lines 221–236).

**Mapping LUT variable → ORAC → McGarragh** (`src/read_sad_lut.F90` lines
979–1103; `src/fm_routines.F90` lines 32–41; `src/set_crp_solar.F90` lines
128–173; `src/set_crp_thermal.F90` lines 102–110):

| LUT (this generator) | ORAC array / index | McGarragh | Interpolated in | Used by |
|---|---|---|---|---|
| R_0v | Rbd, IR_0v | R_bb(θ₀,θ_v,Δφ) | τ, θ_v, θ₀, Δφ, r_e | solar |
| R_dd | Rfd, IR_dd | R_dd | τ, r_e | solar |
| T_00 | Tb, IT_00 | T↓_bb(θ₀) | τ, θ₀, r_e | solar |
| T_00 (same table at θ_v) | Tb, IT_vv | T↑_bb(θ_v) | τ, θ_v, r_e | solar |
| T_0d | Tfbd, IT_0d | T↓_bd(θ₀) | τ, θ₀, r_e | solar |
| T_dv | Td, IT_dv | T↑_db(θ_v) | τ, θ_v, r_e | solar, thermal |
| R_dv | Rd, IR_dv | R_db(θ_v) | τ, θ_v, r_e | thermal (and 2-layer solar) |
| E_md | Em, IEm | ε(θ_v) | τ, θ_v, r_e | thermal |
| R_0d, T_dd | Rfbd, Tfd | – | – | two-layer forms only |

`get_T_dv_from_T_0d` (`read_driver.F90` lines 1278–1297) replaces T_dv by T_0d
only for the legacy 3-character text-LUT classes. For the NetCDF V2 LUTs
(`read_sad.F90` lines 62–80) the LUT's own T_dv is used, as here.

### 9.3 The validation implementation

`validation/v24/disort_options/orac_fm.py` mirrors the routines of §9.2 line by
line, with a citation for each formula.

**Checks** (`check_fm.py`, all pass):
- bicubic node values are exact;
- gradients at nodes equal ORAC's corner differences;
- the bicubic surface is continuous across cells;
- analytic solar and thermal Jacobians equal finite differences of the full
  interpolate-then-forward-model chain (to 8×10⁻¹¹);
- a black surface gives Tac₀ᵥ·R_0v;
- a Lambertian surface reproduces the paper's Eqs. (39)–(40);
- the pre-processor BRDF is reciprocal (ρ_0d(θ) = ρ_dv(θ)).

**Clear-sky terms** (`atmosphere.py`). ORAC obtains these from RTTOV.
- **Sources:** the generator's own MODTRAN mid-latitude-summer gas optical
  depths (`input_files/gas/ModtranGasOpd_A2_*`), the `mls.atm` temperature
  profile and the LUT's band Planck coefficients.
- **Cloud top:** 4 km (liquid) and 9 km (ice).
- **Surface:** emissivity 0.98 at the profile surface temperature.
- **Treatment:** non-scattering Schwarzschild sums along the view.
- They are identical for every DISORT configuration, so only their magnitude
  matters.

**Surfaces:**

| Surface | Definition |
|---|---|
| black | diagnostic only |
| dark | Lambertian, 0.05 (ocean-like albedo; ORAC's Cox–Munk sea BRDF was not ported) |
| vegetated land | ORAC pre-processor Ross-thick/Li-sparse-R with representative MODIS-like kernel weights by band; ρ = 0.03 at 3.7 µm, following the pre-processor's 1 − emissivity |
| bright, snow-like | Lambertian, 0.90 / 0.88 / 0.80 / 0.10 / 0.05 at 0.5 / 0.7 / 0.9 / 1.6 / 2.1 µm |

No ORAC reference surface dataset exists in this repository, so these are
representative, not reference, surfaces.

**States and geometries.**
- Cloud states: all 30 (τ, r_e) nodes and all 20 cell mid-points in
  (log₁₀τ, r_e).
- Geometries: all 112 grid geometries, and 148 in the near-backscatter set.
- Interpolation and Jacobians: ORAC's bicubic method.
- SLSTR: the oblique-view channels are copies of the nadir operators made by
  the generator, so only nadir channels are evaluated.

## 10. Results: the complete ORAC measurement

Source files are in `validation/v24/results/disort_options/`:
- `fm_measurement.csv`: per case, experiment, comparison, surface, channel and
  σ definition;
- `fm_ladder.csv` and `fm_breakdown.csv`: by scattering-angle class;
- `fm_jacobians.csv` and `fm_jacobian_ladder.csv`;
- `fm_cancellation_baseline.csv`;
- `fm_arrays_*.npz`.

In this compact grid, scattering angles are either < 170° or exactly 180°
(SZA = VZA, RAA = 0; 4 of 112 geometries). §10.4 adds the intermediate angles.

### 10.1 NSTR = 60 against the reference

Production σ. Dark, vegetated and snow surfaces (the black surface behaves as
dark). Worst point over channels, states and geometries.

| Case | Solar, Θ < 170° | Solar, Θ = 180° | Mixed 3.7 µm, Θ < 170° | Mixed 3.7 µm, Θ = 180° | Thermal |
|---|---|---|---|---|---|
| MODIS liquid | 0.24 | 4.96 | 0.23 | 1.50 | 0.001 |
| MODIS Baum ice | 0.38 | 380 (14 with signal ≥ 0.01) | 0.61 | 6.39 | 0.002 |
| SLSTR liquid | 0.29 | 4.90 | 0.26 | 1.66 | 0.001 |
| SLSTR ice spheres | 0.30 | 6.47 | 1.10 | 4.23 | 0.002 |

**Share of evaluated points above each level** (`fm_breakdown.csv`):
- Θ < 170°: no point above 0.5σ for any solar or thermal channel, case or
  surface. For the ice-sphere 3.7 µm channel, 1.7% are above 0.5σ and none
  above 1σ.
- Θ = 180°: 10–70% above 0.5σ and up to 39% above 2σ, varying by channel. For
  example, SLSTR ice spheres at 1.6 µm, dark surface: 72% above 0.5σ, 39% above
  2σ.

**Reference uncertainty.** A point counts as "resolved" where the NSTR = 60
error exceeds twice |double NSTR 136 − reference|. The resolved shares are
nearly the same as the raw shares, so the backscatter differences are not
reference noise.

**Instrument noise alone.** With σ from SNR or `rua` only (no 1% and 2%
terms), worst over channels and surfaces:

| Case | Solar, Θ < 170° | Solar, all Θ (signal ≥ 0.01) | Thermal |
|---|---|---|---|
| MODIS liquid | 2.32 (0.86 µm) | 56 | 0.01 |
| MODIS Baum ice | 2.05 | 102 | 0.02 |
| SLSTR liquid | 1.51 (0.87 µm) | 26 | 0.03 |
| SLSTR ice spheres | 1.50 | 32 | 0.04 |

For the mixed 3.7 µm channel with NEBT alone (no 0.5 K and 0.15 K terms, no
1% and 2% radiance terms), NSTR = 60 against the reference reaches:
- away from backscatter: 4.8σ (MODIS liquid), 4.0σ (Baum), 10.9σ (SLSTR
  liquid) and 10.8σ (ice spheres);
- overall: 17–81σ.

ORAC's 1% and 2% terms dominate the solar σ at high SNR, which is why the
production-σ values are 5–10 times smaller. The conclusion "converged away
from backscatter" holds for the uncertainty ORAC uses. It would not hold if
ORAC's solar σ were reduced towards instrument noise.

**Where the extreme Baum values sit.** The values of order 100σ occur at 2.1 µm
at exact backscatter, for r_e = 93 µm. There the reflectance is ≈ 0 or
negative, so the relative σ → 0. At exactly this point the 1000-moment Legendre
expansion of the Baum aggregate phase function that DISORT receives is negative
at 180° (P = −0.12, against +0.08 at 178°, §12.4). These values do not drive
the conclusions; with signal ≥ 0.01 the Baum maximum is 14σ.

### 10.2 Stream-count convergence

`fm_ladder.csv`; `figures/fm_stream_convergence.png`. Worst over channels and
surfaces, production σ.

**Away from backscatter (Θ < 170°).**
- The single-precision ladder is flat from NSTR = 32 upwards. For MODIS liquid
  it is 0.26, 0.20, 0.16, **0.24**, 0.29 and 0.41σ at NSTR = 32, 48, 56,
  **60**, 80 and 100.
- At several other NSTR, single precision shows isolated jumps: NSTR 24 (Baum,
  9σ), 40 (ice spheres, 17σ), 136 (MODIS liquid, 10σ), 152 (ice spheres, 7σ)
  and 168 (SLSTR liquid, 3.6σ).
  - Each is confined to one channel, one solar zenith and one cloud-state
    node (with the interpolated mid-points around it): 8–49 of 5,600 points
    per channel and surface, spread over several view zeniths.
  - So each is a partial failure of a single direct-beam call: some of its
    views are wrong, but too few to move the call's median over views (no
    call differs from the reference by more than 0.01 in median).
  - None occurs at NSTR = 60 on these grids.
  - They are not explained by a view direction lying close to a quadrature
    node (checked); the mechanism was not identified.
- The double-precision runs show nothing of the kind (NSTR 136: ≤ 0.27σ).
- Single-precision DISORT 2.0 is therefore unreliable at the level of
  individual calls. The whole-call failures of §12 show the same weakness at
  NSTR = 60.

**Exact backscatter (Θ = 180°)**: slow, monotonic convergence.

| NSTR | 16 | 32 | 48 | **60** | 80 | 100 | 136 | 168 | 200 | 236 (single) |
|---|---|---|---|---|---|---|---|---|---|---|
| MODIS liquid solar | 8.6 | 6.6 | 5.6 | **5.0** | 4.0 | 3.2 | 2.6 | 1.2 | 0.59 | 0.49 |
| SLSTR ice spheres solar | 10.5 | 8.7 | 7.3 | **6.5** | 5.3 | 4.2 | 3.4 | 1.6 | 2.6 | 0.66 |
| MODIS liquid mixed 3.7 µm | 0.76 | 1.04 | 1.36 | **1.50** | 1.46 | 0.91 | 0.11 | 0.03 | 0.08 | 1.24 |

Double precision gives 4.96 (NSTR 60) and 2.05 (NSTR 136) for MODIS liquid.

**Thermal channels.** ≤ 0.002σ for all NSTR ≥ 32. NSTR = 16 already gives
0.07σ.

### 10.3 Retrieval Jacobians

`fm_jacobian_ladder.csv`. |ΔJ| × half-cell / σ, NSTR = 60 against the
reference, worst over channels and surfaces:

| Case | dlog₁₀τ, Θ < 170° | dlog₁₀τ, all | dr_e, Θ < 170° | dr_e, all |
|---|---|---|---|---|
| MODIS liquid | 0.17 | 3.26 | 0.09 | 0.79 |
| MODIS Baum ice | 0.23 | 119 (signal ≈ 0) | 0.05 | 20 (signal ≈ 0) |
| SLSTR liquid | 0.28 | 3.20 | 0.18 | 0.77 |
| SLSTR ice spheres | 0.24 | 3.57 | 0.14 | 1.83 |

- Mixed channels away from backscatter: ≤ 0.08–0.26σ (dlog₁₀τ) and
  ≤ 0.03–0.33σ (dr_e).
- Thermal: ≤ 0.001σ.
- So the Jacobians are converged wherever the measurement is, and they
  converge together with it at backscatter. No region was found where the
  measurement is accurate but a Jacobian is not.

### 10.4 Near backscatter

`near_backscatter.py` → `near_backscatter/near_backscatter.csv` and
`extent.csv`; `figures/near_backscatter.png`.

**Grid.** Grid B nodes SZA, VZA ∈ {30, 33.75}° and RAA every 5° from 0 to 180°
(148 geometries). Scattering angles cover 150–180°; the closest to exact
backscatter are 176.25° and about 177.5°.

**NSTR = 60 against the reference.** Production σ, worst over channels,
surfaces and states:

| Θ (°) | < 170 | 170–174 | 174–176 | 176–177 | 177–178 | 180 |
|---|---|---|---|---|---|---|
| MODIS liquid, solar | 0.09 | 0.06 | 0.18 | 0.19 | 0.31 | 3.98 |
| SLSTR liquid, solar | 0.05 | 0.07 | 0.19 | 0.17 | 0.29 | 3.90 |
| SLSTR ice spheres, solar | 0.04 | 0.06 | 0.17 | 0.30 | 0.22 | 5.84 |
| MODIS Baum, solar | 0.33 | 0.52 | 0.39 | 0.41 | 0.70 | 69 |
| MODIS liquid, mixed 3.7 µm | 0.02 | 0.06 | 0.25 | 0.26 | 0.35 | 0.75 |
| SLSTR ice spheres, mixed 3.7 µm | 0.08 | 0.19 | 0.70 | 0.82 | 0.98 | 1.94 |
| MODIS Baum, mixed 3.7 µm | 0.16 | 0.20 | 0.36 | 0.75 | 0.56 | 6.39 |

- **Width of the defect.**
  - For droplets and ice spheres in solar channels the NSTR = 60 error stays
    ≤ 0.31σ at every resolved angle up to 177.5°, and exceeds 0.5σ only at
    Θ = 180° itself.
  - The defect is therefore a peak narrower than 2.5°: the glory and
    backscatter peak that 60 streams cannot resolve.
  - The mixed 3.7 µm channel rises earlier, from about 175°, but stays ≤ 1σ
    below 180°.
- **In ORAC.** The LUT holds this peak only at its exact-backscatter nodes
  (SZA = VZA, RAA = 0). ORAC interpolates trilinearly in angle, so a node error
  of 4–6σ enters every measurement within one Grid B cell of exact backscatter
  (5° in RAA, 3.75° in SZA and VZA), falling linearly to zero across the cell.
  The affected region of a retrieval is therefore set by the LUT grid, not by
  the physical width of the peak. A 5° grid cannot represent the real glory
  peak either, at any NSTR.
- **Convergence at 180°.** Liquid: 3.98 (NSTR 60) → 2.57 (100) → 1.68 (double
  136) → 0.57 (200), confirming §10.2 on this denser set. Baum stays large
  (69 → 15 at NSTR 200), where the phase-function expansion is negative at 180°
  (§12.4).
- **Single precision at high NSTR.** The single-precision NSTR 136 runs show
  isolated 2–5σ jumps at all angles (dashed curves in the figure), the partial
  call failures of §10.2; double NSTR 136 does not.

### 10.5 Operator cancellation or amplification

`fm_cancellation_baseline.csv`, via `cancellation.py`. For each point, the
NSTR = 60 operators are substituted one at a time into the reference set and
compared with substituting all of them together:

| Channels | Dominant operator where \|ΔF\| ≥ 10% of its maximum | Median \|ΔF_total\| / Σ_k \|ΔF_k\| | Min | Comment |
|---|---|---|---|---|
| Solar, all cases and surfaces | R_0v (100% of points) | 1.00 | 0.64–1.00 (snow), ≥ 0.94 otherwise | the R_0v error passes unchanged, scaled by Tac₀ᵥ |
| Mixed 3.7 µm | R_0v (100%) | 1.00 | ≥ 0.99 | as solar |
| Thermal 11, 12 µm | E_md (70–88%), else T_dv | 1.00 | 0.14–0.36 | E_md and T_dv errors partly cancel; totals ≤ 0.002 K |

The forward model therefore neither hides nor magnifies the DISORT stream
error in solar channels:
- the surface-coupling operators (T_00, T_0d, T_dv, R_dd) change by
  ≤ 5.4×10⁻⁴ (§13);
- their contribution is ≤ 0.5–3.5% of the R_0v contribution, even over a
  0.8–0.9-albedo surface.

### 10.6 Other DISORT options at NSTR = 60

Worst over all four cases, production σ (`figures/fm_options.png`):

| Option | Against the reference: Θ < 170° | All Θ (signal ≥ 0.01) | Against production: Θ < 170° | Verdict |
|---|---|---|---|---|
| Production (NSTR 60) | 0.38 | 14.4 | – | – |
| ACCUR 10⁻², 10⁻³, 10⁻⁴, 10⁻⁶, 0 | identical to production | identical | 0 (bitwise) | no effect: the azimuthal series is never truncated early at these geometries; no run-time change |
| CORINT off | 6.1 | 167 | 6.2 | far worse; essential |
| DELTAM and CORINT off | 1,219 | 15,854 | 1,219 | unusable |
| Moments truncated to NSTR + 1 | 1,714 | 21,067 | 1,714 | unusable; the full expansion is essential |
| Double precision, NSTR 60 | 0.38 | 14.5 | 0.28 (99th percentile 0.07) | no systematic change; removes the failed calls of §12 |
| CORINT off, NSTR 136 / 236 | 10.0 / 4.8 | 175 / 163 | – | more streams do not substitute for the correction |
| Moments NSTR + 1, NSTR 136 | 2,148 | 36,538 | – | unusable |
| ACCUR 10⁻³ at NSTR 32 / 40 | 1.8 / 17 | 18 / 17 | – | bitwise equal to ACCUR 10⁻⁸ at the same NSTR; the 17σ is the NSTR = 40 single-precision jump of §10.2, not an ACCUR effect |

### 10.7 Nakajima–King space

Figures `figures/nk_fm_*.png` show the complete measurement for NSTR = 60
(solid) and the reference (dashed). They label the τ and r_e grid values, show
the production σ as error bars, and give the displacement in units of σ.

- **Away from backscatter** (MODIS liquid 0.86/2.1 µm over vegetation, SZA =
  VZA = 30°, RAA = 90°): the two grids are indistinguishable, and the largest
  two-channel displacement is 0.014σ. The Baum and SLSTR liquid panels at the
  same geometry look the same.
- **At exact backscatter** (same case, RAA = 0°): the NSTR = 60 grid is
  displaced by up to 3.6σ, mainly along the non-absorbing (0.86 µm) axis, and
  grows with r_e. In a retrieval this would appear mainly as an optical-depth
  bias, with a smaller r_e bias at large r_e.

## 11. The numerical reference and its limits

### 11.1 Choice

Double precision, NSTR = 236:
- the highest usable NSTR on Grid B (§4);
- MXCMU = 256 copy of the production DISORT;
- all other settings as production.

### 11.2 Convergence evidence

Measured against the reference:
- double NSTR 136 is within 0.27σ away from backscatter;
- single NSTR 200 and 236 are within 0.59σ and 0.49σ at exact backscatter
  for liquid;
- the backscatter sequence is still decreasing at NSTR 200–236.

So the reference is converged to well below σ away from backscatter, and to
about 0.5σ at exact backscatter.

### 11.3 Non-converging oscillation (ice spheres, 1.6 and 3.7 µm)

For large ice spheres the 3.7 µm brightness temperature oscillates with NSTR
at a few cloud states:
- range up to 1.8 K over NSTR = 32–236, 99th percentile 0.87 K, median
  0.04 K;
- example: τ = 32, r_e = 47 µm, Θ = 120°, a sphere-rainbow region;
- the 1.6 µm reflectance oscillates by up to 8% in the same way.

Single and double precision agree at fixed NSTR, so this is the discrete-ordinate
solution itself, not rounding, and no NSTR in range converges it. The ~1σ
NSTR = 60 difference there is within this oscillation, not an NSTR = 60
deficiency.

### 11.4 Reference breakdown

The double-precision NSTR 236 run returned NaN at one (channel, τ, r_e) point:
Baum, 0.86 µm, τ = 64, r_e = 5 µm (112 R_0v values and a few diffuse values).
- This is a numerical breakdown for nearly conservative scattering in the
  `-fdefault-real-8` build; DISORT's single-scattering-albedo dither is derived
  from machine ε.
- That point and the bicubic stencils it touches are excluded and counted
  (`n_nonfinite_excluded`).
- The double-precision runs at NSTR 60 and 136 have no non-finite value in any
  case.

### 11.5 Plane-parallel limitation

The reference is plane-parallel DISORT 2.0, like production. Sphericity
(relevant above ~75°) and polarization are outside both.

## 12. Failed single-precision DISORT calls in production products

### 12.1 What fails

One direct-beam DISORT call produces R_0v at every view zenith and relative
azimuth for one (channel, τ, r_e, solar zenith). In a failed call the whole
block is wrong while the neighbouring calls are correct. Four failures were
reproduced on small Grid B neighbourhoods (`run_variant.py` cases `defect_*`):

| Product (V24) | Call | Production R_0v range | Double, NSTR 56, 68, 100 (all agree) | Reproduction of the product |
|---|---|---|---|---|
| sentinel-3a liquid p240 | 0.66 µm (ch 2), τ 64, r_e 37, SZA 60° | −4.82 to 5.71 | 0.38–0.68 | bitwise |
| aqua MODIS Baum agg | 4.47 µm (ch 24), r_e 9, SZA 60°, τ ≥ 16 | −0.10 to 0.28 | 0.010–0.16 | bitwise |
| terra MODIS Baum GHM (= V23) | 1.24 µm (ch 5), τ 11.3, r_e 21, SZA 75° | −3.93 to 4.58 | 0.13–0.81 | 1.5×10⁻⁶ |
| terra MODIS Baum GHM (= V23) | 0.905 µm (ch 17), τ 0.5, r_e 85, SZA 60° | −0.16 to 0.49 (−0.10 to 0.49 at the reproduced views) | 0.028–0.30 | bitwise |

In every case:
- only the production configuration, single precision at NSTR = 60, deviates;
- the other direct-beam outputs of the same call (R_0d, T_0d, T_00) agree
  with double precision to ≤ 10⁻⁴;
- the neighbouring calls agree with double precision to the normal
  single-precision level.

**Effect on the ORAC measurement** (`defect_impact.csv`). Over a black
surface, equation form 3 reduces to Ref = Tac₀ᵥ·R_0v, and ORAC's solar σ is
relative, so r = |ΔR_0v| / (σ_rel·R_0v) exactly:

| Call | Views with r > 2 | Median r | Max r |
|---|---|---|---|
| SLSTR 0.66 µm | 95% | 104 | 602 |
| MODIS GHM 1.24 µm | 63% | 244 | 856 |
| MODIS GHM 0.905 µm | 72% | 13 | 116 |
| MODIS agg 4.47 µm | – (no ORAC σ: SNR = 0, not a retrieval channel) | – | – |

These values are at the LUT node itself. ORAC's bicubic interpolation in
(log₁₀τ, r_e) uses finite-difference derivatives over neighbouring nodes, and
the angular interpolation is trilinear. So a failed node also corrupts the
interpolated measurement and Jacobians in the adjacent τ, r_e and solar-zenith
cells.

### 12.2 Prevalence in V23 and V24 products

`corrupted_calls.py` → `corrupted_calls.csv` (per product) and
`corrupted_calls_isolated.csv` (per flagged cell).

- **Detector.** A cell (channel, τ, r_e, SZA) is flagged when it is
  inconsistent with its two solar-zenith neighbours:
  - score = median over (VZA, RAA) of |R_0v − neighbour mean|, divided by the
    median neighbour value, above 0.2;
  - median absolute anomaly above 0.002;
  - the score is a local maximum in SZA (a failed call's neighbours inherit
    about half its score).
- **Validation of the detector.**
  - All four reproduced failures are found (scores 0.25–2.8).
  - Smooth structure, including the GHM large-r_e ringing of §12.4, scores
    ≤ 0.13.
  - Cells with scores between 0.15 and 0.2 are listed too, to show the
    threshold's sensitivity.
- **Excluded:** SLSTR oblique channels that repeat a nadir channel's
  computation exactly.

**Results** (75 products: all 40 V23, and the 35 of 37 V24 products written
when the census ran):

| | V23 | V24 |
|---|---|---|
| Products with at least one whole-call failure | 17 of 40 | 14 of 35 |
| Failed cells (channel, τ, r_e, SZA) | 44 | 32 |
| Distinct (channel, r_e, SZA); one failure can repeat across saturated τ | 28 | 17 |
| Share of direct-beam cells | 5.1×10⁻⁶ | 4.6×10⁻⁶ |
| Most in one product | 12 (aqua MODIS Baum agg, 4.47 µm) | 12 (the same, unchanged) |
| Partial failures with R_0v outside [−0.1, 1.5], not detected as whole-call | 3 | 2 |
| Cells with minimum in [−0.1, −0.01) (mostly §12.4) | 5,459 | 5,449 |

- **Channels.** Most failures are in strongly absorbing or non-retrieval MODIS
  bands: 4.47 µm (26 cells), 3.96–4.05 µm (14), 1.38 µm and 4.52 µm. Some are
  in ORAC retrieval channels:
  - MODIS 0.86 µm (3 cells);
  - SLSTR 0.66 µm (2), 1.6 µm (4), 2.25 µm (2) and 3.7 µm (2);
  - every particle family (droplets, ice spheres, Baum).
- **A range check alone is not enough.** Of the 76 detected cells, only 33 have
  any R_0v outside [−0.1, 1.5]; the other 43 have physically plausible values
  and would pass it.
- **Partial failures.** Some calls are wrong in only part of their views. In
  aqua MODIS liquid p263 V24, at 1.24 µm, τ 2.8, r_e 23, SZA 67.5°, 94 of 777
  views lie between −2.36 and 3.22 while the rest are correct. The five such
  calls with values outside [−0.1, 1.5] have 9–94 wrong views each.
- **Every R_0v value outside [−0.1, 1.5]** in these products lies in a whole or
  partial failed call.

**Limitation.** The detector compares medians over all views, so it cannot see
partial call failures (the type of §10.2), and they do occur at NSTR = 60 in
the products (above). They are counted only when they leave the physical range.
A view-level detector was tried and rejected: it flagged thousands of cells of
legitimate angular structure. The census is therefore a lower bound on the
number of affected calls.

### 12.3 Interpretation

- The failures are a single-precision weakness of DISORT 2.0. They depend
  sensitively on the inputs (same NSTR, same inputs, double precision: correct)
  and are not a stream-count error (NSTR 56 and 68 are correct).
- The same weakness produces the partial call failures at other NSTR
  (§10.2).
- They are rare per call but catastrophic where they occur, and they recur
  across products: V23 and V24 differ in which calls fail where the inputs
  differ, and agree where the inputs are unchanged (the Baum classes).

### 12.4 A separate effect: negative backscatter from the Baum phase-function expansion

The Baum GHM products contain about 2,700 further cells whose R_0v minimum lies
between −0.1 and −0.01. Examples are MODIS channels 21–25 (3.96–4.52 µm) at
r_e 81–93 µm, at every solar zenith and every τ ≥ 0.5. These are not failed
calls:
- they are smooth in SZA, τ and r_e (detector score ≤ 0.13);
- they sit at exact backscatter (SZA = VZA, RAA = 0) on a signal of about
  0.005;
- they are unchanged in double precision, converge to −0.0127 at NSTR = 136,
  and change only when CORINT is switched off (`ringing_modis_ghm24`).

The cause is the input phase function:
- the 1000-moment Legendre expansion of the Baum GHM phase function at 4.47 µm
  (g ≈ 0.95) gives P(180°) = −0.14 to −0.17 for r_e = 85–93 µm, where the
  tabulated phase function is about +0.03;
- DISORT 2.0's intensity correction evaluates the single-scattering term from
  that expansion;
- the same expansion is negative at 180° for the Baum aggregate at 0.86 µm,
  1.6 µm, 2.1 µm and 3.7 µm at r_e = 93 µm (−0.02 to −0.12), which is where the
  extreme Baum backscatter values of §10.1 sit.

This is a scattering-input truncation effect (nmom = 1000), outside the DISORT
option set. It affects exact backscatter in large-particle Baum classes only.

## 13. Operator-level diagnostics

Retained from the first part of this investigation, but not decisive.

### 13.1 NSTR = 60 against the reference, LUT operators

`lut_operators.csv`, `operator_r0v_emd.csv`:

| Case | R_0v max abs (max rel) | R_0d, R_dd, R_dv, T_dd, T_dv, T_0d max abs | T_00 max abs | E_md max abs (ΔBT bound) |
|---|---|---|---|---|
| MODIS liquid | 2.8×10⁻² (11%) | ≤ 2.2×10⁻⁴ | 1.2×10⁻⁷ | 2.2×10⁻⁵ (0.0006 K) |
| MODIS Baum ice | 5.5×10⁻³ (at R ≈ 0) | ≤ 5.4×10⁻⁴ | 1.8×10⁻⁷ | 3.1×10⁻⁵ (0.002 K) |
| SLSTR liquid | 2.7×10⁻² (12%) | ≤ 2.4×10⁻⁴ | 1.8×10⁻⁷ | 2.0×10⁻⁵ (0.0005 K) |
| SLSTR ice spheres | 3.0×10⁻² (17%) | ≤ 2.1×10⁻⁴ | 1.2×10⁻⁷ | 2.8×10⁻⁵ (0.0014 K) |

The R_0v maxima are all at exact backscatter.

### 13.2 Reproducibility floor

Single-precision DISORT turns last-bit changes of its inputs into up to
3.6×10⁻⁴ relative change in the diffuse operators (§2).

### 13.3 Negative R_0v in current products

Census of the V23 and V24 cloud products available when it ran (65 products):
`negative_r0v.csv`, `negative_radiance.py`.
- 0.2–0.3% of values are negative at all.
- Values below −10⁻⁴ number up to 30,545 per product (terra MODIS Baum GHM,
  V23).
- The minimum is −67, in a failed call (aqua MODIS ice spheres, V23).
- They occur in absorbing channels (1.6, 2.1–2.25, 3.7–4.5 µm). For Baum they
  sit mostly at Θ ≥ 175°; for spheres and droplets, at specific r_e and
  Θ ≈ 55–100°.
- The V24 counts equal the V23 counts for the unchanged Baum classes.

Three sources are now identified:
- failed single-precision calls, whole or partial (§12.1–12.2), which account
  for every value below −0.1 or above 1.5;
- the Baum phase-function expansion at backscatter (§12.4);
- small negative values at other angles in strongly absorbing,
  forward-peaked cases. Their origin was not investigated further; they
  include DISORT 2.0 intensity-correction artefacts.

A TOA radiance cannot be negative. ORAC treats a non-positive *measurement*
specially, but nothing protects against a negative *forward-model* value.

## 14. Runtime

`timing.py` → `runtime.csv`; `figures/runtime_vs_streams.png`.
- **Method:** DISORT CPU time per compact case, median of 3 repeats on the
  same host, one thread, launched together.
- **Kernel:** the production library up to NSTR 100; the MXCMU = 256 copy
  above. The latter's constant array-zeroing overhead (measured as NSTR 100
  with MXCMU 256 minus NSTR 100 with MXCMU 100, 590–850 s) is subtracted.
- **DISORT share:** the direct-beam calls are 96–99% of DISORT time. DISORT
  is almost all of the harness run time; in production the scattering
  calculation comes on top.

| NSTR | 32 | 48 | 56 | **60** | 68 | 80 | 100 | 136 | 200 | 236 |
|---|---|---|---|---|---|---|---|---|---|---|
| Time / NSTR 60 (repeated, 4 cases) | 0.19–0.20 | 0.45–0.48 | 0.73–0.79 | **1** | 1.48–1.65 | 2.48–2.63 | 4.9–5.2 | 15.7–16.7 | ≈ 49–58 | ≈ 84–97 |

- **Scaling.** The repeat spread is ≤ 22%. The cost grows roughly as NSTR³
  (local exponents 2.6–3.5 between NSTR 32 and 136).
- **NSTR 200 and 236:** estimated from the single ladder runs, scaled to the
  repeated NSTR 136 ratio; their single-run timings are noisier because they
  ran concurrently with other work.
- **Options at NSTR 60:** ACCUR 10⁻³, CORINT off and truncated moments cost the
  same as production (0.95–1.04).
- **Double precision:**
  - NSTR = 60: ≤ 3.7–3.9×. This is an upper bound: it was timed with the
    MXCMU = 256 double build, of which only the single-precision overhead could
    be subtracted (a MXCMU = 100 double build corrupts its heap and could not
    be used).
  - At high NSTR the double/single ratio falls to about 1.5–1.7 (NSTR 136) and
    1.2 (NSTR 236), so double precision at NSTR 200 costs about 60–90×.
- **Production scale:** V24 production jobs take 22–42 h. A DISORT cost
  ≤ 3.8× would make them at most 3.8× longer (up to about 80–160 h) before any
  optimisation. NSTR 200 in double precision would take months per job.

## 15. Conclusions by option, and the recommendation

| Setting | Production | Finding | Recommendation |
|---|---|---|---|
| NSTR | 60 | Converged to ≤ 0.38σ (measurement) and ≤ 0.28σ (Jacobians) away from backscatter; 5–6.5σ at exact backscatter, needing NSTR ≳ 200 in double precision (60–90× cost) for < 1σ | Keep 60 |
| ACCUR | 10⁻⁸ | No effect on results or run time at any value 0–10⁻² | Keep; no saving available |
| CORINT (Nakajima–Tanaka) | on | Essential | Keep |
| DELTAM | on | Essential | Keep |
| Moments | full adaptive L (Mie), 1000 (Baum) | Truncation to NSTR + 1 fails by orders of magnitude | Keep (validates the V24 Legendre change); review the Baum expansion at backscatter separately (§12.4) |
| Precision | single | ≤ 0.28σ systematic effect, but sporadic catastrophic failed calls (up to 856σ) in current products, removed by double precision | **Candidate change for a future version** (below) |
| Surface (black, Lambertian), plane-parallel, IBCND, ONLYFL | as in §4 | Design choices, not accuracy levers | Unchanged |

**Recommendation.** For a future production version, protect the LUTs against
failed single-precision DISORT calls. Leave every other DISORT setting
unchanged. There are two ways; the first can serve until the second is ready:

1. **Detect and recompute** (an interim guard, usable now):
   - after each LUT is written, flag cells with the §12.2 detector, plus any
     R_0v outside [−0.1, 1.5];
   - recompute only those direct-beam calls in double precision (or at NSTR 56
     and 68 for a consistency check).

   The cost is negligible (a few calls per product). It catches every
   whole-call failure and the extreme partial ones, but not partial failures
   with plausible values. It is a guard, not a cure.
2. **Double-precision DISORT at NSTR = 60** for all calls (the complete
   remedy). This removes the whole class of failures, including the partial
   ones, at ≤ 3.8× DISORT cost. It first needs:
   - a robust double-precision build: the MXCMU = 100 double build corrupted
     its heap here, and the MXCMU = 256 build returned NaN at one
     near-conservative point at NSTR 236 (§11.4);
   - a full-product comparison against the single-precision product away from
     the failed calls (expected ≤ 0.28σ, §10.6);
   - a production timing test.

**Evidence that should accompany either change:**
- the census of §12.2 before and after, on every V24 product;
- bitwise agreement elsewhere (option 1) or ≤ 0.3σ agreement (option 2);
- for option 1, the list of recomputed calls written into the product
  metadata.

**Not recommended:** raising NSTR. The exact-backscatter error is real (5–6.5σ)
but confined to a few degrees. If it matters for retrievals, the preferred
treatment is targeted: a retrieval-side screen or inflated uncertainty near
Θ = 180°, or a backscatter correction. A global stream increase would cost
60–90×.

**Separately:** the negative backscatter values from the 1000-moment Baum
expansion (§12.4) should be reviewed with the scattering inputs. They are not
a DISORT option.

## 16. Files, tests and final state

**Location.** All of this work is under `validation/`, which is git-ignored.
Nothing was committed.

**Code** (`validation/v24/disort_options/`):

| File | Purpose |
|---|---|
| `build_variants.py` | experimental DISORT builds (MXCMU, DELTAM/CORINT from COMMON, ACCUR, double precision) in `validation/tmp` |
| `run_variant.py`, `launch.sh` | run the production generator for one case and experiment; scattering computed once and reused |
| `check_baseline.py` | harness against the V24 product at matching nodes |
| `orac_fm.py`, `atmosphere.py` | ORAC forward-model mirror (§9.3) and clear-sky terms |
| `check_fm.py` | 8 forward-model checks |
| `propagate.py` | all experiments through the forward model, with Jacobians |
| `analyse_fm.py`, `cancellation.py` | ladder, breakdown, Jacobian and cancellation tables |
| `near_backscatter.py` | §10.4 |
| `timing.py` | §14 |
| `corrupted_calls.py`, `defect_impact.py` | §12 |
| `negative_radiance.py`, `severe_negatives.py`, `analyse.py` | operator-level diagnostics (§13) |
| `make_figures.py`, `nk_fm.py` | figures |
| `compact_disort_*.lut`, `defect_*.lut`, `ringing_modis_ghm24.lut` | Grid B test grids |

**Results** (`validation/v24/results/disort_options/`, 76 MB):
- CSV tables named in each section;
- `fm_arrays_*.npz`;
- `near_backscatter/`;
- `figures/` (convergence, options, runtime, near backscatter and four
  Nakajima–King plots).

**Scratch** (`validation/tmp/disort_options/`, about 650 MB, disposable):
- builds, run LUTs with `record.json`, scattering caches;
- full forward-model arrays and logs.

It is kept so that every number above can be regenerated without rerunning
DISORT.

**Experiments run:**
- 4 compact cases × 38 configurations;
- 3 × 4 × 13 timing repeats;
- 4 near-backscatter cases × 7 configurations;
- 5 reproduction cases: 4 failed calls and the GHM ringing (26 runs).

The operator-level runs of the first part were reused, not repeated.

**Tests:**
- `check_fm.py`: 8 of 8 pass;
- baseline reproduction: bitwise (agg compact case, three of the four failed
  calls; the fourth to 1.5×10⁻⁶);
- `kernel_check`: bitwise;
- repository suite (`pytest`): 291 passed, 14 skipped and 5 failed.
  - All five failures are missing local files: `docs/idl_python_port_status.md`
    (in no commit) and `validation/reference_cases/`, `validation/diagnostics/`
    (git-ignored; absent since before this work, the `validation/` directory
    having last changed on 2026-10-03).
  - None involves DISORT or LUT code.

**Production unchanged.**
- `git status` and `git diff --stat` show exactly the state at the start:
  ` M AGENTS.md`, ` M tests/test_baum.py` (+20 lines in total, pre-existing),
  `?? documents/`, `?? scripts/submit_v23_slstr_cloud_luts.sh`,
  `?? tests/test_nakajima_king.py` and the 40 `?? runs/*_v23.run`.
- HEAD is still `51717b8`.
- No file under `src/`, `create_orac_luts.py`, `create_orac_lut/`, `mie/`,
  `scripts/`, `runs/` or `tests/` has a modification time after the start of
  this work, including the compiled DISORT library.
- The ORAC reference clone stayed in session scratch outside the repository,
  unmodified.
