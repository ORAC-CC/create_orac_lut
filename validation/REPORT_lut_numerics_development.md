# ORAC LUT numerics: integration limits and adaptive Legendre expansion

Date: 2026-10-02. Working tree: `/home/g/grainger/project-oraclut` (`main`).

## 0. Summary and revisions

| Revision | Content | Role in comparisons |
|---|---|---|
| `9d663e9` (start) | Produced the V23 reference LUTs. Fixed 0.001–100 µm Mie integration; fixed NMom = 1000. | "V23" / "reference" |
| `c6ad545` | Liquid-water modified-gamma integration stops at the first legacy-lattice node ≥ 3.5 r_e | "adopted" (Phase 2), "fixed" (Phase 3) |
| `dbcc42c` | Adaptive Legendre expansion of Mie phase functions | "adaptive" (Phase 3) |
| `96bebbe` | Refined (nested) radius integration of liquid-water and ice-sphere Mie distributions (§14–§15) | "adopted" (radius grid) |

**V24** (§16 onwards) is the V23 product set and grids (Grid B) computed with
`c6ad545`, `dbcc42c` and `96bebbe`.

The LUT product/grid version (V23, Grid B) is not changed by either source
change. Both are documented in the `create_orac_luts.py` module docstring under
"Versions and numerical changes", which separates the three items. No V23 LUT
was regenerated, modified or deleted. All compact test products were written
under the ignored `validation/tmp/` and then deleted.

The starting working tree already had changes not made by this work: modified
`AGENTS.md` and `tests/test_baum.py`; untracked `documents/`,
`scripts/submit_v23_slstr_cloud_luts.sh`, `tests/test_nakajima_king.py` and 40
`runs/*_v23.run`. None of these was committed.

## 1. Directory restructuring (Phase 1)

The validation layout is now one directory per investigation:

| Before | After |
|---|---|
| `validation/compare_nakajima_king.py`, `validation/README_nakajima_king.md`, `validation/results/earthcare_msi_liquid_water_stg_v21_vs_v22.png` | `validation/nakajima_king/{compare_nakajima_king.py, README.md, results/…}` |
| `validation/REPORT_legendre_expansion_length.md` | `validation/legendre_expansion/REPORT_legendre_expansion_length.md` |
| `validation/REPORT_repository_development_workflow.md` | `validation/repository_workflow/…` (with a note on the moved paths it cites) |
| caches: `validation/__pycache__`, `size_distribution_limits/{__pycache__,.matplotlib-cache}` | removed; matplotlib cache defaults now `validation/tmp/matplotlib-cache` |

Paths, defaults and help text were updated. The `tests/test_nakajima_king.py`
import is now `validation.nakajima_king.compare_nakajima_king`.

`validation/` is wholly git-ignored (`.gitignore:91`) and the test file is
untracked, so `git mv` could not apply: the moves preserve no Git history,
because there was none to preserve. Narrowing that ignore rule remains the
separate policy decision recommended in the workflow review.

Checks after the move:
- `tests/test_nakajima_king.py` and `size_distribution_limits/test_run_study.py`: 9 passed;
- the moved scripts compile;
- `--help` of both command-line scripts runs.

## 2. Integration-limit methodology (Phase 2)

**Production quadrature.** `create_bwgp.mie_size_dist_new` (IDL
`mie_size_dist_new.pro`) integrates on a linear-radius trapezoid. Its details:

- npts = max(int(2π(r_u − r_l)/(λ·0.4)), 200) nodes, spacing h = (r_u − r_l)/(npts − 1).
- Every node position and h depend on the endpoint.
- The distribution weights and the Mie kernel are evaluated at every node.
- Nothing is cached across effective radii. For fixed limits the grid is
  identical for all r_e at a wavelength, so the identical Mie calculation is
  repeated for every radius.
- At λ > 7.85 µm the 200-node floor gives h = 0.5 µm.

**Why an endpoint change alters results.** The trapezoid samples the narrow Mie
resonance structure at Δx ≈ 0.4. Moving the endpoint moves every node relative
to that structure. Endpoint jitter of ±0.25–0.5% at 100 µm, where the tail is
negligible, reproduces the reference's whole error spread:

- extinction up to 6×10⁻³;
- g up to 3×10⁻³;
- side-scattering phase function up to 25%.

**Invariant comparison.** On a *lattice* variant, the reference nodes truncated
at (or extended to) the first node ≥ k·r_e, the difference from the reference
is exactly the tail. This is also the production rule that was adopted.

**Experiments**, all calling the production routines:
1. `quadrature_endpoint_diagnostic.py`: reference, naive candidate, jitter and
   lattice integrations against a converged truth (xres 0.05, 0.001–150 µm).
   - Main set: liquid stg, r_e 5–40 µm at 0.646, 0.868, 1.63, 2.11, 3.79 and
     11.03 µm.
   - Small-radius thermal set: r_e 1–5 µm at 3.79, 11.03 and 12.02 µm.
   - Coefficients use exact (refined) Gauss rules above each degree bound,
     comparing l = 0–127, 128–999 and ≥ 1000.
2. `compact_lut_comparison.py`: compact MODIS (channels 1, 2, 6, 7, 20, 31) and
   dual-view SLSTR (1, 3, 5, 6, 7, 8) LUTs through `create_orac_cloud_lut`.
   - Grid: τ 0.5–128, r_e 5–40 µm, 3×3×3 geometry.
   - Compares the optics, the 1000 production moments and every LUT variable.
   - Retrieval level (2C): R_0v interpolated at cell mid-points in (log τ, r_e)
     with `oraclut.validation.quadrature`, and its log τ and r_e Jacobians,
     normalised by the SLSTR `.inst` noise (`oldnefr` = 0.005). The MODIS file
     has no noise figures.
3. `benchmark_integration_limits.py`: single-thread CPU over the 24 Grid B
   liquid radii.

## 3. Reference against adaptive-limit results

**Error against the converged truth**, worst over wavelengths:

| r_e (µm) | Variant | Extinction (rel.) | ssa | g | χ_l, l 0–127 | χ_l, l 128–999 | χ_l, l ≥ 1000 |
|---|---|---|---|---|---|---|---|
| 5–30 | reference (V23) | ≤ 5.9e-3 | ≤ 1.5e-3 | ≤ 3.4e-3 | ≤ 3.4e-3 | ≤ 5.3e-5 | ≤ 4e-7 |
| 5–30 | adopted 3.5 r_e | same values | same | same | same | same | same |
| 5–30 | adopted − reference | ≤ 1.5e-6 | ≤ 3.3e-7 | ≤ 4.0e-7 | ≤ 1.6e-6 | ≤ 4.5e-7 | ≤ 2.6e-7 |
| 40 | reference (100 µm = 2.5 r_e) | 1.1e-3 | 3.6e-4 | 1.5e-4 | 3.0e-4 | 1.3e-4 | 8.1e-5 |
| 40 | adopted (≈140 µm) | 7.1e-4 | 3.0e-4 | 1.7e-4 | 2.1e-4 | 8.0e-6 | 8.0e-7 |

**Removed tail, thermal small-radius regime** (r_e 1–5 µm, 3.8–12 µm):

| Factor | Scattering (rel.) | g | χ_l, l 0–127 |
|---|---|---|---|
| 3.0 r_e | 5.8e-4 | 2.5e-4 | 2.5e-4 |
| 3.5 r_e | 2.9e-5 | 1.6e-5 | 1.6e-5 |

**Compact LUTs**, maximum over all variables (MODIS / SLSTR):

| Variant vs reference | max \|ΔR_0v\| | max rel. ΔR_0v | max \|Δg\| | SLSTR max \|ΔR\|/σ | SLSTR max \|ΔJ_logτ\|/σ |
|---|---|---|---|---|---|
| naive 3 r_e (production formula) | 2.6e-2 / 4.9e-2 | 22% / 33% | 3.7e-3 / 5.0e-3 | 5.0 | 3.3 |
| lattice 3 r_e | 1.4e-4 / 2.3e-4 | 0.8% / 1.8% | 1.5e-5 | 0.021 | 0.030 |
| **adopted, lattice 3.5 r_e** | **1.7e-4 / 2.3e-4** | 0.9% / 0.9% | 1.6e-5 | **0.023** | **0.034** |

The remaining adopted differences come from r_e = 40 µm, which is the removal of
the reference truncation error. At r_e 5–30 µm the optics change by ≤ 1.7×10⁻⁶.
The production implementation reproduces the run-time lattice variant
**bitwise** in every LUT variable and every scattering-cache array, for MODIS
and SLSTR.

## 4. Quadrature findings

- **Node placement dominates.** The legacy xres = 0.4 trapezoid is the dominant
  numerical error of the V23 Mie optics, against a converged truth: up to
  6×10⁻³ in extinction, 1.5×10⁻³ in ssa, 3×10⁻³ in g and 10–25% in the
  side/back-scattering phase function.
- **Small radii, thermal channels.** At thermal wavelengths the 200-node floor
  leaves r_e = 1 µm with about 0.9% extinction error (h = 0.5 µm across a
  distribution about 1 µm wide).
- **Naive endpoint changes.** Changing the endpoint naively re-randomises these
  errors, which is why the naive candidate moved LUT reflectances by up to 5σ.
- **Required future work, not done here.** A converged radius quadrature (finer
  xres, or nodes adapted to the distribution) is needed. It would change V23
  products by up to the amounts above and needs its own validation. The present
  change was deliberately constructed to be independent of that question.

## 5. Implementation decision

Implemented (`c6ad545`). The upper limit is the first node of the legacy lattice
at or beyond **3.5 × r_e**, extending past 100 µm when needed (r_e ≥ 29 µm). It
applies to **liquid-water modified-gamma components only**
(`generate_scattering_properties.radius_upper_factor`,
`create_bwgp.lattice_upper_radius`). The lower limit is unchanged at 0.001 µm.

**Excluded:**
- **Log-normal aerosols** are unchanged; they already use quantile limits.
- **Ice spheres** keep 0.001–100 µm, as specified. Their r_e up to 93 µm is
  then heavily truncated, as before.
- **Baum and T-matrix** are unchanged.
- **Effective variance above 0.1111111** is refused as unvalidated.

**Why 3.5 rather than 3.** The tail fraction beyond k·r_e depends only on v_e
and on the weighting, Q(α + 1 + n, k/v_e) with α = (1 − 3v_e)/v_e:

| Weighting | Beyond 3 r_e | Beyond 3.5 r_e |
|---|---|---|
| Area (n = 2, extinction by large particles) | 1.8e-5 | 6.6e-7 |
| r⁶ (scattering by small particles) | 1.0e-3 | 6.5e-5 |

**Why the lattice.** It removes the 5σ artefacts of a naive endpoint change.

## 6. Computational saving of the integration limit

Single-thread Mie step at 1000 angles, 24 Grid B liquid radii
(`results/integration_limit_benchmark.csv`):

| λ (µm) | 0.646 | 1.63 | 2.11 | 3.79 | 11.03 | Total |
|---|---|---|---|---|---|---|
| Reference CPU (s) | 99.1 | 15.1 | 9.2 | 3.1 | 0.6 | 127.2 |
| Adopted CPU (s) | 50.0 | 8.5 | 5.2 | 1.7 | 0.4 | 65.8 (**1.93×**) |

On the six-channel `generate_scattering_properties` benchmark (§9) the step fell
from 409 s (`9d663e9`) to 233 s (`c6ad545`).

## 7. Adaptive Legendre methodology (Phase 3)

The method follows Grainger (1990), §4.4–4.5 (`documents/Grainger 1990.pdf`).
It is implemented in `src/oraclut/idl_mirror/legendre_expansion.py` and
`generate_scattering_properties._adaptive_mie_wavelength`. For every
size-distribution-averaged Mie phase function (wavelength, r_e):

1. **Degree bound.** D_r = 2·NStop(x_max) + 2, from the largest size parameter
   of *its own* integration with the adopted limits (`mie_integration_limits`).
   Coefficients beyond D_r are exactly zero.
2. **Quadrature order.** Nq_r is above D_r plus a noise band, rounded up to a
   multiple of 64. It is per radius, so small radii are not sampled at a large
   radius's order. Nodes and weights are Newton-refined, because numpy/scipy
   weights err by about 2×10⁻⁶ at the forward node near N = 6000.
3. **Expansion length.** All ω_l for l < Nq_r are computed. L is King's
   criterion over the whole tail: |ω_l| < 10⁻⁹ for all l ≥ L, up to D_r.
   - Where the rounding noise measured just above D_r exceeds 10⁻⁹, the
     coefficients down to the noise level are kept, so L ≤ D_r + 1.
   - A shared, largest-radius bound let that noise inflate small-radius L
     (r_e = 10 µm gave 2819 instead of 699), so it was replaced by D_r.
4. **Acceptance.** The double-precision series must reproduce the directly
   calculated averaged phase function to ≤ 5×10⁻⁷ at 459 independent angles,
   and the noise must be ≤ 10⁻⁶. Otherwise Nq is raised by 50%, and the run
   fails after three attempts. A Mie-kernel angle limit (10 000) is checked
   explicitly. No 1000 cap remains.
5. **DISORT.** It receives max(L, NSTR + 1) moments per call, zero-padded.
   `setup_disort` and the Rayleigh `GETMOM` are called per phase function.
6. **Configuration.** `nmom` is obsolete for Mie classes and refused if set. It
   is still required for Baum/T-matrix classes, whose path is bitwise
   unchanged (compact Baum LUT: 37/37 variables identical between `c6ad545`
   and `dbcc42c`).

## 8. Legendre validation results

**Coefficient level** (`adaptive_legendre_validation.py`, 38 phase functions).
Each production result is compared with a run at 1.5–2× the order (the long
reference: every exact coefficient to D_r) and with the directly calculated
averaged phase function at ~1400 dense angles.

**Representative expansion lengths L:**

| Case | r_e (µm) | 0.412 µm | 0.555 µm | 0.646 µm | 2.11 µm | 3.79 µm | 11.03 µm |
|---|---|---|---|---|---|---|---|
| liquid stg | 5 | 550 | 411 | 355 | 117 | 67 | 26 |
| liquid stg | 10 | 1083 | 817 | 699 | 220 | 125 | 47 |
| liquid stg | 20 | 2225 | 1605 | 1439 | 429 | 244 | 88 |
| liquid stg | 40 | 4383 | 3273 | 2739 | 848 | 478 | 169 |
| ice spheres | 41 / 93 | 3151 / 3151 | – | 2031 / 2031 | – | – | – |

For aerosol a76, r_e 0.1 / 1 / 3.36 / 10 µm gives L = 30 / 477 / 817 / 2399 at
0.55 µm and 7 / 31 / 50 / 132 at 10.8 µm.

**Worst values over all 38 phase functions:**

| Check | Worst value |
|---|---|
| Reconstruction error (rel.), 0–1° / 1–30° / 30–180° | 6.1e-9 / 1.2e-8 / 1.3e-7 |
| Low orders, l ≤ 60, \|Δχ_l\| vs long reference | 9.6e-12 |
| High orders, 61 ≤ l < L, \|Δχ_l\| | 3.5e-11 |
| Largest reference χ_l dropped beyond L | 1.3e-11 |
| \|ω₁/3 − g_Mie\| | 2.5e-11 |
| Change of L when Nq is raised | ≤ 2.9% |

The two largest changes of L are noise-limited cases near D_r: liquid r_e
40 µm at 0.646 µm (2739 → 2819) and a76 r_e 10 µm at 0.638 µm (2133 → 2071).

**Radiative transfer, compact LUTs** (`compact_lut_legendre.py`). Each case is
compared with a long-expansion reference LUT (twice the order, every exact
coefficient):

| Case | Adaptive: max \|Δ\| (rel.) | Fixed 1000 (`c6ad545`): max \|Δ\| (rel.) | Max L adaptive / long |
|---|---|---|---|
| MODIS liquid | 2.4e-7 (1.0e-6) | 3.9e-2 (16%) | 2739 / 2819 |
| SLSTR liquid, dual view | 3.1e-6 (9.9e-6) | 5.7e-2 (22%) | 3277 / 3277 |
| SEVIRI aerosol a79 | 4.8e-8 (1.0e-6) | 3.6e-5 | 2095 / 2095 |
| SEVIRI aerosol a76 / a78 | 8.9e-8 (2.5e-7) | 5.3e-4 | 2133 / 2133 |

The figures are over R_0v, R_0d, R_dd, R_dv, T_*, E_md, g and ssa. The moments
passed to DISORT match the long reference to ≤ 1.5×10⁻⁸. The fixed 1000
misses χ_l up to 1.7×10⁻² beyond l = 1000 and aliases low orders up to 3×10⁻⁴.

## 9. Computational effect of adaptive moments

Single-thread `generate_scattering_properties`, liquid stg, 24 Grid B radii,
MODIS channels 3, 1, 6, 7, 20 and 31 (0.47–11.03 µm), in one call including the
0.55 µm reference:

| Revision | CPU (s) | Mean L (max) |
|---|---|---|
| `9d663e9` (V23) | 409 | 1000 (1000) |
| `c6ad545` | 233 | 1000 (1000) |
| `dbcc42c` | 641 | 682 (3789) |

The adaptive expansion costs 2.75× the fixed 1000 on this set, and 1.57× V23
overall.

- **Short wavelengths.** The cost is concentrated at 0.47–0.65 µm (mean L 1700
  at 0.469 µm; Nq up to about 4000), exactly where 1000 was aliased and wrong.
- **Infrared.** From 1.6 µm onwards it is cheaper (mean L 486 at 1.63 µm and
  77 at 11 µm).
- **0.55 µm reference.** It alone costs about 180 s, against 65 s with 1000
  points. Its moments are not passed to DISORT; computing it with the bulk
  quantities only would remove most of the excess (not done; see §13).
- **Radiative transfer.** DISORT cost changed little: compact-LUT RT time
  differed by under 10% between the adaptive and long variants.

The fixed-variant MODIS/SLSTR compact timings were taken without pinned threads
and are not comparable; the aerosol ones are. On those, adaptive took 32–60 s
against 12–23 s for the fixed 1000, because coarse modes need Nq ≈ 2100–2600.

## 10. Files changed

**Committed:**
- `c6ad545`: `src/oraclut/idl_mirror/create_bwgp.py`,
  `src/oraclut/idl_mirror/generate_scattering_properties.py`,
  `create_orac_luts.py`, `tests/test_idl_mirror.py`,
  `tests/test_size_integration_limits.py` (new).
- `dbcc42c`: `src/oraclut/idl_mirror/legendre_expansion.py` (new),
  `create_bwgp.py`, `generate_scattering_properties.py`,
  `create_orac_luts.py`, `src/oraclut/master_run.py`, `runs/template.run`,
  `tests/test_idl_mirror.py`, `tests/test_legendre_expansion.py` (new).

**Not committed** (`validation/` is ignored):
- `validation/size_distribution_limits/`: `quadrature_endpoint_diagnostic.py`,
  `compact_lut_comparison.py`, `benchmark_integration_limits.py`,
  `compact_liquid_grid.lut`, updated `README.md`, and `results/`.
- `validation/legendre_expansion/`: `README.md`,
  `adaptive_legendre_validation.py`, `compact_lut_legendre.py`,
  `compact_aerosol_grid.lut`, `benchmark_legendre.py`, and `results/`.
- `validation/nakajima_king/`, `validation/repository_workflow/` (moved), and
  this report.

`tests/test_nakajima_king.py` (untracked; import path updated) was also not
committed.

## 11. Tests and standalone validations run

**pytest:**
- Full suite with slow tests, at `dbcc42c` content: 276 passed, 14 skipped and
  5 failed.
- The 5 failures fail identically at `9d663e9`: missing `docs/` and
  `validation/diagnostics` files removed by the repository clean-up.
- New tests:
  - lattice nodes and > 100 µm extension;
  - tail-only liquid change;
  - liquid-only rule and variance refusal;
  - legacy arithmetic for other distributions;
  - Gauss–Legendre accuracy;
  - King whole-tail versus single-coefficient stopping;
  - degree bound;
  - L ≈ 220 (< 1000), ≈ 1039 (≈ 1000) and ≈ 2739 (> 1000);
  - the ice-sphere case needing > 3000 moments;
  - stability with increased Nq;
  - bitwise reproducibility;
  - refusal of `nmom` for Mie and requirement for Baum.
- The bitwise port check of the bulk optics against the validated path is kept,
  with the legacy limits imposed.

**Standalone, run directly:**
- all scripts in §2 and §8 (all values finite);
- production-routine consistency of the diagnostic (≤ 4×10⁻¹⁵);
- bitwise lattice/production and Baum equivalence checks;
- all changed modules compiled.

No new plots were produced; the tables above are the evidence.

## 12. Git commit and status

- Start: `9d663e9` on `main`.
- End: `dbcc42c` on `main`, with `c6ad545` between. Nothing was pushed.
- Final `git status --short` matches the starting state: ` M AGENTS.md`,
  ` M tests/test_baum.py`, `?? documents/`,
  `?? scripts/submit_v23_slstr_cloud_luts.sh`, `?? tests/test_nakajima_king.py`
  and 40 `?? runs/*_v23.run`.

## 13. Remaining scientific uncertainties and recommendations

1. **Radius quadrature** (§4). The xres = 0.4 trapezoid limits V23 optics
   accuracy (up to 0.6% extinction, about 20% side-scattering phase function),
   and the r_e = 1 µm thermal case is under-resolved. This is the largest
   remaining numerical error and needs a separate validated change.
   *Resolved for liquid water and ice spheres by `96bebbe` (§14).*
2. **Ice spheres** stay at 0.001–100 µm, truncating r_e > 33 µm distributions.
   Extending the liquid rule to them would raise the Mie and Legendre cost
   sharply (r_u up to about 325 µm).
3. **Baum and T-matrix** Legendre convergence is unsolved; they keep fixed
   nmom.
4. **Cost.** The adaptive expansion is 1.6× V23 on the benchmark set. Computing
   the 0.55 µm reference without moments (unused by DISORT) is the obvious
   saving. The redundant per-radius Mie recomputation on identical grids is
   another (the adopted lattice makes grids nested).
5. **Single-precision DISORT.** DISORT takes single-precision moments, and
   last-digit changes can move single-precision DISORT results (studied
   earlier). This is unchanged and was not part of the L criterion.
6. **Product provenance.** No product-version change was made. Because a LUT
   regenerated from a V23 run file with current code would differ numerically,
   recording the source revision in products, or in the manifests, is
   advisable. The V23 run files still set `nmom = 1000`, which current code
   refuses for Mie classes, so they cannot silently run the new algorithm.
7. **Sibling copy.** The earlier development copy
   `/home/g/grainger/project-oraclut-next` lies outside the repository. It was
   not modified or removed by this task; delete it if no longer wanted.

## 14. Radius-integration grid (V24 radius-grid change)

The full account is in `validation/radius_grid/REPORT_radius_grid.md`. Scripts,
CSV results and figures are in `validation/radius_grid/`.

**1. Legacy algorithm.**
- Linear-radius trapezoid, npts = max(int(2π(100 − 0.001)/(0.4 λ)), 200) on
  the 0.001–100 µm lattice: size-parameter step 0.4, or h = 0.5025 µm at
  λ > 7.85 µm.
- c6ad545 keeps the lattice and ends liquid water at 3.5 r_e.
- The Mie kernel runs at every node and production angle; nothing is reused
  between effective radii.
- Δx = 0.4 is about two nodes per period of the Mie ripple structure (≈ 0.8 in
  x), so the grid aliases it. Liquid r_e = 1 µm at 11 µm has only 8 nodes.

**2. Converged reference.**
- Nested dyadic refinement h₀/2^k, k = 0–6, from a single Mie pass on the
  finest nodes; level 0 reproduces production to ≤ 5×10⁻¹⁴.
- 54 cases: liquid stg r_e 1–40 µm at 0.41–12.0 µm, and ice spheres r_e
  10–93 µm.
- Level 5 against 6 agrees to 10⁻⁷–10⁻¹² in absorbing bands. In the visible it
  agrees to ≤ 3×10⁻⁵ (bulk) and ≤ 4×10⁻³ (back-scattering), where narrow
  resonances make convergence slow and non-monotonic.

**3. V23 grid error** (level 0 against level 6, worst case over all 54 cases):

| Quantity | Error |
|---|---|
| Extinction | 9.1×10⁻³ (r_e 1 µm, 12 µm); 5.9×10⁻³ for r_e ≥ 5 µm |
| ω₀ | 1.5×10⁻³ |
| g | 3.7×10⁻³ |
| Phase function, 5–150° | 5.8×10⁻² |
| Phase function, 150–180° | 2.4×10⁻¹ |
| Phase function, RMS | 3.2×10⁻² |
| χ_l | 3.7×10⁻³ |

The earlier estimates are confirmed for extinction (0.5–0.9%). The "~20% side
scattering" is confirmed in the back-scattering region (13–24%); between 5° and
150° the error is 0.6–6%.

**4. Candidates.**
- Uniform trapezoid, nested refinement — selected.
- Log-r trapezoid and one-panel Gauss–Legendre: 1–2 orders worse at equal node
  count.
- Richardson extrapolation: no gain on the resonance error.
- A-posteriori adaptive refinement: rejected, because its error estimate is
  unreliable when convergence is non-monotonic.

**5. Selected algorithm.**
- k is the smallest level with size-parameter step ≤ Δx_max and spacing
  ≤ r_e √v_e / 3.
- Δx_max = 0.025 for liquid water (k = 4) and 0.05 for ice spheres (k = 3):
  `generate_scattering_properties.refined_xres`,
  `create_bwgp.radius_refinement_level`.
- Limits and legacy nodes are unchanged (nested).
- The Mie sum is blocked (bitwise identical) to bound memory.
- Aerosol, Baum and T-matrix classes are unchanged.

**6. Convergence criterion.** Against level 6:
- R_0v ≤ 1.5×10⁻³ at every compact LUT node, glory included (0.3 σ of the 0.005
  SLSTR noise-equivalent reflectance);
- extinction (rel.) and g ≤ 2×10⁻⁴;
- ω₀ ≤ 10⁻⁵.

The coarsest level meeting all three is chosen per class.

**7–9. Optical properties, phase function and Legendre expansion.** Worst case
over the matrix, against level 6:

| Class | Adopted level | Extinction | g | Phase 150–180° | Phase RMS | χ_l |
|---|---|---|---|---|---|---|
| Liquid water | 4 | 5.6e-5 | 4.7e-5 | 1.3e-2 | 1.0e-3 | 4.7e-5 |
| Ice spheres | 3 | 4.1e-5 | 3.6e-5 | 5.1e-3 | 7.4e-4 | 3.6e-5 |

- The adaptive Legendre acceptance test passed at the first quadrature order
  for every compact phase function.
- L changed by ≤ 4.

**10. Compact LUTs** (production path; MODIS liquid, dual-view SLSTR liquid,
MODIS ice spheres; levels 0–6 plus the production code). Maximum |ΔR_0v|
against level 6:

| Case | V23 grid: glory / elsewhere | Adopted: glory / elsewhere | Reference uncertainty |
|---|---|---|---|
| MODIS liquid | 4.0e-2 / 1.8e-2 | 1.4e-3 / 1.1e-3 | 5.2e-4 |
| SLSTR liquid | 8.1e-2 / 2.3e-2 | 5.7e-4 / 2.9e-4 | 4.2e-4 |
| MODIS ice spheres | 1.8e-2 / 2.8e-3 | 8.8e-4 / 3.8e-4 | 4.0e-4 |

- The production-code products are bitwise identical to the forced levels 4
  and 3.
- Adopted bulk optics are within ≤ 1.6×10⁻⁴ (extinction), ≤ 1.3×10⁻⁴ (g) and
  ≤ 9×10⁻⁶ (ω₀).

**11. Retrieval relevance** (SLSTR, normalised by the 0.005 noise-equivalent
reflectance; reflectance and Jacobians of the bilinear interpolant at cell
mid-points):

| Grid | ΔR/σ | Δ(∂R/∂log τ)/σ | Δ(∂R/∂r_e)/σ per µm |
|---|---|---|---|
| V23 | 15.2 | 10.4 | 14.5 |
| Adopted | 0.09 | 0.07 | 0.06 |
| Level 5 (reference uncertainty) | 0.06 | 0.11 | 0.04 |

**12. Cost.** atmlxint5, one thread, MODIS channels 3, 1, 6, 7, 20, 31 plus
0.55 µm, Grid B radii:

| Case | V23 | dbcc42c | Adopted | Adopted / V23 |
|---|---|---|---|---|
| Liquid | 441 s | 680 s | 9 784 s | 22× |
| Ice spheres | 528 s | 1 316 s | 10 285 s | 19.5× |

- The Mie kernel accounts for ≥ 92% of each time.
- The cost model (within 2% of the measurements) predicts 12.6 h (MODIS liquid)
  and 13.1 h (MODIS ice spheres) of Mie per full product, against V23 job times
  of 22–28 h.
- Expected V24 Mie-class job times are therefore about 35–42 h (MODIS) and
  7–9 h (SLSTR), within the 72 h limit.

**Figures** (`validation/radius_grid/results/figures/`):
- `convergence_by_level.png`
- `compact_lut_r0v_by_level.png`
- `phase_error_re5_0.646um.png`

**Tests:**
- New `tests/test_radius_grid.py`:
  - blocked Mie sum bitwise equal to a single call;
  - refinement levels;
  - nesting and unchanged limits;
  - log-normal unaffected;
  - convergence towards a finer grid;
  - class selection.
- Updated port check: IDL grid imposed; production change bounded by the
  measured V23 error.
- Updated Legendre checks: refined direct phase function.
- Full suite: 291 passed, 14 skipped; the same 5 pre-existing failures (missing
  `docs/` and `validation/diagnostics` files).

## 15. Radius-grid production commit

`96bebbe` ("Refine the radius integration of liquid-water and ice-sphere Mie
distributions") contains:
- `src/oraclut/idl_mirror/create_bwgp.py`
- `src/oraclut/idl_mirror/generate_scattering_properties.py`
- `create_orac_luts.py` (docstring item 4; scattering-cache marker)
- `tests/test_radius_grid.py` (new)
- `tests/test_idl_mirror.py`
- `tests/test_legendre_expansion.py`

`c6ad545` and `dbcc42c` are preserved. The validation material is not
committed (`validation/` is ignored by the repository policy), as for the
earlier commits.

## 16. V24 definition

V24 is the V23 cloud product set, recomputed with the numerical changes
`c6ad545` (3.5 r_e liquid limit), `dbcc42c` (adaptive Legendre) and `96bebbe`
(refined radius grid).

**Unchanged:**
- platforms, channels, microphysical models, atmosphere, gas and Rayleigh
  settings;
- the Grid B sampling: τ, r_e and geometry.

**Coverage:** the complete V23 cloud coverage, determined from the 40 V23 run
files and the V23 products in the ORAC_LUTS directory:
- Aqua and Terra MODIS (channels 1–36);
- Sentinel-3A and -3B SLSTR (channels 1–9, dual view);
- each with liquid water 240, 253, 263, 273, old and stg (liquid-water Grid B)
  and ice sph, agg, ghm and src (ice Grid B).

That is 40 products: 24 liquid-water and 16 ice.

**Naming and documentation:**
- `runs/<platform>_<instrument>_cloud_<model>_v24.run` (40 files). They are
  identical to the V23 run files in every setting except `version = 24`;
  `nmom` is removed for the Mie classes, where it is obsolete, and kept at
  1000 for Baum.
- `runs/V24_PRODUCTION.md` documents the definition.
- The `create_orac_luts.py` docstring states it.
- Products are named `<platform>_<instrument>_m_<substance>_a01_p<model>_v24.nc`
  in the same permanent directory as the `_v23.nc` products. No V23 file is
  renamed or overwritten; the generator refuses to overwrite.

**Provenance:**
- `scripts/submit_v24_cloud_luts.sh` refuses to submit unless the production
  paths are clean and HEAD is contained in the fetched `origin/main`.
- It appends time, host, revision, job ID, run file and product to
  `validation/v24/production_submissions.tsv`.
- `scripts/run_oraclut_lut.slurm` now writes the Git revision and the
  clean/modified state of the source into every job log.
- The NetCDF format itself is unchanged.

## 17. V24 production-preparation commit

`51717b8` ("Prepare V24 cloud LUT production") contains:
- the 40 `runs/*_v24.run` and `runs/V24_PRODUCTION.md`;
- `scripts/submit_v24_cloud_luts.sh`;
- the job-log provenance lines in `scripts/run_oraclut_lut.slurm`;
- the V24 sentence of the `create_orac_luts.py` docstring.

It was kept separate from the scientific change `96bebbe`. The production
source (`src/`, `mie/`, `create_orac_lut/`) is identical in the two commits.

## 18. Smoke tests

`validation/v24/smoke_v24.py`.

**Before the commit**, each representative V24 run file was executed through
`create_orac_luts.run` with three channels, the compact grid and a scratch
output directory, all other settings unchanged:
- MODIS liquid stg, MODIS ice sph, MODIS Baum agg;
- dual-view SLSTR-A liquid 263, SLSTR-B ice sph.

**Results** (all PASS):
- One `_v24.nc` product per case, written only to the scratch directory.
- All floating-point variables finite.
- Variable and dimension names identical to the corresponding V23 products,
  dual view included.
- Liquid water: refined step 0.025, with 13 617 = 851×16 + 1 nodes for r_e
  10 µm at 0.646 µm, ending at 35.036 µm (the lattice node ≥ 3.5 r_e).
- Ice spheres: refined step 0.05, with 19 433 = 2429×8 + 1 nodes on
  0.001–100 µm.
- Baum agg: bitwise identical to the same run from `dbcc42c`.
- Grid B coordinates expanded by `load_lutstr` (τ, r_e, SZA, VZA, RAA) equal
  those of the V23 products (`terra_modis_…_pstg_v23.nc`,
  `sentinel-3a_slstr_…_psph_v23.nc`).

**Wall times** (single thread, atmlxint5):
- Baum: 14 s;
- Mie classes: 15–29 min.

**Peak resident memory:** ≤ 392 MB (Mie), 119 MB (Baum).

The adaptive Legendre expansion was exercised in every Mie run; it raises an
error if its acceptance test fails. In the production-code compact runs of
§14, the scattering caches carry the marker `nested-refinement-1` with
L = 10–3883 (liquid) and 28–2031 (ice spheres).

**After the push:**
- `git diff 96bebbe 51717b8 -- src mie create_orac_lut` is empty, so the
  smoke-tested numerics are those pushed.
- The submit script's dry run confirmed clean production paths, HEAD
  `51717b8` contained in `origin/main`, 40 run files consistent, no V24
  product present and all 40 V23 products present.

## 19. GitHub push

- Before the push: `origin/main` = `9d663e9`, local `main` 4 commits ahead,
  remote `git@github-orac:ORAC-CC/create_orac_lut.git`.
- `git push origin main` (no force): `9d663e9..51717b8  main -> main`.
- `git ls-remote origin refs/heads/main` returned
  `51717b82dece213d0f93ff1476ff1806983bbbe4`.
- Both `96bebbe` and `51717b8` are ancestors of the fetched `origin/main`.

**Production source commit:**
`51717b82dece213d0f93ff1476ff1806983bbbe4` on `origin/main`.

## 20. V24 production submission

- Submitted 2026-10-03 01:59:50–02:00:14 UTC from atmlxint5 with
  `scripts/submit_v24_cloud_luts.sh all`.
- Each job went through `scripts/submit_oraclut_lut.sh`, forwarded to
  atmlxint7, with `--time=72:00:00 --mem=24G` and partitions from
  `SBATCH_PARTITION` (shared, priority-eodg).
- Revision `51717b8` for every job.
- Before submission, the size and mtime of the 40 V23 products were recorded
  in `validation/v24/v23_products_before_v24_submission.txt`.
- Manifest: `validation/v24/production_submissions.tsv`. Submission output:
  `validation/v24/submission_log.txt`.

**Output location:** `/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/`.

**Expected products:**
- 24 liquid-water: 4 platforms × 6 models (job IDs 485440–485445,
  485450–485455, 485460–485465, 485470–485475);
- 16 ice: 4 platforms × 4 models (485446–485449, 485456–485459,
  485466–485469, 485476–485479).

**Job logs:** `validation/slurm/v24_*_<jobid>_*.out`.

| Job ID | Job name | Run file | Product |
|---|---|---|---|
| 485440 | v24_aqua_water_old | `aqua_modis_cloud_liquid-water_old_v24.run` | `aqua_modis_m_liquid-water_a01_pold_v24.nc` |
| 485441 | v24_aqua_water_stg | `aqua_modis_cloud_liquid-water_stg_v24.run` | `aqua_modis_m_liquid-water_a01_pstg_v24.nc` |
| 485442 | v24_aqua_water_240 | `aqua_modis_cloud_liquid-water_240_v24.run` | `aqua_modis_m_liquid-water_a01_p240_v24.nc` |
| 485443 | v24_aqua_water_253 | `aqua_modis_cloud_liquid-water_253_v24.run` | `aqua_modis_m_liquid-water_a01_p253_v24.nc` |
| 485444 | v24_aqua_water_263 | `aqua_modis_cloud_liquid-water_263_v24.run` | `aqua_modis_m_liquid-water_a01_p263_v24.nc` |
| 485445 | v24_aqua_water_273 | `aqua_modis_cloud_liquid-water_273_v24.run` | `aqua_modis_m_liquid-water_a01_p273_v24.nc` |
| 485446 | v24_aqua_ice_sph | `aqua_modis_cloud_water-ice_sph_v24.run` | `aqua_modis_m_water-ice_a01_psph_v24.nc` |
| 485447 | v24_aqua_ice_agg | `aqua_modis_cloud_water-ice_agg_v24.run` | `aqua_modis_m_water-ice_a01_pagg_v24.nc` |
| 485448 | v24_aqua_ice_ghm | `aqua_modis_cloud_water-ice_ghm_v24.run` | `aqua_modis_m_water-ice_a01_pghm_v24.nc` |
| 485449 | v24_aqua_ice_src | `aqua_modis_cloud_water-ice_src_v24.run` | `aqua_modis_m_water-ice_a01_psrc_v24.nc` |
| 485450 | v24_terra_water_old | `terra_modis_cloud_liquid-water_old_v24.run` | `terra_modis_m_liquid-water_a01_pold_v24.nc` |
| 485451 | v24_terra_water_stg | `terra_modis_cloud_liquid-water_stg_v24.run` | `terra_modis_m_liquid-water_a01_pstg_v24.nc` |
| 485452 | v24_terra_water_240 | `terra_modis_cloud_liquid-water_240_v24.run` | `terra_modis_m_liquid-water_a01_p240_v24.nc` |
| 485453 | v24_terra_water_253 | `terra_modis_cloud_liquid-water_253_v24.run` | `terra_modis_m_liquid-water_a01_p253_v24.nc` |
| 485454 | v24_terra_water_263 | `terra_modis_cloud_liquid-water_263_v24.run` | `terra_modis_m_liquid-water_a01_p263_v24.nc` |
| 485455 | v24_terra_water_273 | `terra_modis_cloud_liquid-water_273_v24.run` | `terra_modis_m_liquid-water_a01_p273_v24.nc` |
| 485456 | v24_terra_ice_sph | `terra_modis_cloud_water-ice_sph_v24.run` | `terra_modis_m_water-ice_a01_psph_v24.nc` |
| 485457 | v24_terra_ice_agg | `terra_modis_cloud_water-ice_agg_v24.run` | `terra_modis_m_water-ice_a01_pagg_v24.nc` |
| 485458 | v24_terra_ice_ghm | `terra_modis_cloud_water-ice_ghm_v24.run` | `terra_modis_m_water-ice_a01_pghm_v24.nc` |
| 485459 | v24_terra_ice_src | `terra_modis_cloud_water-ice_src_v24.run` | `terra_modis_m_water-ice_a01_psrc_v24.nc` |
| 485460 | v24_sentinel-3a_water_old | `sentinel-3a_slstr_cloud_liquid-water_old_v24.run` | `sentinel-3a_slstr_m_liquid-water_a01_pold_v24.nc` |
| 485461 | v24_sentinel-3a_water_stg | `sentinel-3a_slstr_cloud_liquid-water_stg_v24.run` | `sentinel-3a_slstr_m_liquid-water_a01_pstg_v24.nc` |
| 485462 | v24_sentinel-3a_water_240 | `sentinel-3a_slstr_cloud_liquid-water_240_v24.run` | `sentinel-3a_slstr_m_liquid-water_a01_p240_v24.nc` |
| 485463 | v24_sentinel-3a_water_253 | `sentinel-3a_slstr_cloud_liquid-water_253_v24.run` | `sentinel-3a_slstr_m_liquid-water_a01_p253_v24.nc` |
| 485464 | v24_sentinel-3a_water_263 | `sentinel-3a_slstr_cloud_liquid-water_263_v24.run` | `sentinel-3a_slstr_m_liquid-water_a01_p263_v24.nc` |
| 485465 | v24_sentinel-3a_water_273 | `sentinel-3a_slstr_cloud_liquid-water_273_v24.run` | `sentinel-3a_slstr_m_liquid-water_a01_p273_v24.nc` |
| 485466 | v24_sentinel-3a_ice_sph | `sentinel-3a_slstr_cloud_water-ice_sph_v24.run` | `sentinel-3a_slstr_m_water-ice_a01_psph_v24.nc` |
| 485467 | v24_sentinel-3a_ice_agg | `sentinel-3a_slstr_cloud_water-ice_agg_v24.run` | `sentinel-3a_slstr_m_water-ice_a01_pagg_v24.nc` |
| 485468 | v24_sentinel-3a_ice_ghm | `sentinel-3a_slstr_cloud_water-ice_ghm_v24.run` | `sentinel-3a_slstr_m_water-ice_a01_pghm_v24.nc` |
| 485469 | v24_sentinel-3a_ice_src | `sentinel-3a_slstr_cloud_water-ice_src_v24.run` | `sentinel-3a_slstr_m_water-ice_a01_psrc_v24.nc` |
| 485470 | v24_sentinel-3b_water_old | `sentinel-3b_slstr_cloud_liquid-water_old_v24.run` | `sentinel-3b_slstr_m_liquid-water_a01_pold_v24.nc` |
| 485471 | v24_sentinel-3b_water_stg | `sentinel-3b_slstr_cloud_liquid-water_stg_v24.run` | `sentinel-3b_slstr_m_liquid-water_a01_pstg_v24.nc` |
| 485472 | v24_sentinel-3b_water_240 | `sentinel-3b_slstr_cloud_liquid-water_240_v24.run` | `sentinel-3b_slstr_m_liquid-water_a01_p240_v24.nc` |
| 485473 | v24_sentinel-3b_water_253 | `sentinel-3b_slstr_cloud_liquid-water_253_v24.run` | `sentinel-3b_slstr_m_liquid-water_a01_p253_v24.nc` |
| 485474 | v24_sentinel-3b_water_263 | `sentinel-3b_slstr_cloud_liquid-water_263_v24.run` | `sentinel-3b_slstr_m_liquid-water_a01_p263_v24.nc` |
| 485475 | v24_sentinel-3b_water_273 | `sentinel-3b_slstr_cloud_liquid-water_273_v24.run` | `sentinel-3b_slstr_m_liquid-water_a01_p273_v24.nc` |
| 485476 | v24_sentinel-3b_ice_sph | `sentinel-3b_slstr_cloud_water-ice_sph_v24.run` | `sentinel-3b_slstr_m_water-ice_a01_psph_v24.nc` |
| 485477 | v24_sentinel-3b_ice_agg | `sentinel-3b_slstr_cloud_water-ice_agg_v24.run` | `sentinel-3b_slstr_m_water-ice_a01_pagg_v24.nc` |
| 485478 | v24_sentinel-3b_ice_ghm | `sentinel-3b_slstr_cloud_water-ice_ghm_v24.run` | `sentinel-3b_slstr_m_water-ice_a01_pghm_v24.nc` |
| 485479 | v24_sentinel-3b_ice_src | `sentinel-3b_slstr_cloud_water-ice_src_v24.run` | `sentinel-3b_slstr_m_water-ice_a01_psrc_v24.nc` |

**Status at submission:** all 40 jobs PENDING (Priority). Expected run times,
from §14 and V23: about 35–42 h (MODIS Mie classes), 18–24 h (MODIS Baum),
7–9 h (SLSTR Mie classes) and 5–6 h (SLSTR Baum).

## 21. Completed V24 jobs

**None had completed when this report was written; all 40 were queued.**

When they finish, check the following with `validation/v24/smoke_v24.py`
(adapted to the product files) or equivalent:
- 40 non-empty `_v24.nc` files;
- dimensions and Grid B coordinates equal to V23;
- finite fields;
- each job log showing `Git revision: 51717b8…` and `Source state: clean`;
- the V23 files unchanged against
  `validation/v24/v23_products_before_v24_submission.txt`.

The expected V24–V23 differences are the deliberate ones:
- up to several 10⁻² in R_0v near the glory;
- 10⁻³–10⁻² elsewhere in the visible and NIR;
- up to ~1% in the thermal extinction for r_e = 1–2 µm;
- the 3.5 r_e tail and adaptive-Legendre changes of `c6ad545` and `dbcc42c`.

Baum products should be identical to V23 apart from the version in the file
name.

## 22. Remaining scientific issues (radius grid and V24)

1. **Visible back-scattering residual.**
   - The adopted grid leaves up to ~1% (max) in the 150–180° phase function at
     r_e ≈ 5–10 µm (0.4–0.9 µm), and 1.4×10⁻³ in R_0v at the exact glory.
   - That is at the tolerance and about 3× the reference uncertainty.
   - Halving it costs a factor of two in Mie time, because the error comes from
     unresolved narrow resonances.
2. **Log-normal aerosol components** keep the legacy Δx = 0.4 grid. The same
   aliasing mechanism applies; a separate, validated change would be needed
   before the next aerosol production.
3. **Narrower size distributions** (e.g. liquid ap02–ap10) rely on the
   untested distribution-resolution term; validate before use.
4. **Ice spheres** keep the 0.001–100 µm limits, which truncate
   r_e > 33 µm distributions. This is unchanged from V23 by design (physical
   model not changed).
5. **Mie cost** dominates the V24 Mie-class jobs (about 13 h per MODIS product).
   Reuse of Mie evaluations across effective radii on the shared lattice is
   possible but needs a common angle set.
6. **Baum / T-matrix** Legendre convergence is still unsolved (fixed `nmom`).
7. **Product metadata.** The NetCDF files do not carry the Git revision; it is
   recorded in the job logs and the submission manifest.
