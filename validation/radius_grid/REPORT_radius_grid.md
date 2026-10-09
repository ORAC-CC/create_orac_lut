# Radius-integration grid of the Mie size distributions (V24)

Written for the repository owner and for reviewers of the V24 ORAC cloud LUTs.
The summary of record is `validation/REPORT_lut_numerics_development.md` §14–§22;
this report holds the detail.

Commits: `96bebbe` (radius grid) and `51717b8` (V24 preparation), both pushed
to `origin/main`.

Scratch products (convergence `.npz`, compact LUTs, revision exports and
smoke tests) were written under `validation/tmp/` and deleted once the
summaries in `results/` had been written. Every script here regenerates them.

## 1. Phase 1: the production radius grid (trace)

Path for a modified-gamma (liquid water, ice spheres) Mie class at source revision
`dbcc42c`:

1. `create_orac_luts.create_orac_cloud_lut` calls `generate_scattering_properties`.
2. `_adaptive_mie_scattering_properties` → `_adaptive_mie_wavelength`, once per
   wavelength (0.55 µm reference and each channel) and effective radius.
3. That calls `create_bwgp`, which calls `mie_integration_limits`. The limits are
   [0.001 µm, the lattice node ≥ 3.5 r_e] for liquid (c6ad545) and
   [0.001, 100] µm for ice spheres.
4. `mie_size_dist_new` then integrates with `quadrature("T", npts)`: a
   linear-radius trapezoid with
   npts = max(int(2π(100 − 0.001)/(λ·0.4)), 200) on the legacy interval, so the
   spacing is h₀ = (100 − 0.001)/(npts − 1).
   - For λ < 7.85 µm the size-parameter step is Δx = 2πh₀/λ ≈ 0.40.
   - Beyond 7.85 µm the 200-node floor gives h₀ = 0.5025 µm (Δx = 0.29 at 11 µm).
5. The modified-gamma number density times the trapezoid weight is applied, then
   the preserved Fortran `mieint` at every node and angle, then sums for
   extinction, scattering, g, F11 and the average volume.

Notes on the grid:
- The grid is uniform in r. There is no adaptivity, no error estimate and no
  reuse between effective radii.
- Typical V23 node counts:
  - r_e = 10 µm: 852 nodes at 0.646 µm, 261 at 2.11 µm and 71 at 11 µm.
  - Ice spheres: 2430 / 742 / 200 at all r_e.
  - Liquid r_e = 1 µm at 11 µm: **8** nodes (h₀/r_e = 0.5).
- Δx = 0.4 is about two nodes per period of the Mie ripple structure of water or
  ice (≈ 0.8 in x for m ≈ 1.33). The grid therefore aliases that structure,
  which is why moving the endpoint by a fraction of h₀ changed the V23-era
  results (`lattice_upper_radius` docstring).

## 2. Phase 2: converged reference

`radius_convergence.py` evaluates the Mie kernel once on the finest of a family of
nested dyadic refinements of the production grid, h_k = h₀/2^k with k = 0 … 6,
on the same endpoints. Every coarser level is a subset, so a single Mie pass gives
the trapezoid rule at every level and the Richardson value (4T_k − T_{k−1})/3.

- Level 0 equals production to ≤ 5×10⁻¹⁴ (relative, phase) and ≤ 10⁻¹⁵ (bulk) in
  all 54 cases (the `check` field).
- Angles: the production Gauss–Legendre nodes for the case plus ~1400 dense
  angles (0–1° every 0.005°, 1–30° every 0.05°, 30–180° every 0.25°).
- Legendre χ_l = ω_l/(2l + 1) is computed per level.

Matrix (54 cases):
- liquid stg, r_e 5, 10, 20, 30, 40 µm at 0.646, 0.868, 1.63, 2.11, 3.79 and
  11.03 µm (0.412 µm for r_e 5, 10, 20);
- r_e 1, 2, 3 µm at 3.79, 11.03 and 12.02 µm;
- ice spheres, r_e 10, 30, 60, 93 µm at 0.646, 2.11 and 11.03 µm.

Results: `results/convergence_matrix.csv` (every level and both rules) and
`results/convergence_matrix.txt`.

**Convergence of the reference, level 5 against 6:**

| Regime | Extinction, g, χ | Side 5–150° | Back 150–180° |
|---|---|---|---|
| Absorbing: 3.79, 11–12 µm; 2.11 µm r_e ≥ 30 | 10⁻⁷–10⁻¹² | ≤ 10⁻⁵ | ≤ 10⁻⁵ |
| 1.63–2.11 µm | ≤ 1.3×10⁻⁵ | ≤ 3×10⁻⁴ | ≤ 5×10⁻⁴ |
| 0.41–0.87 µm (k ≈ 10⁻⁹–10⁻⁷) | ≤ 3×10⁻⁵ | ≤ 1.5×10⁻³ | ≤ 4×10⁻³ |

- In the visible, convergence is close to first order and non-monotonic. Narrow,
  essentially unresolvable Mie resonances are sampled pseudo-randomly; the
  trapezoid error then behaves like a sampling error.
- The level-6 reference is therefore uncertain by ~10⁻³ in the visible
  back-scattering. That is two orders below the V23 error it is used to measure.
- Richardson extrapolation brings no gain there (it assumes smooth h² error) and
  is not used.

**Alternative quadratures at equal node count** (`alternative_grids.py`, r_e 10 µm,
2.11 µm; `results/alternative_grids_re10_2.11um.csv`):
- The uniform-r trapezoid is the most accurate, by 1–2 orders of magnitude
  against a log-r trapezoid and a single-panel Gauss–Legendre rule.
- On an oscillatory integrand sampled on a uniform grid the trapezoid is the
  natural rule. Gauss–Legendre does not converge until it resolves every ripple.

## 3. Phase 3: error of the V23 grid

Level 0 (the V23 / dbcc42c grid) against level 6; maxima over the matrix:

| Quantity | V23 grid error | Worst case |
|---|---|---|
| Extinction (rel.) | **9.1×10⁻³** | liquid r_e 1 µm, 12.02 µm (8 nodes) |
| Extinction (rel.), r_e ≥ 5 µm | 5.9×10⁻³ | r_e 10, 3.79 µm |
| Single-scattering albedo (abs.) | 1.5×10⁻³ | r_e 10, 3.79 µm |
| g (abs.) | 3.7×10⁻³ | r_e 5, 0.868 µm |
| Phase function, forward 0–5° (max rel.) | 1.1×10⁻² | |
| Phase function, side 5–150° (max rel.) | **5.8×10⁻²** | r_e 5, 0.646 µm |
| Phase function, back 150–180° (max rel.) | **2.4×10⁻¹** | r_e 5, 0.868 µm |
| Phase function, RMS over 0–180° (rel.) | 3.2×10⁻² | |
| χ_l, l ≤ 60 (abs.) | 3.7×10⁻³ | |
| χ_l, 61–999 (abs.) | 3.6×10⁻⁴ | |
| χ_l, ≥ 1000 (abs.) | 1.5×10⁻⁶ | |

**Earlier estimates.** The estimates came from the c6ad545 endpoint study and
were "~0.5–0.6% extinction and ~20% side scattering".
- **Extinction: confirmed.** The error is 0.48–0.59% at r_e 5–10 µm in the NIR
  and 3.8 µm, and 0.9% for r_e = 1–2 µm in the thermal infrared, where only 8–15
  nodes cover the distribution.
- **Phase function: confirmed only in the back-scattering region (150–180°),**
  where errors reach 13–24% in the visible and NIR at every r_e. In the side
  region (5–150°) the error is 0.6–6%. The earlier "side" figure included the
  back-scattering angles.

**Where the V23 grid is adequate:**
- Thermal infrared for r_e ≥ 3 µm, and all ice spheres at 11 µm: bulk
  quantities and χ_l are within 7×10⁻⁶ and the phase function within 6×10⁻⁵,
  because absorption damps the ripples.
- The defects are in the weakly absorbing visible–SWIR, at 3.79 µm (aliasing of
  the damped ripple), and for small r_e in the thermal infrared (distribution
  resolution).

Levels needed (`results/convergence_matrix.csv`) — worst case over all 54 cases:

| Level k (Δx at λ < 7.85 µm) | Ext. | g | Side | Back | RMS |
|---|---|---|---|---|---|
| 0 (0.4) | 9.1e-3 | 3.7e-3 | 5.8e-2 | 2.4e-1 | 3.2e-2 |
| 1 (0.2) | 9.9e-4 | 8.5e-4 | 4.3e-2 | 9.7e-2 | 1.9e-2 |
| 2 (0.1) | 8.1e-4 | 7.0e-4 | 2.5e-2 | 6.5e-2 | 9.2e-3 |
| 3 (0.05) | 2.7e-4 | 2.1e-4 | 9.1e-3 | 2.0e-2 | 3.2e-3 |
| 4 (0.025) | 5.6e-5 | 4.7e-5 | 4.8e-3 | 1.3e-2 | 1.0e-3 |
| 5 (0.0125) | 2.9e-5 | 2.5e-5 | 1.5e-3 | 4.0e-3 | 3.6e-4 |

The worst cases at k ≥ 2 are r_e = 5 µm at 0.646–0.868 µm. That distribution
spans the fewest resonances, so their sampling error averages least.

<!-- Phases 4-6 appended below -->

## 4. Phase 4: candidate grids and selection

**Candidates considered:**

1. **Uniform-r trapezoid, nested dyadic refinement of the legacy lattice**
   (h₀/2^k).
   - Every legacy node is kept and the limits are unchanged. This preserves the
     c6ad545 lattice construction, so the 3.5 r_e upper limit and the
     ice-sphere 0.001–100 µm limits are untouched.
   - It is extension-invariant: adding a tail adds nodes without moving any.
   - Its error at level k can be measured against level k + 1 from the same
     nodes.
2. **Log-r trapezoid** at equal node count: 1–2 orders worse (§2).
3. **Single-panel Gauss–Legendre in r** at equal node count: 1–2 orders worse,
   and not nested (§2).
4. **Richardson / Simpson combination of two nested levels**: no gain on the
   resonance-dominated visible error (§2). Not used.
5. **A-posteriori adaptive refinement** (refine until |T_k − T_{k−1}| < tol).
   Rejected, because the visible error is non-monotonic (§2): the
   successive-difference estimator can be accidentally small. A fixed,
   a-priori rule validated against the reference is more predictable and costs
   the same.

**Selected rule** (`create_bwgp.radius_refinement_level`,
`generate_scattering_properties.refined_xres`):

- k is the smallest level for which
  - the size-parameter step is 2π (h₀/2^k)/λ ≤ Δx_max, and
  - the spacing is h₀/2^k ≤ r_e √v_e / 3 (at least three nodes per
    area-weighted standard deviation of the distribution).
- Δx_max = **0.025 for liquid water** and **0.05 for ice spheres**.
  - At λ < 7.85 µm this is k = 4 (liquid) or k = 3 (ice); in the thermal
    infrared (h₀ = 0.5025 µm) the same targets give k = 4 / 3.
  - The bounds carry a 2% margin because the legacy step exceeds its nominal
    0.4 by up to 1%.
- Log-normal (aerosol) components, Baum and T-matrix classes are unchanged.

**The convergence criterion (tolerance) used for the selection:**

| Quantity | Tolerance (vs level 6) |
|---|---|
| R_0v at every compact LUT node, glory included | ≤ 1.5×10⁻³ absolute (0.3 σ of the 0.005 SLSTR SAD noise-equivalent reflectance) |
| Extinction (rel.), g (abs.) | ≤ 2×10⁻⁴ |
| ω₀ (abs.) | ≤ 1×10⁻⁵ |

Each class gets the coarsest level that meets the tolerance. Ice spheres
(m ≈ 1.31) meet it one level earlier than liquid water (m ≈ 1.33). This was
measured both in the matrix (ice r_e 10 µm, 0.646 µm, level 3 back-scattering
error 3.0×10⁻³ against 1.6×10⁻² for liquid) and in the compact LUTs below.
Weaker morphology-dependent resonances at the lower refractive index are the
plausible cause.

**Distribution-resolution term.** It never binds for the V23/V24 grids
(r_e ≥ 1 µm, v_e = 0.111) once Δx ≤ 0.05. It guards narrower or smaller
distributions, which have not been validated here. The level-0 thermal
infrared error for r_e = 1–2 µm (0.9%) is removed by the Δx term.

**Interaction with the earlier changes:**
- c6ad545: the refined grid has the same nodes and end points, so the
  3.5 r_e lattice end is unchanged.
- dbcc42c: the degree bound D_r depends only on the upper limit, so Nq_r is
  unchanged. L changes by ≤ 4 (compact LUTs, against level 6), because the
  better-resolved phase function has slightly different noise-limited tails.
- **Memory.** The Mie kernel is now called on blocks of ≤ 4×10⁶ node × angle
  elements, accumulated in the single-call order. The result is bitwise that
  of one call (`tests/test_radius_grid.py`), so V23-path calculations are
  unaffected while refined grids fit the 24 GB job memory.

## 5. Phase 5: validation

**5.1 Microphysics and phase function** (matrix, against level 6; §3 table and
`results/figures/convergence_by_level.png`):

| Class | Adopted level | Ext. (rel.) | g (abs.) | ω₀ (abs.) | Side 5–150° | Back 150–180° | RMS | χ_l |
|---|---|---|---|---|---|---|---|---|
| Liquid (k = 4) | worst | 5.6e-5 | 4.7e-5 | 2.8e-6 | 4.8e-3 | 1.3e-2 | 1.0e-3 | 4.7e-5 |
| Ice spheres (k = 3, 12 cases) | worst | 4.1e-5 | 3.6e-5 | 7.5e-8 | 6.2e-3 | 5.1e-3 | 7.4e-4 | 3.6e-5 |
| V23 grid (k = 0) | worst | 9.1e-3 | 3.7e-3 | 1.5e-3 | 5.8e-2 | 2.4e-1 | 3.2e-2 | 3.7e-3 |

- The remaining phase-function error is concentrated in the visible
  back-scattering and is within 1–4× the uncertainty of the reference itself
  (level 5 vs 6, §2).
- `results/figures/phase_error_re5_0.646um.png` shows the angular structure for
  the worst case (r_e 5 µm, 0.646 µm).

**5.2 Legendre expansion.** The adaptive expansion (dbcc42c) is applied to the
refined phase functions unchanged. Its acceptance test (series against the
directly calculated phase function at 459 angles, ≤ 5×10⁻⁷) passed at the
first quadrature order for every phase function of every compact run: none
of the 24 logs contains a retry message.
Expansion length changes are ≤ 4 against level 6.

**5.3 Compact LUTs through the production path** (`compact_lut_radius.py`; summaries
in `results/compact_lut/`).

Cases:
- MODIS liquid (channels 3, 1, 2, 6, 7, 20, 31);
- dual-view SLSTR liquid (channels 1, 3, 5, 6, 7, 8);
- MODIS ice spheres (channels 1, 2, 6, 7, 20, 31).

Grid: τ 0.5–128 (5 values), SZA 0/30/60°, VZA 0/30/55°, RAA 0/90/180°;
r_e 1, 2, 3, 5, 10, 20, 30, 40 µm (liquid) and 5, 10, 30, 60, 93 µm (ice).

Levels 0–6 were forced at run time, and "adopted" is the production code. The
adopted products are bitwise identical to the forced level 4 (liquid) and
level 3 (ice) in every variable.

Maximum |ΔR_0v| against level 6 (`k6_r0v_by_geometry.csv`,
`results/figures/compact_lut_r0v_by_level.png`):

| Case | V23 grid Θ ≥ 170° | V23 grid Θ < 170° | Adopted Θ ≥ 170° | Adopted Θ < 170° | Level 5 vs 6 (reference uncertainty) |
|---|---|---|---|---|---|
| MODIS liquid | 4.0e-2 (17%) | 1.8e-2 | 1.4e-3 | 1.1e-3 | 5.2e-4 |
| SLSTR liquid | 8.1e-2 (34%; 16 σ) | 2.3e-2 (4.5 σ) | 5.7e-4 (0.11 σ) | 2.9e-4 (0.06 σ) | 4.2e-4 |
| MODIS ice spheres | 1.8e-2 (12%) | 2.8e-3 | 8.8e-4 | 3.8e-4 | 4.0e-4 |

- The largest V23 error is at SZA = VZA = 0 (exact backscatter, the glory):
  SLSTR 0.555 µm, r_e 2 µm, τ = 8.
- Off the glory the V23 error is still 2–4.5 σ for liquid water.

**Other LUT variables** (max abs. against level 6, adopted / V23 grid):

| Case | R_0d | R_dd | T_dd | E_md | Extinction ratio |
|---|---|---|---|---|---|
| MODIS liquid | 6.2e-4 / 1.6e-2 | 1.9e-4 / 1.3e-2 | 2.2e-4 / 1.3e-2 | 1.2e-5 / 7.9e-3 | 1.5e-4 / 2.9e-2 |
| SLSTR liquid | 3.8e-4 / 1.9e-2 | 2.1e-4 / 1.5e-2 | 1.1e-4 / 1.5e-2 | 1.2e-5 / 7.9e-3 | 5.2e-5 / 2.9e-2 |
| MODIS ice spheres | 2.6e-4 / 1.8e-3 | 1.9e-4 / 1.5e-3 | 2.1e-4 / 3.9e-3 | 1.1e-5 / 3.6e-3 | 1.0e-4 / 5.9e-3 |

**Bulk optics in the LUT** (scattering cache, adopted against level 6, worst
channel and radius):

| Case | Extinction (rel.) | g | ω₀ |
|---|---|---|---|
| MODIS liquid | 1.6e-4 | 1.3e-4 | 4.4e-6 |
| SLSTR liquid | 3.3e-5 | 3.0e-5 | 9.1e-6 |
| MODIS ice spheres | 1.3e-4 | 1.1e-4 | 1.0e-6 |

With the V23 grid the corresponding errors are up to 2.0×10⁻², 1.4×10⁻² and
1.1×10⁻³. Thermal channels improve from up to 8.6×10⁻³ (extinction, r_e 1 µm)
to ≤ 10⁻⁷.

**5.4 Retrieval relevance** (`k6_retrieval_level.csv`):
- R_0v is interpolated bilinearly to the (log τ, r_e) cell mid-points, and the
  Jacobians ∂R/∂log₁₀τ and ∂R/∂r_e of the interpolant are formed. Both are
  normalised by the SLSTR SAD noise-equivalent reflectance (0.005). The MODIS
  instrument file carries none, so MODIS is given in reflectance units.

| Case | Variant | max \|ΔR\| | ΔR/σ | Δ(∂R/∂log τ)/σ | Δ(∂R/∂r_e)/σ (per µm) |
|---|---|---|---|---|---|
| SLSTR liquid | V23 grid | 7.6e-2 | 15.2 | 10.4 | 14.5 |
| SLSTR liquid | adopted | 4.7e-4 | 0.09 | 0.07 | 0.06 |
| SLSTR liquid | level 5 (ref. unc.) | 2.9e-4 | 0.06 | 0.11 | 0.04 |
| MODIS liquid | V23 grid / adopted | 3.0e-2 / 1.2e-3 | – | – | – |
| MODIS ice spheres | V23 grid / adopted | 1.1e-2 / 6.9e-4 | – | – | – |

- With the adopted grid, the radius-quadrature error in reflectances and
  Jacobians is about a tenth of the measurement noise, at the level of the
  reference's own uncertainty.
- With the V23 grid it reached 10–15 σ at the glory and several σ elsewhere,
  in values and in the r_e Jacobian.

## 6. Phase 6: cost

**Benchmark** (`benchmark_radius_grid.py`, `results/benchmark/`):
- Host atmlxint5, Intel Xeon Gold 6240R @ 2.40 GHz.
- One thread (OMP/OPENBLAS/MKL_NUM_THREADS = 1), one process per variant, run
  with `nice`, other jobs on the host.
- Workload: `generate_scattering_properties` for Terra MODIS channels 3, 1, 6,
  7, 20, 31 plus the 0.55 µm reference, on the V23 Grid B radii.
- The Mie kernel was wrapped to count nodes and node × angle evaluations and to
  time it separately.

| Case | Variant | Wall (s) | CPU (s) | Mie CPU (s) | Other CPU (s) | Radius nodes | Node × angle |
|---|---|---|---|---|---|---|---|
| liquid (24 r_e) | V23 (`9d663e9`) | 441 | 441 | 438 | 3 | 263 448 | 2.63e8 |
| liquid | `dbcc42c` | 680 | 680 | 624 | 56 | 161 504 | 3.93e8 |
| liquid | adopted (k = 4) | 9 784 | 9 783 | 9 688 | 95 | 2 581 544 | 6.29e9 |
| ice spheres (29 r_e) | V23 (`9d663e9`) | 529 | 528 | 525 | 4 | 318 333 | 3.18e8 |
| ice spheres | `dbcc42c` | 1 316 | 1 316 | 1 301 | 15 | 318 333 | 8.42e8 |
| ice spheres | adopted (k = 3) | 10 285 | 10 284 | 10 230 | 54 | 2 545 243 | 6.73e9 |

- **Liquid:** the adopted scattering step costs 22× V23 and 14× `dbcc42c`.
- **Ice spheres:** it costs 19.5× V23 and 7.8× `dbcc42c`.
- **Reference (level 6):** a further 4× (liquid) or 8× (ice) on the same
  workload. It was not run on this full workload. On the compact runs it cost
  12 842 s (MODIS liquid) against 3 288 s adopted and 249 s for level 0.

**Full products.**
- `cost_model.py` sums NStop(x) × angles over every node of every Mie
  evaluation. Calibrated on the `dbcc42c` benchmark, it predicts the adopted
  benchmark to within 2% and the node counts exactly.
- Predicted single-thread Mie CPU per product on the benchmark host:

| Product | dbcc42c | V24 (adopted) |
|---|---|---|
| MODIS liquid (36 channels) | 0.8 h | 12.6 h |
| MODIS ice spheres | 1.6 h | 13.1 h |
| SLSTR liquid (9 channels) | 0.16 h | 2.5 h |
| SLSTR ice spheres | 0.33 h | 2.7 h |

- The V23 jobs took 22–24 h (MODIS liquid), 25–28 h (MODIS ice spheres),
  18–24 h (MODIS Baum), 4.6–4.8 h (SLSTR liquid) and 5.4–5.6 h (SLSTR ice), the
  RT dominating.
- The V24 Mie-class jobs are therefore expected to take about 35–42 h (MODIS)
  and 7–9 h (SLSTR) on comparable nodes, inside the 72 h limit with a margin of
  about 1.7.
- Level 4 for ice spheres (26 h of Mie) would have reduced that margin to about
  1.3 for no measurable LUT benefit (§5.3).
- Baum classes are unchanged.

**Memory.** The blocked Mie summation bounds each per-call phase array to
32 MB, independent of refinement level. Production keeps `--mem=24G`; the
peak resident memory of the V24 smoke tests is recorded in
`validation/REPORT_lut_numerics_development.md` §18.

## 7. Remaining issues

- **Visible back-scattering.** It converges slowly under any uniform grid,
  because narrow resonances are sampled pseudo-randomly. The adopted grid leaves
  ~1% (max) in the 150–180° phase function at r_e ≈ 5–10 µm in the visible and
  NIR. That is 1.4×10⁻³ in R_0v at the glory, at the tolerance and about 3× the
  reference uncertainty. An analytic or resonance-aware treatment would be a
  separate study.
- **Log-normal aerosol components** keep the legacy Δx = 0.4 grid. The same
  aliasing mechanism applies to them, but they were not part of this task's
  validation, and no aerosol product is regenerated.
- **Distribution-resolution term** (three nodes per area-weighted standard
  deviation). It is inactive for v_e = 0.111 and has not been validated for
  narrower distributions (e.g. the liquid-water ap02–ap10 models).
- **Cost.** Mie evaluations are not shared between effective radii, although on
  the liquid lattice the nodes are common. Reuse would need a common angle set
  across radii (adaptive Nq_r differs) and is left as an optimisation.
