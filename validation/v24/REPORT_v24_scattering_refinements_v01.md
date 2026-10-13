# V24: three refinements of the Mie scattering properties on the V23 baseline

Date: 2026-10-13.  Release documentation: `V24_RELEASE.md`.  Scripts and
results: `validation/v24/scattering_refinements/` (`results/*.json`; the
NetCDF and npz work files are git-ignored).

**Status:** implemented and validated on compact calculations.  The full V24
products have **not** been generated; V24 is not validated as a product set
until they are produced and checked (section 12).

## 1. The three methods, recovered from Git

The historical implementation was read with `git show` only
(`c6ad545`, `dbcc42c`, `96bebbe`; production commit `51717b8`) together with
`validation/REPORT_lut_numerics_development.md` sections 2-9 and 14-16.  It was
written for the pre-V22 code (`9d663e9` lineage), before the IDL-compatible
arithmetic of V22, so it was re-implemented, not restored.

| Method | Historical commit | Definition recovered |
|---|---|---|
| Liquid-water integration limit | `c6ad545` | Liquid-water modified-gamma components only: the upper limit is the first node of the legacy lattice (the 0.001-100 um trapezoid with npts = max(int(2 pi (100 - 0.001) nu / 0.4), 200) nodes) at or beyond 3.5 r_e, beyond 100 um when needed; npts = intervals + 1; lower limit 0.001 um; v_e > 0.1111111 refused.  Applied at 0.55 um and every channel.  Why 3.5: r^6-weighted tail 6.5e-5, area-weighted 6.6e-7 at v_e = 1/9. |
| Adaptive Legendre expansion | `dbcc42c` | Classes whose components are all Mie (`opac`/`user`, not T-matrix).  Per wavelength and r_e: degree bound D = 2 NStop(x_max) with NStop = int(x + 4.05 x^(1/3) + 2) + 1 (mieint.f series length, conservative), x_max from that radius's own integration limits; Gauss-Legendre order Nq >= D + 1 + (32 + D // 32), rounded up to a multiple of 64; nodes from numpy leggauss Newton-refined three times, weights 2 / ((1 - x)(1 + x) P_N'^2); omega_l = (2l + 1)/2 sum(p P_l w) for all l < Nq in double precision; L = 1 + the highest l <= D with |omega_l| >= 1e-9 (King's criterion on the whole tail; coefficients beyond D are zero); noise = max |omega_l| for D < l <= D + band; acceptance: reconstruction error at 459 independent angles (0-1 deg every 0.01, 1-180 every 0.5) <= 5e-7 relative and noise <= 1e-6, otherwise Nq x 1.5 (at most three orders), then failure; Mie kernel limit 10000 angles enforced; chi_l = omega_l / (2l + 1).  Baum / T-matrix classes keep the fixed nmom expansion. |
| Nested radius refinement | `96bebbe` | Liquid-water and water-ice (sphere) modified-gamma Mie components: every interval of the legacy lattice (with the c6ad545 limit for liquid water) divided into 2^k parts, k the smallest level with spacing / 2^k <= 1.02 min(xres / (2 pi nu), r_e sqrt(v_e) / 3), xres = 0.025 (liquid) or 0.05 (ice spheres); npts = (n - 1) 2^k + 1; limits and legacy nodes kept; Mie kernel called on blocks of at most 4e6 nodes x angles with the sum accumulated in single-call order. |

No numerical defect was found in the historical methods (section 9 lists the
deliberate differences; section 5 one documented limitation).

## 2. The V23 baseline (protected)

* Git: `main` at `209e145`, equal to `origin/main`; uncommitted at the start:
  ` M create_orac_lut/makerunfile.pro` (owner), ` M validation/v23_terra_modis/README.md`
  and the untracked V22/V23 consistency check (previous task), `documents/`, `plots/`.
* V23 release `releases/lut_v23.release` (105 entries, unchanged), scientific
  code `lut-code-v22.4` with Grid B; pinned scientific files:
  `create_orac_luts_v04.py`, `create_bwgp_v03.py`,
  `generate_scattering_properties_v03.py`, `legacy_disort_v04.py`, ...;
  DISORT `DISORTtoIDL.so` sha256 `59307e42604ca0f7...`.
* V23 products (`staging/lut_v23/`, unchanged): pold `bc2450ca26ca3bf4...`,
  p240 `e6d4548d8b0a1fcb...`, pagg `d99fbfeb0792ea94...`.
* The working `create_bwgp.py` (= `_v04`) differs from the pinned `_v03` only in
  a docstring; `legacy_disort.py` (= `_v06`, not changed here) only in its
  default precision text.  V24 starts from these files.
* Routines that needed changes: `create_bwgp.py` (methods 1 and 3),
  `generate_scattering_properties.py` (methods 1-3: class rules, adaptive
  path), a new `legendre_expansion.py` (method 2), and the driver
  `create_orac_luts.py` (DISORT moment count per phase function).  Nothing in
  DISORT, the Mie kernel, Baum, the readers, the writer or the inputs.

## 3. Implementation

Only the four files below differ from V23 (`releases/lut_v24.release` against
`releases/lut_v23.release`: CODE_VERSION, driver, two modules, one new module).
Default arguments keep every routine's IDL behaviour; the new behaviour is
selected per class by `generate_scattering_properties`.

| File | SourceVault | Changes |
|---|---|---|
| `src/oraclut/idl_mirror/legendre_expansion.py` | new, `_v01` | Method 2, written for V24 (not copied): `gauss_legendre` (Newton-refined), `mie_nstop`, `degree_bound`, `noise_band`, `quadrature_order`, `legendre_coefficients`, `rounding_noise`, `king_length`, `series`, `expand`; constants of section 1. |
| `src/oraclut/idl_mirror/create_bwgp.py` | `_v04` -> `_v05` | `size_integration_limits` and `legacy_npts` (the IDL inline code as routines, identical arithmetic); `mie_size_dist_new(npts=None)` and blocked Mie sums (sequential IDL TOTAL carried between blocks); `lattice_upper_radius` (method 1), `radius_refinement_level` and `mie_integration_limits` (methods 1 and 3); `create_bwgp(radius_upper_factor=None, refined_xres=None, phase_dtype=float32)`. |
| `src/oraclut/idl_mirror/generate_scattering_properties.py` | `_v03` -> `_v04` | `LIQUID_UPPER_RADIUS_FACTOR`, `LIQUID_REFINED_XRES`, `ICE_SPHERE_REFINED_XRES`, `radius_upper_factor`, `refined_xres`, `uses_adaptive_legendre`; Mie classes branch to `_adaptive_mie_scattering_properties` / `_adaptive_mie_wavelength`; the first return value is `lmom[m, l, r]` (nmom everywhere for Baum / T-matrix).  The class combination keeps the V23 FLTARR storage and IDL TOTAL order. |
| `create_orac_luts.py` | `_v04` -> `_v05` | `_disort_moments` (NMOM = max(nmom, L), zero padding, refusal of any non-zero moment beyond those passed); `setup_disort` and Rayleigh `GETMOM` once per distinct NMOM; scattering-cache scheme marker; configuration print and provenance `numerics` describing the V24 methods (V23 text kept for Baum / T-matrix classes). |

Also: `CODE_VERSION` (`_v05` holds V23's; V24 content: `lut-code-v24.1`,
`lut_version` 24, no compatible versions), `releases/lut_v24.release`
(`scripts/make_lut_v24_manifest.py`), the three `runs/*_v24.run`,
`scripts/{run,submit,monitor}_v24_release.*`, `releases/lut_v23.dependencies`
(the V23 record its manifest refers to, written from the V23 products'
provenance by the new `scripts/make_release_dependencies.py`), tests (section 11) and
documentation (`V24_RELEASE.md`, `LUT_CODE_VERSIONS.md` section 11,
`releases/README.md`, `AGENTS.md`).

## 4. DISORT receives max(1000, L) moments (`disort_nmom.py`)

Source facts (`create_orac_lut/disort2/src/DISORTfunctions.f`, the compiled
file): `PMOM(0:MAXMOM, *)` is dimensioned by `MAXMOM`, which the wrapper and
`legacy_disort_v04` set equal to NMOM; the only fixed limit, `MXSQT = 1000`, is
checked against 2 max(100, NSTR) = 200 (lines 536-538, 5417-5418); CHEKIN
requires NMOM >= NSTR and every PMOM in [-1, 1] (5233-5244).  DELTAM and
CORINT are always true (622-623): moments beyond NSTR enter only the
intensity corrections, `INTCOR` (sum K = 1 ... NMOM, 2474-2484) and `SECSCA`
(K = NSTR ... NMOM, 2684).

Tests with the original Intel kernel, NSTR 60, Grid B view geometry, single
cloud layer, tau 1 / 8 / 64, three source cases each:

| Case (um, r_e) | L | max(1000,L) = L bitwise | L+512 padded bitwise | max(L,61) bitwise (L<1000) | fluxes bitwise under truncation | max dR at N = 1000 / 1500 / 2000 / 2500 | s/call at 1000 / production |
|---|---|---|---|---|---|---|---|
| 0.469 21 | 2057 | True | True | n/a | True | 1.0e-01 / 5.0e-04 / 0.0e+00 / - | 0.087 / 0.100 |
| 0.646 21 | 1509 | True | True | n/a | True | 1.4e-03 / 0.0e+00 / - / - | 0.086 / 0.092 |
| 0.646 39 | 2751 | True | True | n/a | True | 1.4e+00 / 4.7e-02 / 8.3e-04 / 3.6e-06 | 0.089 / 0.108 |
| 0.858 11 | 581 | True | True | True | True | - / - / - / - | 0.090 / 0.095 |
| 3.79 21 | 255 | True | True | True | True | - / - / - / - | 0.100 / 0.097 |

* **Accepted and used:** no DISORT error up to NMOM = 2751; truncating the
  exact expansion changes R_0v and the change vanishes monotonically as
  N -> L (exactly 0 at L).
* **Fluxes** are bitwise independent of every moment beyond NSTR (all
  truncations), as the source requires.
* **Padding** is bitwise neutral; the production rule equals NMOM = L bitwise,
  and equals the historical max(L, 61) bitwise when L < 1000.
* **Beyond King's L** (IR/NIR cases, L < D): adding the coefficients below
  1e-9 (|chi| <= 1e-12) changes R by exactly 0.
* **Cost:** +7% to +21% per DISORT call at L = 1509-2751 (mean over the three
  optical depths); none when L <= 1000.
* The truncation rows are a diagnostic of the moments alone, not of V23: V23's
  1000 moments were aliased (a 1000-point rule), which makes its series
  interpolate the phase function at the nodes; section 7 shows what V23 lost
  (the glory).

## 5. Convergence of the radius integration and the Legendre expansion (`convergence.py`)

| Case (model, um, r_e) | nodes V24 / reference / V23 | V24 vs 4x finer: ext, g, phase 5-150, 150-180 | V23 grid vs 4x finer: ext, g, phase 150-180 | 3.5 vs 5 r_e: ext | L (2 Nq) | max dchi at 2 Nq | reconstruction | noise |
|---|---|---|---|---|---|---|---|---|
| liquid 240 0.469 5 | 9377 / 37505 / 3349 | 1.4e-05, 1.2e-05, 2.2e-03, 6.9e-03 | 1.4e-03, 1.1e-03, 1.2e-01 | 5.9e-07 | 487 (487) | 3e-13 | 3e-09 | 2e-10 |
| liquid 240 0.646 21 | 28593 / 114369 / 2431 | 3.0e-06, 2.7e-06, 1.4e-03, 1.3e-03 | 1.9e-05, 1.7e-05, 1.4e-01 | 6.2e-07 | 1509 (1449) | 3e-12 | 4e-09 | 5e-09 |
| liquid 240 0.646 39 | 53073 / 212289 / 2431 | 6.5e-07, 5.4e-07, 1.2e-03, 6.2e-04 | 2.7e-04, 6.1e-06, 5.3e-02 | 6.5e-07 | 2751 (2751) | 6e-12 | 9e-09 | 1e-08 |
| liquid 240 0.646 5 | 6817 / 27265 / 2431 | 1.7e-05, 1.4e-05, 2.4e-03, 2.9e-03 | 2.4e-03, 1.8e-03, 1.9e-01 | 5.8e-07 | 359 (359) | 2e-13 | 9e-10 | 9e-11 |
| liquid 240 0.858 11 | 11281 / 45121 / 1830 | 1.2e-05, 9.7e-06, 1.7e-03, 2.3e-03 | 7.6e-04, 6.3e-04, 1.4e-01 | 5.5e-07 | 581 (581) | 2e-13 | 3e-09 | 3e-10 |
| liquid 240 11.03 1 | 113 / 449 / 200 | 0.0e+00, 2.2e-08, 7.6e-08, 4.5e-08 | 8.9e-03, 2.7e-04, 8.1e-04 | 2.1e-06 | 10 (10) | 2e-16 | 4e-11 | 1e-14 |
| liquid 240 11.03 21 | 2353 / 9409 / 200 | 0.0e+00, 0.0e+00, 5.9e-12, 2.7e-12 | 5.7e-07, 0.0e+00, 1.0e-07 | 5.7e-07 | 92 (92) | 2e-14 | 6e-08 | 1e-12 |
| liquid 240 2.13 21 | 8657 / 34625 / 737 | 1.2e-07, 0.0e+00, 8.7e-05, 2.5e-05 | 4.4e-05, 1.0e-04, 7.9e-02 | 6.0e-07 | 447 (447) | 3e-13 | 1e-09 | 1e-10 |
| liquid 240 3.79 1 | 241 / 961 / 414 | 0.0e+00, 0.0e+00, 4.7e-09, 3.8e-10 | 1.7e-04, 2.9e-05, 2.2e-04 | 1.4e-06 | 17 (17) | 4e-16 | 3e-09 | 6e-14 |
| liquid 240 3.79 21 | 4865 / 19457 / 414 | 0.0e+00, 0.0e+00, 1.1e-07, 5.8e-08 | 2.8e-04, 2.2e-04, 6.2e-02 | 5.8e-07 | 255 (255) | 3e-14 | 1e-09 | 2e-11 |
| ice sph 0.646 41 | 19441 / 77761 / 2431 | 1.8e-06, 1.7e-06, 3.2e-03, 3.0e-03 | 6.9e-06, 6.3e-06, 9.6e-02 | - | 2031 (2031) | 6e-12 | 2e-08 | 2e-08 |
| ice sph 0.646 9 | 19441 / 77761 / 2431 | 3.2e-05, 2.7e-05, 4.2e-03, 2.9e-03 | 1.4e-04, 1.1e-04, 6.6e-02 | - | 2031 (839) | 5e-13 | 7e-10 | 2e-09 |
| ice sph 11.03 93 | 1593 / 6369 / 200 | 7.9e-08, 0.0e+00, 8.8e-07, 6.7e-10 | 6.4e-06, 0.0e+00, 4.8e-08 | - | 134 (134) | 1e-13 | 2e-08 | 2e-11 |
| ice sph 2.13 41 | 5889 / 23553 / 737 | 6.4e-08, 0.0e+00, 9.3e-05, 3.6e-05 | 3.8e-07, 6.0e-08, 5.8e-02 | - | 615 (615) | 5e-13 | 4e-09 | 7e-10 |

* **Radius integration:** against two further nested levels (4x the nodes),
  the production grid is within <= 3.2e-5 in extinction, <= 2.7e-5 in g,
  <= 6e-8 (liquid) in w, <= 4.2e-3 in the 5-150 deg and <= 6.9e-3 in the
  150-180 deg phase function.  The V23 grid against the same reference: up to
  8.9e-3 in extinction (r_e 1 um, 11 um), 1.8e-3 in g and 19% at back-scatter.
* **Integration limit:** 3.5 r_e against 5 r_e changes extinction by
  <= 2.1e-6 and the phase function by <= 1.5e-5 (liquid).
* **Legendre expansion:** doubling Nq changes the moments by <= 6e-12;
  reconstruction errors <= 5.6e-8 (tolerance 5e-7); noise <= 2.1e-8
  (limit 1e-6); |omega_0 - 1| <= 4e-12; omega_1 / 3 = g to the float32 rounding of g.
* **Noise-limited L** (documented property of the historical rule): when the
  rounding noise just above D exceeds 1e-9, L extends to D + 1 with pure-noise
  coefficients (|chi| ~ 1e-12).  For liquid water L changes by <= 4% when Nq is
  doubled (0.646 um, r_e 21: 1509 -> 1449); for ice spheres at r_e 9 um, whose
  integration still runs to 100 um, 2031 -> 839.  These coefficients change
  DISORT results by exactly 0 (section 4) and cost only DISORT time.

## 6. Scattering properties, V24 against V23 (`optics.py`, `compare_optics.py`)

Terra MODIS channels 1, 2, 3, 6, 7, 20, 31, 32, r_e 1, 5, 11, 21, 39 um
(ice spheres 5, 41, 93 um), each code from its own release (V23 sealed,
V24 working tree = `lut_v24.release`):

| Model | max abs rel. extinction | max abs rel. extinction ratio | max abs dw | max abs dg | L range | max abs dchi (l < 1000) | consistency |
|---|---|---|---|---|---|---|---|
| liquid-water_240 | 9.1e-03 | 1.1e-02 | 8.3e-04 | 3.7e-03 | 9-3789 | 3.4e-02 | w in [0,1] True, chi_0 = 1 True, max abs chi 1, abs(chi_1 - g) <= 0e+00 |
| liquid-water_old | 9.0e-03 | 9.8e-03 | 8.4e-04 | 2.7e-03 | 9-3789 | 3.4e-02 | w in [0,1] True, chi_0 = 1 True, max abs chi 1, abs(chi_1 - g) <= 0e+00 |
| water-ice_sph | 2.3e-03 | 2.1e-03 | 3.5e-04 | 7.7e-04 | 26-2793 | 1.9e-01 | w in [0,1] True, chi_0 = 1 True, max abs chi 1, abs(chi_1 - g) <= 0e+00 |

* Largest bulk changes where the V23 grid was least resolved: r_e 1 um in the
  thermal infrared (extinction +0.9%, ratio +1.1%), r_e 5-11 um at 3.79 um
  (+0.4 to +0.8%); at r_e 39 um every extinction rises by 2.8e-4 (the tail
  beyond the IDL's 100 um now integrated).
* The first 1000 moments change by up to 3.4e-2 (liquid) and 0.19 (ice
  spheres, r_e 93 um, 0.466 um): V23's 1000-point quadrature aliased them.
* Consistency: 0 <= w <= 1, chi_0 = 1, |chi_l| <= 1, chi_1 = g to float32.

## 7. Compact LUTs, V24 against V23 (`compact_lut.py`, `compare_compact.py`)

Identical compact Grid B subset (`compact_liquid_grid_b.lut`: tau 0.125, 1, 8,
64, 256; r_e 1, 5, 11, 21, 39 um; SZA and VZA 0, 30, 60, 75; RAA 0-180 by 30;
ice: tau 0.125-64, r_e 5, 41, 93) and settings (atmosphere 2, no gas,
Rayleigh, srf_quad 1, NSTR 60, nmom 1000), channels 1, 2, 3, 6, 7, 20, 31, 32:

**liquid-water_240**: variables True, dimensions True, coordinates True, V24 finite True, bitwise identical variables 38, isolated R_0v call outliers 2

| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (liquid 1, 5, 11, 21, 39; ice 5, 41, 93) |
|---|---|---|---|---|---|---|
| extinction_coefficient | 1.2e-01 | 3.8e+00 | 5.0e+00 | {'effective_radius': 11.0, 'channels': 20} | 9.1e-03 | 2e-02, 6e-01, 5e+00, 5e-01, 2e+00 |
| extinction_coefficient_ratio | 3.1e-04 | 6.7e-03 | 9.1e-03 | {'effective_radius': 11.0, 'channels': 20} | 1.1e-02 | 2e-03, 3e-03, 9e-03, 3e-04, 3e-04 |
| single_scatter_albedo | 7.2e-07 | 6.5e-04 | 8.3e-04 | {'effective_radius': 21.0, 'channels': 20} | 1.1e-02 | 4e-04, 3e-04, 3e-04, 8e-04, 2e-04 |
| asymmetry_parameter | 4.6e-05 | 3.4e-03 | 3.7e-03 | {'effective_radius': 11.0, 'channels': 20} | 4.6e-03 | 2e-03, 3e-03, 4e-03, 5e-04, 3e-04 |
| R_0v | 1.2e-04 | 1.6e-02 | 1.1e-01 | {'relative_azimuth': 0.0, 'satellite_zenith': 30.0, 'solar_zenith': 30.0, 'optical_depth': 256.0, 'effective_radius': 39.0, 'solar_channels': 3} | 4.1e+00 | 6e-03, 3e-02, 9e-02, 1e-02, 1e-01 |
| R_0d | 1.6e-04 | 3.9e-03 | 6.0e-03 | {'solar_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 5.0, 'solar_channels': 2} | 3.1e-02 | 2e-03, 6e-03, 5e-03, 8e-04, 3e-03 |
| T_0d | 7.8e-05 | 3.7e-03 | 6.2e-03 | {'solar_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'solar_channels': 20} | 7.5e-02 | 2e-03, 6e-03, 6e-03, 2e-03, 3e-03 |
| T_00 | 3.5e-09 | 9.9e-04 | 3.0e-03 | {'solar_zenith': 0.0, 'optical_depth': 1.0, 'effective_radius': 11.0, 'solar_channels': 20} | 3.4e-02 | 7e-04, 9e-04, 3e-03, 9e-05, 1e-04 |
| R_dv | 6.0e-05 | 3.7e-03 | 6.0e-03 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 5.0, 'channels': 2} | 3.2e-02 | 2e-03, 6e-03, 5e-03, 8e-04, 3e-03 |
| T_dv | 3.9e-05 | 3.7e-03 | 6.2e-03 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'channels': 20} | 7.5e-02 | 4e-03, 6e-03, 6e-03, 2e-03, 3e-03 |
| R_dd | 5.9e-05 | 3.4e-03 | 4.9e-03 | {'optical_depth': 8.0, 'effective_radius': 5.0, 'channels': 2} | 1.7e-02 | 2e-03, 5e-03, 4e-03, 5e-04, 3e-03 |
| T_dd | 3.8e-05 | 3.9e-03 | 4.9e-03 | {'optical_depth': 8.0, 'effective_radius': 5.0, 'channels': 2} | 7.4e-02 | 4e-03, 5e-03, 4e-03, 1e-03, 3e-03 |
| E_md | 2.1e-05 | 3.7e-03 | 3.9e-03 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 1.0, 'thermal_channels': 32} | 1.4e-02 | 4e-03, 3e-03, 4e-03, 2e-03, 3e-04 |

**liquid-water_old**: variables True, dimensions True, coordinates True, V24 finite True, bitwise identical variables 38, isolated R_0v call outliers 0

| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (liquid 1, 5, 11, 21, 39; ice 5, 41, 93) |
|---|---|---|---|---|---|---|
| extinction_coefficient | 4.6e-02 | 3.4e+00 | 4.0e+00 | {'effective_radius': 11.0, 'channels': 20} | 9.0e-03 | 1e-02, 3e-01, 4e+00, 2e+00, 2e+00 |
| extinction_coefficient_ratio | 3.9e-04 | 6.0e-03 | 7.1e-03 | {'effective_radius': 11.0, 'channels': 20} | 9.8e-03 | 1e-03, 4e-03, 7e-03, 8e-04, 1e-04 |
| single_scatter_albedo | 1.4e-06 | 7.3e-04 | 8.4e-04 | {'effective_radius': 11.0, 'channels': 20} | 1.1e-02 | 4e-04, 2e-04, 8e-04, 3e-04, 2e-05 |
| asymmetry_parameter | 4.4e-05 | 2.4e-03 | 2.7e-03 | {'effective_radius': 11.0, 'channels': 6} | 3.2e-03 | 6e-04, 2e-03, 3e-03, 5e-04, 9e-05 |
| R_0v | 1.1e-04 | 1.2e-02 | 7.3e-02 | {'relative_azimuth': 0.0, 'satellite_zenith': 30.0, 'solar_zenith': 30.0, 'optical_depth': 256.0, 'effective_radius': 39.0, 'solar_channels': 3} | 3.6e-01 | 1e-02, 2e-02, 3e-02, 4e-03, 7e-02 |
| R_0d | 1.1e-04 | 7.1e-03 | 9.2e-03 | {'solar_zenith': 0.0, 'optical_depth': 256.0, 'effective_radius': 11.0, 'solar_channels': 6} | 2.7e-02 | 9e-04, 7e-03, 9e-03, 2e-03, 3e-03 |
| T_0d | 7.3e-05 | 3.7e-03 | 1.1e-02 | {'solar_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'solar_channels': 6} | 2.0e-01 | 9e-04, 4e-03, 1e-02, 2e-03, 3e-03 |
| T_00 | 6.4e-09 | 1.4e-03 | 2.3e-03 | {'solar_zenith': 0.0, 'optical_depth': 1.0, 'effective_radius': 11.0, 'solar_channels': 20} | 2.8e-02 | 4e-04, 2e-03, 2e-03, 3e-04, 5e-05 |
| R_dv | 3.6e-05 | 6.9e-03 | 9.2e-03 | {'satellite_zenith': 0.0, 'optical_depth': 256.0, 'effective_radius': 11.0, 'channels': 6} | 2.8e-02 | 9e-04, 7e-03, 9e-03, 2e-03, 3e-03 |
| T_dv | 4.5e-05 | 3.6e-03 | 1.1e-02 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'channels': 6} | 2.0e-01 | 4e-03, 4e-03, 1e-02, 2e-03, 3e-03 |
| R_dd | 4.8e-05 | 5.9e-03 | 8.5e-03 | {'optical_depth': 256.0, 'effective_radius': 11.0, 'channels': 6} | 1.4e-02 | 7e-04, 6e-03, 9e-03, 2e-03, 3e-03 |
| T_dd | 5.1e-05 | 3.2e-03 | 9.0e-03 | {'optical_depth': 8.0, 'effective_radius': 11.0, 'channels': 6} | 2.0e-01 | 3e-03, 3e-03, 9e-03, 2e-03, 3e-03 |
| E_md | 2.2e-05 | 3.5e-03 | 3.6e-03 | {'satellite_zenith': 30.0, 'optical_depth': 8.0, 'effective_radius': 1.0, 'thermal_channels': 32} | 1.4e-02 | 4e-03, 8e-04, 4e-03, 4e-04, 5e-05 |

**water-ice_agg**: variables True, dimensions True, coordinates True, V24 finite True, bitwise identical variables 52, isolated R_0v call outliers 0

| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (liquid 1, 5, 11, 21, 39; ice 5, 41, 93) |
|---|---|---|---|---|---|---|
| extinction_coefficient | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| extinction_coefficient_ratio | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| single_scatter_albedo | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| asymmetry_parameter | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_0v | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'relative_azimuth': 0.0, 'satellite_zenith': 0.0, 'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_0d | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_0d | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_dv | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'satellite_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_dv | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'satellite_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_dd | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_dd | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| E_md | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'satellite_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'thermal_channels': 20} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |

Attribution of the largest differences (`validation/v24/scattering_refinements/`):

* **pagg (Baum):** all 52 variables bitwise identical: Baum classes are
  calculated as in V23.
* **Fluxes, E_md, transmissions** (<= 1.1e-2): follow the optics changes of
  section 6 (e.g. pold T_dd at tau 64, 1.63 um, r_e 11 um changes by 20% of a small
  value because w rises by 5.6e-4, i.e. the co-albedo falls by 7%).
* **R_0v at exact back-scatter** (Theta = 180 deg: SZA = VZA, RAA 0), r_e 39 um,
  0.469 um: +0.07 to +0.11.  This is the glory.  The directly calculated Mie
  phase function at 180 deg is 1.092 (0.466 um, r_e 39 um); V24's 3789-term
  series gives 1.092, V23's 1000 aliased moments 0.484 (0.646 um: 0.930 /
  0.930 / 0.756; at 179 deg all agree within 0.02).  V24 is correct; V23
  under-represented the narrow glory peak.
* **R_0v, p240, 3.79 um, tau 1, r_e 11 um, SZA 30 deg:** the V23 values
  oscillate with azimuth and are negative (-0.010, -0.026), while V23's other
  solar zeniths agree with V24 to 2e-3 and V24 is smooth and positive; L is
  139 there, so V23's moments were complete and its optics differ by <= 0.8%.
  This is an isolated single-precision DISORT failure of the V23 call (the
  known, accepted V23 limitation), absent in V24: not a V24 change.
* **The two R_0v call "outliers" of p240** (1.63 um tau 0.125 r_e 21; 0.646 um
  tau 256 r_e 21; both SZA 30) are geometric: only the SZA = VZA = 30 deg row
  of the compact grid samples the glory (180 deg) and the rainbow (~138 deg,
  RAA 90), so the phase-function change shows there alone; V24 also removes a
  negative V23 value at Theta = 120 deg.  No isolated V24 DISORT failure was
  found in the compact products.

No unexplained discrepancy or discontinuity remains.

**Release launcher.**  `launcher_smoke.run` (the V24 p240 run file with the
compact grid and channels 1 and 31; SLURM 495970) ran through
`create_orac_luts_release_v04.py` with `releases/lut_v24.release`: 106 files
verified, 34 Python modules loaded from them (V23's 33 plus
`legendre_expansion_v01`), code release `lut-code-v24.1`, Intel single-precision
kernel, "DISORT NMOM = max(1000, L)"; its product is bitwise identical, in
every RT operator and optical property of channels 1 and 31, to the
working-tree compact product above, and its embedded provenance records LUT
release V24, `lut-code-v24.1`, Grid B and the module files
`create_bwgp_v05`, `generate_scattering_properties_v04`, `legendre_expansion_v01`.

## 8. Runtime (`timing.py`, `runtime_estimate.py`; one CPU, atmnode010)

Mie step for all 24 Grid B liquid-water radii at one wavelength, each
refinement added in turn (wall time; the 1000-point rule for the
non-adaptive variants):

| um | V23 | + 3.5 r_e | + 3.5 r_e + refinement | + 3.5 r_e + adaptive Legendre | V24 (all three) | V24 / V23 |
|---|---|---|---|---|---|---|
| 0.469 | 241 s | 125 s | 2009 s | 369 s | 5990 s | 24.9 |
| 0.646 | 131 s | 69 s | 1089 s | 156 s | 2512 s | 19.2 |
| 1.64 | 24 s | 13 s | 213 s | 18 s | 284 s | 11.7 |
| 3.79 | 5 s | 3 s | 46 s | 3 s | 43 s | 8.3 |
| 11.03 | 1 s | 1 s | 11 s | 1 s | 8 s | 6.5 |

* **3.5 r_e limit (method 1)** halves the Mie cost: the IDL integrated every
  radius to 100 um.
* **Nested refinement (method 3)** is responsible for the slowdown: 2^4 = 16
  times the radius nodes for liquid water, about 16x the cost.
* **Adaptive Legendre (method 2)** costs 1.3-3x the fixed 1000 angles at
  visible wavelengths (orders up to ~3800 + 459 check angles for large
  radii) and less in the infrared; it raises the DISORT time per call by up
  to 25% where L > 1000 (section 4), about 3% over a whole product.
* **Compact runs** (8 channels, 5 radii, 0.55 um reference included): class
  optics 155 s (V23) / 2992-3359 s (V24); compact LUTs 227 s / 3035 s.
* **Full products** (`results/runtime_estimate.json`): the V23 jobs took 17.1
  (p240), 18.6 (pold) and 18.7 h (pagg).  V24 adds about 15 h of Mie and 3%
  of DISORT: **about 32-34 h** for each liquid-water product, **about 19 h**
  for pagg (unchanged).  The SLURM limit is 72 h
  (`scripts/run_v24_release.slurm`); priority-eodg has no partition limit.
  The historical V24 estimate was 35-42 h.
* The cost is that of a converged radius quadrature; no optimisation was
  attempted.  Reusing Mie evaluations across effective radii on the shared
  nested lattice would need a common angle set per wavelength (the historical
  report's recommendation); it is not done here.

## 9. Deviations from the historical methods

1. **DISORT moments:** NMOM = max(1000, L) (owner decision, 2026-10-13) instead
   of max(L, NSTR + 1); bitwise equivalent (section 4).  `nmom = 1000` stays
   in every run file (minimum DISORT moments for Mie classes; the fixed
   expansion for Baum / T-matrix) instead of being refused for Mie classes.
2. **IDL-compatible arithmetic:** the lattice and refinement use V23's FLOAT
   limits (0.001, 100), xres (0.4) and wavenumber, so the lattice nodes are
   V23's legacy nodes; the historical code (pre-V22) used doubles.
3. **Double-precision phase function** for the expansion (`phase_dtype`): the
   V23 `create_bwgp` stores Phi as FLTARR, whose rounding (~6e-8 of a forward
   peak of 1e4-1e5) would exceed King's 1e-9 by orders of magnitude; the bulk
   optics keep the V23 float32 storage and IDL TOTAL order.
4. **Not restored:** the historical `master_run.py`, `runs/template.run` and
   submission-script changes, the 40-product run set and the historical
   cache marker (V24 uses `v24-scattering-refinements-1`).
5. The modules were rewritten in the V23 style; the Gauss-Legendre rules,
   degree bounds, quadrature orders and King lengths are identical to the
   historical module's (checked bitwise up to 6016 points).

No numerical defect was identified in the historical methods.

## 10. All other scientific settings are those of V23

`releases/lut_v24.release` differs from `lut_v23.release` in exactly five
entries (CODE_VERSION, the driver, `create_bwgp`, `generate_scattering_properties`,
the new `legendre_expansion`); every other entry, including DISORT
(`DISORTtoIDL.so` `59307e42...`, its sources and the Intel run time), the Mie
kernel, Baum, readers, writer, NetCDF schema, provenance, and the pinned inputs
(Grid B `liquid-water-cloud.lut` / `ice-cloud.lut` `_v02`, instrument, SRFs,
solar spectrum, microphysics, refractive indices, atmosphere), is copied with
its digest.  The run files equal V23's except `version` and `out_path`.
NSTR 60, nmom 1000, atmosphere 2, no gas, Rayleigh on, srf_quad 1, 36 channels.

## 11. Tests

* New `tests/test_v24_scattering_refinements.py` (11 tests): Gauss-Legendre
  exactness; degree bound against mieint.f; Henyey-Greenstein coefficients and
  King length; default `create_bwgp` bitwise equal to the sealed V23
  `create_bwgp_v03`, including blocked Mie sums; lattice and refinement
  nesting and bounds; class selection and variance refusal; an adaptive
  expansion of a real averaged phase function; `_disort_moments`
  (max rule, padding, refusal of truncation).
* V21/V22 reproduction tests now run the V22 code from its sealed release
  (`tests/sealed_release.py`, `releases/lut_v23.release`):
  `test_compact_case_a_reproduces_v21_idl_bitwise` (passes: V21 IDL products
  and scattering arrays bitwise), the optics port check and the two compact
  product checks of `test_idl_mirror.py`.
* `test_idl_mirror` exercises the run files of this release; the historical
  `lut-code-v24` run files unchanged since `51717b8` are recognised by their
  Git blobs and not exercised.  `test_version` takes the declared version from
  CODE_VERSION.
* Full working suite (non-archived test files) with the V24 code: 400 passed,
  16 skipped, 6 failed; the 6 failures are those of the pre-V24 baseline on
  the same files (388 passed, 17 skipped, 6 failed: files removed in the
  earlier repository clean-up, the port-status inventory and one
  manifest-parser message), none related to V24.

## 12. Production (not submitted)

    bash scripts/submit_v24_release.sh        # from atmlxint1/2/3/5/6 or atmlxint7; three jobs on priority-eodg
    bash scripts/monitor_v24_release.sh --compare

Products: `staging/lut_v24/terra_modis_m_{liquid-water_a01_pold,liquid-water_a01_p240,water-ice_a01_pagg}_v24.nc`.
After production: `python validation/v24/scattering_refinements/compare_full_products.py`
(structure, Grid B, NaN, provenance, pagg identical to V23, V24 - V23
statistics), then `scripts/make_release_dependencies.py` for
`releases/lut_v24.dependencies`, then the release tag.
