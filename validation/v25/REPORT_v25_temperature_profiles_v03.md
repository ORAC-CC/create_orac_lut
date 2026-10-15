# ORAC LUT V25: vertically varying cloud temperature in the thermal emission

Written for the owner of the ORAC LUT generator.  2026-10-15/16.

V25 is the definitive V24 calculation (`releases/lut_v24.release`) with one
addition: a vertically varying cloud temperature in the thermal-emission
DISORT call.  This report covers the recovery of the historical method, its
re-implementation, the correction of two numerical limitations found in
compact validation, and the validation of the final treatment.

**Owner decisions (2026-10-15):**
* the cloud top stays at 250 K, the V21-V24 emission temperature;
* liquid water uses the historical adiabat as it stands;
* the thin-cloud failure and the thick-cloud layering error are corrected,
  both approved.

**Status:** see section 12.

Numerical improvements (sections 5-7) and the prescribed physical profiles
(sections 2-4) are kept separate below.  The profiles are exactly those of
the historical method, with the 250 K top.

## 1. Historical V25 method and source commits

| Commit | Date | Content | Status |
|---|---|---|---|
| `7f6b4e2` | 2026-10-06 08:36 | saturated adiabats for both phases, constant beta, depth caps, 240 K / 628 hPa | superseded (jobs 488174-488213 cancelled) |
| `13b9df3` | 09:56 | caps and altitude removed, z = tau / beta | ice branch superseded (1630 K at tau 256) |
| `5722280` | 18:18 | ice: the supplied cirrostratus profile | final ice method |
| `990d970` | 21:45 | liquid: the `13b9df3` wet adiabat | final liquid method; historical products |
| `f066fdd` | 10-09 | preservation (tag `lut-code-v25`) | source of record |

The final historical method (`validation/v25/results/REPORT_v25_cloud_temperature_profile.md`):
* the emission call splits the cloud layer into 40 equal-optical-depth
  sub-layers;
* the boundary temperatures are the profile at the cumulative particle
  optical-depth fraction F of the node's tau_055;
* E_md is normalised by B(T_top) and may exceed 1;
* the diffuse and direct-beam calls are unchanged.

The historical cloud top was 240 K, commented "as the existing calculation",
whereas every earlier version uses 250 K.  The "unexplained anomaly" of the
historical V25 tables (`LUT_CODE_VERSIONS.md` section 6) was never written
down.  The thin-cloud failure of section 5 is present in the historical ice
products and is the probable cause.

## 2. Cirrostratus profile data (water ice)

* **File:** `references/data/ocalut_cloudprofile_Cirrostratus.dat`, committed
  in `5722280`; 92 819 bytes.
  * md5 `7d7b472f3ff9681901c03afe85c3202b`.
  * sha256 `8e9e68bba7f65ef0231dd309d1c129a0292fed64b9a49f53ec70458d23042aa6`.
  * Pinned in `releases/lut_v25.release`.
* **Header:** `Origin: /tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf`,
  `Cloud type: Cirrostratus`.
* **Supplier:** P. Watts' vertically inhomogeneous LUT tests for OCA
  (EUMETSAT).  No document names the Met Office; provenance beyond the header
  is not recorded.
* **Content:** 13 profiles (COT 0.0625-256) of 99 rows.
* **Checks** (`read_cloud_profile` on every read): one extinction shape for
  all COTs, the layer and cumulative optical depths consistent, rows equally
  spaced in depth, ΔT = 8.000 K/km × depth.
* **Interpretation:** only the extinction shape and the depth H(COT)
  (2-11 km) are data; the temperature is a constant 8 K/km lapse.
* **Mapping:** through the cumulative optical-depth fraction F (midpoint rule
  for the rows):
  * T(F; tau) = 250 K + 8 K/km × H(tau) × zeta(F);
  * H is linear in COT, held at 2 km below COT 0.0625, never extrapolated.
  * No extinction coefficient enters.

## 3. Liquid-water saturated adiabat

* **Equations:** dT/dz = Γ_w(T, p) (pseudo-adiabat over liquid water) and
  dp/dz = p g / (R_d T).
* **Inputs:** e_sw from Murphy and Koop (2005) Eq. 10; L_v from Rogers and
  Yau (1989).
* **Integration:** RK4 in 10 m steps from (250 K, 628 hPa) along
  z = F tau_055 / 20 km^-1, with no cap.
* **Reproduction:** from 240 K the code gives the historical 317.49 K and
  2957.2 hPa at 12.8 km (test).
* **V25 values:** tau 16: 256.45 K; tau 64: 273.37 K; tau 256: 320.13 K and
  2867 hPa.
* **Limitation:** beyond tau_055 = 74 (z = 3.7 km) the pressure exceeds
  1013 hPa, so the deep states are not a physical cloud.  This was accepted by
  the owner as the historical method.

## 4. Conventions

* **Temperature and normalisation:** the cloud top is at 250 K, the
  normalising temperature of E_md, as in V24 (owner decision).
* **Depth convention:**
  * depth increases downward from the top (F = 0);
  * DISORT's TEMPER(0) is the cloud top;
  * one temperature profile per node (tau_055), the same for every channel;
  * each channel's sub-layers follow its own optical depth D = tau_055 ×
    extinction ratio (+ Rayleigh), as V24.
* **Interpretation of E_md:** unchanged, E_md = I / B(250 K).  It may exceed
  1, so its NetCDF valid range is [0, FLT_MAX] in V25.

## 5. Thin clouds: cause of the historical failure and the correction

**Cause.**
* The cirrostratus profile keeps a 2 km, 16 K span as tau → 0.  The 40 equal
  sub-layers then have optical depths of order tau/40 but finite temperature
  steps.
* DISORT's thermal source in a layer is B(t) = XR0 + XR1 t with
  XR1 = ΔB/Δτ.  It is expanded about the column top (XR0 = B − XR1 TAUCPR).
* In single precision the particular solution then loses about
  2×10⁻⁸ (ΔB/B)/Δτ² of the emitted radiance
  (`temperature_profiles/conditioning.py`; a single gradient layer against the
  exact radiance: 1.9×10⁻² at Δτ 10⁻³, 3.5×10⁻⁴ at 10⁻², 3×10⁻⁵ at
  3×10⁻²).
* At tau = 1e-10 the result is garbage: E_md = −1.45 (ch 20), −0.41 (ch 31)
  and −0.59 (ch 32) in the compact case.  The historical V25 product has
  −1.66…+0.55 at tau = 1e-10 and negative values at 0.0078.
* Liquid water was not affected, because the adiabat's span vanishes with
  tau.

**Correction (`cloud_temperature.emission_layers`).**
* The layering is chosen by the numerical resolution DISORT can achieve, not
  by cloud type.
* A linear-source sub-layer is accepted only if
  (|ΔB_k| / B_mean) / Δτ_k² ≤ `PLANCK_SLOPE_LIMIT` = 1000 (round-off
  ≲ 2×10⁻⁵).
* When no linear layering of at least 23 sub-layers satisfies this (an
  optically thin cloud whose temperature spans the profile), the cloud
  becomes an **isothermal staircase**:
  * K = min(23, floor(D / 10⁻³)) equal cells, each isothermal at the
    temperature whose window Planck radiance is the cell mean of B[T(F)];
  * the cells are joined by layers of zero optical depth that carry the
    temperature steps (DISORT gives them no source slope);
  * isothermal cells avoid the cancellation entirely; cells thinner than
    10⁻³ only add DISORT's ordinary thin-layer round-off, so below that the
    cloud is a single layer;
  * the profile is never switched off.

**Optically thin limit.**
* For a layer of optical depth D and single-scattering albedo w, emitting
  B(t), the emergent radiance is

      I(mu) = (1 - w) / mu ∫_0^D B(t) dt + O(D²) = (1 - w) D B_mean / mu + O(D²),
      B_mean = (1/D) ∫_0^D B_nu[T(t)] dt.

* The O(D²) terms are the attenuation of the deeper emission and the
  scattering of the cloud's own emission.  To first order the source enters
  only through the optical-depth mean of the Planck radiance.
* The mean is taken over the **radiance**, not the temperature: B_nu(T) is
  convex, and its curvature differs between windows, so the effective
  temperature does too (thin cirrostratus: 259.03 K at 3.7 µm, 258.70 K at
  11 µm).
* The mean is computed in the channel's DISORT window (PLKAVG of the centre
  wavenumber ± 0.5 %, the V24 convention) by Gauss-Legendre quadrature over
  the profile's kinks.
* The staircase preserves ∫B over every cell, so it is exact to first order
  in D and contains the second-order effects through its cells.  As
  D → 0, E_md → (1 − w) D B_mean / (mu B(250 K)) → 0.
* Absorption and scattering are DISORT's: the staircase only changes the
  source.

## 6. Thick clouds: non-uniform emission layering

* **Layering:** n = min(46, floor(D/h)) linear sub-layers, MXCLY = 46 (the
  production cloud is one layer).
  * Equal sub-layers while D ≤ n h; otherwise a geometric progression from a
    top sub-layer h = `TOP_SUBLAYER_TAU` = 0.01, its ratio found by
    bisection.
  * h is doubled only if the resolution criterion of section 5 fails.
  * The emitting top of the cloud is resolved (about 20 sub-layers in the
    top optical depth of a tau 256 cloud), and every sub-layer keeps a
    similar depth-to-thickness ratio, which also limits DISORT's
    cancellation.
* **Source:** the L2 projection of B[T(F)] onto continuous piecewise-linear
  functions on the sub-layers (mass-matrix solve).  This preserves the
  integral of B, so the linear and staircase representations share the same
  thin limit.
* **DISORT inputs:**
  * DISORT receives at each boundary the temperature whose window Planck
    radiance is the projected value (Newton on PLKAVG; ≤ 1e-5 in radiance,
    the resolution of PLKAVG with a float32 temperature);
  * these boundary temperatures follow the profile to ≤ 0.9 K (deep, coarse
    sub-layers of tau 256 ice) and ≤ 0.1 K at the top;
  * the sub-layer optical depths add up to the layer's to ≤ 1 float32 ulp;
  * none is negative;
  * the single-scattering albedo and phase moments are the layer's.

## 7. Validation of the numerical treatment

**Independent references** (`temperature_profiles/layering.py`):
* the emission-call inputs of real cases, recorded from LUT runs
  (`emission_inputs.py`): ice agg r_e 5/33/93 µm, liquid water p240 r_e
  3/11/35 µm and pold r_e 11 µm; channels 20, 31, 32; all Grid B optical
  depths; view zeniths 0/60/75°;
* exact radiance without scattering (albedo set to 0, analytic integral on
  4000 nodes refined at the top);
* convergence with scattering (46 against 34 sub-layers; Richardson estimate
  E_inf);
* the first-order thin-cloud radiance.

A superposition of hat-function responses could not serve as a scattering
reference: single-precision DISORT is wrong by up to 0.56 for such steep
sources (recorded in `results/layering.json`).

**Thick clouds** (tau_055 = 1-256; maximum |ΔE_md|):

| | historical layering | V25 layering |
|---|---|---|
| no scattering vs exact, ice ch 20 / 31 / 32 | 1.2×10⁻² / 3.3×10⁻³ / 3.1×10⁻³ | 3.0×10⁻⁶ / 3.0×10⁻⁶ / 2.4×10⁻⁶ |
| no scattering vs exact, liquid ch 20 / 31 / 32 | 1.0×10⁻³ / 4.1×10⁻⁵ / 2.0×10⁻⁵ | ≤ 2.0×10⁻⁶ |
| with scattering vs E_inf, ice ch 20 / 31 / 32 | 1.1×10⁻² / 2.8×10⁻³ / 2.6×10⁻³ | 8.1×10⁻⁵ / 4.8×10⁻⁶ / 1.4×10⁻⁵ |
| with scattering vs E_inf, liquid | ≤ 1.5×10⁻³ | ≤ 4.6×10⁻⁶ |

* The largest V25 value (ice ch 20, tau 16) is the Richardson estimate's own
  uncertainty: 46 against 34 layers converges monotonically (23 → 28 → 34 →
  46).
* The result is insensitive to h = 0.003-0.03 and to a slope limit of 300.

**Thin clouds** (tau_055 = 1e-10 … 0.25, 147 cases):
* At tau = 1e-10, |E_md| ≤ 1.5×10⁻¹⁰ (ice; the V24 isothermal call gives
  7×10⁻¹¹); historical ice −1.45.
* No negative E_md from tau 0.0078.
* Without scattering, the error is ≤ 8×10⁻⁵ relative.
* E_md is 0.998-1.025 times the first-order radiance at tau = 0.0078 (the
  O(D/mu) terms).

**Continuity** (sweeps of the channel optical depth 10⁻⁵-5 for ice r_e 5/33
and liquid r_e 11, channels 20/31; every change of representation located by
bisection):
* changes of staircase cell count: ≤ 4.8×10⁻⁶ in E_md (≤ 9.5×10⁻⁴
  relative; DISORT's thin-layer round-off, the same for the nearly isothermal
  liquid case);
* staircase → linear and the added linear sub-layers: ≤ 5.3×10⁻⁵
  (≤ 4×10⁻⁴ relative);
* no negative values;
* largest error without scattering 1.2×10⁻⁴ (ice ch 20, in the staircase
  range near D = 0.2).

**Maximum remaining numerical error** (E_md):
* about 1×10⁻⁴ for thick clouds and in the staircase-linear transition
  (channel optical depth 0.1-0.5, ice, 3.7 µm);
* ≤ 10⁻⁵ elsewhere;
* below 10⁻³ relative for thin clouds;
* compared with the historical 1.2×10⁻² (thick) and −1.45 (thin).

## 8. Compact V24/V25 comparisons (Terra MODIS)

**Setup:**
* Grids `compact_{liquid,ice}_grid_b.lut`: all 24 Grid B optical depths,
  r_e 3/11/35 µm and 5/33/93 µm, SZA 0/30/60°, VZA 0/30/60/75°, three
  azimuths.
* Channels 1, 6, 7 (solar), 20 (3.7 µm), 31 and 32.
* Settings as in production; models liquid-water_240, liquid-water_old and
  water-ice_agg.
* The V24 baseline is the sealed `releases/lut_v24.release`.
* Round 2 (final code, SLURM 496973-496990, `results/compact_*.json`;
  round 1 = the candidate, in `results/candidate1/`).

**Liquid water** (p240 / pold):
* E_md − V24 from 0 to +0.116 / +0.117 (ch 20), ≤ +0.033 / +0.042 (ch 31)
  and ≤ +0.018 / +0.025 (ch 32), growing with tau and saturating from
  tau ≈ 16;
* mirrored profile: −0.088 / −0.091 … 0;
* thin clouds: ratio = thin limit (1.000-1.003);
* 34 against 46 sub-layers: ≤ 3.8×10⁻⁶;
* E_md B(250) ≤ B(T_base) everywhere.

Three thin liquid nodes (tau ≤ 0.031, r_e 35 µm, ch 32), where the physical
change is about 10⁻⁷ (ΔT ≈ 0.002 K), have |ΔE_md| ≤ 2.9×10⁻⁶ of the wrong
sign.  This is below the single-precision noise of the isothermal limit
(8×10⁻⁶).  Along tau the saturated plateau shows the same float32 ripples
(about 10⁻⁵) as the candidate.


**Test A (regression):** without a profile V25 equals sealed V24 **bitwise in
every variable**, in all three models.

**Test F (channels):** with a profile only E_md differs.  Every solar R/T
operator, including the reflected solar part of channel 20, and every
optical property is bitwise V24 in every variant and model.

**Test C (thin), ice agg** (V25/V24 ratio over views and r_e; thin limit
B_mean/B(250)):

| tau_055 | ch 20 | limit | ch 31 | limit |
|---|---|---|---|---|
| 1e-10 | E_md ±2×10⁻¹⁰ | – | 0 | – |
| 0.0078 | 1.696-1.698 | 1.691 | 1.193 | 1.192 |
| 0.0156 | 1.694-1.697 | 1.691 | 1.193 | 1.192 |
| 0.0312 | 1.691-1.698 | 1.691 | 1.192-1.193 | 1.192 |
| 0.0625 | 1.684-1.697 | 1.691 | 1.190-1.193 | 1.192 |
| 0.125 | 1.705-1.731 | 1.726 | 1.195-1.200 | 1.199 |
| 0.25 | 1.744-1.801 | 1.798 | 1.203-1.214 | 1.215 |

* The profile acts at every optical depth; the ratio is smooth in tau.
* Round 1 (historical layering) had ratios 1.69-2.42 at tau = 0.0078 and
  E_md = −1.47 at 1e-10.

**Test E (direction), ice:**
* warming profile: E_md − V24 from 0 to +2.90 (ch 20), +0.42 (ch 31), +0.36
  (ch 32), never negative;
* mirrored (cooling) profile: −0.58…0, never positive;
* E_md B(250) ≤ B(T_base) everywhere;
* one maximum along tau (tau ≈ 8-16, the supplied H(COT));
* finite everywhere.
* Emitted brightness temperature, ch 31, nadir, r_e 33:

| tau | 0.0625 | 0.5 | 2.8 | 8 | 22.6 | 64 | 181 |
|---|---|---|---|---|---|---|---|
| V24 (K) | 148.8 | 192.3 | 236.2 | 248.9 | 249.9 | 249.9 | 249.9 |
| V25 (K) | 151.8 | 198.6 | 251.5 | 266.1 | 262.4 | 256.3 | 252.8 |

**Test D (convergence), ice:** 34 against 46 sub-layers: ≤ 6.9×10⁻⁵ (ch 20),
≤ 1.2×10⁻⁵ (ch 31/32).

## 9. Isothermal limit and regression

| Test | Ice agg | Liquid p240 | Liquid pold |
|---|---|---|---|
| V25 without profile against sealed V24 | bitwise, every variable | bitwise, every variable | bitwise, every variable |
| constant profile through the V25 layering: max abs / max rel (E_md > 10⁻³) | 1.3×10⁻⁵ / 7.9×10⁻⁴ | 8.2×10⁻⁶ / 3.9×10⁻⁴ | 5.7×10⁻⁶ / 4.7×10⁻⁴ |

* At unit level (`tests/test_v25_cloud_temperature.py`) a constant profile
  gives boundary temperatures of exactly 250 K at every optical depth.
* **Release launcher:** `temperature_profiles/launcher_smoke.run` (the V25
  pagg run file on the compact ice grid, channels 1/6/7/20/31/32), run through
  `releases/lut_v25.release`.  108 files were verified and 35 modules loaded
  from the sealed copies.  The product equals the working-tree compact profile
  run in all 52 variables.
* The test suite gives the 6 pre-existing failures only (section 12).

## 10. Historical V25 comparisons

* **Recovered method:** at the historical 240 K top with the historical
  layering (round 1), the code reproduced the historical V25 Terra products to
  ≤ 9×10⁻⁵ (ice, tau ≥ 0.0625) and ≤ 1.2×10⁻⁵ (p240).  It reproduced the
  historical compact profile effect to 6.8×10⁻⁵.
* **Final V25 at 240 K against the historical products:**
  * ice, tau ≥ 0.0625: ≤ 1.3×10⁻² (ch 20) and ≤ 3.1×10⁻³ (ch 31/32);
  * against the historical compact profile effect: ≤ 2.6×10⁻³.
  * These are the historical layering errors of section 7 (1.2×10⁻² /
    3.3×10⁻³).
* **Liquid water, final V25 at 240 K against the historical p240 / pold
  products:**
  * tau ≥ 0.0625: ≤ 2.4×10⁻³ (ch 20) and ≤ 1.2×10⁻⁴ (ch 31/32), the
    historical liquid layering error (round 1 against round 2:
    1.5×10⁻³ / 4×10⁻⁵);
  * thin: ≤ 1.7×10⁻⁵.
* **Thin ice, tau < 0.0625:** the historical failure (−1.66…+0.55 at
  1e-10) is not reproduced.
* **250 K against 240 K:** E_md is lower by up to 0.58 (ch 20) and 0.044
  (ch 31).

## 11. Runtime

* **Compact (ice):** an emission call costs 35 ms with the V25 layering
  against 4 ms (V24).  The run took 162 s against 155 s (isothermal).
* **Compact (liquid):** emission takes 7.5-8.9 s against 1.1 s per 216 calls (about +36 ms per call).  The run totals of 1250-1300 s are dominated by the Mie calculation and vary more than that between jobs.
* **Production:** from the V24 logs, 9 216 (pold, p240) and 11 136 (pagg)
  emission calls at about +40 ms each (DISORT plus layering, measured), i.e. +0.10 h (liquid) and +0.12 h (ice).
* **Measured V24 times:** pold 31.4 h, p240 33.2 h, pagg 18.8 h.
* **V25 estimates:** pold ≈ 31.5 h, p240 ≈ 33.3 h, pagg ≈ 18.9 h (the Mie and the other DISORT calls are as in V24).

## 12. Assumptions, limitations, status

**Physical:**
* the 250 K reference top (E_md depends on it through the nonlinearity of B);
* one cirrostratus profile for all ice habits;
* a constant 8 K/km lapse;
* the liquid adiabat's prescribed beta and its unphysical deep states;
* the interpretive choices: midpoint rule, H linear in COT, H held below
  COT 0.0625.

**Numerical:**
* single-precision DISORT round-off (about 10⁻⁴ in the staircase-linear
  transition, 10⁻³ relative for thin clouds; the V24 isothermal call has the
  same order);
* `emission_layers` supports one in-cloud layer, as in every production
  microphysics file (refused otherwise).

**Status:**
* compact validation passed;
* `releases/lut_v25.release` sealed;
* committed and pushed;
* full production not submitted: see `V25_RELEASE.md`.
