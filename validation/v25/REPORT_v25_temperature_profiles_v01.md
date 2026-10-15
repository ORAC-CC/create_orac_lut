# ORAC LUT V25: vertically varying cloud temperature in the thermal emission

Written for the owner of the ORAC LUT generator.  2026-10-15.

**Status: implementation and compact validation done; STOPPED for an owner
decision before the release is prepared (section 9).**  Nothing has been
committed, no manifest has been sealed and no production job has been
submitted.  The re-implemented historical method is reproduced exactly, but it
carries a numerical defect for optically thin ice cloud.  The historical V25
ice products have the same defect.

## 1. Historical V25 method and source commits

| Commit | Date | Content | Status |
|---|---|---|---|
| `7f6b4e2` | 2026-10-06 08:36 | saturated adiabats for both phases, constant beta (liquid 20, ice 1 km^-1), depth caps, cloud top 240 K / 628 hPa from the production atmosphere | superseded (jobs 488174-488213 cancelled) |
| `13b9df3` | 09:56 | caps and altitude removed: z = tau / beta, no limit | ice branch superseded (1630 K at tau 256) |
| `5722280` | 18:18 | ice: the supplied cirrostratus profile | final ice method |
| `990d970` | 21:45 | liquid: the `13b9df3` wet adiabat restored | final liquid method; products made here |
| `f066fdd` | 10-09 | preservation commit (tag `lut-code-v25`) | historical source of record |

Final historical method (`990d970`; `validation/v25/results/REPORT_v25_cloud_temperature_profile.md`):

* E_md from one emission DISORT call in which each in-cloud layer is split
  into 40 equal-optical-depth sub-layers with the layer's single-scattering
  albedo and moments.  The boundary temperatures follow the cumulative particle
  optical-depth fraction F of the node's tau_055.  E_md is normalised by
  B(T_top) and may exceed 1.  The diffuse and direct-beam calls are unchanged.
* **Ice:** T = T_top + 8 K/km × H(tau_055) × zeta(F).
  * zeta(F) comes from the supplied extinction shape (midpoint rule).
  * H is linear in COT between the supplied profiles: 2 km at COT 0.0625,
    held at 2 km below it, 11 km at COT 128-256.
  * The mapping is **through cumulative optical depth**; no extinction
    coefficient or altitude enters.
* **Liquid:** a saturated liquid-water pseudo-adiabat, integrated by RK4 in
  10 m steps, from (T_top, 628 hPa) along z = F tau_055 / 20 km^-1.
  * The mapping is **through a prescribed physical path length**.
  * There is no cap, no altitude, no surface and no freezing limit.
* **Cloud top:** 240 K, commented "as the existing calculation"; in fact the
  IDL and every Python version up to V24 use 250 K.
* **Removed during the history:** H_max, depth caps, beta for ice, the ice
  adiabat, altitude and atmosphere dependence.

The "unexplained anomaly in the V25 tables" (`LUT_CODE_VERSIONS.md` section 6)
was never described.  The thin-cloud defect of section 8.3 is present in the
historical ice products and is a probable candidate.

## 2. Cirrostratus profile data and provenance

* **File:** `references/data/ocalut_cloudprofile_Cirrostratus.dat`, committed
  in `5722280` and unchanged since.
  * 92 819 bytes, md5 `7d7b472f3ff9681901c03afe85c3202b`.
  * sha256 `8e9e68bba7f65ef0231dd309d1c129a0292fed64b9a49f53ec70458d23042aa6`.
* **Header:** `Origin: /tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf`,
  `Cloud type: Cirrostratus`.
* **Supplier:** P. Watts' vertically inhomogeneous LUT tests for OCA
  (EUMETSAT), according to the historical report.  Neither the file nor the
  repository names the Met Office; the provenance chain beyond the header is
  not recorded.
* **Content:** 13 profiles, COT 0.0625-256 (doubling), each of 99 rows.  The
  columns are % extinction, layer dCOT, cumulative COT, km below the cloud top
  and K above the cloud top.
* **Checks** (made by `read_cloud_profile` on every read):
  * one extinction shape for all COTs, summing to 100;
  * dCOT = COT f / 100, and the cumulative column is its running sum;
  * rows equally spaced in depth from 0 to H(COT);
  * ΔT = 8.000 K/km × depth in every row.
* **Interpretation:** the temperature is a constant 8 K/km lapse rate.  Only
  the extinction shape and H(COT) are data.

## 3. Liquid-water saturated adiabat

* **Equations:**
  * dT/dz = Γ_w(T, p), with
    Γ_w = g (1 + L_v r_s / (R_d T)) / (c_pd + L_v² r_s ε / (R_d T²)) and
    r_s = ε e_sw / (p − e_sw);
  * dp/dz = p g / (R_d T).
* **Inputs:** e_sw from Murphy and Koop (2005) Eq. 10; L_v from Rogers and Yau
  (1989); g 9.80665, R_d 287.04, c_pd 1005.7, ε = R_d / 461.5.
* **Integration:** RK4 in 10 m steps from (T_top, 628 hPa) to max(tau)/20 km^-1.
* **Reproduction:** from 240 K the new code gives the historical values
  317.49 K and 2957.2 hPa at 12.8 km (`tests/test_v25_cloud_temperature.py`).
* **V25** (from 250 K):

| tau_055 | z (km) | T (K) | p (hPa) |
|---|---|---|---|
| 1 | 0.05 | 250.42 | 632 |
| 16 | 0.8 | 256.45 | 700 |
| 64 | 3.2 | 273.37 | 953 |
| 128 | 6.4 | 291.46 | 1403 |
| 256 | 12.8 | 320.13 | 2867 |

The deep states are not a physical cloud: the pressure exceeds 1013 hPa
beyond z = 3.7 km, which is tau_055 = 74, at 276.5 K.  This was accepted by the owner (2026-10-15) as
the final historical method.

## 4. Cloud-top temperature and depth conventions

* **Owner decision 2026-10-15:** the cloud top is at **250 K**, the V21-V24
  emission temperature (`create_orac_cloud_lut.pro` line 561).
  * E_md is normalised by B(250 K).
  * ORAC's retrieved cloud temperature keeps its V24 meaning, the temperature
    of the cloud top.
  * With a constant profile, V25 is V24.
  * At 240 K (the historical value) E_md would be up to 0.58 larger at 3.7 µm
    and 0.044 larger at 11 µm (section 8.6).
* **Depth convention:** depth increases downward from the cloud top (z = 0, F = 0).
  * TEMPER(0) of DISORT is the cloud top.
  * The analytic test of section 8.2 verifies the order.

## 5. Source changes

| File | SourceVault | Change |
|---|---|---|
| `src/oraclut/idl_mirror/cloud_temperature.py` | new, `_v01` | profile reader and check, H(COT), zeta(F), wet adiabat, `setup_cloud_profile`, `cloud_temperature`, `emission_layers`; plain functions and structures (no class hierarchy) |
| `create_orac_luts.py` | `_v05` → `_v06` | run-file key `cloud_vertical_profile` (absent = `isothermal` = V24); substance pairing; profile set up once per LUT; in the emission block the in-cloud layer is replaced by `emission_layers` when a profile is set; E_md valid range; provenance text `cloud_emission` |
| `src/oraclut/idl_mirror/write_v2_lut.py` | `_v03` → `_v04` | optional `em_valid_range` (default [0, 1], as V24); V25 writes [0, FLT_MAX] because E_md may exceed 1 |
| `CODE_VERSION` | `_v06` → `_v07` | `lut_version = 25`, `code_release = lut-code-v25.2` (format required by `oraclut.version`; no tag; part of the deferred naming review) |
| `tests/test_idl_mirror.py` | `_v04` → `_v05` | historical V25 run files (`f066fdd`) excluded like the historical V24 ones |
| `tests/test_v25_cloud_temperature.py` | new | 9 tests (section 8.2) |
| `runs/terra_modis_cloud_{liquid-water_old,liquid-water_240,water-ice_agg}_v25.run` | Git (historical content at `f066fdd`) | V24 run files with `version = 25`, `staging/lut_v25`, `cloud_vertical_profile` |

The diffuse and direct-beam calls are untouched.  So are the scattering
properties, the Mie, Baum and DISORT code, the grids, the NetCDF variables and
dimensions, and the provenance infrastructure.

## 6. V24 baseline and scattering unchanged

* `releases/lut_v24.release`: 106 entries, 0 hash mismatches.
  * Driver `create_orac_luts_v05`; DISORT `DISORTtoIDL.so` sha256 `59307e42…`.
  * NSTR 60; NMOM = max(1000, L); Grid B.
* V24 products: pold `e45035bd…`, p240 `1287c5db…`, pagg `ffc61063…`.  The
  copies in the shared ORAC LUT directory have these checksums.
* Compact runs, against the sealed V24 release:
  * **without a profile, the V25 code is bitwise V24** in every variable of
    every case;
  * with any profile variant, **E_md is the only variable that changes**: all
    R/T operators (including channel 20), extinction, SSA and g are bitwise V24.

## 7. Compact validation setup

* **Grids:** `validation/v25/temperature_profiles/compact_{liquid,ice}_grid_b.lut`.
  * All 24 Grid B optical depths (1e-10 … 256).
  * r_e 3/11/35 µm (liquid) and 5/33/93 µm (ice).
  * SZA 0/30/60°, VZA 0/30/60/75°, three azimuths.
* **Channels:** Terra MODIS 1, 6, 7 (solar), 20 (3.7 µm, mixed), 31 and 32
  (thermal).
* **Settings:** production (atmosphere 2, Rayleigh, no gas, NSTR 60,
  single-precision Intel DISORT).
* **Models:** liquid-water_240 (p240), liquid-water_old (pold) and
  water-ice_agg (Baum, pagg).
* **Variants** (`compact_lut.py`; SLURM 496906-496923):
  * v24: the sealed release;
  * isothermal;
  * profile;
  * constant: profile path with all temperatures 250 K;
  * mirror: T = 500 K − T(F);
  * ttop240: the historical cloud top;
  * layers20 and layers46.

## 8. Results

### 8.1 Identity

All cases: the isothermal variant equals V24 in every variable; every profile
variant differs from V24 in E_md only (`results/compact_*.json`).

### 8.2 Planck source and the DISORT interface (`planck_check.py`, tests)

* **Test case:** a non-scattering cloud on DISORT's cold black lower boundary,
  which has an exact top radiance.
* **Agreement:** with the 40 sub-layers made by `emission_layers`, DISORT
  agrees with the exact integral of the linear-in-τ source to ≤ 4×10⁻⁷
  (cirrostratus τ 1-256 and wet adiabat τ 1-256; channels 20, 31, 32; VZA
  0-75°).
* **Order check:** reversing the boundary order changes the result by
  7 %-4700 %.
* **Conclusion:** each boundary temperature enters at its depth, top down.

### 8.3 Optically thin ice: defect inherited from the historical method

* **Cause:** below COT 0.0625 the cirrostratus cloud keeps H = 2 km (ΔT_base
  16 K) as τ → 0.  The 40 sub-layers then become optically tiny while their
  boundary temperatures differ.  Single-precision DISORT's linear-in-τ Planck
  source (slope ΔB/Δτ) loses all precision in that case.
* **Non-scattering, DISORT against exact:**

| tau_055 | 1e-6 | 1e-4 | 1e-3 | 1e-2 | 3e-2 | 0.1 |
|---|---|---|---|---|---|---|
| ch 20 relative error | 7×10⁴ | 2.5 | 7×10⁻² | 1.1×10⁻³ | 8×10⁻⁵ | 7×10⁻⁶ |

  * At τ ≤ 1e-8 DISORT returns 0.
* **Compact ice LUT (with scattering):**
  * τ = 1e-10 (the Grid B "no cloud" node): E_md = −1.47…−1.42 (ch 20),
    −0.42…−0.20 (ch 31), −0.59…−0.26 (ch 32), where it should be ≈ 2×10⁻¹⁰.
  * τ = 0.0078: the V25/V24 ratio is 1.69-2.42 (ch 20) and 1.08-1.20 (ch 31),
    against the thin limits 1.691 and 1.192.
  * τ = 0.0156-0.031: within about 6 % and 2 % of the limits.
  * τ ≥ 0.0625: consistent.
* **Historical V25 Terra pagg product (2026-10-06):** E_md −1.66…−1.60
  (ch 20) and −0.46…−0.22 (ch 31) at τ = 1e-10, and negative values down to
  −0.0072 at τ = 0.0078.
* **Liquid water is not affected:** the adiabat's span vanishes with τ, so
  dT/dτ = Γ/β stays finite.

### 8.4 Ice cloud, tau ≥ 0.0625 (water-ice_agg)

* **Change in E_md** (V25 − V24, all nodes):
  * ch 20: +0.0055 … +2.90 (V25 max 3.75);
  * ch 31: +0.0048 … +0.42;
  * ch 32: +0.0049 … +0.36.
  * The change is never negative.
  * Each curve has one peak, at τ ≈ 11-16, as the supplied H(COT) prescribes.
  * E_md is finite everywhere, and E_md B(250) ≤ B(T_base) holds with margin
    ≥ 0.56.
* **Mirror** (cooling with depth): −0.58 … −0.003, never positive.
* **Emitted brightness temperature, ch 31, nadir, r_e 33:**

| tau | 0.25 | 2.8 | 11.3 | 45 | 181 |
|---|---|---|---|---|---|
| V24 (K) | 175.7 | 236.2 | 249.7 | 249.9 | 249.9 |
| V25 (K) | 180.5 | 251.5 | 266.0 | 258.0 | 252.6 |

* **Isothermal limit** (constant variant against V24): |ΔE_md| ≤ 2.5×10⁻⁵.
  * This is single-precision DISORT's sensitivity to splitting one layer into
    40.
  * It is at most 2.2×10⁻³ relative where E_md > 10⁻³.
* **Layer count** (thermal channels):
  * 46 against 40 layers: ≤ 2.7×10⁻³ (ch 20) and ≤ 1.3×10⁻³ (ch 31/32);
  * 20 against 40: ≤ 2.0×10⁻².
  * Against a converged 4000-layer reference (non-scattering, section 8.2
    script), 40 equal-τ layers leave errors up to 1.2×10⁻² B(250) at 3.7 µm
    and 3.5×10⁻³ at 11-12 µm for thick ice (τ ≥ 64), about 0.2 K.  The
    sub-layers there are 1.6-6.4 optical depths thick.  Layers graded toward
    the top (fractions (k/40)³) would reduce this to ≤ 5×10⁻⁴.

### 8.5 Liquid water (liquid-water_240, liquid-water_old)

* **Identity:** as in 8.1.
* **Change in E_md** (V25 − V24, every node from τ = 1e-10):
  * ch 20: 0 … +0.117 (p240) and +0.118 (pold);
  * ch 31: ≤ +0.033 / +0.042;
  * ch 32: ≤ +0.018 / +0.025.
  * The change is never negative, and E_md B(250) ≤ B(T_base) holds.
  * The change grows with τ and saturates from τ ≈ 16.  The saturated plateau
    carries float32 ripples of about 10⁻⁵ and a slow rise of about 2×10⁻⁵ to
    τ = 256 (the deep warm adiabat).
* **Mirror:** −0.088 … 0, never positive.
* **Thin cloud:** well conditioned.  At τ ≤ 0.031 the V25/V24 ratio is
  1.000 ± 0.003 (the thin limit is 1.000, since ΔT → 0 with τ), and at
  τ = 1e-10 |E_md| ≤ 3×10⁻¹⁰.
* **Isothermal limit:** |ΔE_md| ≤ 1.6×10⁻⁵.
* **Layer count** (p240): 46 against 40 ≤ 3.6×10⁻⁴ (ch 20) and ≤ 2×10⁻⁵
  (ch 31/32); 20 against 40 ≤ 4×10⁻³.
* **Emitted brightness temperature, nadir, r_e 11 µm, p240:**

| tau | 0.0625 | 0.5 | 2.8 | 8 | 22.6 | 64 | 181 |
|---|---|---|---|---|---|---|---|
| ch 31, V24 (K) | 149.6 | 193.5 | 237.4 | 249.2 | 249.9 | 249.9 | 249.9 |
| ch 31, V25 (K) | 149.6 | 193.6 | 237.8 | 249.9 | 250.8 | 250.8 | 250.8 |
| ch 20, V24 (K) | 186.0 | 207.7 | 230.8 | 242.8 | 246.5 | 246.6 | 246.6 |
| ch 20, V25 (K) | 186.0 | 207.8 | 231.2 | 244.1 | 248.4 | 248.6 | 248.6 |

* **Historical V25 p240 product:** the ttop240 variant agrees to ≤ 1.2×10⁻⁵
  (ch 20/31/32) at every compact node, including τ = 1e-10.  The V24
  scattering refinements therefore leave E_md of the historical liquid V25
  practically unchanged at these nodes.  The historical liquid products have
  no thin-cloud defect.
* **Runtime:** compact runs take 1244 s (V24) against 1246 s (V25); emission
  takes 0.8 s against 7.6 s.

### 8.6 Comparison with historical V25

* **Profile effect:** at the historical compact nodes (63: ch 20/31/32,
  τ 0.0625-256, r_e 5/33/93, nadir), the ttop240 variant reproduces the
  historical V25 − V24 change to 6.8×10⁻⁵ (1.1 %).
  * The new V24 E_md equals the historical V24 E_md to 8×10⁻⁷ at these nodes:
    the Baum path is unchanged.
* **Full historical V25 Terra pagg product:** ttop240 agrees to ≤ 9.4×10⁻⁵
  (ch 20) and ≤ 2×10⁻⁵ (ch 31/32) at all compact nodes with τ ≥ 0.0625.
  * At τ < 0.0625 the two differ by up to 1×10⁻² (ch 20): the defect of
    section 8.3, present in both.
* **Effect of the 250 K top** (owner decision): E_md is lower than at 240 K by
  up to 0.58 (ch 20) and 0.044 (ch 31).
* **Liquid water** (section 8.5): the historical p240 product is reproduced to
  ≤ 1.2×10⁻⁵.
* **Conclusion:** the recovered method is implemented exactly.  The differences
  from historical V25 come from the 250 K top (owner decision).  The
  historical ice products share the thin-cloud defect (section 8.3).

## 9. Decisions needed before the release

1. **Thin ice cloud (required).** The historical behaviour is not acceptable:
   E_md is about −1.5 at the τ = 1e-10 node that ORAC interpolates from.  The
   options are below.  Each changes only the cirrostratus path at small τ, and
   each would be revalidated on the compact cases.
   * **(a) Recommended.** Keep every sub-layer at or above a minimum optical
     depth: N = min(40, floor(τ_layer / δ)), with δ ≈ 2.5×10⁻³ (the
     non-scattering error is then ≤ 10⁻⁵).  When even one layer is thinner
     than δ, use a single layer at the window Planck-mean temperature of the
     profile; this is the optically thin limit, exact to first order in τ.
   * **(b)** Use only the thin-limit scaling below a threshold τ.
   * **(c)** Give no profile below COT 0.0625 (isothermal).  This creates a
     step in E_md at 0.0625.
2. **Thick ice discretisation (optional).** Keep the historical 40 equal-τ
   sub-layers, at about 0.2 K for thick ice at 3.7 µm, or grade them toward the
   top (≤ 5×10⁻⁴), respecting δ.

## 10. Runtime

* **Compact:** an emission call costs about 33 ms with 40 sub-layers against
  3.2 ms (V24).  The ice compact run takes 137 s against 131 s.
* **Production:** from the V24 logs, 9 216 (pold, p240) and 11 136 (pagg)
  emission calls, which adds about 0.1 h.
* **Full V25 estimates** from the measured V24 SLURM times: pold ≈ 31.5 h,
  p240 ≈ 33.3 h, pagg ≈ 18.9 h.

## 11. Remaining scientific limitations

* E_md depends on the 250 K reference top through the nonlinearity of B.
* There is one ice profile (cirrostratus) for all ice habits.  The 8 K/km
  lapse rate is a constant; it is not measured.
* The liquid adiabat's depth scale beta = 20 km^-1 is prescribed, and its
  deep states are unphysical (section 3).
* The interpretive choices without the supplier's DISORT example are the
  midpoint rule, H linear in COT, and H held below COT 0.0625.

## 12. Production

Not prepared: it waits for the decision of section 9.  The scripts are
`scripts/submit_v25_release.sh`, `run_v25_release.slurm` and
`monitor_v25_release.sh`.  They refuse to run until the release files are
committed and `releases/lut_v25.release` verifies.
