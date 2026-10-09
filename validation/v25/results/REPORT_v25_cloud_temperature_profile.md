# ORAC cloud LUT Version 25: vertically varying cloud temperature in the thermal emission

Written for the owner of the ORAC LUT generator and for reviewers of the V25
production LUTs. Dates are 2026-10-06 unless stated. This revision records
the definitive V25 formulation: the supplied cirrostratus profile for ice
and the saturated liquid-water (wet) adiabat for liquid water.

## 1. Executive summary

- **Definitive formulation (user decisions).** Two backends of one
  mechanism — a vertical distribution of the cloud optical depth paired with
  a temperature at every emission-layer boundary, at the reference cloud-top
  temperature 240 K, `E_md` normalised by B(240 K), reflection and
  transmission unchanged:
  - **ice LUTs** (`cloud_vertical_profile = 'cirrostratus'`): the supplied
    vertically inhomogeneous cirrostratus profile of P. Watts (OCA /
    EUMETSAT), with its distribution of optical depth with depth below the
    cloud top, ΔT = 8 K/km × depth, and a cloud depth H(COT) of 2 km (COT
    0.0625) to 11 km (COT ≥ 128), never deeper;
  - **liquid-water LUTs** (`cloud_vertical_profile = 'wet_adiabat'`): the
    saturated liquid-water adiabat followed continuously from the reference
    state (240 K, 628 hPa) over the path length dz = dτ₀.₅₅ / β with
    β_ext,055 = 20 km⁻¹, i.e. z = τ₀.₅₅ / 20 km⁻¹, with no H, no H_max, no
    altitude or surface constraint, no freezing transition and no ceiling;
    the pressure is evolved hydrostatically along the path from the
    reference state.
  - The Grid B optical-depth endpoint τ₀.₅₅ = 256 is retained (accepted by
    the user); no larger grid or asymptotic study was made.
- **History.** (1) A constant extinction coefficient with a depth cap H_max
  and a saturated adiabat for both phases (`7f6b4e2`; jobs 488174–488213
  cancelled by the user) — superseded. (2) The uncapped z = τ/β adiabat for
  both phases (`13b9df3`): its ice branch gave 1630 K reference states for
  thick ice and was superseded by the supplied profile (`5722280`); its
  liquid branch is the definitive liquid treatment, restored bitwise.
- **What changes in the LUT.** Only `E_md` of the thermal and mixed
  channels. Reflection and transmission operators are bitwise V24 (§13).
- **Validation.** 26 focused tests pass; the repository suite gives 317 passed, 14 skipped, 5 failed, the five failures being the pre-existing ones (missing local files: `docs/idl_python_port_status.md`, `validation/reference_cases/`, `validation/diagnostics/`), identical to every earlier run; no new failure. On the compact cases
  (four ice, two liquid): the isothermal path is bitwise identical to the
  pre-V25 generator on the established compact grids; only `E_md` changes;
  the ice results are bitwise identical to the `5722280` runs; `E_md` is
  finite, smooth and bounded (ice ≤ 1.46 in the window channels and ≤ 4.3
  at 3.7 µm, peaking at τ ≈ 16; liquid ≤ 1.05 / 1.07); the emission layering
  converges (40 layers within 1.4×10⁻³ (ice) and 3×10⁻⁵ (liquid, window) of
  46); TOA effects are bounded and monotonic in τ (ice up to +16–18 K at
  11 µm for τ = 16; liquid ≤ +2.4 K at τ = 256).
- **Liquid extreme.** τ₀.₅₅ = 256 → z = 12.8 km; the reference adiabat
  reaches 317.49 K at 2957 hPa (recomputed; finite, monotonic, converged).
- **Reference-temperature dependence.** The normalised `E_md` depends on the
  240 K reference (7–8 % in the window channels and 24–43 % at 3.7 µm for
  thick ice between 220/260 K and 240 K, §8): an accepted limitation of V25.
- **Production.** the definitive formulation is validated, committed (`990d970`), pushed and
  passes the provenance gate; the submission of the 40 jobs (24 liquid-water,
  16 ice) was attempted with the established machinery and is blocked by an
  expired Kerberos ticket (the ssh hop to atmlxint7 needs it); a watcher
  submits the set automatically once a ticket is renewed (§18).

## 2. Repository state and the superseded V25 jobs

- Start of this stage: branch `main`, HEAD = `origin/main` = `5722280`, with
  the pre-existing uncommitted ` M AGENTS.md`, ` M tests/test_baum.py`,
  `?? documents/`, `?? scripts/submit_v23_slstr_cloud_luts.sh`,
  `?? tests/test_nakajima_king.py` and the 40 `?? runs/*_v23.run`, all
  untouched.
- V24 is preserved exactly: the V24 run files, scripts and products are
  unchanged, and the legacy isothermal path reproduces the unmodified
  generator bitwise (§13).
- Jobs **488174–488213** (submitted 2026-10-06 07:36 UTC from `7f6b4e2`)
  were cancelled by the user at 08:12:40 UTC in their scattering phase; no
  product was written. They are superseded and not authoritative; their rows
  stay in `validation/v25/production_submissions.tsv` and their logs in
  `validation/slurm/v25_*_4881xx_*.out` for provenance.

## 3. The supplied cirrostratus profile (ice)

**File.** `references/data/ocalut_cloudprofile_Cirrostratus.dat` (copied
unchanged from `documents/`; md5 `7d7b472f3ff9681901c03afe85c3202b`, 1317
lines, 92819 bytes). Header:

    13  # COTs
    99  # layers
    Origin: /tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf
    Cloud type: Cirrostratus

supplied by Alessio in connection with Phil Watts' vertically inhomogeneous
LUT tests for OCA / EUMETSAT. Each of the 13 blocks starts with
`<COT> LUT Cot value` and the column header
`% Extinction: layer dCOT: Cumul COT: Km-below CT: Temp above CT:`, followed
by 99 rows.

**Structure, established numerically** (`load_cloud_profiles` validates all
of this on every load; `figures/profile_supplied.png`):

| Property | Result |
|---|---|
| COTs | 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256 (13) |
| Rows per profile | 99 in every block; 13 × 99 × 5 values, all finite |
| `% Extinction` | sums to 100.000 in every profile; **identical in all 13 profiles** (max difference 0): one extinction shape f_i, minimum 0.26 % (top), maximum 1.22 % (row 75, z/H = 0.77), 0.81 % at the base |
| `layer dCOT` | = COT × f_i / 100 (relative deviation ≤ 2.7×10⁻³, the file's six-decimal rounding); non-negative; sums to the stated COT (relative error ≤ 3.2×10⁻⁵) |
| `Cumul COT` | the running sum of `layer dCOT` (≤ 5×10⁻⁶ absolute), strictly increasing, ending at the stated COT (≤ 3×10⁻⁵ relative) |
| `Km-below CT` | 0 in row 1, strictly increasing, uniformly spaced: z_i = H(COT) i/98 (max deviation 4×10⁻⁷ of H) |
| `Temp above CT` | 0 in row 1, strictly increasing; **ΔT = 8 K/km × z in every row of every profile** (max deviation 4×10⁻⁶ K) |

Supplied cloud-base values (all reproduced by the parser and the model):

| COT | depth below CT (km) | ΔT at base (K) |
|---|---|---|
| 0.0625 | 2.0000 | 16.0000 |
| 0.125 | 2.0774 | 16.6194 |
| 0.25 | 2.2323 | 17.8581 |
| 0.5 | 2.5419 | 20.3355 |
| 1 | 3.1613 | 25.2903 |
| 2 | 4.4000 | 35.2000 |
| 4 | 5.2000 | 41.6000 |
| 8 | 7.0000 | 56.0000 |
| 16 | 9.7200 | 77.7600 |
| 32 | 10.8067 | 86.4533 |
| 64 | 10.9100 | 87.2800 |
| 128 | 11.0000 | 88.0000 |
| 256 | 11.0000 | 88.0000 |

The cloud does not become 256 km deep at COT 256: H(COT) grows exactly
linearly in COT from 0.0625 to 1 (1.2387 km per unit COT), then 4.4, 5.2,
7.0, 9.72, 10.81, 10.91 km and 11.0 km at COT 128 and 256; above COT ≈ 32 the
depth is essentially constant and the extinction per unit depth grows.

**Other profiles / Phil's example.** The repository holds no other cloud-type
profile and no example of Phil's DISORT implementation (searched for Phil,
Watts, cirrostratus, cloudprofiles, METimCTP, OCA, EUMETSAT, Carbajal,
Henken, "vertically inhomogeneous", "inhomogeneous LUT", dCOT, "Temp above
CT"; only the MODTRAN and DISORT documentation matched). The interpretative
choices below are made from the data and documented.

## 4. Mathematical interpretation of the supplied model

With f_i the common extinction shape (Σ f_i = 100), for a cloud of total
0.55 µm optical depth COT:

    dτ_i(COT) = COT f_i / 100,            i = 0 … 98
    z_i(COT)  = H(COT) i / 98
    ΔT_i      = 8 K/km × z_i(COT)

so the only COT-dependent quantity is the cloud depth H(COT). The 99 rows
span the cloud from z = 0 (top) to z = H (base) inclusive and each carries an
optical depth dτ_i: they are sampled profile points with optical-depth
weights. The relation between cumulative optical depth and temperature is
taken at the row midpoints, F_i = (Σ_{j ≤ i} f_j − f_i/2)/100, with F = 0 at
the top and F = 1 at the base, which gives a COT-independent relation ζ(F)
between the cumulative optical-depth fraction and the depth fraction, and

    ΔT(F; COT) = 8 K/km × H(COT) × ζ(F).

The alternative readings (row i = the bottom or the top of layer i) differ
by half a row, at most ΔT_base / 196 ≈ 0.45 K locally.

## 5. Interpolation to the ORAC optical-depth grid (ice)

Grid B (`ice-cloud-grid-b.lut`): τ = 10⁻¹⁰, 2⁻⁷ … 2⁻¹, 1, 2, 2√2, …, 256 (24
nodes), maximum 256 = the largest supplied COT, so no extrapolation above the
data is needed. f_i and ζ(F) are COT-independent. H(COT) is interpolated
**linearly in COT** between the 13 supplied values (the data are exactly
linear in COT from 0.0625 to 1 and nearly so to 8; linear-in-log(COT)
interpolation differs by ≤ 0.23 km / 1.9 K at the intermediate Grid B nodes,
`profile_grid_b.csv`), held at 2 km below COT 0.0625; every supplied profile
is reproduced exactly at its own COT. As τ → 0 the emission vanishes with τ.

## 6. Reduction to the DISORT layer representation (both backends)

The production DISORT accepts MXCLY = 46 layers; the emission call holds the
in-cloud layers only (one atmosphere layer in every production microphysics
file). Each in-cloud layer is divided into N = 40 equal-optical-depth
sub-layers (boundaries F_k = k/N) with the layer's single-scattering albedo
and phase moments and optical depth dtau/N (the channel's particle optical
depth, scaled by the spectral extinction ratio exactly as V24); the boundary
temperatures are T_k = 240 K + ΔT(F_k; τ₀.₅₅) from the backend; DISORT's
Planck source is linear in optical depth within each sub-layer. The total
optical depth, SSA and phase moments are preserved exactly; only the thermal
source distribution changes.

**Layer count** (`compact_layers.csv`; maximum |ΔE_md| against 46 layers):

| Layers | Ice 11/12 µm | Ice 3.7 µm | Liquid 11/12 µm | Liquid 3.7 µm |
|---|---|---|---|---|
| 12 | 2.4×10⁻² (≈ 0.9 K) | 4.1×10⁻² | 1.0×10⁻³ | 1.9×10⁻² |
| 24 | 8.3×10⁻³ (≈ 0.3 K) | 1.6×10⁻² | 2.6×10⁻⁴ | 4.6×10⁻³ |
| 36 | 2.5×10⁻³ (≈ 0.09 K) | 5.1×10⁻³ | 6×10⁻⁵ | 1.2×10⁻³ |
| **40** | **1.4×10⁻³ (≈ 0.05 K)** | **2.8×10⁻³** | **3×10⁻⁵** | **6×10⁻⁴** |

The ice ladder converges at about second order; the residual at 40 layers
is ≤ 0.1–0.2 K at 11 µm (below the 0.52 K production uncertainty). The
liquid ladder converges faster (the temperature span is 1–78 K over the
cloud) and is converged at any count ≥ 24. One common production value
**EMISSION_LAYERS = 40** satisfies both and stays within MXCLY = 46.

## 7. The liquid-water wet adiabat

For a liquid-water LUT the backend `WetAdiabatProfile` tabulates once per
LUT, from z = 0 to z_max = max(τ₀.₅₅)/β, the coupled system

    dT/dz = Γ_w(T, p),    dp/dz = p g / (R_d T)   (hydrostatic, dry air)

from the reference state (240 K, 628 hPa) by fourth-order Runge–Kutta in
10 m steps, with the saturated (pseudo-adiabatic) lapse rate over liquid
water

    Γ_w = g (1 + L_v r_s / (R_d T)) / (c_pd + L_v² r_s ε / (R_d T²)),   r_s = ε e_sw / (p − e_sw),

e_sw of Murphy and Koop (2005, Eq. 10, valid for supercooled water), L_v of
Rogers and Yau (1989), g = 9.80665 m s⁻², R_d = 287.04, c_pd = 1005.7
J kg⁻¹ K⁻¹, ε = 0.622. The cumulative optical-depth fraction F of a node
maps to z = F τ₀.₅₅ / 20 km⁻¹ and T(z) is interpolated from the table. There
is no H, no H_max, no cap, no altitude, no surface, no freezing transition
and no ceiling; the liquid adiabat is followed through and above 273.15 K.
The reference pressure 628 hPa is a convention of the reference state (the
value the earlier V25 formulation used); it enters the lapse rate only
through r_s (4 % of the lapse rate at 240 K). The restored backend is
bitwise identical to the validated liquid branch of `13b9df3`.

**Production grid** (`liquid-water-cloud-grid-b.lut`, τ_max = 256): the
adiabat is finite and strictly monotonic in T and p; a 1 m step changes T at
12.8 km by 1.5×10⁻¹² K. Representative nodes:

| τ₀.₅₅ | z (km) | T (K) | p (hPa) | 40-layer max ΔT step (K) |
|---|---|---|---|---|
| 0.0625 | 0.0031 | 240.03 | 628.3 | 0.00 |
| 1 | 0.05 | 240.45 | 632.5 | 0.01 |
| 4 | 0.2 | 241.80 | 646.1 | 0.05 |
| 16 | 0.8 | 247.09 | 702.6 | 0.18 |
| 64 | 3.2 | 266.22 | 966.8 | 0.72 |
| 128 | 6.4 | 286.60 | 1435.3 | 1.44 |
| **256** | **12.8** | **317.49** | **2957.2** | **2.87** |

## 8. Reference-temperature approximation

Both backends give a temperature departure from the cloud top; V25
generates `E_md` at the reference T_top = 240 K and ORAC multiplies by
B(T_retrieved). Because B is nonlinear, the normalised `E_md` depends on the
reference temperature. The diagnostic made in the previous stage for ice
(`compact_ttop.csv`: at τ = 16, |ΔE_md| = 0.10–0.11 / 0.07–0.08 (7.5–7.8 % /
5.4–5.6 %) at 11 µm, 0.08 / 0.06 (6.5 % / 4.8 %) at 12 µm and 1.6–1.9 /
0.9–1.0 (42–43 % / 23–24 %) at 3.7 µm for T_top = 220 / 260 K against 240 K)
is **accepted as a limitation of V25** by the user: 240 K is retained for
both treatments, no cloud-temperature dimension is added, and the
investigation is not repeated.

## 9. Phase coverage and the run-file selector

`cloud_vertical_profile = 'cirrostratus'` in the 16 ice run files (4 models ×
Aqua, Terra, Sentinel-3A, Sentinel-3B) and `'wet_adiabat'` in the 24
liquid-water run files (6 models × 4 platforms); every other setting is the
V24 one (asserted by a test). Absent or `'isothermal'` reproduces V24. The
generator refuses a backend for the other phase (cirrostratus + liquid
water, wet_adiabat + water ice: `PROFILE_SUBSTANCES`).

## 10. Code and API

`src/oraclut/cloud_temperature.py`:

| Component | Description |
|---|---|
| `load_cloud_profiles(path)` → `CloudProfileSet` | parser and validator of the supplied file (every check of §3); depth_at (H linear in COT), depth_fraction ζ(F), delta_t_at, layer_boundaries |
| `CloudVerticalProfile` | the ice model of one LUT: profile set + T_top; `boundary_temperatures_k`, `cloud_base_temperature_k`, `describe`, `global_attributes` |
| `saturation_vapour_pressure_liquid`, `latent_heat_vaporisation`, `wet_adiabatic_lapse_rate`, `wet_adiabat` | the liquid thermodynamics (restored from `13b9df3`, liquid part only) |
| `WetAdiabatProfile` | the liquid model of one LUT: tabulated adiabat, z = x/β; the same interface as `CloudVerticalProfile` |
| `cloud_vertical_profile_model(name, max_tau_055, t_top_k=240)` | dispatch: `'cirrostratus'` (PROFILE_FILES) or `'wet_adiabat'` |
| `emission_layers(dtau, ssalb, pmom, scatreltau, tau_055, model, nlayers)` | the sub-layered DISORT inputs for either model (refuses more than 46 layers) |
| Constants | `T_TOP_K = 240`, `LAPSE_K_PER_KM = 8` (verified), `BETA_EXT_055_LIQUID_PER_KM = 20`, `P_TOP_HPA = 628`, `EMISSION_LAYERS = 40`, `DISORT_MAX_LAYERS = 46`, `PROFILE_SUBSTANCES` |

Removed / absent from the production path: H_max and any bounded-depth
remapping, β_ice and the ice adiabat, any absolute cloud altitude or
atmosphere lookup, and any freezing / melting logic (a test asserts this).
`create_orac_luts.py`: `cloud_vertical_profile` ∈ {isothermal, cirrostratus,
wet_adiabat}; the model is built once per LUT from the setting and
`max(lutstr.opd)`; the phase pairing is enforced; the V25 global attributes
record the treatment. `scripts/submit_v25_cloud_luts.sh` checks the setting
per family and takes an optional `ice | liquid | both` selector.

## 11. Tests

`tests/test_cloud_temperature.py`, 26 tests, all pass. Ice (as in the
previous stage): file provenance and shape; layer sums and cumulative;
common normalised extinction; equally spaced depth and ΔT = 8 z; supplied
base values; native profiles reproduced at their COT; parser rejection of
corrupted files; Grid B interpolation (no extrapolation, monotonic, top
ΔT = 0, held below 0.0625); linear-in-COT interpolation; τ → 0; layer
conservation and the MXCLY guard; one profile for all channels; absence of
β / cap / altitude / ice-adiabat / melting logic in the ice path. Liquid
(new): reference state and constants (240 K, 628 hPa, β = 20 km⁻¹);
z = τ/β exactly (10 → 0.5 km) with no saturation; thermodynamics (Murphy–Koop
and Rogers–Yau values, lapse-rate limits, every tabulated step equals the
liquid lapse rate, hydrostatic pressure in the adiabat's own temperature);
finite and monotonic over the production grid with 300 K < T_base(256) <
330 K; τ → 0; layer conservation and channel independence; the phase
pairing. End to end (slow): isothermal = default bitwise and cirrostratus
changes only `E_md` (ice test grid); wet adiabat changes only `E_md` (liquid
test grid), with the V25 attributes and no depth-cap attribute. Run files:
40, with 16 ice (cirrostratus) and 24 liquid (wet_adiabat).

Repository suite: the repository suite gives 317 passed, 14 skipped, 5 failed, the five failures being the pre-existing ones (missing local files: `docs/idl_python_port_status.md`, `validation/reference_cases/`, `validation/diagnostics/`), identical to every earlier run; no new failure.

## 12. Compact radiative-transfer validation

Six cases on the Grid B subsets `validation/v25/compact_v25_ice.lut` and
`compact_v25_liquid.lut` (τ = 0.0625, 0.25, 1, 4, 16, 64, 256; r_e = 5–93 µm
ice, 3–35 µm liquid; SZA/VZA 0, 30, 60, 75°; 7 azimuths), all channels of the
DISORT-options cases (MODIS 8, 3, 1, 2, 6, 7, 20, 31, 32; SLSTR 1, 3, 5, 6, 7,
8, 9), production kernel, NSTR 60 (`validation/v25/run_compact.py`,
`compare_compact.py`):

| Case | Instrument | Microphysics | V25 backend |
|---|---|---|---|
| modis_ice_agg | Terra MODIS | Baum severely rough aggregate | cirrostratus |
| modis_ice_sph | Terra MODIS | ice spheres (Mie) | cirrostratus |
| slstr_ice_sph | Sentinel-3A SLSTR | ice spheres | cirrostratus |
| slstr_ice_agg | Sentinel-3A SLSTR | Baum aggregate | cirrostratus |
| modis_liquid | Terra MODIS | liquid water stg | wet_adiabat |
| slstr_liquid | Sentinel-3A SLSTR | liquid water stg | wet_adiabat |

**E_md** (`compact_e_md.csv`; `figures/e_md_toa.png`, upper row): V25 − V24
at view zenith 0, worst over τ and r_e, with the V25 maximum:

| Case | 3.7 µm (mixed) | 11 µm | 12 µm |
|---|---|---|---|
| MODIS Baum aggregate | +3.48 at τ 16, r_e 5 (max 4.33) | +0.42 at τ 16 (1.42) | +0.35 at τ 16 (1.34) |
| MODIS ice spheres | +3.00 at τ 16 (3.82) | +0.38 at τ 4 (1.37) | +0.34 at τ 4 (1.32) |
| SLSTR ice spheres | +3.09 at τ 16 (3.95) | +0.42 at τ 16 (1.42) | +0.34 at τ 4 (1.32) |
| SLSTR Baum aggregate | +3.40 at τ 16 (4.27) | +0.47 at τ 16 (1.46) | +0.35 at τ 16 (1.34) |
| MODIS liquid water | +0.156 at τ 256, r_e 3 (1.07) | +0.052 at τ 256 (1.05) | +0.030 at τ 256 (1.02) |
| SLSTR liquid water | +0.154 at τ 256 (1.07) | +0.057 at τ 256 (1.05) | +0.030 at τ 256 (1.02) |

Ice: the change is positive everywhere, grows with τ while the cloud
deepens, peaks at τ ≈ 16 (H 9.7 km, ΔT_base 78 K) and falls for thicker
cloud (E_md 1.04 at τ = 256), as the supplied profile prescribes; the
superseded constant-β adiabat rose without bound to 15.7. Liquid: the change
is positive, grows with τ and saturates from τ ≈ 16 (0.896 → 1.006 at 3.7 µm,
0.999 → 1.019 at 11 µm for r_e 15); the deep part of a thick liquid cloud
(317 K at 12.8 km for τ = 256) contributes little because the emission comes
from the top few optical depths. Both are finite and smooth across τ, r_e and
view zenith. The isothermal runs are bitwise identical to the pre-V25
generator (`51717b8`) on the established compact grids (§13) and the ice V25
runs are bitwise identical to the `5722280` runs.

## 13. Reflection, transmission and the Nakajima–King control

`compact_operator_changes.csv`: V25 against isothermal, every reflection and
transmission operator (R_0v, R_0d, R_dd, R_dv, T_00, T_0d, T_dd, T_dv) and
every optical property are bitwise identical in all six cases (33–35
variables); `E_md` alone changes. Intended: both backends prescribe
temperature structure, not microphysics, every sub-layer has the LUT
particle model's single-scattering albedo and phase function, and in
plane-parallel radiative transfer a homogeneous layer's reflection and
transmission depend on its optical depth, SSA and phase function only; the
redistribution enters the emission call alone, and the V24 diffuse and
direct-beam calls on the V24 atmosphere are kept.

`compact_v24_reproduction.csv`: the isothermal runs match the V24 products
at the matching Grid B nodes bitwise (Baum cases) or within the known
host-dependent single-precision floor (≤ 1.3×10⁻⁴; MODIS liquid and ice
spheres, SLSTR ice spheres); the definitive code's isothermal default is
bitwise identical to the pre-V25 generator's output on the established
compact grids (modis_liquid, modis_ice_agg; `validation/tmp/v25/legacy_check`).
The reflectance Nakajima–King control (`figures/nk_r0v_*.png`) shows
identical grids (max |ΔR| = 0).

## 14. Top-of-atmosphere brightness temperature

`compact_toa.csv`; `figures/e_md_toa.png`, lower row; `figures/nk_bt.png`.
The ORAC thermal forward model of `fm_thermal.F90` (McGarragh et al. 2018
Eq. 47; `validation/v24/disort_options/orac_fm.py`, `atmosphere.py`:
mid-latitude-summer clear-sky terms, surface emissivity 0.98, cloud top
9 km (ice) / 4 km (liquid); B(T_c) at the reference 240 K; mixed channels
with the daytime solar term at SZA 30°, RAA 90°, black surface, as `fm.F90`):

| Case | 3.7 µm (mixed) ΔBT max | 11 µm ΔBT max | 12 µm ΔBT max | r = ΔBT/σ (11 µm): max, median |
|---|---|---|---|---|
| MODIS Baum aggregate | +14.2 K (τ 16, r_e 93) | +16.5 K (τ 16, r_e 5) | +15.2 K | 32, 6.6 |
| MODIS ice spheres | +13.3 K | +14.8 K (τ 16, r_e 5) | +14.3 K | 28, 6.7 |
| SLSTR ice spheres | +12.7 K | +16.1 K (τ 16, r_e 5) | +14.2 K | 31, 6.8 |
| SLSTR Baum aggregate | +13.8 K | +17.7 K (τ 16, r_e 5) | +15.1 K | 34, 6.4 |
| MODIS liquid water | +0.32 K (τ 256, r_e 35) | +2.20 K (τ 256, r_e 3) | +1.36 K | 4.2, 0.43 |
| SLSTR liquid water | +0.29 K | +2.36 K (τ 256, r_e 3) | +1.34 K | 4.5, 0.44 |

σ is the production ORAC measurement uncertainty (0.52 K at 11 µm). Ice:
along τ (MODIS Baum, r_e 33, VZA 0, 11 µm) +0.15, +0.65, +3.4, +13.3, +15.3,
+6.3, +2.0 K for τ = 0.0625 … 256; against view zenith at τ = 16: +15.3,
+13.8, +9.3, +5.8 K at 0, 30, 60, 75°. Liquid (MODIS, r_e 15): +0.00, +0.00,
+0.05, +0.47, +0.87, +0.88, +0.88 K at 11 µm; at τ = 64: +0.88, +0.76, +0.43,
+0.22 K at 0, 30, 60, 75°. All differences are positive and smooth (one
−0.001 K at the float32 noise floor), the brightness temperature decreases
monotonically with τ in every case (`figures/nk_bt.png`: the dashed V25 τ
isolines keep their order; the split-window grids move along the diagonal),
so the LUTs remain retrieval-compatible. Retrieved cloud-top temperatures of
intermediate-τ ice cloud will be colder than with V24 by up to ≈ 15 K, of
thick liquid cloud by ≤ 2.4 K.

## 15. Runtime

`compact_timing.csv`, single-threaded, all channels: the emission call costs
3.3–4.1 s per 90 calls with 40 layers against 0.17–0.45 s with one (≈ 40 ms
instead of 3 ms per call), below 3 % of the DISORT time (the direct-beam
calls dominate; the totals of concurrent runs vary by more than that).
Building the models costs 6 ms (profile file) and 50 ms (adiabat to 12.8 km)
once per LUT. For a production job (≈ 10⁴ emission calls) the addition is of
order 400 s on a 22–42 h job.

## 16. Figures

- `figures/profile_supplied.png`: the supplied extinction profile; cumulative
  optical depth against normalised depth; H against COT (supplied points,
  the linear-in-COT interpolation used and the log alternative); base ΔT
  against COT; ΔT against depth; the remapped 40-layer boundaries.
- `figures/e_md_toa.png`: E_md(V25) − E_md(V24) and ΔBT against τ by case and
  channel (six cases; the superseded constant-β result dotted for ice).
- `figures/layers_ttop.png`: layer-count convergence (ice); the T_top
  diagnostic of §8.
- `figures/nk_bt.png`: brightness-temperature Nakajima–King diagrams (11 vs
  12 µm; 3.7 vs 11 µm), V24 solid, V25 dashed, six cases.
- `figures/nk_r0v_*.png`: reflectance Nakajima–King control (identical).

## 17. Git

- Start: `main` = `origin/main` = `5722280`.
- Definitive commit: **`990d970`** "Add the saturated liquid-water adiabat as
  the V25 liquid-water backend" (`create_orac_luts.py`,
  `src/oraclut/cloud_temperature.py`, `tests/test_cloud_temperature.py`, the
  24 liquid-water `runs/*_v25.run`, `runs/V25_PRODUCTION.md`,
  `runs/template.run`, `scripts/submit_v25_cloud_luts.sh`), on top of
  `5722280` (the cirrostratus ice implementation and the profile file).
  Pushed; `git fetch origin main` then `git merge-base --is-ancestor HEAD
  origin/main` succeeds (HEAD = origin/main = `990d970`);
  `scripts/submit_v25_cloud_luts.sh --dry-run` reports "clean production
  paths; contained in origin/main" and 40 products to submit.
- Not committed (pre-existing, untouched): ` M AGENTS.md`,
  ` M tests/test_baum.py`, `?? documents/`, `?? scripts/submit_v23_slstr_cloud_luts.sh`,
  `?? tests/test_nakajima_king.py`, the 40 `?? runs/*_v23.run`. The validation
  area (`validation/`) is git-ignored by repository policy.
- V24 run files, scripts and products and the IDL sources are unchanged.

## 18. Production submission and status

**Attempt.** `scripts/submit_v25_cloud_luts.sh all` was run at 20:4x UTC from
revision `990d970` after the gate passed. The wrapper's ssh hop to
atmlxint7 failed (exit 255: `Permission denied (publickey,keyboard-interactive)`):
the user's Kerberos ticket had expired at 17:58 BST, and without it sshd on
atmlxint7 cannot read `~/.ssh/authorized_keys` on the Kerberos-protected
home, so the agent's keys are refused (this morning's hops succeeded right
after a `kinit`). Nothing was submitted; the manifest has no rows at
`990d970`.

**Automatic submission.** `validation/v25/autosubmit_when_authenticated.sh`
is running (`nohup`, polling `klist -s` every two minutes for up to seven
days; log `validation/tmp/v25/logs/autosubmit.log`). When a valid ticket
appears it checks HEAD = `990d970`, confirms the ssh hop, runs
`scripts/submit_v25_cloud_luts.sh all` once, counts the manifest rows and
runs `validation/v25/monitor_jobs.sh`. To cancel it:
`pkill -f autosubmit_when_authenticated`. The manual alternative is
`kinit` followed by `scripts/submit_v25_cloud_luts.sh all`.

**Expected set:** 40 jobs — MODIS Aqua (6 liquid-water + 4 ice), Terra
(6 + 4), SLSTR Sentinel-3A (6 + 4), Sentinel-3B (6 + 4) — each recorded in
`validation/v25/production_submissions.tsv` with time, host, revision, job
ID, job name, run file and product. After submission verify: 40 valid job
IDs, no immediate failure, every log at revision `990d970` with source state
clean, the configuration lines `Cloud profile: wet_adiabat` (24 jobs) and
`Cloud profile: cirrostratus` (16 jobs); then `validation/v25/check_products.py`
when products appear. Jobs 488174–488213 remain the cancelled, superseded
pre-final runs (§2).

**Post-submission verification:** to be recorded here (job IDs, immediate SLURM states, log revision and `Cloud profile` lines, 24 + 16 count) once the watcher or a manual `scripts/submit_v25_cloud_luts.sh all` has submitted the set; `validation/tmp/v25/logs/autosubmit.log` holds the watcher's record.

## 19. Remaining limitations

- The 240 K reference temperature: the normalised `E_md` carries the
  temperature contrast of a 240 K top (§8); accepted for V25.
- The ice profile is one cirrostratus case; the liquid treatment is a
  reference adiabat with a prescribed β = 20 km⁻¹ (a 317 K, 3 bar reference
  state at τ = 256 contributes little to the emission but is unphysical as a
  cloud).
- Interpretative choices made without Phil's example: rows as optical-depth-
  weighted samples (midpoint rule), linear-in-COT interpolation of H,
  equal-optical-depth DISORT layers (§4–§6).
- The infrared asymptote beyond τ = 256 was not studied; the Grid B endpoint
  is accepted. A larger-τ check is a possible future refinement only.
- `python -m oraclut.generate` (`src/oraclut/pipeline.py`) keeps the V24
  isothermal emission; V25 is implemented in the production path only.
- No compiled ORAC retrieval was available on this system; the reader check
  is the static trace (`read_sad_lut.F90`, `ncdf_read_template.inc`: no range
  clipping of `E_md`) and the repository's Python reader.
