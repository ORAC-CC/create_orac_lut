# Aerosol channel-1 layer-input diagnostic (r_eff = 0.01 µm)

Diagnostic only — no scientific code changed. Date: 2026-09-14.

## State

Meteosat-10 SEVIRI channel 1 (0.63818 µm), `aerosol_a79.mm`, r_eff = 0.01 µm,
τ₅₅₀ = 1.0, p_surf = 950 hPa, atmosphere code 2 (`mls.atm`), gas on, Rayleigh on.
Chosen as the maximum relative `R_dd` discrepancy in the channel-1 aerosol
comparison (legacy 0.038124, Python 0.037788, −0.88 %). `R_dd` is
diffuse→diffuse, so SZA/VZA/RAA do not enter. Particle optics at this state
(identical in both codes): SSA 0.026926, b_ext ratio 0.845771, g 0.021705.

Both implementations build **45 layers** from 46 levels (100 → 0 km, reversed
MODTRAN order) and **1000 moments**; layer ordering (top → surface) and moment
ordering (0 → 999) agree.

## Strict array comparison (Python − legacy)

| Array | max abs | RMS | mean abs | differing > 1e-6 | note |
| --- | --- | --- | --- | --- | --- |
| level heights, pressures, layer heights | 0 | 0 | 0 | 0 | bitwise identical |
| profile heights / amounts (from `.mm`) | 0 | 0 | 0 | 0 | bitwise identical |
| gas levels, `tau_gas` | 0 | 0 | 0 | 0 | bitwise identical |
| `column_tau_rayleigh`, `rayleigh_level`, `tau_rayleigh` | 3.0e-8 / 7.5e-9 | — | — | 0 | float32 rounding only |
| molecular moments (GETMOM) | 0 | 0 | 0 | 0 | bitwise identical |
| particle moments (1000) | 1.9e-9 | 7.5e-11 | — | 0 | identical to float precision |
| **`relative_tau_raw`** | **3.064e+2** | 4.57e+1 | — | **1 (layer 44)** | **first material difference** |
| `relative_tau` (normalised) | 8.08e-2 | 1.49e-2 | — | 3 (layers 42–44) | consequence |
| `tauscat` | 6.84e-2 | 1.26e-2 | — | 3 | consequence |
| `dtau` | 6.84e-2 | 1.26e-2 | 3.04e-3 | 3 | max rel 0.17 at layer 44 |
| `ssa` | 3.58e-3 | 7.96e-4 | 2.03e-4 | 3 | max rel 0.074 at layer 42 |
| `pmom` | 9.95e-4 | 6.9e-6 | — | 12 | all in layers 42–44 |

`pmom` by moment: m=0 identical; m=1 max 9.9e-4 (layer 44: 0.013534 vs
0.012540); m=2 1.7e-5; m=10 1.5e-11; m=100, 500, 999 identical (zero).

## Where the difference enters

**Before vertical mixing**, in the interpolation of the aerosol relative-amount
profile onto layer mid-heights:

| layer | height | legacy weight | Python weight |
| --- | --- | --- | --- |
| 42 | 2.5 km | 472.37 | 472.37 |
| 43 | 1.5 km | 778.80 | 778.80 |
| 44 | 0.5 km | **1085.23** | **778.80** |

The `.mm` profile is defined at 1.5–5.5 km. Legacy uses IDL `INTERPOL`
(`create_orac_aerosol_lut.pro` line 342), which **extrapolates linearly**
below 1.5 km (slope from the two lowest nodes → 1085.23 at 0.5 km). Python's
`_interpolate_profile` uses `np.interp`, which **clamps** to the end value
(778.80). After normalisation the legacy column puts 46.4 % / 33.3 % / 20.2 %
of the aerosol in the 0.5 / 1.5 / 2.5 km layers; Python puts 38.4 % / 38.4 % /
23.3 %. Column totals are identical, which is why `T_dd` agrees to 1e-5.

The Rayleigh/gas combination, the SSA weighting and the moment mixing are
arithmetically identical: every downstream difference is inherited from the
three affected `tauscat` values.

## Physical decomposition (layers 42–44, legacy → Python)

| layer | aerosol scattering τ | aerosol absorption τ | Rayleigh τ | gas τ |
| --- | --- | --- | --- | --- |
| 42 (2.5 km) | 0.00460 → 0.00530 | 0.16640 → 0.19151 | 0.00519 | 0.00154 |
| 43 (1.5 km) | 0.00759 → 0.00874 | 0.27433 → 0.31574 | 0.00577 | 0.00228 |
| 44 (0.5 km) | 0.01058 → 0.00874 | 0.38227 → 0.31574 | 0.00639 | 0.00323 |

Rayleigh and gas terms are identical; only the aerosol terms move between
layers. The earlier hypothesis — a different Rayleigh/absorber *mixing rule* —
is **not supported**: the mixing is the same. What differs is *where* the
near-pure absorber (SSA 0.027) sits relative to the Rayleigh scattering, and
with Rayleigh dominating `R_dd` for this particle that placement matters; for
r_eff = 0.1 µm (SSA 0.85) the aerosol's own scattering dominates and the same
weight shift changes `R_dd` by only ~5e-4.

## DISORT-input equivalence

**NO.** First divergence: `relative_tau_raw`, layer 44 (0.5 km), before any
optical mixing.

## Controlled DISORT check (production kernel called from Python)

| Inputs | R_dd | T_dd |
| --- | --- | --- |
| Python captured arrays | 0.0377878 | 0.2499366 |
| legacy captured arrays | 0.0381495 | 0.2499354 |
| Python arrays with legacy profile weights substituted | 0.0381610 | 0.2499354 |
| legacy product (IDL end to end) | 0.038124 | — |
| Python product | 0.037788 | — |

Python-kernel-on-Python-inputs reproduces the Python product exactly; swapping
in only the legacy profile weights recovers the legacy value to 1e-4 (the
remaining 2.6e-5–3.7e-5 between 0.03815/0.03816 and 0.038124 is the same
kernel-process-level residual seen in the cloud case). The profile
interpolation convention therefore accounts for the whole channel-1 aerosol
discrepancy.

## Implication (not acted on here)

The difference is a boundary-extrapolation convention in
`_interpolate_profile` (`src/oraclut/pipeline.py`) versus IDL `INTERPOL`: to
reproduce legacy, the Python profile interpolation would need linear
extrapolation beyond the profile's end nodes rather than clamping. Whether the
legacy extrapolation is the *intended* physics (it invents aerosol below the
lowest defined node) is a scientific decision for the user; reproduction first,
then any deliberate change with its own validation. No production code was
modified in this diagnostic.

## Files

`aerosol_layer_capture.pro` (IDL, reuses legacy loaders + `legacy_scatfile_ch01.sav`),
`capture_python_state.py`, `compare_layer_inputs.py`, `legacy_state.sav`,
`legacy_state.npz`, `python_state.npz`, `comparison.json`, this file.
