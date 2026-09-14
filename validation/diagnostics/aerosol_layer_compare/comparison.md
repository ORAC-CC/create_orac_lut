# Aerosol channel-1 layer-input diagnostic (r_eff = 0.01 µm) — after the profile-interpolation correction

State: Meteosat-10 SEVIRI ch 1, `aerosol_a79.mm`, r_eff 0.01 µm, τ 1.0, 950 hPa, atmosphere 2, gas on, Rayleigh on.
The pre-correction diagnostic (which located the divergence) is preserved in `before_profile_fix/`.

## Correction applied

`_interpolate_profile` now reproduces IDL `INTERPOL`: linear inside the tabulated profile and
linear extrapolation from the two nearest end nodes on both sides (segment index clamped to the
end segments, `lib/interpol.pro` line 206), single precision, ascending or descending nodes.

## Post-correction array comparison (Python − legacy)

| Array | max abs | note |
| --- | --- | --- |
| `relative_tau_raw`, `relative_tau` | 0 | **bitwise identical** (0.5 km layer: 1085.23 both) |
| `tauscat`, `dtau` | 3.0e-8 | float32 rounding |
| `ssa` | 3.3e-7 | float32 rounding |
| `pmom` | 8.4e-9 | float32 rounding |
| gas, Rayleigh, levels, molecular/particle moments | as before: identical | |

First material difference: **none** (`first_material_difference: null`).

Lowest three normalised layer weights, legacy / Python: 0.20218 / 0.20218, 0.33333 / 0.33333, 0.46449 / 0.46449.

## Controlled DISORT check

| Inputs | R_dd |
| --- | --- |
| Python captured arrays (post-correction) | 0.0381610 |
| legacy captured arrays | 0.0381495 |
| legacy LUT | 0.038124 |

The remaining 1e-5–4e-5 is the kernel-process-level residual also seen in the cloud case.
DISORT inputs are now effectively identical: **YES**.
