# V23 Terra MODIS validation

* `smoke.py`: compact pre-production smoke test and benchmark through the
  release launchers (V23 on Grid B angles, V22 on the coarse angles; p240 and
  pagg, channels 1 and 20, optical depths 1 and 8, radii 5 and 9 um); results
  in `smoke_results.json`.
* `compare_with_v22.py`: validation of the complete V23 products (structure,
  coordinates against the Grid B definitions, NaN/Inf, fill as in V22,
  provenance) and comparison with the definitive V22 products at exactly
  shared nodes; results in `results/`.

## Intensities at shared nodes depend on the number of user angles

Compact smoke test (2026-10-12): at shared nodes the flux quantities R_dd,
T_dd, R_0d, T_0d, T_00 and T_dv are bitwise equal to V22, while R_0v, R_dv
(up to 7e-6) and E_md (up to 2e-7) differ.  The cause was isolated on
captured V21 IDL DISORT calls (`tests/data/v21_idl_single/`), passed through
the production Intel backend with identical inputs and only the user angles
changed: the upwelling nadir intensity UU at azimuths 0/90/180 changes when the
user polar angles are those of Grid B instead of V22 (diffuse call: no change;
direct-beam calls: 6e-4, and 0.41 at the ill-conditioned tau 256 / r_e 29 /
SZA 80 call; emission: 6e-11), whereas changing only the azimuths, or
comparing the fluxes, gives no difference.  The intensity arithmetic of the
original single-precision DISORT kernel depends on NUMU (vectorised loops over
user angles; the alignment study of the V21/V22 work showed the same
sensitivity).  The V23 scientific code is the V22 code.
