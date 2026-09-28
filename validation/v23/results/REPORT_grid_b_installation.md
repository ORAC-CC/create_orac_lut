# Grid B operational LUT installation

Grid B is now the recommended operational LUT sampling specification for liquid-water cloud, ice cloud, and generic PA76 aerosol. It is independent of ORAC software, processing, generator, and product version numbers. No RT or NetCDF generation was performed.

## Selected operational files

| family | installed file | counts (tau × r_eff × SOZ × SAZ × RAA) |
|---|---|---:|
| water | `create_orac_lut/input_files/lut/liquid-water-cloud-grid-b.lut` | 24 × 24 × 21 × 21 × 37 |
| ice | `create_orac_lut/input_files/lut/ice-cloud-grid-b.lut` | 24 × 29 × 21 × 21 × 37 |
| aerosol_pa76 | `create_orac_lut/input_files/lut/aerosol-pa76-grid-b.lut` | 29 × 27 × 16 × 21 × 13 |

## Source and exact grid

The files were populated directly from the Stage 18 JSON selections:
- `validation/v23/results/stage18_final_grid_water.json`
- `validation/v23/results/stage18_final_grid_ice.json`
- `validation/v23/results/stage18_final_grid_aerosol_pa76.json`

All three use explicit `uneven_linear` r_eff and SAZ arrays, `uneven_logarithmic` tau, family-specific SOZ, and parser-required linear RAA endpoints. SAZ is exactly 0–75° for every family.

## Validation

- Parser round-trip: **PASS** for tau, r_eff, SOZ, SAZ, and RAA.
- Representative instruments: **PASS** for EarthCARE 20°, AATSR 56°, VIIRS 57°, MODIS 65°, and geostationary SEVIRI advertising 90°.
- SETDIS/DISORT: **PASS** with NSTR=60, NMOM=1000, NUMU/NPHI within compiled limits.
- Generator preflight: **PASS** for the three validation runs under `validation/v23/grid_b/runs/`.

## Documentation

- `docs/grid_b_lut_definitions.md` identifies Grid B as the current recommended choice.
- `validation/v23/results/GRID_B.md` gives exact nodes and five-dimensional point counts.
- Stage 18 JSON and scientific reports remain the provenance/evidence layer.

## Git

Commit and push identifiers are recorded here after the controlled commit/push step.
