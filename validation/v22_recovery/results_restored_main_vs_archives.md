# Restored main tree (lut-code-v22.1) against the archived V22 and IDL V21 campaign products

Regenerated products: `validation/tmp/v22_restored/out`.  Archives: `/network/scratch/grainger/quota-relief/project-oraclut-bloated-backup` (read-only).

| case | instrument | model | qm | regen vs archived V22 | regen vs IDL V21: worst abs (family) | archived V22 vs IDL V21: worst abs | envelopes | attributes vs archived V22 |
|---|---|---|---|---|---|---|---|---|
| 000 | aqua_modis | aerosol aerosol_a79.mm | 1 | DIFFERS (51/52 variables bitwise) | 0.00293 (metadata), 2.38e-07 (optics), 8.81e-05 (operators) | 0.00293, 2.38e-07, 8.81e-05 | within; same as archive | names/values identical; 87 NC_CHAR->NC_STRING |
| 004 | earthcare_msi | aerosol aerosol_a79.mm | 1 | IDENTICAL (52/52 variables bitwise) | 0.000671 (metadata), 1.79e-07 (optics), 8.72e-05 (operators) | 0.000671, 1.79e-07, 8.72e-05 | within; same as archive | names/values identical; 87 NC_CHAR->NC_STRING |
| 034 | meteosat-10_seviri | cloud liquid-water_stg.mm | 1 | IDENTICAL (51/51 variables bitwise) | 0.000977 (metadata), 0.000397 (optics), 0.000126 (operators) | 0.000977, 0.000397, 0.000126 | within; same as archive | names/values identical; 84 NC_CHAR->NC_STRING |
| 035 | meteosat-10_seviri | cloud liquid-water_stg.mm | 2 | IDENTICAL (51/51 variables bitwise) | 0.000977 (metadata), 4.58e-05 (optics), 0.000136 (operators) | 0.000977, 4.58e-05, 0.000136 | within; same as archive | names/values identical; 84 NC_CHAR->NC_STRING |
| 062 | noaa-20_viirs | cloud liquid-water_stg.mm | 1 | IDENTICAL (51/51 variables bitwise) | 0.000977 (metadata), 0.00562 (optics), 0.000192 (operators) | 0.000977, 0.00562, 0.000192 | within; same as archive | names/values identical; 84 NC_CHAR->NC_STRING |
| 086 | sentinel-3a_slstr | cloud liquid-water_stg.mm | 1 | IDENTICAL (53/53 variables bitwise) | 0.00488 (metadata), 0.000381 (optics), 0.000121 (operators) | 0.00488, 0.000381, 0.000121 | within; same as archive | names/values identical; 88 NC_CHAR->NC_STRING |

All regenerated products bitwise identical to the archived V22 products in every variable, with attributes identical in name and value: **False**.

Envelopes: metadata 0.0400390625, optics 0.0088958740234375, operators 0.0008418560028076172 (V22_FREEZE.md).  'other' variables (channel ids, flags, grids, strings) are compared for identity only.
