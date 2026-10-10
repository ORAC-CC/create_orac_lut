# V22 Terra MODIS cloud LUT validation

## terra_modis_m_liquid-water_a01_pold_v22.nc

SLURM job 490632 (v22_terra_water_old); run file `runs/terra_modis_cloud_liquid-water_old_v22.run`; submitted from revision `4fee1b4`.

```
PASS  opens; V2 reader accepts it; format NETCDF4; size 41330203 bytes
PASS  dimensions {optical_depth: 17, effective_radius: 20, solar_zenith: 10, satellite_zenith: 10, relative_azimuth: 11, channels: 36}
PASS  optical_depth equals liquid-water-cloud.lut (17 values, 1e-10..256)
PASS  effective_radius equals liquid-water-cloud.lut (20 values, 1..39)
PASS  solar_zenith equals liquid-water-cloud.lut (10 values, 0..90)
PASS  satellite_zenith equals liquid-water-cloud.lut (10 values, 0..65)
PASS  relative_azimuth equals liquid-water-cloud.lut (11 values, 0..180)
PASS  channel_id 1..36 (36 channels)
PASS  solar_channel_flag matches the instrument file
PASS  thermal_channel_flag matches the instrument file
      solar channels 26, thermal 16, mixed 6
PASS  T_dv                           finite; range -2.555e-15..1; outside valid_range [0.0, 1.0]: 1041 of 122400
PASS  T_dd                           finite; range 0..1; outside valid_range [0.0, 1.0]: 222 of 12240
PASS  R_dv                           finite; range -5.34e-06..0.983; outside valid_range [0.0, 1.0]: 1714 of 122400
PASS  R_dd                           finite; range -7.296e-05..0.9787; outside valid_range [0.0, 1.0]: 130 of 12240
PASS  R_0v                           finite; range -0.0001018..1.22; outside valid_range [0.0, 1.0]: 107787 of 9724000
PASS  R_0d                           finite; range -0.0001813..0.9912; outside valid_range [0.0, 1.0]: 215 of 88400
PASS  T_0d                           finite; range -6.284e-05..0.8183; outside valid_range [0.0, 1.0]: 183 of 88400
PASS  T_00                           finite; range 0..1; outside valid_range [0.0, 1.0]: 0 of 88400
PASS  E_md                           finite; range -9.112e-11..0.9989; outside valid_range [0.0, 1.0]: 939 of 54400
PASS  extinction_coefficient         finite; range 0.2044..7434; outside valid_range [0.0, 3.4028234663852886e+38]: 0 of 720
PASS  extinction_coefficient_ratio   finite; range 0.03711..1.445; outside valid_range [0.0, 3.4028234663852886e+38]: 0 of 720
PASS  single_scatter_albedo          finite; range 0.02356..1; outside valid_range [0.0, 1.0]: 0 of 720
PASS  asymmetry_parameter            finite; range 0.07353..0.9754; outside valid_range [-1.0, 1.0]: 0 of 720
PASS  NC_STRING contract: {'NC_STRING': 86}; ORAC spacing reads {'optical_depth': 'uneven_logarithmic', 'effective_radius': 'linear', 'solar_zenith': 'linear', 'satellite_zenith': 'linear', 'relative_azimuth': 'linear'}
PASS  provenance sidecar terra_modis_m_liquid-water_a01_pold_v22.provenance.json
PASS  provenance lut_version 22, code_release lut-code-v22.1, scientific_config v22
PASS  provenance commit aa62864 = submitted revision 4fee1b4 or differs from it in no production path
PASS  commit in repository; product SHA-256 matches
PASS  run-file SHA-256 equals the committed run file
PASS  numerics: fixed, nmom = 1000 (IDL generate_scattering_properties); legacy IDL linear-radius trapezoid, 0.001-100 um, size-param...; DISORT 60 streams; isothermal cloud (legacy IDL)
      working tree at generation: MODIFIED; tracked changes ['create_orac_lut/initpath.bash', 'create_orac_lut/makerunfile.pro', 'create_orac_lut/terra_modis_run', 'validation/v22_terra_modis/validate_products.py']; generated 2026-10-10T04:33:38Z on atmnode010.atm.ox.ac.uk
```

Cross-check with `terra_modis_m_liquid-water_a01_pold_v23.nc` (same scientific code, Grid B) at shared nodes:

```
  shared optical_depth: 17 values [0.0, 0.007799999788403511, 0.015599999576807022, 0.031199999153614044, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
  shared effective_radius: 20 values [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 37.0, 39.0]
  shared solar_zenith: 3 values [0.0, 30.0, 60.0]
  shared satellite_zenith: 1 values [0.0]
  shared relative_azimuth: 3 values [0.0, 90.0, 180.0]
  T_dv                           nodes     12240  BITWISE  max|d| 0  max rel 0
  T_dd                           nodes     12240  BITWISE  max|d| 0  max rel 0
  R_dv                           nodes     12240  BITWISE  max|d| 0  max rel 0
  R_dd                           nodes     12240  BITWISE  max|d| 0  max rel 0
  R_0v                           nodes     79560  BITWISE  max|d| 0  max rel 0
  R_0d                           nodes     26520  BITWISE  max|d| 0  max rel 0
  T_0d                           nodes     26520  BITWISE  max|d| 0  max rel 0
  T_00                           nodes     26520  BITWISE  max|d| 0  max rel 0
  E_md                           nodes      5440  BITWISE  max|d| 0  max rel 0
  extinction_coefficient         nodes       720  BITWISE  max|d| 0  max rel 0
  extinction_coefficient_ratio   nodes       720  BITWISE  max|d| 0  max rel 0
  single_scatter_albedo          nodes       720  BITWISE  max|d| 0  max rel 0
  asymmetry_parameter            nodes       720  BITWISE  max|d| 0  max rel 0
  channel_id                     nodes        36  BITWISE  max|d| 0  max rel 0
  central_wavelength             nodes        36  BITWISE  max|d| 0  max rel 0
  summary: all shared values bitwise identical
```

IDL V21 reference: none exists for Terra MODIS with this model (the V21/V22 campaign's Terra MODIS cases are compact liquid-water stg products); no V21 comparison is possible.

## terra_modis_m_liquid-water_a01_p240_v22.nc

SLURM job 490633 (v22_terra_water_240); run file `runs/terra_modis_cloud_liquid-water_240_v22.run`; submitted from revision `4fee1b4`.

```
PASS  opens; V2 reader accepts it; format NETCDF4; size 41330203 bytes
PASS  dimensions {optical_depth: 17, effective_radius: 20, solar_zenith: 10, satellite_zenith: 10, relative_azimuth: 11, channels: 36}
PASS  optical_depth equals liquid-water-cloud.lut (17 values, 1e-10..256)
PASS  effective_radius equals liquid-water-cloud.lut (20 values, 1..39)
PASS  solar_zenith equals liquid-water-cloud.lut (10 values, 0..90)
PASS  satellite_zenith equals liquid-water-cloud.lut (10 values, 0..65)
PASS  relative_azimuth equals liquid-water-cloud.lut (11 values, 0..180)
PASS  channel_id 1..36 (36 channels)
PASS  solar_channel_flag matches the instrument file
PASS  thermal_channel_flag matches the instrument file
      solar channels 26, thermal 16, mixed 6
PASS  T_dv                           finite; range -2.582e-15..1; outside valid_range [0.0, 1.0]: 1111 of 122400
PASS  T_dd                           finite; range 0..1; outside valid_range [0.0, 1.0]: 212 of 12240
PASS  R_dv                           finite; range -5.34e-06..0.9827; outside valid_range [0.0, 1.0]: 1713 of 122400
PASS  R_dd                           finite; range -0.0005941..0.9783; outside valid_range [0.0, 1.0]: 132 of 12240
PASS  R_0v                           finite; range -0.4739..1.191; outside valid_range [0.0, 1.0]: 103898 of 9724000
PASS  R_0d                           finite; range -0.001416..1.628; outside valid_range [0.0, 1.0]: 218 of 88400
PASS  T_0d                           finite; range -5.139e-05..0.8195; outside valid_range [0.0, 1.0]: 166 of 88400
PASS  T_00                           finite; range 0..1; outside valid_range [0.0, 1.0]: 0 of 88400
PASS  E_md                           finite; range -9.533e-11..0.9989; outside valid_range [0.0, 1.0]: 944 of 54400
PASS  extinction_coefficient         finite; range 0.192..7443; outside valid_range [0.0, 3.4028234663852886e+38]: 0 of 720
PASS  extinction_coefficient_ratio   finite; range 0.03467..1.484; outside valid_range [0.0, 3.4028234663852886e+38]: 0 of 720
PASS  single_scatter_albedo          finite; range 0.02382..1; outside valid_range [0.0, 1.0]: 0 of 720
PASS  asymmetry_parameter            finite; range 0.07662..0.9767; outside valid_range [-1.0, 1.0]: 0 of 720
PASS  NC_STRING contract: {'NC_STRING': 86}; ORAC spacing reads {'optical_depth': 'uneven_logarithmic', 'effective_radius': 'linear', 'solar_zenith': 'linear', 'satellite_zenith': 'linear', 'relative_azimuth': 'linear'}
PASS  provenance sidecar terra_modis_m_liquid-water_a01_p240_v22.provenance.json
PASS  provenance lut_version 22, code_release lut-code-v22.1, scientific_config v22
PASS  provenance commit aa62864 = submitted revision 4fee1b4 or differs from it in no production path
PASS  commit in repository; product SHA-256 matches
PASS  run-file SHA-256 equals the committed run file
PASS  numerics: fixed, nmom = 1000 (IDL generate_scattering_properties); legacy IDL linear-radius trapezoid, 0.001-100 um, size-param...; DISORT 60 streams; isothermal cloud (legacy IDL)
      working tree at generation: MODIFIED; tracked changes ['create_orac_lut/initpath.bash', 'create_orac_lut/makerunfile.pro', 'create_orac_lut/terra_modis_run', 'validation/v22_terra_modis/validate_products.py']; generated 2026-10-10T04:32:00Z on atmnode010.atm.ox.ac.uk
```

Cross-check with `terra_modis_m_liquid-water_a01_p240_v23.nc` (same scientific code, Grid B) at shared nodes:

```
  shared optical_depth: 17 values [0.0, 0.007799999788403511, 0.015599999576807022, 0.031199999153614044, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
  shared effective_radius: 20 values [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 37.0, 39.0]
  shared solar_zenith: 3 values [0.0, 30.0, 60.0]
  shared satellite_zenith: 1 values [0.0]
  shared relative_azimuth: 3 values [0.0, 90.0, 180.0]
  T_dv                           nodes     12240  BITWISE  max|d| 0  max rel 0
  T_dd                           nodes     12240  BITWISE  max|d| 0  max rel 0
  R_dv                           nodes     12240  BITWISE  max|d| 0  max rel 0
  R_dd                           nodes     12240  BITWISE  max|d| 0  max rel 0
  R_0v                           nodes     79560  BITWISE  max|d| 0  max rel 0
  R_0d                           nodes     26520  BITWISE  max|d| 0  max rel 0
  T_0d                           nodes     26520  BITWISE  max|d| 0  max rel 0
  T_00                           nodes     26520  BITWISE  max|d| 0  max rel 0
  E_md                           nodes      5440  BITWISE  max|d| 0  max rel 0
  extinction_coefficient         nodes       720  BITWISE  max|d| 0  max rel 0
  extinction_coefficient_ratio   nodes       720  BITWISE  max|d| 0  max rel 0
  single_scatter_albedo          nodes       720  BITWISE  max|d| 0  max rel 0
  asymmetry_parameter            nodes       720  BITWISE  max|d| 0  max rel 0
  channel_id                     nodes        36  BITWISE  max|d| 0  max rel 0
  central_wavelength             nodes        36  BITWISE  max|d| 0  max rel 0
  summary: all shared values bitwise identical
```

IDL V21 reference: none exists for Terra MODIS with this model (the V21/V22 campaign's Terra MODIS cases are compact liquid-water stg products); no V21 comparison is possible.

## terra_modis_m_water-ice_a01_pagg_v22.nc

SLURM job 490634 (v22_terra_ice_agg); run file `runs/terra_modis_cloud_water-ice_agg_v22.run`; submitted from revision `4fee1b4`.

```
PASS  opens; V2 reader accepts it; format NETCDF4; size 49582811 bytes
PASS  dimensions {optical_depth: 17, effective_radius: 24, solar_zenith: 10, satellite_zenith: 10, relative_azimuth: 11, channels: 36}
PASS  optical_depth equals ice-cloud.lut (17 values, 1e-10..256)
PASS  effective_radius equals ice-cloud.lut (24 values, 1..93)
PASS  solar_zenith equals ice-cloud.lut (10 values, 0..90)
PASS  satellite_zenith equals ice-cloud.lut (10 values, 0..65)
PASS  relative_azimuth equals ice-cloud.lut (11 values, 0..180)
PASS  channel_id 1..36 (36 channels)
PASS  solar_channel_flag matches the instrument file
PASS  thermal_channel_flag matches the instrument file
      solar channels 26, thermal 16, mixed 6
PASS  T_dv                           finite; range -8.604e-15..1; outside valid_range [0.0, 1.0]: 2534 of 146880
PASS  T_dd                           finite; range 0..1; outside valid_range [0.0, 1.0]: 254 of 14688
PASS  R_dv                           finite; range -5.34e-06..0.9901; outside valid_range [0.0, 1.0]: 2043 of 146880
PASS  R_dd                           finite; range -0.0001777..0.9877; outside valid_range [0.0, 1.0]: 140 of 14688
PASS  R_0v                           finite; range -6.125e-05..1.08; outside valid_range [0.0, 1.0]: 111845 of 11668800
PASS  R_0d                           finite; range -0.0009975..0.9948; outside valid_range [0.0, 1.0]: 350 of 106080
PASS  T_0d                           finite; range -0.0006489..0.6948; outside valid_range [0.0, 1.0]: 250 of 106080
PASS  T_00                           finite; range 0..1; outside valid_range [0.0, 1.0]: 0 of 106080
PASS  E_md                           finite; range -9.49e-11..0.9986; outside valid_range [0.0, 1.0]: 1152 of 65280
PASS  extinction_coefficient         finite; range 1.23..2.342; outside valid_range [0.0, 3.4028234663852886e+38]: 0 of 864
PASS  extinction_coefficient_ratio   finite; range 0.5949..1.133; outside valid_range [0.0, 3.4028234663852886e+38]: 0 of 864
PASS  single_scatter_albedo          finite; range 0.3279..1; outside valid_range [0.0, 1.0]: 0 of 864
PASS  asymmetry_parameter            finite; range 0.7426..0.983; outside valid_range [-1.0, 1.0]: 0 of 864
PASS  NC_STRING contract: {'NC_STRING': 86}; ORAC spacing reads {'optical_depth': 'uneven_logarithmic', 'effective_radius': 'linear', 'solar_zenith': 'linear', 'satellite_zenith': 'linear', 'relative_azimuth': 'linear'}
PASS  provenance sidecar terra_modis_m_water-ice_a01_pagg_v22.provenance.json
PASS  provenance lut_version 22, code_release lut-code-v22.1, scientific_config v22
PASS  provenance commit aa62864 = submitted revision 4fee1b4 or differs from it in no production path
PASS  commit in repository; product SHA-256 matches
PASS  run-file SHA-256 equals the committed run file
PASS  numerics: fixed, nmom = 1000 (IDL generate_scattering_properties); legacy IDL linear-radius trapezoid, 0.001-100 um, size-param...; DISORT 60 streams; isothermal cloud (legacy IDL)
      working tree at generation: MODIFIED; tracked changes ['create_orac_lut/initpath.bash', 'create_orac_lut/makerunfile.pro', 'create_orac_lut/terra_modis_run', 'validation/v22_terra_modis/validate_products.py']; generated 2026-10-10T04:55:24Z on atmnode010.atm.ox.ac.uk
```

Cross-check with `terra_modis_m_water-ice_a01_pagg_v23.nc` (same scientific code, Grid B) at shared nodes:

```
  shared optical_depth: 17 values [0.0, 0.007799999788403511, 0.015599999576807022, 0.031199999153614044, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
  shared effective_radius: 24 values [1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0, 29.0, 33.0, 37.0, 41.0, 45.0, 49.0, 53.0, 57.0, 61.0, 65.0, 69.0, 73.0, 77.0, 81.0, 85.0, 89.0, 93.0]
  shared solar_zenith: 3 values [0.0, 30.0, 60.0]
  shared satellite_zenith: 1 values [0.0]
  shared relative_azimuth: 3 values [0.0, 90.0, 180.0]
  T_dv                           nodes     14688  BITWISE  max|d| 0  max rel 0
  T_dd                           nodes     14688  BITWISE  max|d| 0  max rel 0
  R_dv                           nodes     14688  BITWISE  max|d| 0  max rel 0
  R_dd                           nodes     14688  BITWISE  max|d| 0  max rel 0
  R_0v                           nodes     95472  BITWISE  max|d| 0  max rel 0
  R_0d                           nodes     31824  BITWISE  max|d| 0  max rel 0
  T_0d                           nodes     31824  BITWISE  max|d| 0  max rel 0
  T_00                           nodes     31824  BITWISE  max|d| 0  max rel 0
  E_md                           nodes      6528  BITWISE  max|d| 0  max rel 0
  extinction_coefficient         nodes       864  BITWISE  max|d| 0  max rel 0
  extinction_coefficient_ratio   nodes       864  BITWISE  max|d| 0  max rel 0
  single_scatter_albedo          nodes       864  BITWISE  max|d| 0  max rel 0
  asymmetry_parameter            nodes       864  BITWISE  max|d| 0  max rel 0
  channel_id                     nodes        36  BITWISE  max|d| 0  max rel 0
  central_wavelength             nodes        36  BITWISE  max|d| 0  max rel 0
  summary: all shared values bitwise identical
```

IDL V21 reference: none exists for Terra MODIS with this model (the V21/V22 campaign's Terra MODIS cases are compact liquid-water stg products); no V21 comparison is possible.

Overall: PASS

## Notes

Runtimes (SLURM jobs on atmnode010, partition priority-eodg, 1 CPU, 8 GB,
24 h limit; started 2026-10-09 22:49:50 BST):

| job | product | end | wall time |
|---|---|---|---|
| 490632 | terra_modis_m_liquid-water_a01_pold_v22.nc | 2026-10-10 05:33:38 | 6 h 44 min |
| 490633 | terra_modis_m_liquid-water_a01_p240_v22.nc | 2026-10-10 05:32:02 | 6 h 42 min |
| 490634 | terra_modis_m_water-ice_a01_pagg_v22.nc | 2026-10-10 05:55:24 | 7 h 06 min |

The jobs were submitted to `shared` with 16 GB; every shared node had all its
memory allocated, so they were changed in place (`scontrol update
Partition=shared,priority-eodg MinMemoryNode=8192`) and started at once.

**Corrupted DISORT call (240 K product).**  One direct-beam DISORT call is
wrong: channel 6 (1.64 um), tau = 256, r_e = 29 um, solar zenith 80 deg.  It
gives `R_0d` = 1.628 and `R_0v` down to -0.474 at 52 view geometries of that
call; every other value of the three products lies within [-0.01, 1.3].  This
is the known sporadic whole-call failure of single-precision DISORT 2.0 at
NSTR 60 (memory note / validation/v24/results/REPORT_disort_options.md:
17 of 40 V23 and 14 of 35 V24 products affected), inherent in the legacy V22
solver; solar zenith 80 deg is not a Grid B node, so the V23 cross-check
cannot see it.  It was not corrected, because a fix (recomputation or double
precision) is a change of the scientific algorithm.  ORAC users of this table
should be aware of the defective node; the small negative values elsewhere
(down to -1.8e-4) are the usual legacy near-zero noise.

DISORT "SGECO says matrix near singular" warnings: 5 (old), 28 (240) and
6 (agg), against 71, 61 and 55 in the corresponding V23 runs on the larger
Grid B.

**Provenance timing.**  The sidecar is written when a job finishes and
records the Git state at that moment (here identical to the start-time
banner: `aa62864`, MODIFIED by the owner's uncommitted
`create_orac_lut/{initpath.bash,makerunfile.pro,terra_modis_run}` edits, which
the Python generator does not read, and the then-uncommitted
`validation/v22_terra_modis/validate_products.py`).  Committing during a run
would make the record name a later commit; capturing the state at start-up
would remove that hazard (an infrastructure change not made here).
