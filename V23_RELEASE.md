# ORAC LUT V23: the V22 scientific code on the approved Grid B

V23 is a **grid-only successor of the definitive V22** (`V21_V22_REFERENCE.md`).

| | |
|---|---|
| LUT release | **V23** — `releases/lut_v23.release` |
| Scientific code | **`lut-code-v22.4`, unchanged**: every code, library and Intel run-time entry of the manifest is copied from `releases/lut_v22_single-r2.release` (driver `create_orac_luts_v04.py`, `legacy_disort_v04`, `idl_math_v02`, `baum_v03`, ...); original single-precision Intel DISORT `DISORTtoIDL.so`, NSTR 60, NMom 1000 |
| Grid | **Grid B**, installed as `liquid-water-cloud.lut` and `ice-cloud.lut` (SourceVault `liquid-water-cloud_v02.lut`, `ice-cloud_v02.lut`; the V21/V22 grids are the `_v01` copies) |
| CODE_VERSION | `code_release = lut-code-v22.4`, `scientific_config = v22-idl-single` (unchanged), `lut_version = 23`, compatible 21 and 22 |
| Run files | `runs/terra_modis_cloud_{liquid-water_old,liquid-water_240,water-ice_agg}_v23.run` (identical to the V22 run files except `version = 23` and the output directory) |
| Products | `staging/lut_v23/terra_modis_m_{liquid-water_a01_pold,liquid-water_a01_p240,water-ice_a01_pagg}_v23.nc` |

## The approved grid

Grid B is the outcome of the V23 value-and-Jacobian grid study (September
2026; the acceptance criterion was the interpolated forward model and its
optical-depth and effective-radius Jacobians over the ORAC operational domain,
solar and satellite zenith 0-75 deg, with liquid water, ice and aerosol
optimised separately).  Evidence:

* commit `4bb0ac7` "Add definitive Grid B ORAC LUT definitions": the files
  were populated from the study's Stage 18 final selections
  (`validation/v23/results/stage18_final_grid_{water,ice}.json`, kept in the
  study's backup tree);
* `4bb0ac7:validation/v23/results/grid_b_manifest.json` records the SHA-256 and
  the exact axes of each file: `b640697cf268936666cf...` (liquid water) and
  `022aba67de745e3716b3...` (ice), which the installed definitions match byte
  for byte;
* `GRID_B.md`, `docs/grid_b_lut_definitions.md` (at `4bb0ac7`) and
  `LUT_CODE_VERSIONS.md` (V23 = V22 code with Grid B); the files were unchanged
  until `d38af77` moved them to `create_orac_lut/input_files/lut/archive/`.

| Axis | V22 liquid | **V23 liquid** | V22 ice | **V23 ice** |
|---|---|---|---|---|
| optical depth | 17: 1e-10, 2^-7 ... 2^8 | **24**: V22 nodes plus 2^1.5 ... 2^7.5 | 17 | **24** (as liquid) |
| effective radius (um) | 20: 1, 3, ..., 39 | **24**: 1-9 by 1, 11-39 by 2 | 24: 1, 5, ..., 93 | **29**: 1-21 by 2, 25-93 by 4 |
| solar zenith (deg) | 10: 0-90 by 10 | **21**: 0-75 by 3.75 | 10 | **21** |
| satellite zenith (deg) | 10: 0-65 (MODIS maximum) | **21**: 0-75 by 3.75 | 10 | **21** |
| relative azimuth (deg) | 11: 0-180 by 18 | **37**: 0-180 by 5 | 11 | **37** |
| grid points per channel | 374 000 | **9 398 592** | 448 800 | **11 356 632** |

Grid B stops at 75 deg (no 80/90 deg solar-zenith nodes) and samples satellite
zenith to 75 deg for every instrument.

## Equivalence with V22

At the nodes common to both grids (optical depths of V22, shared radii, solar
zenith 0/30/60, satellite zenith 0, azimuth 0/90/180) V23 reproduces V22
**bitwise** in every quantity that does not depend on the user viewing angles
(R_dd, T_dd, R_0d, T_0d, T_00) and in T_dv.  The upwelling intensities R_0v,
R_dv and the emissivity E_md differ slightly (up to ~1e-5 in the compact
check): with identical inputs the original DISORT kernel's intensities at a
given view depend on the set of user polar angles requested (NUMU 42 instead
of 20 changes its arithmetic order; demonstrated directly on captured calls,
`validation/v23_terra_modis/README.md`).  This is the unchanged V22 code, not a
code difference.

## Reproduce

    bash scripts/submit_v23_release.sh              # SLURM, three jobs, 48 h limit each
    bash scripts/monitor_v23_release.sh [--compare]
    # or directly:
    ORAC_LUT_RELEASE_MANIFEST=lut_v23.release ORAC_LUT_DISORT_PRECISION=single \
    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/terra_modis_cloud_liquid-water_240_v23.run

The launcher refuses to run unless the sealed sources and the pinned inputs
(including the Grid B definitions) are unchanged.

## Earlier "V23" products

The shared ORAC LUT directory still holds V23 products made on 2026-09-29 by
`lut-code-v23` (the uncorrected V22 code of `9d663e9`, gfortran DISORT, with
Grid B).  They carry the same file names; the new products stay in
`staging/lut_v23/` until the owner decides how to promote them.
