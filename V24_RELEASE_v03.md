# ORAC LUT V24: the V23 calculation with refined Mie scattering properties

V24 is the **V23 release (`V23_RELEASE.md`) with three refinements of the
scattering properties of Mie classes, and nothing else**.  The methods are
those of the historical commits `c6ad545`, `dbcc42c` and `96bebbe`
(`validation/REPORT_lut_numerics_development.md`), re-implemented in the V23
code; the historical V24 source was not restored.

| | |
|---|---|
| LUT release | **V24**: `releases/lut_v24.release` |
| Scientific code | **`lut-code-v24.1`** = `lut-code-v22.4` (V22/V23) with `create_orac_luts_v05`, `create_bwgp_v05`, `generate_scattering_properties_v04` and the new `legendre_expansion_v01`; every other file is that of V23 |
| Unchanged | original single-precision Intel DISORT `DISORTtoIDL.so` (binary, wrapper, precision, calling convention), NSTR 60, Mie kernel `mieint.f`, Baum interpolation, IDL-compatible arithmetic, atmosphere, gas, Rayleigh, SRF integration, instrument, microphysics and refractive-index inputs, isothermal cloud emission, NetCDF schema, provenance scheme |
| Grid | **Grid B, as V23** (`liquid-water-cloud.lut`, `ice-cloud.lut`: SourceVault `_v02`), pinned by the manifest |
| CODE_VERSION | `code_release = lut-code-v24.1` (the historical tag `lut-code-v24` is `51717b8`; whether this extra code-release number is tagged is part of the deferred naming review, `LUT_CODE_VERSIONS.md` section 12), `lut_version = 24`, no compatible versions (earlier versions come from their sealed releases), `scientific_config = v24-scattering-refinements` |
| Run files | `runs/terra_modis_cloud_{liquid-water_old,liquid-water_240,water-ice_agg}_v24.run`: the V23 run files with `version = 24` and `staging/lut_v24` (the historical content of these names is in Git at `51717b8`) |
| Products | `staging/lut_v24/terra_modis_m_{liquid-water_a01_pold,liquid-water_a01_p240,water-ice_a01_pagg}_v24.nc` (the shared ORAC LUT directory holds historical `lut-code-v24` products with the same names; by owner decision of 2026-10-13 the new V24 supersedes them, and once the three new products have completed and passed validation they may replace those files; nothing is overwritten before that: `LUT_CODE_VERSIONS.md` section 12) |

## The three refinements (Mie classes only)

1. **Liquid-water integration to 3.5 r_e** (`c6ad545`).  A liquid-water
   modified-gamma integration ends at the first node of the legacy radius
   lattice (the 0.001-100 um trapezoid, size-parameter step 0.4) at or beyond
   3.5 r_e, beyond 100 um when needed, instead of at 100 um.  Effective
   variances above 0.1111111 are refused.  Ice spheres keep 0.001-100 um.
2. **Nested radius refinement** (`96bebbe`).  For liquid-water and ice-sphere
   modified-gamma components every interval of that lattice is divided into
   2^k equal parts, k the smallest level giving a size-parameter step
   <= 0.025 (liquid water) or 0.05 (ice spheres) and at least three nodes per
   r_e sqrt(v_e); legacy nodes and limits are kept.  The Mie kernel is called
   on blocks of radius nodes with the IDL TOTAL summation order preserved.
3. **Adaptive Legendre expansion** (`dbcc42c`).  Every size-averaged phase
   function (one per wavelength and r_e) is sampled on a Newton-refined
   Gauss-Legendre rule above its polynomial-degree bound D = 2 NStop(x_max)
   (no aliasing) and keeps L coefficients by King's criterion
   |omega_l| < 1e-9 for all l >= L (Grainger 1990), accepted only if the
   series reproduces the phase function to 5e-7 and the rounding noise beyond
   D is below 1e-6.

**DISORT moments.**  Each call receives **NMOM = max(1000, L)** moments:
the established 1000-moment array, zero-padded beyond L, whenever L <= 1000,
and every moment of the expansion when L > 1000 (owner decision,
2026-10-13; the historical method used max(L, NSTR + 1)).  No moment is ever
truncated (`create_orac_luts._disort_moments` refuses a non-zero moment
beyond those passed).  On the original binary, padding is bitwise neutral,
max(1000, L) is bitwise max(L, 61) for L < 1000, the fluxes are bitwise
independent of the moments beyond NSTR, and the reflectances converge as
NMOM approaches L.

**Baum and T-matrix classes** (here pagg) are calculated exactly as in V23:
the compact V24 and V23 pagg products are identical in every variable.

## Validation

`validation/v24/REPORT_v24_scattering_refinements.md`: method recovery,
scattering-property convergence, DISORT moment verification, compact V23/V24
LUT comparison, runtime, deviations from the historical methods.  The full V24
products are validated after production with
`validation/v24/scattering_refinements/compare_full_products.py`.

## Reproduce

    bash scripts/submit_v24_release.sh               # SLURM, three jobs (priority-eodg; about 32-34 h liquid water, 19 h pagg; limit 72 h)
    bash scripts/monitor_v24_release.sh [--compare]
    # or directly:
    ORAC_LUT_RELEASE_MANIFEST=lut_v24.release ORAC_LUT_DISORT_PRECISION=single \
    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/terra_modis_cloud_liquid-water_240_v24.run

The launcher refuses to run unless every sealed source and pinned input
(including the Grid B definitions) is unchanged.  V21, V22 and V23 remain
reproducible from their own sealed releases; the V21/V22 regression tests run
the V22 code from `releases/lut_v23.release` (`tests/sealed_release.py`).

## Status (2026-10-15)

The three full products were produced (SLURM 495982-495984, `e2c666a`,
18.8-33.2 h) and **passed validation**
(`validation/v24/REPORT_v24_scattering_refinements.md` section 13):

| Product (`staging/lut_v24/`) | SHA-256 |
|---|---|
| `terra_modis_m_liquid-water_a01_pold_v24.nc` | `e45035bda89116305ebdee4ab42ffc5ef252016b5b2e918ed1265a00b6fcbee0` |
| `terra_modis_m_liquid-water_a01_p240_v24.nc` | `1287c5dbdbc96811ca4f74e6b93e4b481e335d413ec0754cabaf85e1d12c5c68` |
| `terra_modis_m_water-ice_a01_pagg_v24.nc` | `ffc61063f542e217936c71275c2e764dae2dc522ef344f203725e64aef779fa1` |

`releases/lut_v24.dependencies` records their inputs and external shared
objects.  Known limitation, as in V23: isolated single-precision DISORT
failures at a few nodes (p240: channel 17, tau 8, r_e 17 um, SZA 37.5 deg and
channel 16, tau 1/64, r_e 17 um, SZA 67.5 deg; pold: a partial failure at
channel 12, tau 0.5, r_e 19 um, SZA 33.75 deg).

**Promotion (pending, owner action).**  The historical V24 files are moved
unchanged to a superseded directory, as for V22, and the validated products
installed:

    L=/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS
    S=$L/superseded_v24_lut-code-v24_20261015
    mkdir "$S"
    for f in terra_modis_m_liquid-water_a01_pold_v24 terra_modis_m_liquid-water_a01_p240_v24 terra_modis_m_water-ice_a01_pagg_v24; do
       mv -n "$L/$f.nc" "$S/"
       cp -n staging/lut_v24/$f.nc staging/lut_v24/$f.provenance.json "$L/"
       chmod 644 "$L/$f.nc" "$L/$f.provenance.json"
    done
    sha256sum $L/terra_modis_m_*_v24.nc        # must equal the table above
