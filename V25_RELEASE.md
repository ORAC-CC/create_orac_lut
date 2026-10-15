# ORAC LUT V25: the V24 calculation with vertically varying cloud temperature

V25 is the **V24 release (`V24_RELEASE.md`) with one addition: a vertically
varying cloud temperature in the thermal-emission calculation**.  The
cloud-temperature method is that of the historical V25 (commits `7f6b4e2`,
`13b9df3`, `5722280`, `990d970`; tag `lut-code-v25` = `f066fdd`),
re-implemented in the V24 code.  The historical V25 source was not restored.
Two numerical limitations of the historical emission layering were corrected
with the owner's approval.  The new V25 supersedes the historical V25.

| | |
|---|---|
| LUT release | **V25**: `releases/lut_v25.release` |
| Scientific code | the V24 code with `create_orac_luts_v07`, `write_v2_lut_v04` and the new `cloud_temperature_v02`; every other file is that of V24.  `CODE_VERSION` declares `lut-code-v25.2` (the format `oraclut.version` requires; no tag; part of the deferred naming review, `LUT_CODE_VERSIONS.md` section 12) |
| Unchanged | everything that is not thermal emission: Grid B, the V24 scattering refinements, the Mie kernel, Baum and T-matrix classes, the original single-precision Intel DISORT `DISORTtoIDL.so` (binary, wrapper, precision, calling convention), NSTR 60, NMOM = max(1000, L), atmosphere, gas, Rayleigh, SRF integration, instrument and microphysics inputs, the diffuse and direct-beam DISORT calls, the NetCDF variables and dimensions, the provenance scheme |
| Input added | `references/data/ocalut_cloudprofile_Cirrostratus.dat` (sha256 `8e9e68bb…`), pinned |
| Run files | `runs/terra_modis_cloud_{liquid-water_old,liquid-water_240,water-ice_agg}_v25.run`: the V24 run files with `version = 25`, `staging/lut_v25` and `cloud_vertical_profile` (`wet_adiabat`, `wet_adiabat`, `cirrostratus`).  The historical content of these names is in Git at `f066fdd` |
| Products | `staging/lut_v25/terra_modis_m_{liquid-water_a01_pold,liquid-water_a01_p240,water-ice_a01_pagg}_v25.nc`.  The shared ORAC LUT directory holds the historical V25 products under these names.  Once the new products have passed validation they may replace them; nothing is overwritten before that |

## The thermal emission (thermal channels, E_md only)

**Profiles** (physics, as the historical method):

* water ice (`cirrostratus`): the supplied cirrostratus profile, mapped
  through the cumulative 0.55 µm optical depth: T = 250 K + 8 K/km × H(tau)
  × zeta(F), with the depth H of 2-11 km taken from the file;
* liquid water (`wet_adiabat`): the saturated liquid-water adiabat from 250 K
  and 628 hPa along z = tau / 20 km^-1;
* the cloud top is at the V24 emission temperature, 250 K, and E_md is
  normalised by B(250 K), as in V24;
* absent `cloud_vertical_profile` means V24.

**Numerics** (corrected, `cloud_temperature.emission_layers`):

* **Thick clouds:** up to 46 geometric sub-layers from a 0.01 top sub-layer.
  Their source is the integral-preserving (L2) projection of the profile's
  Planck radiance in the channel's window.
* **Thin clouds:** where single-precision DISORT cannot resolve a source
  gradient, isothermal cells at their mean Planck radiance are joined by
  zero-thickness layers.  Below a channel optical depth of 10⁻³ the cloud is
  one layer with the optical-depth mean Planck radiance, the exact optically
  thin limit.
* **Effect:** this removes the historical failure (E_md ≈ −1.5 at
  tau = 1e-10 in ice) and the historical layering error (up to 1.2×10⁻² at
  3.7 µm).

## Validation

`validation/v25/REPORT_v25_temperature_profiles.md`:

* method recovery;
* the profile data;
* the corrected layering against independent references (exact without
  scattering, convergence with it, the optically thin limit, continuity
  sweeps);
* compact V24/V25 LUTs for the three Terra MODIS models (bitwise V24 without
  profile; only E_md changes; isothermal limit; direction; thin clouds;
  historical V25);
* runtime.

The full products are validated after production.

## Reproduce

    bash scripts/submit_v25_release.sh               # SLURM, three jobs (priority-eodg; estimates in the report section 11; limit 72 h)
    bash scripts/monitor_v25_release.sh
    # or directly:
    ORAC_LUT_RELEASE_MANIFEST=lut_v25.release ORAC_LUT_DISORT_PRECISION=single \
    ORAC_LUT_DISORT_SP_KERNEL=intel_dlm python create_orac_luts_release_v04.py runs/terra_modis_cloud_water-ice_agg_v25.run

The launcher refuses to run unless every sealed source and pinned input,
including the profile data and the Grid B definitions, is unchanged.  V21-V24
remain reproducible from their own sealed releases.

## Status

* compact validation passed;
* release sealed and committed;
* **full production not submitted** (owner review first);
* no V25 tag (deferred naming review).
