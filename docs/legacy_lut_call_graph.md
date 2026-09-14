# Legacy ORAC LUT call graph

## Primary V2-era call graph

```text
makerunfile_v2.pro
  ├─ selects final platform/instrument assignment
  ├─ selects active microphysical model list
  ├─ maps material prefix to driver/forward-model label
  └─ writes <platform>_<instrument>_run
       ├─ source initpath.bash
       ├─ export CREATE_ORAC_LUT_DRIVER=...
       └─ idl -e "create_orac_lut_wrapper_v2,..."
            └─ create_orac_lut_wrapper_v2
                 ├─ choose driver from keyword/env
                 ├─ parse mandatory and optional driver lines
                 ├─ derive out_path = luts/<driver stem>
                 └─ create_orac_lut_v2(...)
                      ├─ load_inststr
                      │    └─ instrument/channel/SRF filenames and flags
                      ├─ load_lutstr
                      │    └─ optical-depth/size/geometry grids
                      ├─ load_srfstrarr
                      │    ├─ SRF files
                      │    ├─ Gueymard2018.sssi
                      │    └─ bbconstants
                      ├─ load_mmdat
                      │    └─ particle profile, components, sizes, RI, code
                      ├─ load_atmstr
                      │    └─ layer heights, pressure, temperature
                      ├─ optional load_gasstr
                      │    └─ per-channel gas optical depth by layer
                      ├─ interpolate particle profile to RT layers
                      ├─ generate_scattering_properties
                      │    ├─ read_ri
                      │    ├─ create_range
                      │    ├─ create_bwgp
                      │    │    ├─ Mie DLM/backend
                      │    │    └─ Dubovik T-matrix backend
                      │    ├─ read_baran / read_baum
                      │    ├─ phase normalisation and Legendre expansion
                      │    └─ save/reload scatfile.sav
                      ├─ setup_disort
                      │    └─ initialise DISORT common block
                      ├─ per spectral point/grid/geometry:
                      │    ├─ construct gas + Rayleigh + particle tau
                      │    ├─ mix SSA and phase moments
                      │    └─ call_disort
                      │         └─ DISORT DLM / Fortran solver
                      ├─ SRF-weight and normalise operators
                      ├─ write_v2_lut
                      │    ├─ NetCDF4 dimensions/metadata/operators
                      │    └─ driver copy
                      └─ timestamp.txt
```

The historical V2 source also includes `@create_orac_lut_scat.pro` and
`@create_orac_lut_disort.pro` from `create_orac_lut_v2.pro`. Those source files
are deleted in the current working tree, while the current split files use
`generate_scattering_properties.pro`, `setup_disort.pro`, and
`call_disort.pro`. This is a source-state discrepancy to resolve before
attempting an IDL build; it was not repaired here.

## Driver and override precedence

The effective configuration is assembled in layers:

```text
driver file mandatory values
  + driver file optional keyword values
  + explicit keywords in generated IDL command
  = effective call to create_orac_lut_v2 / split generator
```

The generated command explicitly supplies `srf_quad`, `mmfile`, `lutfile`,
`tmatrix_path`, and `version`. Some later aerosol scripts additionally supply
`atmospheres=2,gas=1`. The command therefore has scientific precedence over
same-named driver values. The wrapper computes its output directory from the
driver filename, not from the model or LUT filename.

## Later split path

The current working tree also contains this distinct path:

```text
makerunfile.pro
  -> <platform>_<instrument>_run
  -> source initpath.bash; load Intel compiler module
  -> CREATE_ORAC_AEROSOL_LUT_DRIVER / CREATE_ORAC_CLOUD_LUT_DRIVER
  -> create_orac_aerosol_lut_wrapper / create_orac_cloud_lut_wrapper
  -> create_orac_aerosol_lut / create_orac_cloud_lut
  -> shared input/scattering/DISORT helpers
  -> write_v2_lut
```

The aerosol and cloud wrappers parse drivers in the same general way but call
different top-level procedures. Their current outputs differ structurally:

| path | pressure dimension | Rayleigh reference | writer call |
|---|---:|---|---|
| aerosol split | yes, when `/include_pressure` is used | LUT surface pressure | `write_v2_lut,/include_pressure` |
| cloud split | no | fixed 1013 hPa | `write_v2_lut` |

Both top-level procedures contain solar, mixed, and thermal channel handling.
The aerosol/cloud label is therefore not a reliable proxy for “particle versus
cloud particle”; the microphysical model is selected independently by `mmfile`.

## Scientific data flow

```text
instrument file ───────┐
SRF files + solar ─────┼─> channel/spectral state ─────────┐
LUT grid ──────────────┤                                    │
atmosphere ────────────┤                                    v
gas tables (optional) ─┤                         spectral RT loop
microphysics + RI ─────┘                                    │
       │                                                     │
       v                                                     v
particle Bext, SSA, g, phase/moments ──> layer optical state ──> DISORT
                                                               │
                                                               v
                                           SRF-integrated LUT operators + NetCDF
```

The particle-property boundary is the output of
`generate_scattering_properties`: wavelength-dependent extinction and ratios,
single-scattering albedo, asymmetry parameter, phase function/Legendre moments,
and average volume per particle over the effective-radius grid. The RT code
then combines those properties with atmospheric and Rayleigh terms. This is the
natural future interchange point for legacy optics and POM.

## Git-history call-graph evidence

Commit `e5758f6` is the clearest historical V2 snapshot: it adds the V2 top-level
routine/wrapper and the V2 run-file generator. Commit `b3c8a91` later removes
the separate top-level generators and wrappers and changes `makerunfile.pro` to
emit calls named `create_orac_lut_aerosol` and `create_orac_lut_cloud`. Those
entry points are not present in the inspected `b3c8a91` tree. That branch is
therefore evidence of an intended V2.1 redesign, not a verified operational
call graph.

## Known anomalies to preserve in diagnostics

- `makerunfile_v2.pro` overwrites its own model and platform selections.
- Older generated scripts use the unified wrapper; newer generated scripts use
  split wrappers.
- Some newer scripts export a split environment variable but echo the old
  `$CREATE_ORAC_LUT_DRIVER` variable.
- Wrapper comments retain obsolete driver/environment descriptions.
- `create_orac_lut_v2.pro` includes source files deleted from the current
  working tree.
- The current V2 source contains a stray diagnostic `print` in the thermal
  calculation area.
- Driver copies in output directories do not always contain every command-line
  override used to generate the associated NetCDF product.

These are provenance and reproducibility findings, not fixes to apply during
the first port.
