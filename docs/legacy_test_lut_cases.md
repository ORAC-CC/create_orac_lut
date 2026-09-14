# Existing legacy test LUT definitions

The repository already contains small LUT-grid definitions intended for fast
legacy development runs. They are preferable to inventing a new reduced grid:
the current loaders and generators consume them unchanged, while the
particle, atmosphere, SRF, Mie, and DISORT inputs remain the normal legacy
inputs.

## First fast development reference

The first case is:

```text
Meteosat-10 SEVIRI
current create_orac_cloud_lut_wrapper -> create_orac_cloud_lut path
liquid-water_stg.mm
liquid-water-cloud_test.lut
SRF mode 1, Rayleigh enabled, gas disabled, 60 DISORT streams
```

The exact existing grid file is
`create_orac_lut/input_files/lut/liquid-water-cloud_test.lut`:

| coordinate | spacing | values | count |
|---|---|---|---:|
| optical depth | linear | `0.5`, `10` | 2 |
| effective radius | linear | `1`, `5.5`, `10` | 3 |
| solar zenith | linear | `0`, `90` degrees | 2 |
| satellite zenith | linear | `0`, `90` degrees | 2 |
| relative azimuth | linear | `0`, `180` degrees | 2 |

There is no pressure coordinate in a cloud LUT. `load_lutstr.pro` expands the
two-endpoint linear definitions using `findgen(N)/(N-1)`; this is why this
existing test file uses two or three points rather than one-point invented
grids. The values propagate directly into `lutstr.opd`, `lutstr.efr`,
`lutstr.soz`, `lutstr.saz`, and `lutstr.raa`, which size the current generator
loops and output arrays.

The model remains the production STG definition
`create_orac_lut/input_files/microphysics/liquid-water_stg.mm`: one modified-
gamma water component, `Rm=12`, `S=0.1111111`, Mie scattering, and
`H2O_Segelstein_1981.ri`. Only the sampling grid is reduced.

The prepared runner selects channel 1 through the existing wrapper
`channelid` keyword. This is a visible solar SEVIRI channel, so the first
reference exercises STG -> Mie -> phase/moment processing -> production
DISORT -> V2 solar operators without thermal calculations. The source
instrument still has 11 available channels; the explicit selection reduces
the resolved calculation to one channel. The test case therefore has
`1 * 2 * 3 * 2 * 2 * 2 = 48` particle/geometry combinations before the
current solar operator branches, with one SRF point because `srf_quad=1`.

The full production acceptance target remains separate:

```text
Meteosat-10 SEVIRI liquid-water STG, V2, revision v21
liquid-water-cloud.lut, all 11 channels and the complete production grid
```

Passing the test grid does not establish full production equivalence.

## Other existing test grids

These files were inventoried but are not run by the first case:

| file | dimensions | pressure |
|---|---|---|
| `aerosol_test.lut` | 2 optical-depth, 2 effective-radius, 2 solar-zenith, 2 satellite-zenith, 2 relative-azimuth | 3 surface-pressure |
| `ash-plume_test.lut` | 2 optical-depth, 2 effective-radius, 2 solar-zenith, 2 satellite-zenith, 2 relative-azimuth | none |
| `biomass-plume_test.lut` | 2 optical-depth, 2 effective-radius, 2 solar-zenith, 2 satellite-zenith, 2 relative-azimuth | none |
| `ice-cloud_test.lut` | 2 optical-depth, 2 effective-radius, 2 solar-zenith, 2 satellite-zenith, 2 relative-azimuth | none |
| `sulphuric-acid-cloud_test.lut` | 2 optical-depth, 2 effective-radius, 2 solar-zenith, 2 satellite-zenith, 2 relative-azimuth | none |

`aerosol_test.lut` is the only listed test grid with a pressure section. The
other files are compact cloud/plume grids and retain the same five cloud-grid
coordinates. They are future fast validation assets, not alternative models
for the first STG case.

## Reproduction and capture

Source `scripts/legacy_idl_env.sh` in a shell with the AOPP module system. The
safe developer runner is:

```bash
scripts/generate_liquid_water_cloud_test_reference.sh
```

It checks the resolved dimensions and all required input/DLM files before IDL
starts, resolves the current cloud-path routines, and writes below
`validation/generated/meteosat10_seviri_liquid_water_stg_cloud_test/`.
The legacy IDL generator assumes its working directory is
`create_orac_lut`: several loaders construct relative names such as
`input_files/srf/...`, and the wrapper constructs `luts/<driver-stem>`. The
outer runner hides this requirement by changing to `create_orac_lut` for the
IDL process. The current cloud generator determines its canonical output
directory from the active driver and writes the NetCDF under
`create_orac_lut/luts/meteosat-10_seviri_cloud/`; the wrapper's nominal
`out_path` does not provide a reliable relocation interface. The runner
fingerprints that canonical file before and after the run, rejects an
unchanged file, then copies the verified result into a hidden repository-local
staging directory and promotes it to the validation directory. Failed runs
remove only staging state, and a completed captured reference is never
overwritten automatically.

The relevant path audit is:

- `create_orac_cloud_lut_wrapper.pro` opens the driver path supplied through
  `CREATE_ORAC_CLOUD_LUT_DRIVER` and constructs `out_path` as
  `luts/<driver-stem>`.
- `create_orac_cloud_lut.pro` constructs instrument, microphysics, LUT, and
  atmosphere paths below the driver `in_path` (the test uses `input_files`),
  and constructs the solar path as `input_files/sun/Gueymard2018.sssi`.
- `load_srfstrarr.pro` directly opens `input_files/srf/<SRF filename>`.
- `load_mmdat.pro` calls `read_ri` with `input_files/ri/<RI filename>`.
- `load_inststr.pro`, `load_lutstr.pro`, `load_atmstr.pro`, and
  `read_srfstr.pro` open the paths passed to them; those paths inherit the
  relative-root assumptions above.
- The cloud generator writes `scatfile.sav`, the NetCDF, and `timestamp.txt`
  below `out_path`; `write_v2_lut.pro` receives the already-constructed NetCDF
  path and does not relocate it.

For the current canonical Meteosat-10 cloud driver, the expected output is:

```text
create_orac_lut/luts/meteosat-10_seviri_cloud/
  meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc
```

The validation copy uses the unmistakable filename
`meteosat10_seviri_liquid_water_stg_cloud_test_legacy_reference_v21.nc` and
records the source path, source hash, source size, capture time, IDL version,
and production DLM paths in `capture_manifest.json`.

The current generator already saves `scatfile.sav`. After a successful run,
the validation-only capture procedure writes text arrays for the saved
reference extinction, SSA, asymmetry, phase, Legendre moments, extinction
ratio, and average-volume quantities. The generated V2 NetCDF contains the
operator quantities written by `write_v2_lut`. The unmodified production
routine does not expose each per-call DISORT layer array, so those inputs are
explicitly recorded as a remaining diagnostic boundary rather than recreated
by a guessed hook.

The Python comparison command is prepared as:

```bash
PYTHONPATH=src python3 -m oraclut.validation.compare \
  legacy-test.nc python-candidate.nc \
  --json-out validation/comparisons/legacy-test-summary.json
```

It compares shapes, dtypes, units, finite residuals, NaN/Inf masks, and exact
stored dimension ordering. Phase-function normalization and Legendre ordering
are reported as unresolved until the captured arrays are compared under the
legacy convention; the tool does not invent a tolerance or silently transpose
them.

## Environment provenance

The verified user-terminal environment is IDL 8.9.0 with
`intel-compilers/2022` and `idl/890` loaded. The current production DLM paths
are `mie/dlm-code` for `MIE_DLM_SINGLE` and
`create_orac_lut/disort2/src` for `DISORT`, `GETMOM`, and `PLKAVG`. The
`setup_disort`/`call_disort` path is current production; DISORT4 is not used.
Recursive project-wide DLM discovery is deliberately avoided because it can
expose conflicting DISORT interfaces.

The user has verified in that licensed terminal that the Mie entry point is
registered and that `CREATE_BWGP`, `SETUP_DISORT`, and `CALL_DISORT` resolve.

The user's ordinary VS Code terminal acquires the Oxford IDL licence. Codex's
isolated shell can load the modules but has previously failed to acquire that
licence; this is an environment limitation, not a scientific result. No
licence or shell-startup configuration is changed by the runner.
