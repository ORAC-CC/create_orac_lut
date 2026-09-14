# Legacy ORAC LUT workflow

## Scope and evidence

This is a forensic description of the current IDL workflow present in the
repository, retained as input to reproducing the selected current V2 product.
Historical versions are mentioned only where needed to distinguish current
entry points. The legacy source tree and generated drivers were treated as
read-only. The working tree is already substantially dirty, so observations
below distinguish current files from committed historical versions where that
matters.

The strongest evidence is:

1. `create_orac_lut/makerunfile_v2.pro` and the generated `*_run` scripts;
2. the driver, instrument, microphysics, LUT, atmosphere, SRF, gas, and optical
   property files under `create_orac_lut/input_files/`;
3. the IDL procedures called by those scripts;
4. existing NetCDF products under `create_orac_lut/luts/`;
5. current source and output provenance markers; Git history is secondary and
   used only to resolve current-path ambiguity.

The exact papers named in the project brief (`2009Thomas2.pdf`,
`2018Sus1.pdf`, and `2018McGarragh1.pdf`) are not present under the repository
root. Related PDFs are present, but they have not been substituted as if they
were the requested references.

## Current operational workflow

For the selected Meteosat-10 SEVIRI liquid-water product, the current source
path is the cloud split path:

```text
current cloud driver
  -> create_orac_cloud_lut_wrapper
  -> create_orac_cloud_lut
  -> current input/scattering/RT routines
  -> write_v2_lut
```

The exact cloud shell script is not retained, but the archived cloud driver,
current cloud wrapper/generator, output metadata, and `git_revision.txt`
provide the configuration and schema evidence. See
`docs/current_reference_case.md`.

## V2-era path retained for source comparison

The primary V2-era path is:

```text
makerunfile_v2.pro
  -> <platform>_<instrument>_run
  -> source initpath.bash
  -> CREATE_ORAC_LUT_DRIVER=<instrument/forward-model driver>
  -> idl -e create_orac_lut_wrapper_v2(...)
  -> create_orac_lut_wrapper_v2
  -> create_orac_lut_v2
  -> input loaders + generate_scattering_properties
  -> setup_disort / call_disort / DISORT DLM
  -> write_v2_lut
  -> luts/<driver stem>/<platform>_<instrument>_<m|b>_..._v<version>.nc
```

`makerunfile_v2.pro` is a generator, not the LUT calculation itself. It first
defines a broad master list of microphysical model names and then overwrites it
with an active list. It also performs a sequence of platform/instrument
assignments; the last active assignment in the inspected file is
`himawari-8`/`ahi`. With `Test = 0`, `Versions = '21'`, and `srf_quad = 1`, the
generated file is `himawari-8_ahi_run`, and its active model list is the final
ten-entry list in the generator. This “last assignment wins” pattern means the
`.pro` file should not be read as a declarative list of all intended products.

For each model, the generator derives `material` from the part of the model
name before the first underscore. In the current V2 generator it maps:

| material prefix | selected forward-model label | LUT filename prefix |
|---|---|---|
| `aerosol` | `aerosol` | `aerosol` |
| `biomass` | `cloud` in the active code | `biomass-plume` |
| `liquid-water` | `cloud` | `liquid-water-cloud` |
| `sulphuric-acid` | `cloud` | `sulphuric-acid-cloud` |
| `volcanic-ash` | `cloud` | `ash-plume` |
| `water-ice` | `cloud` | `ice-cloud` |

The model name is passed as `mmfile`, the selected LUT definition as
`lutfile`, and the driver as an environment variable. For liquid water, and
for the special `volcanic-ash_htha` name, the generator emits a second command
with `no_rayleigh=1,reuse_scat=1`. The loop only emits `srf_quad=1`; the band
branch exists in the generator but is not active in the inspected file.

The generated V2 script uses the old unified environment variable
`CREATE_ORAC_LUT_DRIVER` and calls `create_orac_lut_wrapper_v2`. For example,
`create_orac_lut/himawari-8_ahi_run` selects
`input_files/driver/himawari-8_ahi_cloud.driver`, passes the model and LUT names
on the IDL command line, supplies a hard-coded Dubovik T-matrix path, and
passes `version=21`.

## Driver interpretation

The wrapper reads the driver named by its keyword or by the relevant
environment variable. It strips blank/comment text, removes trailing comments,
joins lines ending in `$`, and interprets the first five logical values as:

1. input root, normally `input_files`;
2. instrument definition filename;
3. microphysical model filename (usually overridden by the run command);
4. LUT-grid definition filename (usually overridden by the run command);
5. atmosphere selector.

The remaining lines are optional keyword assignments. Observed keywords include
`channelid`, `gas`, `no_rayleigh`, `no_screen`, `reuse_scat`, `scat_only`,
`srf_quad`, `n_theta`, `srfdat`, `tmatrix_path`, `opt_prop_luts`, and `version`.
Explicit wrapper keywords overwrite values read from the driver. The output
directory is derived from the driver basename, for example:

```text
input_files/driver/meteosat-10_seviri_aerosol.driver
  -> luts/meteosat-10_seviri_aerosol/
```

This has two practical consequences for reproducibility:

- the driver remains part of the scientific input even when `mmfile`, `lutfile`,
  `version`, or other options are passed explicitly by the shell script;
- the output directory identifies the driver, but not necessarily every command
  line override used to make a product.

The current aerosol and cloud wrappers use separate environment variables,
`CREATE_ORAC_AEROSOL_LUT_DRIVER` and `CREATE_ORAC_CLOUD_LUT_DRIVER`, and call
`create_orac_aerosol_lut` or `create_orac_cloud_lut`. The comments at the top of
those wrappers still describe the older unified environment name, so the code,
not the stale comment, is authoritative.

## Inputs and calculation stages

`create_orac_lut_v2.pro` loads the instrument, LUT grid, SRFs, microphysics,
atmosphere, and optional gas tables. It maps atmosphere codes `0`–`6` to the
available atmospheric profiles (`midsatm.dat`, `tro.atm`, `mls.atm`, `mlw.atm`,
`sas.atm`, `saw.atm`, and `std.atm`). The particle profile in the microphysics
file is interpolated to the RT layers and normalised to a total relative
scattering optical depth.

Scattering properties are generated once for the requested channel spectral
points and microphysical grid. The result is saved as `scatfile.sav` in the
output directory and then reloaded by the V2 routine. For each channel, size,
particle optical depth, pressure where enabled, and geometry, the code:

- combines particle extinction with Rayleigh and optional gas optical depth;
- forms layer single-scattering albedo and phase-function moments;
- calls the DISORT boundary through `call_disort`;
- repeats over SRF quadrature points;
- normalises the accumulated operators by the SRF weights;
- writes solar, diffuse, direct, thermal, and emissivity operators to NetCDF.

The V2 routine uses 60 DISORT streams. For a solar calculation it uses a fixed
diffuse source and a 50-degree solar zenith for the diffuse operators, then
separate direct-beam solar geometries. Thermal calculations use the instrument
thermal flag, a zero-illumination Planck calculation, and a 250 K source in the
historical V2 file. The current aerosol source contains a working-tree change
to 270 K; this is one reason not to silently mix current aerosol and historical
V2 results.

The Rayleigh column optical depth in historical V2 is calculated from the
surface pressure relative to 1013 hPa and a wavelength-dependent expression.
The current aerosol split path adds a pressure dimension and scales Rayleigh
optical depth with the LUT pressure; the current cloud path retains the fixed
1013 hPa treatment and no pressure dimension. These are real downstream
implementation differences, not merely names for particle types.

## V2 meaning

“V2” is best treated as a product/workflow generation family, not as a unique
particle class or a single immutable entry-point name. Evidence for this is:

- commit `e5758f6` introduces `create_orac_lut_v2.pro`, its wrapper, the V2
  scattering/DISORT helpers, V2 drivers, and `makerunfile_v2.pro`;
- the V2 output writer is `write_v2_lut`, which creates NetCDF4 products and
  uses the `m`/`b`, atmosphere, pressure, particle, and version fields in the
  filename;
- later current scripts call separate aerosol/cloud routines that also write
  the V2 NetCDF schema;
- commit `b3c8a91` labels a different, incomplete job-sheet/path-handling
  redesign as V2.1 and removes the separate generator/wrapper sources.

Therefore, a Python reproduction must record at least the exact source path,
wrapper, driver, overrides, version, and output schema. A `v21` filename alone
does not uniquely identify the IDL call graph.

## Dependencies and execution assumptions

The source requires an IDL runtime and IDL NetCDF support. It also expects
compiled or externally supplied Fortran/C components for Mie/T-matrix
calculations and DISORT. The verified user-terminal environment loads
`intel-compilers/2022` and `idl/890`, with IDL 8.9.0 licensed by Oxford. The
current production DLM directories are the repository Mie DLM and
`create_orac_lut/disort2/src`; `setup_disort` and `call_disort` invoke the
production `DISORT` DLM. DISORT4 is not part of this path.

Codex's isolated shell can locate the IDL executable and load the modules but
has previously failed to acquire the Oxford licence. The repository-local
`scripts/legacy_idl_env.sh` records the temporary environment without editing
startup or licence configuration. The user's ordinary VS Code terminal is the
execution environment for the explicit fast test runner when Codex cannot
obtain a licence. No production driver is run by that runner.

## Current reproduction sequence

Use the existing products in this order:

1. Read and structurally validate the external Meteosat-10 SEVIRI cloud,
   `liquid-water_stg`, V2 reference;
2. reproduce its current cloud-path intermediate inputs and operators;
3. validate the current pressure-aware aerosol path with the Meteosat-10
   `a79` V2 product;
4. progress to the current V2 Meteosat-12 FCI volcanic-ash products in the
   external archive.

The first target should preserve the legacy array ordering and NetCDF variable
dimensions exactly. Compare intermediate scattering properties before comparing
the full RT operators.
