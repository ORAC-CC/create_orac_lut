# Legacy ORAC reconciliation

The root repository preserves the GitHub ORAC LUT V2.1 history imported from
commit `b3c8a9145880cafa7b4293aeb2d38078531d510c` under `create_orac_lut/`.
That history is the definitive historical baseline. Root commit `17591c9` is
the authoritative reference tree. The current local tree is treated as
evidence of later development, not as an authoritative replacement. The
Python migration must initially reproduce V2.1; deviations require an explicit
scientific reason, isolated validation, and user approval.

## Current implementation retained

The later split implementations are the current operational legacy path:

```text
create_orac_cloud_lut_wrapper -> create_orac_cloud_lut
create_orac_aerosol_lut_wrapper -> create_orac_aerosol_lut
                         -> setup_disort -> call_disort -> DISORT2
```

The split cloud and aerosol routines and wrappers are retained. They both use
the shared scattering-property generation and the production DISORT2 path.
The monolithic `create_orac_lut.pro` and its V2 wrapper remain historical or
diagnostic alternatives; their DISORT4 branch is not the current production
path.

## Modified tracked files

- `create_bwgp.pro`: **KEEP LOCAL**. The local `tmatrix_path` interface is
  required by current callers and replaces the historical `tmatrix_dir` name.
- `generate_scattering_properties.pro`: **KEEP LOCAL**. It adds the current
  refractive-index overrides, configurable moment count, Mie forcing, and the
  `tmatrix_path` interface.
- `dubovik/dubovik_lognormal_multiple_eps.pro`: **KEEP LOCAL**. Removes an
  unintended diagnostic print of the base path.
- `input_files/driver/sentinel-3a_slstr_aerosol.driver`: **KEEP LOCAL**. It
  reflects the current aerosol driver format and channel selection.
- `load_atmstr.pro`: **NEEDS USER DECISION**. The local numeric-atmosphere
  branch differs materially from the historical filename branch.
- `load_inststr.pro`: **KEEP LOCAL**. The local changes correct IDL `WHERE`
  no-match handling and normalize parsed channel identifiers.
- `load_lutstr.pro`: **NEEDS USER DECISION**. The local version removes the
  pressure-grid `Invalid_N_FLAG` check; this may be a regression and has not
  been overwritten automatically.
- `load_srfstrarr.pro`: **NEEDS USER DECISION**. The local QM=0/QM=1 behavior
  is reversed relative to V2.1 and affects spectral integration semantics.
- `makerunfile.pro`: **NEEDS USER DECISION**. The local file contains many
  commented/active experiment lists and is not a clean canonical driver.
- `read_ri.pro`: **KEEP LOCAL**, subject to a later syntax/runtime check. The
  local compressed-file branch uses the function argument consistently.
- `sentinel-3a_slstr_run`: **KEEP LOCAL only as historical evidence**. It
  invokes the later split aerosol wrapper but contains machine-specific
  paths and should not be treated as a production-safe runner.
- `create_orac_lut/.gitignore`: **NEEDS USER DECISION**. The local version
  drops the useful historical generated-output rules; root-level ignore rules
  now carry the broader policy, but the file was not overwritten here.

## Focused review of four unresolved files

### `load_atmstr.pro` — V2.1 authoritative; local variant forensic only

V2.1 selects the three-column `midsatm.dat` reader only when the selector is
the string `midsatm.dat`; other values use the RFM/MODTRAN array format. The
local version selects the three-column reader when `atmospheres EQ 0` and
otherwise retains the RFM path. Current split cloud and aerosol drivers use
numeric selectors, including `2` for `mls.atm`, and the Python reader maps
numeric code `0` to `midsatm.dat` and code `2` to `mls.atm`. Therefore the local
numeric behavior matches the exercised modern interface, but it is not part of
the authoritative V2.1 baseline. Do not merge it without a separate approved
change and validation.

### `load_lutstr.pro` — KEEP V2.1

V2.1 passes `Invalid_N_FLAG` while constructing the optional pressure grid and
checks that the declared pressure count matches the supplied values. The local
version omits that argument only for pressure, silently weakening validation.
The pressure grid is used by the aerosol formulation; cloud LUT definitions
do not request it. The Python parser validates the declared count for every
grid, including aerosol surface pressure, and the tests assert the expected
`[950, 1013, 1050]` pressure grid. There is no evidence of a legitimate
production case requiring malformed pressure input. Restore the V2.1 pressure
validation.

### `load_srfstrarr.pro` — V2.1 authoritative; local variant forensic only

The V2.1 comments define QM=1 as weighted centre-wavelength/monochromatic
calculation, QM=0 as raw SRF values, and QM=2 as segmented SRF quadrature, but
the V2.1 code implements QM=0 as centre wavelength and QM=1 as raw SRF values.
The local version reverses QM=0 and QM=1 to match the comments and the later
operational interface. Both current split generators default to QM=1, and the
validated shell runners explicitly use `srf_quad=1`. The Python implementation
also requires `srf_quad=1` and computes centre wavelength/wavenumber with unit
SRF weight, matching the local routine. The tiny/full validation documentation
records the same convention. This is an identified Python/local-versus-V2.1
discrepancy, not permission to change the V2.1 baseline. Retain the local
routine only as forensic evidence pending separate validation and approval.

### `makerunfile.pro` — HISTORICAL ONLY

This is an operational run-file generator, not a radiative-transfer or LUT
scientific kernel. The local version contains many successive active `MM`
assignments and platform/instrument selections, so its behavior depends on
which later assignment remains active. Current validation uses explicit,
repository-local shell runners and driver files instead. Do not make the local
experiment selections canonical. Keep the file only as historical evidence;
any useful instrument additions should later be merged into a separately
reviewed, declarative tool.

## Local additions

Retain as later source or operational evidence:

- `create_orac_cloud_lut.pro` and its wrapper;
- `create_orac_aerosol_lut.pro` and its wrapper;
- `generate_scattering_properties_simple.pro`;
- `call_disort4.pro`, `setup_disort4.pro`, and `disort4dlm/` as historical or
  diagnostic DISORT4 material only;
- current run scripts and input additions only where tied to supported drivers
  or documented validation cases.

The monolithic `create_orac_lut.pro`, `create_orac_lut_v2.pro`, and wrapper
variants are duplicate historical formulations. They should remain available
for provenance and comparison, but must not replace the split production path.
Run scripts, backups, `files.txt`, notes, and experiment-specific material
require individual review before staging.

Compiled objects and shared libraries under `disort4dlm/` are generated
artefacts and should not be retained in the baseline. Its source and Makefile
may be retained only as diagnostic/build provenance.

## Deleted V2.1 material

The deleted V2.1 attic/configuration text under `input_files/attic/` and the
legacy helper files needed by the source call graph were restored from the
imported baseline where their absence was not an intentional design decision.
The restoration list is recorded in the reconciliation work log.

The embedded `create_orac_lut/mie/` copy was not restored as an operational
Mie backend. The current validated environment uses the independent top-level
`mie/` repository. Historical embedded Mie material remains recoverable from
the imported V2.1 Git history.

## Mie provenance

Top-level `mie/` is independent upstream material:

- upstream: `https://github.com/eodg-code/mie`
- preserved nested HEAD: `67856daf66fa58a2aa95d391dfc9a6aca13d1ca9`

It should be tracked as current ordinary source files in the root repository,
with compiled products excluded, while its independent Git metadata remains in
the verified preservation record.

## Input data and generated material

The `create_orac_lut/input_files/` tree is scientifically important. Drivers,
instrument definitions, gas data, SRFs, refractive-index data, atmospheres,
microphysics, and LUT grid definitions should be retained when used by the
current supported cases or the Python validation suite. Large generated LUT
collections, compiler products, temporary files, editor backups, and compiled
DLMs should remain outside the baseline through the root ignore policy.

No IDL, Mie, DISORT, or LUT calculation was run during this review.
