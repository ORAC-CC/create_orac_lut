# Workspace state and quota cleanup — 2026-09-13

Scope: workspace stabilisation, quota recovery and read-only Git forensics.
No scientific reconciliation, no source modification, no Git write operations.

## 1. Root Git baseline

The canonical root repository was created manually by the user.

| Commit | Meaning |
| --- | --- |
| `d97bbc1` | Initialize canonical ORAC LUT repository |
| `17591c9` | Import ORAC LUT V2.1 history under `create_orac_lut` (current `main`) |
| `b3c8a91` | Historical V2.1 tip, 2026-05-11, "Update ORAC LUT v2.1 workflow and path handling" |

`refs/heads/main` is the only ref. The V2.1 history is preserved in the root
history and must not be rewritten or re-imported.

## 2. Cloud/aerosol historical finding (read-only investigation)

The four files in question:

- `create_orac_aerosol_lut.pro`
- `create_orac_aerosol_lut_wrapper.pro`
- `create_orac_cloud_lut.pro`
- `create_orac_cloud_lut_wrapper.pro`

History:

- `aed975f` (2023-09-03, "1st version of aerosol and cloud separate") **added**
  all four (1968 insertions), splitting the previously combined
  `create_orac_lut.pro` workflow.
- `e5758f6` (2026-05-07, "Create ORAC LUT version 2") **modified** all four.
- `bf29736` and `012efca` did not touch them.
- `b3c8a91` (2026-05-11, V2.1 tip) **deleted** all four.

The deletion is a plain `D` with no rename or copy detected even at
`--find-renames=30% --find-copies-harder`; the blobs were not moved elsewhere
in the tree. The commit message gives no reason, and the commit otherwise only
adjusts environment/path handling (`setup_oraclut_env.sh` gains
`ORAC_LUT_INPUT_ROOT_DIR`, `ORAC_LUT_OUTPUT_ROOT_DIR`, `ORAC_LUT_SOLAR_FILE`,
`ORAC_LUT_TMATRIX_DIR`) plus edits to `makerunfile.pro`, `load_atmstr.pro`,
`load_srfstrarr.pro`, `create_bwgp.pro`, `generate_scattering_properties.pro`
and one driver file.

Decisive evidence that V2.1 is **not** a complete standalone working tree:

- `makerunfile.pro` at `b3c8a91` writes run scripts that call
  `create_orac_lut_aerosol, ...` and `create_orac_lut_cloud, ...`
  (note the reversed word order relative to the deleted files).
- `sentinel-3a_slstr_run` at `b3c8a91` correspondingly contains
  `idl -e "create_orac_lut_aerosol, 0, 'aerosol_a70', 'aerosol_test', '21'"`.
- No file named `create_orac_lut_aerosol.pro` or `create_orac_lut_cloud.pro`
  exists anywhere in the `b3c8a91` tree, and an exhaustive blob grep of that
  tree finds **no** reference to `create_orac_aerosol_lut` /
  `create_orac_cloud_lut` either.

Interpretation: `b3c8a91` should be read as **V2.1 repository infrastructure**
(environment, drivers, readers, scattering, DISORT, output) whose top-level LUT
entry points were renamed in the user's working area but never committed. The
V2.1 tree as committed cannot execute a LUT generation on its own.

The user's local working files (`create_orac_aerosol_lut.pro` etc., dated
2025-11-14) survive in the preserved local tree on scratch and still define the
older `create_orac_aerosol_lut` / `create_orac_cloud_lut` names, driven by a
local `makerunfile.pro` that exports `CREATE_ORAC_AEROSOL_LUT_DRIVER` /
`CREATE_ORAC_CLOUD_LUT_DRIVER` — i.e. the *pre*-V2.1 calling convention.

This provenance question remains **unresolved** and was deliberately not acted
upon. No source was modified.

## 3. Current `create_orac_lut/` state

`create_orac_lut/` is **incomplete**: an earlier restoration ran out of /home
quota part-way. `input_files/inst/` is absent, among other material. It was
deliberately left unrepaired in this task; nothing was deleted from it.

## 4. Preservation locations on scratch

Root: `/network/scratch/grainger/oraclut-reconciliation/`

| Path | Contents |
| --- | --- |
| `preservation/v21_preserved_current_local_create_orac_lut.tar.gz` | 1.1 GB archive of the current local `create_orac_lut` (moved from `docs/`) |
| `preservation/git_preservation/reconciliation/20260913T135357Z/` | `build/`, `create_orac_lut/`, `mie/`, `create_orac_lut.git.tar`, `mie.git.tar` (moved from `docs/git_preservation/…`) |
| `cleanup-hold/` | `disort.dat` (2.1 GB) and nine `CRE*.tmp` LaTeX autosave copies |
| `imported-v21-tree/` | pre-existing V2.1 checkout (1.1 GB) |
| `failed-v21-checkout/` | pre-existing failed checkout remnants (17 MB) |

Verification hashes (identical before and after the move):

```
b45526bc7648fae29dae0571c51c9fbc93142d534eabde1fd17bdffe434638f7  v21_preserved_current_local_create_orac_lut.tar.gz
730c4227ae6f0d415fb0642c78ce41300dfc0c6596a70125045b7b144a0eee46  create_orac_lut.git.tar
5cd83e376a10e3cfb1106b4c5ef6dd5430c92b2ef87b59b750d1fd0354a3f913  mie.git.tar
4bfe92997d68009a0d6a308a8ef4767200c1b5bb1423707717183b927997bd4e  disort.dat
```

Moved directories verified by file count: `build` 18454 files, `create_orac_lut`
5945 files, `mie` 33 files — identical before and after.

Small textual manifests were kept in the repository at
`docs/git_preservation/reconciliation/20260913T135357Z/`:
`backup.sha256`, `filesystem_inventory.txt`, `root_git_state.txt`,
`source_config_sha256.txt`, `timestamp.txt`.

## 5. Quota-cleanup actions

Project size: **10 GB → 1.3 GB** (≈ 8.7 GB of /home recovered).

Deleted as unequivocally reproducible intermediates (root level only):

- LaTeX intermediates: `*.aux`, `*.bbl`, `*.blg`, `*.dvi`, `*.log`, `*.toc`,
  `*.tps`, `missfont.log`;
- Python caches: all `__pycache__/` under `src/` and `tests/`, and
  `.pytest_cache/`;
- compiler objects `build/oraclut/*.o`.

Retained deliberately: all `*.tex`, `*.tcp`, `*.pro`, `*.bak`, `*.eps`, `*.pdf`,
`*.png`, `*.jpg`, `*.pptx`, `*.ps`, driver/input files, `std.out`,
`validation/` (14 MB, includes regression references), `RFM/`, `mie/`,
and the production DLM shared libraries in `build/oraclut/`
(`liboraclut_mie*.so`, `liboraclut_disort*.so`, gfortran/ifort variants) which
the validated legacy IDL system needs.

`create_orac_lut/` was not touched at all.

## 6. `.gitignore`

Narrowed conservatively (edited, **not staged or committed**): removed the
blanket `*.eps` and `*.pptx` patterns, which were hiding manually created
scientific figures with no reproducible source; added an explanatory note and
`missfont.log`. Generated-material patterns (Python caches, objects, LaTeX
intermediates, `luts/`, `disort.dat`, `*.sav`) are unchanged.

## 7. Verification

- Write/delete test: a 200 MB file was created and removed successfully in the
  project root — ordinary file creation is again possible.
- Fast Python tests: **22 passed, 5 failed, 2 deselected**. All five failures
  are `FileNotFoundError` on
  `create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst` and similar
  missing legacy inputs — environmental, not scientific. The previous
  quota-driven temporary-output failure no longer occurs (6 failures → 5).

## 8. Next step (not performed)

Reconstruct `create_orac_lut/` **once**, cleanly, from
`preservation/…/20260913T135357Z/` on scratch, after the cloud/aerosol
provenance decision is settled — specifically whether the canonical tree should
carry the V2.1 infrastructure plus the local 2025-11-14
`create_orac_aerosol_lut.pro` / `create_orac_cloud_lut.pro` implementations, and
under which procedure names (`create_orac_lut_aerosol` as V2.1's
`makerunfile.pro` expects, or the older `create_orac_aerosol_lut`).

## 9. Update (same day, later session): input hierarchy restored, Python workflow usable

- `create_orac_lut/input_files/` was restored from the verified preservation
  copy `preservation/git_preservation/reconciliation/20260913T135357Z/build/current_local_create_orac_lut/input_files/`
  on scratch using `rsync --ignore-existing` (no active file overwritten):
  `inst/` (24 files), `gas/` (3003 files), 19 top-level helper `.pro`/notes
  files, and `baum/` (uncompressed `.nc` tables and small imager tables; the
  `.nc.gz` byte-duplicates were not copied). `modtran3.5-v1.1/` was verified
  identical by checksum and left alone. Result: 6212 files, 1.5 GB.
- One content conflict was found and **not** replaced: the active
  `input_files/driver/sentinel-3a_slstr_aerosol.driver` is the V2.1 (`b3c8a91`)
  version; the preserved local copy is the older pre-V2.1 layout.
- The restored top-level `input_files/*.pro` helpers and
  `input_files/baum/` are not tracked at `b3c8a91` (`baum/` is Git-ignored).
- `python -m oraclut.generate` now accepts `--config <driver>`, resolves
  `--platform/--instrument` and bare `.mm`/`.lut` names from the input
  hierarchy, derives formulation and LUT family from the material as
  `makerunfile_v2.pro` did, and writes legacy-named products under
  `create_orac_lut/luts/`. Representative configurations: `configs/`.
  User guide: `docs/python_lut_generation.md`.
- The four legacy IDL entry points (`create_orac_{cloud,aerosol}_lut[_wrapper].pro`)
  remain absent from the active tree (see §2). `scripts/generate_aerosol_test_reference.sh`
  therefore cannot run until they are restored from the preserved local tree —
  a deliberate user decision, not made here.

## 10. Update: aerosol legacy reference preparation

The §2 provenance finding was confirmed empirically while preparing the legacy
aerosol reference. Running the preserved 2025-11-14 `create_orac_aerosol_lut`
against the V2.1 routines in `create_orac_lut/` halts at once with
`% Keyword TMATRIX_PATH not allowed in call to: GENERATE_SCATTERING_PROPERTIES`
(V2.1 renamed that keyword to `tmatrix_dir`). Five routines in the aerosol call
graph differ between the generations; most consequentially, V2.1's
`load_srfstrarr` **swaps the meaning of `srf_quad=0` and `srf_quad=1`**, so the
same driver would select band integration rather than the monochromatic
calculation the product filename records as `m`. V2.1 also adds an
`Invalid_N_FLAG` correction to the surface-pressure grid in `load_lutstr` and
fixes channel matching in `load_inststr`; `segment.pro` is absent from V2.1.

Note also that IDL searches the current working directory before `IDL_PATH`, so
running from `create_orac_lut/` silently selects the V2.1 routines whatever the
path says. The legacy aerosol runner therefore executes from a clean directory
and verifies the resolved source file of every routine before calculating.

This remains a provenance question for the user to settle; no source was
changed, and nothing was copied out of the preserved tree.
