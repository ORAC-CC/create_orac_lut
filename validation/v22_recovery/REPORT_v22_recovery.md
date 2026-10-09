# V25 preservation, historical version recovery, code versioning and V22 validation

Date: 2026-10-09.  Repository: `/home/g/grainger/project-oraclut`, branch
`main`, remote `origin` = `git@github-orac:ORAC-CC/create_orac_lut.git`.
Host: atmlxint5 (interactive, `nice -n 19` for every calculation; no SLURM
job was submitted).  No V21-V25 lookup table was modified, regenerated or
deleted; the only generated products are compact validation cases under the
session scratch directory and `validation/tmp/`.

This report is the checkpoint required before the main working tree is
changed to the V22 implementation.  **The restoration (section 14) has not
been performed.**

## 1. Original repository state and preservation actions

State at the start (16:1x BST):

- `main` at `1d6aa74` "Fix ORAC LUT NetCDF string attributes and add
  compatibility validation", equal to `origin/main`; `master` (`c182ef6`) is
  the historical IDL branch; no stash, no other worktree, no dangling commits.
- Tags: only `v0.0.0`, `v0.0.1` (IDL era, both at `c182ef6`).
- Working tree: tracked modifications `AGENTS.md` (the owner's "IDL reference
  source - immutable" rule) and `tests/test_baum.py` (one added test);
  untracked `runs/*_v23.run` (40), `scripts/submit_v23_slstr_cloud_luts.sh`,
  `tests/test_nakajima_king.py`, `documents/` (two reference papers and a
  duplicate of `references/data/ocalut_cloudprofile_Cirrostratus.dat`).
  Ignored (by repository policy): `validation/`, `build/`,
  `create_orac_lut/input_files/baum/`, `scripts/submit_v23_modis_cloud_luts.sh`,
  caches.
- SLURM (queried on atmlxint7): no ORAC LUT job queued or running; no
  process writing into the tree.  The NetCDF-attribute repair of the same
  morning was complete, committed (`1d6aa74`) and re-verified (120 of 120
  products compliant).

Preservation commit `f066fdd` "Preserve the complete V25 LUT generator state
before versioning work": the 40 V23 run files, the SLSTR V23 submission
script, the `test_baum.py` test and the AGENTS.md rule (43 files, 1342
insertions).  Left untracked deliberately: `documents/` (reference papers,
35 MB, not source) and `tests/test_nakajima_king.py` (imports the
git-ignored `validation/nakajima_king` module and would fail on a clean
clone); both are unchanged on disk.  Pre-commit verification: focused tests
123 passed; the uncommitted tests 16 passed; full suite 329 passed / 5
pre-existing failures (missing historical files removed in `9d663e9`).

## 2-3. V25 preservation tag and push

- Annotated tag **`lut-code-v25`** -> `f066fdd55530ed326ad476dc50403305fa79032e`
  (tag object `7a37688`), message describing the V25 scientific content and
  that the 40 V25 products were generated at `990d970` (job logs 489345-489384).
- Pushed: `git push origin main` (`1d6aa74..f066fdd`) and `git push origin
  lut-code-v25`; `git ls-remote origin` shows `refs/heads/main` and
  `refs/tags/lut-code-v25^{}` at `f066fdd`.
- Later the same day (section 8): commit `e202fb8` (versioning) and tag
  `lut-code-v25.1`, both pushed and verified remotely.

## 4. Recovered V22 commit and the evidence for it

**`lut-code-v22` = `9d663e9` (9d663e94e6de364d344b848f99acd983c4066fb7, 2026-09-29 17:09 BST).**

Investigation:

- The Python history is short: `d97bbc1` (init, 2026-09-13), `17591c9`
  (IDL V2.1 import), `0a68fad` / `5e5cc0f` (2026-09-14, configuration-driven
  `oraclut.generate` pipeline validated against legacy compact cases),
  `4bb0ac7`-`33f7b19` (2026-09-28, Grid B definitions), `9d663e9` (2026-09-29,
  "Clean repository and consolidate ORAC LUT generators"), then the V24
  (`c6ad545`, `dbcc42c`, `96bebbe`, `51717b8`) and V25 (`7f6b4e2`…`990d970`)
  changes.  Reflogs (local clone of 2026-09-29 and the 2026-09-29 backup tree
  `/network/scratch/grainger/quota-relief/project-oraclut-bloated-backup`,
  itself at `9d663e9`) contain no further commits; `git fsck` finds no
  dangling commits; `/network/scratch/grainger/oraclut-reconciliation`
  holds only V21 IDL material (2026-09-13/14).
- The V22 campaign products (88 compact Python products, 2026-09-20/21, in
  the backup tree's `luts/`; 90 IDL V21 counterparts in
  `create_orac_lut/luts/`) and the V22 EarthCARE MSI products in `ORAC_LUTS`
  (2026-09-20) were generated through `oraclut.master_run` manifests
  (`/network/scratch/grainger/oraclut/logs/manifests/oraclut_master_*.json`)
  by `create_orac_luts.py`, i.e. the IDL-structured generator.  That generator
  (`create_orac_luts.py`, `src/oraclut/idl_mirror/*`, `makerunfile`) first
  appears in Git in `9d663e9`; no earlier commit contains it.  The V22 freeze
  record (`docs/validation_history/python_v22/V22_FREEZE.md`, 2026-09-21) and
  campaign report (`REPORT.md`, `campaign_summary.json`) exist only in the
  backup tree (`docs/` is git-ignored).
- The generator's own docstring and `validation/REPORT_lut_numerics_development.md`
  state that `9d663e9` kept the IDL numerics (fixed 0.001-100 um Mie
  integration, fixed NMom = 1000) and produced the V23 reference LUTs; V23
  differs from V22 only in the grid files selected by run files.

Numerical proof (`compare_regenerated.py`, `results_regenerated_vs_archives.md`):
a detached worktree at `9d663e9` (Fortran kernels rebuilt with GNU Fortran
12.2.0, the compiler recorded for V22) regenerated 13 campaign cases chosen
to cover cloud and aerosol forward models, SRF quadrature 1 and 2, a
dual-view instrument and five instruments:

| case | instrument | model / qm | regenerated vs archived V22 | worst |diff| vs IDL V21 (metadata / optics / operators) |
|---|---|---|---|---|
| 000 | aqua_modis | aerosol a79 / 1 | 51 of 52 variables bitwise; `R_dv` at 2 points differs by 1 float32 ULP (5.7e-14 at 3.3e-7) | 0.00293 / 2.4e-7 / 8.8e-5 (identical to archive) |
| 002 | aqua_modis | cloud stg / 1 | 51/51 bitwise | 0.00293 / 4.3e-4 / 8.4e-4 |
| 004 | earthcare_msi | aerosol a79 / 1 | 52/52 bitwise | 6.7e-4 / 1.8e-7 / 8.7e-5 |
| 006 | earthcare_msi | cloud stg / 1 | 51/51 bitwise | 6.7e-4 / 1.1e-3 / 1.7e-4 |
| 007 | earthcare_msi | cloud stg / 2 | 51/51 bitwise | 6.7e-4 / 3.1e-5 / 1.2e-4 |
| 032 | meteosat-10_seviri | aerosol a79 / 1 | 52/52 bitwise | 9.8e-4 / 3.0e-7 / 1.2e-4 |
| 033 | meteosat-10_seviri | aerosol a79 / 2 | 52/52 bitwise | 9.8e-4 / 6.0e-8 / 1.0e-4 |
| 034 | meteosat-10_seviri | cloud stg / 1 | 51/51 bitwise | 9.8e-4 / 4.0e-4 / 1.3e-4 |
| 035 | meteosat-10_seviri | cloud stg / 2 | 51/51 bitwise | 9.8e-4 / 4.6e-5 / 1.4e-4 |
| 062 | noaa-20_viirs | cloud stg / 1 | 51/51 bitwise | 9.8e-4 / 5.6e-3 / 1.9e-4 |
| 063 | noaa-20_viirs | cloud stg / 2 | 51/51 bitwise | 9.8e-4 / 2.1e-4 / 9.9e-5 |
| 086 | sentinel-3a_slstr | cloud stg / 1 (dual view) | 53/53 bitwise | 4.9e-3 / 3.8e-4 / 1.2e-4 |
| 087 | sentinel-3a_slstr | cloud stg / 2 (dual view) | 53/53 bitwise | 4.9e-3 / 1.5e-4 / 1.1e-4 |

In every case the worst differences against IDL V21 are identical to those
of the archived V22 product, and all are inside the V22 freeze envelopes
(metadata 0.0400390625, optics 0.0088958740234375, operators
0.0008418560028076172; the largest operator difference, 8.4e-4 in case 002,
is the envelope's own value).  The single last-bit difference (case 000,
archived product written by a SLURM compute node on 2026-09-20) is the
documented DISORT compiler/solver sensitivity.

Classification (the four requirements of the task, section 7):

1. *Exact historical source recovery*: the generator, inputs and legacy grids
   of V22 are exactly `9d663e9`; the 2026-09-21 working tree as a whole was
   never committed and is not separately recoverable.  Hence `lut-code-v22`
   and `lut-code-v23` name the same commit.
2. *Execution of the historical code*: a bare worktree at `9d663e9` builds
   its Mie and DISORT kernels and runs; its own test suite gives 227 passed,
   13 failed, 13 skipped, all 13 failures environmental (7 need the
   git-ignored Baum tables, 1 the ignored `validation/generated` directory,
   5 the historical files already missing on `main`).
3. *Agreement with original V22 outputs*: bitwise in 12 of 13 cases, one ULP
   in two values of the thirteenth.
4. *Agreement with IDL V21*: identical to the frozen campaign result; the
   Nakajima-King comparison of the archived full EarthCARE MSI liquid-water
   `stg` V22 product against IDL V21 (0.668 um vs 2.211 um, nadir) gives a
   median displacement 1.49e-5 and maximum 9.83e-5 in reflectance
   (`results/nakajima_king_earthcare_msi_stg_v21_vs_archived_v22.png`).

## 5-6. V23, V24 and the version-to-commit mapping

| product version | tag | commit | evidence |
|---|---|---|---|
| V21 (IDL) | (none; `17591c9` import, historical `b3c8a91`, branch `origin/create-orac-lut-v2.1`) | 17591c901572a7ef9f2db8434dae7ff17a535377 | `AGENTS.md`, import commit |
| V22 | `lut-code-v22` | 9d663e94e6de364d344b848f99acd983c4066fb7 | section 4 |
| V23 | `lut-code-v23` | 9d663e94e6de364d344b848f99acd983c4066fb7 | generator docstring; numerics report; all 40 V23 jobs started 2026-09-30 10:18 to 2026-10-01 22:01 BST, before `c6ad545` (2026-10-02 15:21); the first two products (2026-09-29 15:18) came from the same tree two hours before the commit |
| V24 | `lut-code-v24` | 51717b82dece213d0f93ff1476ff1806983bbbe4 | 40 job logs: "Git revision: 51717b8…", "Source state: clean" |
| V25 | `lut-code-v25` | f066fdd55530ed326ad476dc50403305fa79032e | products generated at `990d970` (40 job logs); `1d6aa74` NetCDF repair; `f066fdd` preservation |
| V25 maintenance | `lut-code-v25.1` | e202fb82a9d2d1005ba5a69761f39cc90b4be574 | this work |

All five tags are annotated and present on `origin` (`git ls-remote origin
'refs/tags/lut-code-v*'`).  No tag was overwritten; none existed before.

## 7. V21/V22 numerical validation summary

See section 4.  Dimensions, optical-depth and effective-radius grids, channel
metadata, particle optical properties (extinction, extinction ratio, single
scattering albedo, asymmetry parameter), size-distribution integration and
quadrature (through the optical properties), DISORT configuration (60
streams, through the operators), radiative-transfer operators, NetCDF
structure, top-of-atmosphere reflectances and the Nakajima-King diagram were
all exercised: the regenerated products are the archived ones to the bit,
and the archived ones are the frozen campaign's.  Tolerances are the original
envelopes; no reference datum was changed.

## 8-9. Version-identification mechanism and banner

Implemented in commit `e202fb8` (tag `lut-code-v25.1`): `CODE_VERSION`,
`src/oraclut/version.py`, `src/oraclut/provenance.py`, hooks in
`create_orac_luts.py` and `src/oraclut/generate.py`, `tests/test_version.py`,
`scripts/tag_lut_code_release.sh`, `LUT_CODE_VERSIONS.md`, an AGENTS.md rule.
Design: see `LUT_CODE_VERSIONS.md` sections 1-2.  Banner printed by a real
run of the current tree (`validation/tmp/version_smoke/smoke.log`, before the
commit, hence MODIFIED):

```
ORAC LUT Generator
------------------------------------------------------------
Generator:         create_orac_luts.py (run file)
LUT version:       25 (requested 25)
Code release:      lut-code-v25.1
Git commit:        f066fdd (f066fdd55530ed326ad476dc50403305fa79032e)
Git describe:      lut-code-v25-0-gf066fdd
Release tag:       tag lut-code-v25.1 not present in this repository
Working tree:      MODIFIED
    modified:      create_orac_luts.py
    modified:      runs/template.run
    modified:      src/oraclut/generate.py
    untracked:     CODE_VERSION
    untracked:     src/oraclut/provenance.py
    untracked:     src/oraclut/version.py
Scientific config: v25
                   V24 numerics (...) plus the V25 cloud vertical temperature profiles ...
Compatible LUTs:   24
------------------------------------------------------------
```

and after the commit and tag (`scripts/tag_lut_code_release.sh` output):

```
Code release:      lut-code-v25.1
Git commit:        e202fb8 (e202fb82a9d2d1005ba5a69761f39cc90b4be574)
Git describe:      lut-code-v25.1-0-ge202fb8
Release tag:       HEAD is the tagged release lut-code-v25.1
Working tree:      CLEAN
```

The banner appears in every job log because the generator prints it; the
SLURM batch script already prints the Git revision and source state as well.

## 10. Example provenance record

`validation/tmp/version_smoke/out/meteosat-10_seviri_m_liquid-water_a01_pstg_v25.provenance.json`
(compact one-channel smoke product; abridged):

```json
{
 "provenance_format": 1,
 "generated_utc": "2026-10-09T15:48:57Z",
 "generator": "create_orac_luts.py",
 "lut_version": 25,
 "output": {"file": "meteosat-10_seviri_m_liquid-water_a01_pstg_v25.nc", "size_bytes": 53811,
            "sha256": "78dc4a4f5775e0c574a7897bbdca03d736e2494e6b9809b1e2f2d6dc8152ed46"},
 "source": {"code_release": "lut-code-v25.1", "lut_version": 25, "compatible_lut_versions": [24],
            "scientific_config": "v25",
            "git": {"commit": "f066fdd55530ed326ad476dc50403305fa79032e", "branch": "main",
                    "describe": "lut-code-v25-0-gf066fdd", "working_tree": "MODIFIED",
                    "tracked_changes": ["create_orac_luts.py", "runs/template.run", "src/oraclut/generate.py"],
                    "untracked_production_files": ["CODE_VERSION", "src/oraclut/provenance.py", "src/oraclut/version.py"]},
            "software": {"python": "3.11.13", "numpy": "2.2.6", "netCDF4": "1.7.2", "libnetcdf": "4.9.2",
                         "hdf5": "1.14.6", "fortran_compiler_version": "GNU Fortran (Debian 12.2.0-14+deb12u1) 12.2.0",
                         "hostname": "atmlxint5.atm.ox.ac.uk"}},
 "configuration_file": {"path": "validation/tmp/version_smoke/smoke_v25.run", "sha256": "466ee5f6…", "text": "..."},
 "configuration": {"forward_model": "cloud", "platform": "meteosat-10", "instrument": "seviri",
                   "instrument_file": "meteosat-10_seviri_v1.inst", "microphysics_file": "liquid-water_stg.mm",
                   "substance": "liquid-water", "particle_model_shortname": "stg",
                   "components": {"types": ["user"], "names": ["water"]}, "channels": [1], "srf_files": ["..."]},
 "input_files": {"instrument": {"path": "create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst", "sha256": "c023a0cc…"},
                 "microphysics": {"path": "create_orac_lut/input_files/microphysics/liquid-water_stg.mm", "sha256": "40f6983e…"},
                 "lut_definition": {"path": "create_orac_lut/input_files/lut/liquid-water-cloud_test.lut", "sha256": "e0233540…"},
                 "atmosphere": {"path": "create_orac_lut/input_files/atm/mls.atm", "sha256": "45244e3c…"}},
 "grid": {"definition_file": "liquid-water-cloud_test.lut",
          "opd": {"n": 2, "spacing": "linear", "values": [...]}, "efr": {"n": 3, ...}, "saz": ..., "soz": ..., "raa": ...},
 "numerics": {"radius_grid": "nested-refinement-1", "mie_radius_limits_um": [0.001, 100.0],
              "mie_size_parameter_step_legacy": 0.4,
              "legendre_expansion": "adaptive (King's criterion on each averaged Mie phase function)", "nmom": null,
              "srf_quadrature": 1, "scattering_cache_reused": false},
 "radiative_transfer": {"solver": "DISORT (create_orac_lut/disort2, oraclut.radiative_transfer.legacy_disort)",
                        "nstreams": 60, "rayleigh": true, "gas_absorption": false, "atmosphere_code": "2",
                        "atmosphere_file": "mls.atm", "cloud_vertical_profile": "isothermal", "cloud_profile_attributes": null},
 "job": {"slurm_job_id": null, "user": "grainger", "working_directory": "/home/g/grainger/project-oraclut"}
}
```

`python -m oraclut.version --verify <file>` on it reports the commit present
in the repository, the release tags at and containing it, `lut_sha256_matches:
true` and `consistent: false` because the tree was MODIFIED at generation time.

## 11. Tests demonstrating that version mismatches are detected

`tests/test_version.py` (19 tests, all passing; full suite 348 passed, 13
skipped, 5 pre-existing failures, 7 slow deselected):

- `test_run_file_with_a_foreign_lut_version_is_refused_before_any_calculation`:
  a run file requesting V22 on the V25 release returns exit 2 with the
  banner, "requests LUT version 22" on stderr, the cloud generator never
  called and no output directory created;
- `test_configuration_entry_point_refuses_a_foreign_lut_version`:
  `oraclut.generate --version 2` exits 2 before writing;
- `test_check_lut_version_accepts_declared_and_compatible_versions_only`;
- `test_modified_tree_is_never_reported_clean`, `test_commits_after_the_release_tag_are_reported`,
  `test_missing_git_metadata_is_reported_not_invented`;
- `test_provenance_verification_detects_a_changed_product` (SHA-256 mismatch
  after the product is edited).

Manual confirmation: `python create_orac_luts.py runs/template.run` with the
old `version = 22` was refused (section 9 banner followed by the mismatch
message, exit 2); `python -m oraclut.version --check-lut-version 22` exits 1.

## 12. Files that would change when restoring V22

`git diff --name-status lut-code-v22 main` (scientific and infrastructure
paths; the 120 `runs/*_v2{3,4,5}.run` and the git-ignored `validation/` are
not listed):

| file | relation to V22 | proposed action |
|---|---|---|
| `create_orac_luts.py` | V24/V25 numerics switches, cloud profile, nmom deprecation, plus the versioning hooks of `e202fb8` | restore the `9d663e9` content, then re-apply only the version banner / check / provenance hooks (adapted: no `RADIUS_GRID`, no cloud profile) |
| `src/oraclut/idl_mirror/create_bwgp.py` | V24 radius refinement, 3.5 r_e limit | restore `9d663e9` |
| `src/oraclut/idl_mirror/generate_scattering_properties.py` | V24 integration limits, adaptive Legendre, radius refinement | restore `9d663e9` |
| `src/oraclut/idl_mirror/legendre_expansion.py` | V24 (new) | remove |
| `src/oraclut/cloud_temperature.py` | V25 (new) | remove |
| `src/oraclut/idl_mirror/write_v2_lut.py` | V25 global attributes / E_md valid range parameters | restore `9d663e9` (the NC_STRING fix lives in `src/oraclut/io/v2.py` and stays) |
| `src/oraclut/master_run.py` | nmom made optional for the adaptive expansion | restore `9d663e9` |
| `tests/test_idl_mirror.py` | adapted to V24 nmom semantics | restore `9d663e9` |
| `tests/test_cloud_temperature.py`, `tests/test_legendre_expansion.py`, `tests/test_radius_grid.py`, `tests/test_size_integration_limits.py` | test V24/V25 science | remove |
| `CODE_VERSION` | — | `code_release = lut-code-v22.1`, `lut_version = 22`, `compatible_lut_versions =` (empty), `scientific_config = v22` |
| `runs/template.run` | version example | `version = 22` |
| `LUT_CODE_VERSIONS.md`, `AGENTS.md`, `runs/V24_PRODUCTION.md`, `runs/V25_PRODUCTION.md`, `LUT_FORMAT_CONTRACT.md` | documentation | keep; add the restoration record |
| `src/oraclut/io/v2.py`, `io/netcdf_c.py`, `io/string_attributes.py`, `repair_string_attributes.py`, `tests/test_netcdf_string_attributes.py` | NetCDF NC_STRING contract (metadata only) | keep |
| `src/oraclut/version.py`, `src/oraclut/provenance.py`, `src/oraclut/generate.py` hooks, `tests/test_version.py`, `scripts/tag_lut_code_release.sh` | versioning infrastructure | keep |
| `scripts/run_oraclut_lut.slurm` | prints Git revision | keep |
| `scripts/submit_v23_slstr_cloud_luts.sh`, `submit_v24_cloud_luts.sh`, `submit_v25_cloud_luts.sh`, `runs/*_v23/24/25.run` | historical production configuration | keep (a V22 release refuses to run them: the version check makes them inert) |
| `references/data/ocalut_cloudprofile_Cirrostratus.dat` | V25 input data | keep as reference data (no code reads it after restoration) or remove; owner's choice |
| `tests/test_baum.py` | owner's added test (independent of version) | keep |
| `create_orac_lut/`, `mie/` | unchanged since `9d663e9` | nothing to do (V21 reference material preserved) |

Not a hybrid: every scientific module returns to its `9d663e9` content; the
retained later files are I/O typing, repair, versioning, provenance and
documentation, none of which enters a calculation.

## 13. Remaining uncertainties and validation limits

- The exact 2026-09-21 V22 working tree is not in Git; identity rests on the
  commit history (no other candidate contains the generator) and on the
  bitwise regeneration of 13 of the 88 campaign cases.  The remaining 75
  cases were not regenerated (the archived set is available and the
  comparison tool accepts any subset).
- One last-bit difference in a near-zero `R_dv` value (case 000) between a
  compute-node product and an atmlxint5 regeneration; this is the known
  DISORT compiler/node sensitivity, not a code difference.
- Ice Baum models (`agg`, `ghm`, `src`) were not regenerated: the Baum
  tables are git-ignored and absent from a bare worktree (they are present in
  the production tree).
- The full-size V22 products (EarthCARE MSI in `ORAC_LUTS`, 2026-09-20) were
  not regenerated (hours per product); the archived full product was used
  for the Nakajima-King comparison.
- `lut-code-v25.1` is one commit after `lut-code-v25` and declares the same
  science; the documentation commit following it (this report's companion
  update of `LUT_CODE_VERSIONS.md`) leaves HEAD one commit after the tag,
  which the banner reports as such.
- Pre-existing test failures (5) and the two untracked files
  (`documents/`, `tests/test_nakajima_king.py`) are unchanged.
- The detached inspection worktree at `9d663e9` remains registered
  (`git worktree list`) in the session scratch directory
  `/tmp/user/27004/claude-27004/.../scratchpad/worktree-v22`; remove it with
  `git worktree remove --force <path>` or `git worktree prune` once it has
  gone.

## 14. Proposed Git operations for the final restoration (not executed)

All on `main`, as ordinary commits; no reset, no force-push, no branch.

```
cd /home/g/grainger/project-oraclut && git switch main && git status      # must be CLEAN
# 1. scientific files back to their V22 content
git checkout lut-code-v22 -- create_orac_luts.py \
    src/oraclut/idl_mirror/create_bwgp.py src/oraclut/idl_mirror/generate_scattering_properties.py \
    src/oraclut/idl_mirror/write_v2_lut.py src/oraclut/master_run.py tests/test_idl_mirror.py
git rm src/oraclut/idl_mirror/legendre_expansion.py src/oraclut/cloud_temperature.py \
    tests/test_cloud_temperature.py tests/test_legendre_expansion.py tests/test_radius_grid.py \
    tests/test_size_integration_limits.py
# 2. re-apply the version-identification hooks to the restored create_orac_luts.py
#    (banner at run(), check_lut_version before "making ...", _record_provenance after
#    _publish_lut; adapted: no RADIUS_GRID / cloud-profile fields) and set
#    CODE_VERSION: code_release = lut-code-v22.1, lut_version = 22,
#    compatible_lut_versions = (empty), scientific_config = v22; runs/template.run: version = 22;
#    record the restoration in LUT_CODE_VERSIONS.md section 3.
# 3. verify before committing
PYTHONPATH=src python -m pytest tests -q -m "not slow"
python create_orac_luts.py validation/v22_recovery/results/case034_meteosat-10_seviri_cloud_stg_qm1.run   # (out_path edited)
python validation/v22_recovery/compare_regenerated.py <out> --cases 034,035,000,086   # bitwise vs archived V22
# 4. commit, tag, push
git add <the files above> && git commit -m "Restore the V22 scientific implementation as the production generator"
scripts/tag_lut_code_release.sh -m <message> --push          # creates and pushes lut-code-v22.1
git ls-remote origin refs/heads/main refs/tags/lut-code-v22.1
```

Expected end state: one branch `main`; V22 science; versioning retained;
V21 reference material untouched; V23-V25 recoverable from their tags and
preserved on GitHub; no LUT modified; no permanent experimental branch.

## 15. Restoration performed (2026-10-09, approved by the owner)

Preparatory commit `05b1978` "Track the validation source, reports and
tabulated results": the blanket `/validation/` ignore was narrowed so that
validation source, reports, grid/run definitions and small tables (173 files,
3.6 MB) are versioned, while NetCDF products, NumPy arrays, figures, job logs
(`validation/slurm` is a symlink to scratch), caches, pytest scratch and the
large forward-model arrays of `validation/v24/results/disort_options/` stay on
disk; `tests/test_nakajima_king.py` is now tracked and passes (6 tests).  The
V25 release tag was not touched.

Restoration, exactly as in section 14, on `main` by ordinary commits:

- `git checkout lut-code-v22 -- create_orac_luts.py src/oraclut/idl_mirror/create_bwgp.py
  src/oraclut/idl_mirror/generate_scattering_properties.py src/oraclut/idl_mirror/write_v2_lut.py
  src/oraclut/master_run.py tests/test_idl_mirror.py`;
- `git rm src/oraclut/idl_mirror/legendre_expansion.py src/oraclut/cloud_temperature.py
  tests/test_cloud_temperature.py tests/test_legendre_expansion.py tests/test_radius_grid.py
  tests/test_size_integration_limits.py`;
- the version banner, the LUT-version check and the provenance hook were
  re-applied to the restored `create_orac_luts.py`: `git diff lut-code-v22 --
  create_orac_luts.py` adds 90 lines (two imports, `_record_provenance`, its
  two calls after `_publish_lut`, the banner and check in `run()`) and changes
  only the three lines of the `except` clause in `main()`; no scientific
  statement differs from `9d663e9`.  The provenance record describes the V22
  numerics (legacy lattice, fixed `nmom`, isothermal cloud);
- `CODE_VERSION`: `code_release = lut-code-v22.1`, `lut_version = 22`,
  `compatible_lut_versions = 21` (legacy-comparison run files only; justified
  in `LUT_CODE_VERSIONS.md` section 2), `scientific_config = v22`;
- `runs/template.run` restored to its V22 form (`nmom = 1000`, `version = 22`)
  with the version note; `tests/test_idl_mirror.py` (V22 version) iterates
  only the run files whose version the release accepts, because the V24/V25
  run files remain in the tree as history and are not readable by the V22
  reader (`nmom` missing, `cloud_vertical_profile` unknown); one assertion in
  `tests/test_version.py` chose a foreign version dynamically instead of the
  literal 22;
- kept unchanged: `references/data/ocalut_cloudprofile_Cirrostratus.dat`
  (historical V25 reference material), `runs/*_v23/24/25.run`, the V23-V25
  submission scripts, `runs/V24_PRODUCTION.md`, `runs/V25_PRODUCTION.md`, the
  NetCDF NC_STRING writer/repair (`src/oraclut/io/`), the versioning and
  provenance modules, `src/oraclut/generate.py` hooks, `create_orac_lut/`, `mie/`.

Verification of the restored tree (before the commit, hence `Working tree:
MODIFIED` in the banners and provenance):

| check | result |
|---|---|
| banner | `LUT version: 22`, `Code release: lut-code-v22.1`, `Scientific config: v22` |
| version check | V23, V24 and V25 run files refused: V23 by version ("requests LUT version 23 ... produces LUT version 22 (compatible: 21)"), V24/V25 already by the V22 run-file reader (`nmom` missing / `cloud_vertical_profile` unknown) |
| campaign cases regenerated from the restored tree (034, 035, 000, 004, 062, 086) vs archived V22 | bitwise identical in every variable (5 cases); case 000: 51/52, the same two one-ULP `R_dv` values as before; worst differences vs IDL V21 identical to the archived products; `results_restored_main_vs_archives.md` |
| restored tree vs `9d663e9` worktree regenerations | bitwise identical data and dimensions in all six cases |
| NC_STRING repair effect (requirement 2) | attribute names, order and values identical to the archived V22 products; the only differences are 84-88 NC_CHAR -> NC_STRING type changes per file; dimensions, coordinate grids and every data value unchanged; `python -m oraclut.repair_string_attributes --check`: 6 compliant |
| provenance | one `.provenance.json` per product; `--verify` reports the commit, `lut_sha256_matches: true`, V22 numerics recorded |
| tests | full non-slow suite: 284 passed, 13 skipped, 5 failed (the pre-existing missing historical files: `docs/idl_python_port_status.md`, `validation/reference_cases/...`, `validation/diagnostics/aerosol_legacy_provenance.pro`); slow V21-equivalence tests: 2 passed, 3 skipped (their reference products are not present on disk) |

The V22.1 commit, tag and push are recorded at the end of this report.

Restoration commit: `d9a7e8b` "Restore the V22 scientific implementation as
the production generator" on `main`.  Release tag `lut-code-v22.1` is placed
on the follow-up commit that records this hash (so that the tagged state
contains its own documentation); it was created with `git tag -a` rather than
`scripts/tag_lut_code_release.sh` because the owner's concurrent, unrelated
edits to `create_orac_lut/makerunfile.pro` and `create_orac_lut/terra_modis_run`
(20:13 and 20:19 BST, not part of this work, left untouched and uncommitted)
make the script's clean-tree check refuse; the committed content itself is
complete and the banner reports those two files as the only modifications.
