# NetCDF text-attribute repair of the V23, V24 and V25 ORAC LUTs

Date: 2026-10-09.  Host: atmlxint5 (repair ran interactively at `nice -n 19`;
it is file I/O only, about 22 minutes for 120 files).  No radiative-transfer,
scattering or lookup-table numerical calculation was performed or repeated;
only attribute metadata was changed.

## 1. Root cause and exact code locations

**Symptom.** Production ORAC stopped on the new tables with
`ERROR: ncdf_get_string_att(): NetCDF: Attempt to convert between text & numbers`.

**ORAC reader contract** (ORAC-CC/orac `master`, commit
`eb9a1233c1136ef95fa0c37f650fe519e9134ed4`, cloned read-only into the session
scratch directory on 2026-10-09; `git status` clean; nothing copied into this
repository):

- `src/read_sad_lut.F90`, `Read_NCDF_SAD_LUT`, lines 979-983 call
  `sad_dimension_read_nc` for `optical_depth`, `effective_radius`,
  `satellite_zenith`, `solar_zenith`, `relative_azimuth`;
  `sad_dimension_read_nc` (lines 848-903) reads `name:spacing` with
  `ncdf_get_string_att` (line 864) and expects `logarithmic`,
  `uneven_logarithmic`, `linear`, `uneven_linear` or `unknown`.
- `common/orac_ncdf.F90`, `ncdf_get_string_att` (lines 257-295) calls the C
  routine `common/nc_get_string_att.c` (`nc_open`, `nc_inq_varid`,
  `nc_inq_attlen`, `nc_get_att_string` at line 42, first element only); on any
  error `nc_check` (lines 17-23) prints the message above and `exit(1)`.
- `nc_get_att_string` succeeds only for an **NC_STRING** attribute; for an
  **NC_CHAR** attribute netcdf-c returns NC_ECHAR (-56), "Attempt to convert
  between text & numbers".
- No other text attribute is read: the `units` reads in
  `common/ncdf_read_template.inc` (line 43), `ncdf_read_field.inc` (line 69)
  and `ncdf_read_packed_field.inc` (line 45) are under `#ifdef DEBUG` and
  status-checked; `_FillValue`, `scale_factor`, `add_offset`, `valid_min`,
  `valid_max` are read numerically (`common/ncdf_open_field.inc` lines 39-43).

**Attribute types in the files** (`nc_inq_att` through
`src/oraclut/io/netcdf_c.py`):

| file | text attributes | numeric attributes |
|---|---|---|
| `aqua_modis_m_liquid-water_a01_p240_v11.nc` (IDL, working) | 86 NC_STRING, 0 NC_CHAR | 30 NC_FLOAT, 5 NC_INT, 4 NC_SHORT |
| `aqua_modis_m_liquid-water_a01_pmodis_v19.nc` (IDL) | 86 NC_STRING | same |
| `earthcare_msi_m_liquid-water_a01_p240_v21.nc` (IDL) | 84 NC_STRING | 30 NC_FLOAT, 4 NC_INT, 4 NC_SHORT |
| `earthcare_msi_m_liquid-water_a01_p240_v22.nc` (Python) | **84 NC_CHAR** | same as V21 |
| `aqua_modis_m_water-ice_a01_pagg_v23.nc` (Python) | **86 NC_CHAR** | 30 NC_FLOAT, 5 NC_INT, 4 NC_SHORT |
| `aqua_modis_m_liquid-water_a01_p240_v24.nc` (Python) | **86 NC_CHAR** | same |
| `aqua_modis_m_liquid-water_a01_p240_v25.nc` (Python) | **88 NC_CHAR** (2 global) | 33 NC_FLOAT, 6 NC_INT, 4 NC_SHORT |

The hypothesis that V23 and V24 are affected as well as V25 is confirmed: every
Python-written generation (V22 to V25) has every text attribute as NC_CHAR,
including the five axis `spacing` attributes ORAC reads.  The IDL-written
tables are NC_STRING throughout.  Attribute names, values and order are
otherwise identical between V11 and V23-V25 for the shared attributes.

**Cause in the generator.** `src/oraclut/io/v2.py`, `write_v2_lut`, wrote all
attributes with the netCDF4 package's `setncattr` (formerly lines 89,
variable attributes, and 91, global attributes).  `setncattr` stores a Python
`str` as NC_CHAR unless `Dataset.set_ncstring_attrs(True)` was called; IDL's
`NCDF_ATTPUT` stores strings as NC_STRING.  Both production writers reach this
primitive: `src/oraclut/idl_mirror/write_v2_lut.py` line 309 (used by
`create_orac_luts.py`, the V23-V25 path) and `src/oraclut/pipeline.py` line
1006.  The regression entered with the Python writer (first product V22,
EarthCARE MSI, 2026-09-20) and was not caught because the Python reader and
`ncdump` show NC_CHAR and NC_STRING text identically.

## 2. Inventory of V23, V24 and V25 LUTs

Location: `/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS` (gluster
volume, no filesystem snapshots).  120 files, 98.6 GB, all NetCDF-4 (HDF5),
every variable contiguous and uncompressed.  No other copy of these products
exists under the project's directories or scratch.  SLURM (queried on
atmlxint7) showed no ORAC LUT job running or pending; every V25 job log in
`validation/slurm/` (jobs 489345-489384) ends with the completion line, and the
youngest product was last modified 06:09 on 2026-10-09, six hours before the
repair.  `lsof` showed no process holding any file.  Two V23 files carry the
group `subdept-aopp-sags` instead of `aopp-eodg`; the repair preserves mode
and group.

| version | platform / instrument | liquid-water (6 models) | water-ice (4 models) | status before | status after |
|---|---|---|---|---|---|
| V23 | Aqua MODIS | 6 x 985.7 MB | 4 x 1191.1 MB | 10 non-compliant | 10 repaired |
| V23 | Terra MODIS | 6 x 985.7 MB | 4 x 1191.1 MB | 10 non-compliant | 10 repaired |
| V23 | Sentinel-3A SLSTR | 6 x 530.5 MB | 4 x 641.1 MB | 10 non-compliant | 10 repaired |
| V23 | Sentinel-3B SLSTR | 6 x 530.5 MB | 4 x 641.1 MB | 10 non-compliant | 10 repaired |
| V24 | Aqua MODIS | 6 x 985.7 MB | 4 x 1191.1 MB | 10 non-compliant | 10 repaired |
| V24 | Terra MODIS | 6 x 985.7 MB | 4 x 1191.1 MB | 10 non-compliant | 10 repaired |
| V24 | Sentinel-3A SLSTR | 6 x 530.5 MB | 4 x 641.1 MB | 10 non-compliant | 10 repaired |
| V24 | Sentinel-3B SLSTR | 6 x 530.5 MB | 4 x 641.1 MB | 10 non-compliant | 10 repaired |
| V25 | Aqua MODIS | 6 x 985.7 MB | 4 x 1191.1 MB | 10 non-compliant | 10 repaired |
| V25 | Terra MODIS | 6 x 985.7 MB | 4 x 1191.1 MB | 10 non-compliant | 10 repaired |
| V25 | Sentinel-3A SLSTR | 6 x 530.5 MB | 4 x 641.1 MB | 10 non-compliant | 10 repaired (1 as the representative) |
| V25 | Sentinel-3B SLSTR | 6 x 530.5 MB | 4 x 641.1 MB | 10 non-compliant | 10 repaired |

Per-file inventory (path, version, size, original mtime, status, number of
attributes converted, number of axis `spacing` conversions, text-attribute
types before repair): `inventory_v23_v24_v25.tsv` (state as audited at the
start of the batch; the representative file was already repaired by then) and
`inventory_v23_v24_v25_after_repair.tsv` (release check after the batch).  The
full per-file table is appended in section 9.

Other generations in the same directory, audited read-only
(`inventory_other_versions_dry_run.tsv`): 846 files of V00-V21 and V99 are
compliant; the **29 V22 EarthCARE MSI products are non-compliant** in the same
way (NC_CHAR) and were not repaired because they are outside the requested
scope.  They can be repaired with the same command.

## 3. Files repaired, already compliant, pending or inaccessible

- Repaired: 120 of 120 (40 V23, 40 V24, 40 V25); re-checked read-only after the commit on 2026-10-09: 120 compliant (`logs/release_check_2026-10-09_recheck.log`).  The representative file
  `sentinel-3a_slstr_m_liquid-water_a01_p240_v25.nc` was repaired first on its
  own, then the remaining 119 in one batch (`logs/batch_repair.log`).
- Already compliant before the campaign: none of the 120.
- Skipped, failed, pending or inaccessible: none.  No file was being written.
- Not in scope and left unchanged: 29 V22 EarthCARE MSI products (non-compliant).
- Nathalie's files and directories were not located and not touched; only the
  files in the shared `ORAC_LUTS` directory were modified.

## 4. Metadata changes made

For every file, every NC_CHAR text attribute was rewritten as an NC_STRING
attribute of length 1 with the identical UTF-8 value: 86 (V23/V24 MODIS and
SLSTR liquid, 20 files each), 90 (V23/V24 ice, 20 files each), 88, 92 or 96
(V25, which added global `cloud_vertical_profile*` attributes) per file;
10,704 attributes in all, of which 600 are the axis `spacing` attributes ORAC
reads (5 per file).  Attribute names, values, order, and all numeric
attributes (`valid_range`, `_FillValue`, the V25 global numeric attributes)
are unchanged.  No attribute was added or removed, no history or provenance
attribute was written.  File sizes grew by 0 to 8192 bytes (HDF5 heap space
for the variable-length strings); the file modification time is the repair
time (the original mtime is recorded in the inventory and in each JSON log).
Ownership, mode and group are those of the original.  `ncdump -h` of a
repaired file differs from the original only by the `string` keyword before
each text attribute (checked on the representative file; headers otherwise
identical).

Recovery copies: each original was hard-linked (no data copy) into
`ORAC_LUTS/originals_before_nc_string_repair_20261009/` before the atomic
replace; 120 files, 92 GiB (98.6 GB).  They are byte-identical to the
pre-repair files (same inodes).  They should be deleted once Nathalie has run
ORAC successfully with the repaired tables; until then they are the fallback.

## 5. Data-preservation verification results

Method (`src/oraclut/io/string_attributes.py`, `verify_preserved`, run on the
working copy against the untouched original before the replace, with the
NetCDF C library for attributes and raw-byte comparison for data in slabs of
at most 64 MiB along the leading dimension; no tolerance anywhere):

| check | result over 120 files |
|---|---|
| file format, data model | identical (NETCDF4 / HDF5) |
| dimensions: names, order, lengths, unlimited flag | identical, 1,680 dimensions (14 per file) |
| variables: names, order, dtype, dimensions, shape, chunking, filters, byte order | identical, 6,360 variables (53 or 54 per file) |
| attributes: name, order, external type, length, raw value | identical except the listed conversions, 15,452 attributes |
| converted attributes | 10,704, each NC_CHAR -> NC_STRING[1] with equal UTF-8 value |
| data values, raw bytes (NaN payloads, signed zeros, fill values included) | identical, 98,549,551,020 bytes compared, no difference |
| ORAC-style `nc_get_att_string` read of each axis `spacing` in the copy | 600 reads, all equal to the original values |
| post-replace audit | 120 compliant |

(Totals: 119 batch files from `results_repair_summary.tsv` plus the
representative file from `logs/representative_repair.log`: 54 variables, 134
attributes, 530,477,906 data bytes.)  Independent checks on the
representative file: `ncdump -h` headers identical apart from the `string`
keyword; `h5diff -c` between the retained original and the repaired file
reports no dataset differences (exit 0; the 92 converted attributes are listed
as "not comparable" because their HDF5 types differ, as expected).  The unit
tests confirm on a synthetic file that a NaN with payload `0x7fc12345`, a
negative zero and the `-999` fill value survive the repair bit for bit.

## 6. ORAC compatibility test results

No ORAC retrieval build or executable is accessible from this account (search
of the home, group, scratch and `/network/aopp` trees), so a full ORAC
integration run was not possible.  The closest test was done instead: ORAC's
own, unmodified `common/nc_get_string_att.c` (the routine that failed in
production) was compiled against the project's NetCDF library with a 10-line
`main` (`orac_c_reader/main.c`, `orac_c_reader/build_and_run.sh`) and run over
all five axes of every file:

| files | reads | result |
|---|---|---|
| 120 repaired V23/V24/V25 products | 600 | all succeed; values `uneven_logarithmic`, `uneven_linear` x3, `linear` (`orac_c_reader/results_v23_v24_v25_repaired.tsv`) |
| 120 retained originals | 600 | all fail with `ERROR: ncdf_get_string_att(): NetCDF: Attempt to convert between text & numbers` (`results_v23_v24_v25_originals.tsv`) |
| V11 reference `aqua_modis_m_liquid-water_a01_p240_v11.nc` | 1 | succeeds, `uneven_logarithmic` |

The same read is reproduced through ctypes in
`tests/test_netcdf_string_attributes.py` (`read_string_attribute_like_orac`).
No other LUT attribute is read as text by ORAC (section 1), and the numeric
attributes ORAC reads are unchanged, so no further metadata incompatibility is
expected; an end-to-end ORAC run on one repaired table remains to be done by
Nathalie (section 8).

## 7. Writer changes and regression tests

Writer (`src/oraclut/io/v2.py`): new `put_attribute` writes `str` (and
sequences or NumPy arrays of `str`) with `setncattr_string`
(`nc_put_att_string`), other values with `setncattr` and their NumPy dtype,
and refuses `bytes`; both attribute loops use it; after closing the file the
writer calls `check_text_attributes_are_strings` (C-library audit) and raises
if any text attribute is NC_CHAR.  Nothing in the scientific calculation,
grids, variables, values or version numbering changed; the two production
callers (`idl_mirror/write_v2_lut.py`, `pipeline.py`) are unchanged and pick
up the fix through the primitive.  Existing tests of the writer and generator
pass (`tests/test_lut_reader.py`, `tests/test_generator.py`,
`tests/test_port_defect_regressions.py`, `tests/test_master_run.py`).

New modules: `src/oraclut/io/netcdf_c.py` (ctypes binding: `nc_inq_att`,
`nc_get_att_text`, `nc_get_att_string`, `nc_get_att`, test-only writers),
`src/oraclut/io/string_attributes.py` (audit, conversion, verification,
repair, CLI), `src/oraclut/repair_string_attributes.py` (entry point).

Tests (`tests/test_netcdf_string_attributes.py`, 21 tests, all passing):
writer emits NC_STRING for every text attribute and ORAC-style reads succeed;
a plain `setncattr` file reproduces the production failure (NC_ECHAR);
the writer's self-check fails the write when text becomes NC_CHAR (the
"must fail if NC_STRING becomes NC_CHAR" requirement); explicit typing of
`put_attribute`; audit contents; repair conversion with preservation of
attribute order, storage layout, fill value, NaN payload and signed zero;
idempotence; already-correct file untouched; `--attributes` restriction;
multiple affected attributes across variables and the global scope;
unexpected attribute types (invalid UTF-8 NC_CHAR, user-defined vlen type)
block the repair without modification; failed conversion and failed
verification leave the original untouched and remove the working copy;
tampered data, numeric attribute and converted value are detected; recent
files and fresh working copies skipped, stale working copies removed;
original retained as a hard link and never overwritten; bounded slab
iteration and verification in many slabs; CLI dry-run, check, inventory,
log records and repair; the IDL-written V11 reference is compliant.

Full suite (2026-10-09, after the Kerberos ticket was renewed): 329 passed,
13 skipped, 5 failed, 7 deselected (slow).  The 5 failures pre-date this work
and are unrelated: they reference files removed in commit 9d663e9
(`docs/idl_python_port_status.md`, `validation/reference_cases/...`,
`validation/diagnostics/aerosol_legacy_provenance.pro`).  An earlier run with
an expired Kerberos ticket also showed 17 `tests/test_submit_wrapper.py` and
`tests/test_legacy_runner.py` failures ("Key has expired"); they pass with a
valid ticket.  The focused tests (`test_netcdf_string_attributes.py`,
`test_lut_reader.py`, `test_generator.py`, `test_port_defect_regressions.py`,
`test_master_run.py`, `test_submit_wrapper.py`) give 123 passed.

Documentation: `LUT_FORMAT_CONTRACT.md` (contract, root cause, writer fix,
utility usage, tests, release check), a "LUT product format contract" rule in
`AGENTS.md`, and notes in `runs/V24_PRODUCTION.md` and `runs/V25_PRODUCTION.md`.

## 8. Remaining actions before Nathalie can use the tables

1. Source changes committed and pushed: commit `1d6aa74` ("Fix ORAC LUT
   NetCDF string attributes and add compatibility validation") on `main`,
   pushed to `origin/main` on 2026-10-09.  Files: `AGENTS.md` (contract
   rule only), `LUT_FORMAT_CONTRACT.md`, `runs/V24_PRODUCTION.md`,
   `runs/V25_PRODUCTION.md`, `src/oraclut/io/v2.py`,
   `src/oraclut/io/netcdf_c.py`, `src/oraclut/io/string_attributes.py`,
   `src/oraclut/repair_string_attributes.py`,
   `tests/test_netcdf_string_attributes.py`.  This directory
   (`validation/`) is excluded from Git by `.gitignore`, as all validation
   material is.  Pre-existing uncommitted items in the tree (the IDL-reference
   section of `AGENTS.md`, `tests/test_baum.py`, `tests/test_nakajima_king.py`,
   `scripts/submit_v23_slstr_cloud_luts.sh`, `runs/*_v23.run`, `documents/`)
   were left as they were.  Future LUT jobs submitted from a clean tree at or
   after this commit use the corrected writer.
2. **Run production ORAC on one repaired table** (for example the Aqua MODIS
   V25 liquid-water `p240` table).  The attribute read that failed now
   succeeds with ORAC's own C routine, but only a retrieval run proves the
   whole chain.
3. **Delete the retained originals** once step 2 succeeds:
   `rm -r /network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/originals_before_nc_string_repair_20261009`
   (92 GiB).  If ORAC still fails for a reason traced to the repair, the
   originals can be put back by renaming.
4. **Optionally repair the 29 V22 EarthCARE MSI products** with
   `PYTHONPATH=src python -m oraclut.repair_string_attributes --original-dir ... ORAC_LUTS/earthcare_msi_*_v22.nc`.
5. If Nathalie has her own copies of V23-V25 tables (for example in an ORAC
   `sad_dir`), they must be replaced by the repaired files or repaired with the
   same utility; copies made before 2026-10-09 12:19 are the old NC_CHAR files.
6. Before any future release, run
   `PYTHONPATH=src python -m oraclut.repair_string_attributes --check <products>`.

## 9. Per-file inventory

Columns: file, version, size (bytes), original mtime, status at the start of
the batch, attributes converted, axis `spacing` conversions, verification
(variables / attributes / data bytes compared), retained original.

| file | ver | size | original mtime | status | conv | spacing | verified vars/attrs/bytes | retained |
|---|---|---|---|---|---|---|---|---|
| `aqua_modis_m_liquid-water_a01_p240_v23.nc` | 23 | 985,726,652 | 2026-10-01T09:19:38 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p240_v24.nc` | 24 | 985,726,652 | 2026-10-05T03:21:48 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p240_v25.nc` | 25 | 985,727,736 | 2026-10-09T06:01:22 | repaired | 88 | 5 | 52 / 131 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p253_v23.nc` | 23 | 985,726,652 | 2026-10-01T02:51:43 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p253_v24.nc` | 24 | 985,726,652 | 2026-10-05T02:31:26 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p253_v25.nc` | 25 | 985,727,736 | 2026-10-09T05:55:05 | repaired | 88 | 5 | 52 / 131 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p263_v23.nc` | 23 | 985,726,652 | 2026-10-01T12:01:53 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p263_v24.nc` | 24 | 985,726,652 | 2026-10-05T03:53:06 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p263_v25.nc` | 25 | 985,727,736 | 2026-10-09T05:55:33 | repaired | 88 | 5 | 52 / 131 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p273_v23.nc` | 23 | 985,726,652 | 2026-10-01T12:11:50 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p273_v24.nc` | 24 | 985,726,652 | 2026-10-05T02:45:22 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_p273_v25.nc` | 25 | 985,727,736 | 2026-10-09T05:38:49 | repaired | 88 | 5 | 52 / 131 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_pold_v23.nc` | 23 | 985,726,652 | 2026-09-29T15:19:52 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_pold_v24.nc` | 24 | 985,726,652 | 2026-10-05T01:26:15 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_pold_v25.nc` | 25 | 985,727,736 | 2026-10-08T23:29:53 | repaired | 88 | 5 | 52 / 131 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_pstg_v23.nc` | 23 | 985,726,652 | 2026-09-29T15:18:33 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_pstg_v24.nc` | 24 | 985,726,652 | 2026-10-05T03:26:46 | repaired | 86 | 5 | 52 / 125 / 985,669,554 | yes |
| `aqua_modis_m_liquid-water_a01_pstg_v25.nc` | 25 | 985,727,736 | 2026-10-09T06:09:19 | repaired | 88 | 5 | 52 / 131 / 985,669,554 | yes |
| `aqua_modis_m_water-ice_a01_pagg_v23.nc` | 23 | 1,191,073,532 | 2026-10-02T12:44:36 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_pagg_v24.nc` | 24 | 1,191,073,532 | 2026-10-04T15:51:32 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_pagg_v25.nc` | 25 | 1,191,077,260 | 2026-10-08T13:00:56 | repaired | 92 | 5 | 52 / 133 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_pghm_v23.nc` | 23 | 1,191,073,532 | 2026-10-02T12:29:54 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_pghm_v24.nc` | 24 | 1,191,073,532 | 2026-10-04T12:57:33 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_pghm_v25.nc` | 25 | 1,191,077,260 | 2026-10-08T13:03:11 | repaired | 92 | 5 | 52 / 133 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_psph_v23.nc` | 23 | 1,191,073,532 | 2026-10-01T16:34:24 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_psph_v24.nc` | 24 | 1,191,073,532 | 2026-10-05T08:16:06 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_psph_v25.nc` | 25 | 1,191,077,260 | 2026-10-09T00:37:38 | repaired | 92 | 5 | 52 / 133 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_psrc_v23.nc` | 23 | 1,191,073,532 | 2026-10-02T12:27:58 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_psrc_v24.nc` | 24 | 1,191,073,532 | 2026-10-04T17:19:18 | repaired | 86 | 5 | 52 / 125 / 1,191,016,474 | yes |
| `aqua_modis_m_water-ice_a01_psrc_v25.nc` | 25 | 1,191,077,260 | 2026-10-08T11:26:17 | repaired | 92 | 5 | 52 / 133 / 1,191,016,474 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p240_v23.nc` | 23 | 530,534,768 | 2026-10-02T02:37:11 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p240_v24.nc` | 24 | 530,534,768 | 2026-10-04T02:41:41 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p240_v25.nc` | 25 | 530,543,290 | 2026-10-07T22:24:53 | repaired (representative, 12:19) | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p253_v23.nc` | 23 | 530,534,768 | 2026-10-02T02:39:41 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p253_v24.nc` | 24 | 530,534,768 | 2026-10-04T07:06:01 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p253_v25.nc` | 25 | 530,535,779 | 2026-10-08T00:43:41 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p263_v23.nc` | 23 | 530,534,768 | 2026-10-02T02:39:02 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p263_v24.nc` | 24 | 530,534,768 | 2026-10-04T06:57:12 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p263_v25.nc` | 25 | 530,535,779 | 2026-10-08T00:47:14 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p273_v23.nc` | 23 | 530,534,768 | 2026-10-02T02:29:37 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p273_v24.nc` | 24 | 530,534,768 | 2026-10-04T00:19:49 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_p273_v25.nc` | 25 | 530,535,779 | 2026-10-07T22:39:09 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_pold_v23.nc` | 23 | 530,534,768 | 2026-10-02T02:42:44 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_pold_v24.nc` | 24 | 530,534,768 | 2026-10-04T02:44:11 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_pold_v25.nc` | 25 | 530,535,779 | 2026-10-07T22:07:16 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_pstg_v23.nc` | 23 | 530,534,768 | 2026-10-02T02:40:36 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_pstg_v24.nc` | 24 | 530,534,768 | 2026-10-04T02:47:21 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_liquid-water_a01_pstg_v25.nc` | 25 | 530,535,779 | 2026-10-07T21:48:44 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_pagg_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:24:48 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_pagg_v24.nc` | 24 | 641,051,375 | 2026-10-03T22:56:05 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_pagg_v25.nc` | 25 | 641,055,103 | 2026-10-07T22:14:54 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_pghm_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:25:49 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_pghm_v24.nc` | 24 | 641,051,375 | 2026-10-04T03:04:12 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_pghm_v25.nc` | 25 | 641,055,103 | 2026-10-07T22:28:23 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_psph_v23.nc` | 23 | 641,051,375 | 2026-10-02T03:23:41 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_psph_v24.nc` | 24 | 641,051,375 | 2026-10-04T07:11:37 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_psph_v25.nc` | 25 | 641,055,103 | 2026-10-08T02:33:46 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_psrc_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:22:50 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_psrc_v24.nc` | 24 | 641,051,375 | 2026-10-03T22:46:16 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3a_slstr_m_water-ice_a01_psrc_v25.nc` | 25 | 641,055,103 | 2026-10-07T22:28:23 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p240_v23.nc` | 23 | 530,534,768 | 2026-10-01T22:34:51 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p240_v24.nc` | 24 | 530,534,768 | 2026-10-04T05:34:48 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p240_v25.nc` | 25 | 530,535,779 | 2026-10-08T05:55:56 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p253_v23.nc` | 23 | 530,534,768 | 2026-10-01T22:32:59 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p253_v24.nc` | 24 | 530,534,768 | 2026-10-04T07:10:10 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p253_v25.nc` | 25 | 530,535,779 | 2026-10-08T06:44:10 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p263_v23.nc` | 23 | 530,534,768 | 2026-10-01T22:29:27 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p263_v24.nc` | 24 | 530,534,768 | 2026-10-04T07:11:04 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p263_v25.nc` | 25 | 530,535,779 | 2026-10-08T06:44:35 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p273_v23.nc` | 23 | 530,534,768 | 2026-10-01T22:38:42 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p273_v24.nc` | 24 | 530,534,768 | 2026-10-04T07:11:14 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_p273_v25.nc` | 25 | 530,535,779 | 2026-10-08T06:51:05 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_pold_v23.nc` | 23 | 530,534,768 | 2026-10-01T22:41:18 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_pold_v24.nc` | 24 | 530,534,768 | 2026-10-04T01:05:30 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_pold_v25.nc` | 25 | 530,535,779 | 2026-10-08T04:36:03 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_pstg_v23.nc` | 23 | 530,534,768 | 2026-10-01T22:34:12 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_pstg_v24.nc` | 24 | 530,534,768 | 2026-10-04T01:44:17 | repaired | 90 | 5 | 54 / 128 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_liquid-water_a01_pstg_v25.nc` | 25 | 530,535,779 | 2026-10-08T06:10:31 | repaired | 92 | 5 | 54 / 134 / 530,477,906 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_pagg_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:20:17 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_pagg_v24.nc` | 24 | 641,051,375 | 2026-10-03T23:48:43 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_pagg_v25.nc` | 25 | 641,055,103 | 2026-10-08T05:15:32 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_pghm_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:24:29 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_pghm_v24.nc` | 24 | 641,051,375 | 2026-10-04T04:03:07 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_pghm_v25.nc` | 25 | 641,055,103 | 2026-10-08T05:13:45 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_psph_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:31:32 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_psph_v24.nc` | 24 | 641,051,375 | 2026-10-04T02:58:32 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_psph_v25.nc` | 25 | 641,055,103 | 2026-10-08T07:56:55 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_psrc_v23.nc` | 23 | 641,051,375 | 2026-10-01T23:21:59 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_psrc_v24.nc` | 24 | 641,051,375 | 2026-10-04T04:29:03 | repaired | 90 | 5 | 54 / 128 / 640,993,626 | yes |
| `sentinel-3b_slstr_m_water-ice_a01_psrc_v25.nc` | 25 | 641,055,103 | 2026-10-08T05:11:29 | repaired | 96 | 5 | 54 / 136 / 640,993,626 | yes |
| `terra_modis_m_liquid-water_a01_p240_v23.nc` | 23 | 985,726,654 | 2026-10-02T05:44:45 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p240_v24.nc` | 24 | 985,726,654 | 2026-10-05T05:46:04 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p240_v25.nc` | 25 | 985,727,738 | 2026-10-08T21:19:19 | repaired | 88 | 5 | 52 / 131 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p253_v23.nc` | 23 | 985,726,654 | 2026-10-02T05:44:22 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p253_v24.nc` | 24 | 985,726,654 | 2026-10-05T03:00:50 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p253_v25.nc` | 25 | 985,727,738 | 2026-10-08T19:14:58 | repaired | 88 | 5 | 52 / 131 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p263_v23.nc` | 23 | 985,726,654 | 2026-10-02T05:41:17 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p263_v24.nc` | 24 | 985,726,654 | 2026-10-05T03:04:21 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p263_v25.nc` | 25 | 985,727,738 | 2026-10-08T19:25:39 | repaired | 88 | 5 | 52 / 131 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p273_v23.nc` | 23 | 985,726,654 | 2026-10-02T05:27:48 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p273_v24.nc` | 24 | 985,726,654 | 2026-10-05T04:50:52 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_p273_v25.nc` | 25 | 985,727,738 | 2026-10-08T19:11:40 | repaired | 88 | 5 | 52 / 131 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_pold_v23.nc` | 23 | 985,726,654 | 2026-10-01T12:13:29 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_pold_v24.nc` | 24 | 985,726,654 | 2026-10-05T02:50:23 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_pold_v25.nc` | 25 | 985,727,738 | 2026-10-08T22:36:06 | repaired | 88 | 5 | 52 / 131 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_pstg_v23.nc` | 23 | 985,726,654 | 2026-10-02T01:04:43 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_pstg_v24.nc` | 24 | 985,726,654 | 2026-10-05T04:12:07 | repaired | 86 | 5 | 52 / 125 / 985,669,556 | yes |
| `terra_modis_m_liquid-water_a01_pstg_v25.nc` | 25 | 985,727,738 | 2026-10-08T22:39:32 | repaired | 88 | 5 | 52 / 131 / 985,669,556 | yes |
| `terra_modis_m_water-ice_a01_pagg_v23.nc` | 23 | 1,191,073,534 | 2026-10-02T06:22:26 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_pagg_v24.nc` | 24 | 1,191,073,534 | 2026-10-04T19:11:46 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_pagg_v25.nc` | 25 | 1,191,077,262 | 2026-10-08T00:40:51 | repaired | 92 | 5 | 52 / 133 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_pghm_v23.nc` | 23 | 1,191,073,534 | 2026-10-02T06:16:49 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_pghm_v24.nc` | 24 | 1,191,073,534 | 2026-10-05T00:30:04 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_pghm_v25.nc` | 25 | 1,191,077,262 | 2026-10-08T00:13:16 | repaired | 92 | 5 | 52 / 133 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_psph_v23.nc` | 23 | 1,191,073,534 | 2026-10-02T10:27:40 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_psph_v24.nc` | 24 | 1,191,073,534 | 2026-10-05T12:03:07 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_psph_v25.nc` | 25 | 1,191,077,262 | 2026-10-08T23:46:49 | repaired | 92 | 5 | 52 / 133 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_psrc_v23.nc` | 23 | 1,191,073,534 | 2026-10-02T06:41:12 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_psrc_v24.nc` | 24 | 1,191,073,534 | 2026-10-04T16:20:47 | repaired | 86 | 5 | 52 / 125 / 1,191,016,476 | yes |
| `terra_modis_m_water-ice_a01_psrc_v25.nc` | 25 | 1,191,077,262 | 2026-10-08T12:38:37 | repaired | 92 | 5 | 52 / 133 / 1,191,016,476 | yes |
