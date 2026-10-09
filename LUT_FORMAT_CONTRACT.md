# ORAC LUT NetCDF format contract: text attributes

This document records the NetCDF metadata contract that the ORAC retrieval
imposes on the lookup-table (LUT) files written by this repository, the
regression that violated it (V22 to V25 products), the permanent writer
correction, the metadata-only repair utility, the regression tests, and the
release check.  The repair campaign itself is reported in
`validation/netcdf_string_attributes/REPORT_netcdf_string_attribute_repair.md`.

## 1. What ORAC requires

Source: ORAC-CC/orac, `master` at commit `eb9a1233c1136ef95fa0c37f650fe519e9134ed4`
(cloned read-only into a scratch directory on 2026-10-09; nothing from it is
copied into this repository).

- `src/read_sad_lut.F90`, `sad_dimension_read_nc` (lines 848-903), called from
  `Read_NCDF_SAD_LUT` (lines 979-983) for the five axis variables
  `optical_depth`, `effective_radius`, `satellite_zenith`, `solar_zenith` and
  `relative_azimuth`, reads each axis's `spacing` attribute with
  `ncdf_get_string_att` and interprets the values `logarithmic`,
  `uneven_logarithmic`, `linear`, `uneven_linear` and `unknown`.
- `common/orac_ncdf.F90`, `ncdf_get_string_att` (lines 257-295), is a Fortran
  interface to `common/nc_get_string_att.c`, which calls the NetCDF C function
  `nc_get_att_string` (line 42) and exits the program on any error (line 17-23):
  `ERROR: ncdf_get_string_att(): <nc_strerror>`.
- `nc_get_att_string` succeeds only on an attribute of external type
  **NC_STRING**.  On an **NC_CHAR** attribute it returns `NC_ECHAR`, whose
  message is `NetCDF: Attempt to convert between text & numbers`.  This is the
  failure Nathalie observed with the V25 tables.
- No other text attribute is read by ORAC.  The `units` reads in
  `common/ncdf_read_template.inc`, `ncdf_read_field.inc` and
  `ncdf_read_packed_field.inc` are inside `#ifdef DEBUG` and status-checked;
  `_FillValue`, `scale_factor`, `add_offset`, `valid_min` and `valid_max` are
  read numerically (`common/ncdf_open_field.inc`).

The LUTs that ORAC is known to read (IDL-written V11, V19 and V21) store
**every** text attribute, global and per variable, as NC_STRING of length 1
(verified with `nc_inq_att`).  The repository therefore adopts that
representation as the contract, not just the five `spacing` attributes:

> Every text attribute of an ORAC LUT product is NC_STRING.  The only NC_CHAR
> attribute permitted is a `_FillValue` of a character variable, which is a
> fill character rather than text.  Numeric attributes keep their numeric
> types.  The file format is NetCDF-4 (HDF5), the only format with NC_STRING.

The axis `spacing` attributes are the hard requirement; the remaining text
attributes are brought to the same type for consistency with the working
tables and so that any future ORAC string read behaves as it does on V11.

## 2. Root cause of the regression

`src/oraclut/io/v2.py`, `write_v2_lut`, wrote every attribute with the netCDF4
package's `setncattr` (formerly lines 89 and 91).  For a Python `str`,
`setncattr` writes **NC_CHAR** unless the dataset was created with
`set_ncstring_attrs(True)`; the IDL `NCDF_ATTPUT` used by the legacy writer
writes NC_STRING.  Both production entry points route through this primitive:
`src/oraclut/idl_mirror/write_v2_lut.py` (used by `create_orac_luts.py`, the
V23-V25 production path) and `src/oraclut/pipeline.py`.  Every text
attribute of every Python-written product (V22 EarthCARE MSI, V23, V24, V25)
was therefore NC_CHAR, and ORAC stopped at the first axis read.

## 3. Permanent writer correction

`src/oraclut/io/v2.py` now writes attributes through `put_attribute`, which
types text explicitly: `str` (and sequences or NumPy arrays of `str`) go
through `setncattr_string`, i.e. `nc_put_att_string`; other values go through
`setncattr` with their NumPy dtype; `bytes` is refused so that no text reaches
the file with an implicit type.  After the dataset is closed the writer
re-opens the file with the NetCDF C library
(`oraclut.io.string_attributes.check_text_attributes_are_strings`) and raises
if any text attribute is NC_CHAR.  A future library default can therefore not
reintroduce the problem silently: the write fails.

No scientific calculation, grid, variable, numerical value or version number
is affected; only the external type of the text attributes changes.

## 4. Repair utility for existing files

`src/oraclut/io/string_attributes.py` (command line:
`python -m oraclut.repair_string_attributes`) audits and repairs existing
products without recomputing or rewriting any numerical content:

```
PYTHONPATH=src python -m oraclut.repair_string_attributes --dry-run FILE...
PYTHONPATH=src python -m oraclut.repair_string_attributes --check FILE...      # release check, exit 3 if non-compliant
PYTHONPATH=src python -m oraclut.repair_string_attributes \
    [--original-dir DIR] [--attributes spacing,units] [--min-age-seconds N] \
    [--slab-mib N] [--log-dir DIR] [--inventory FILE.tsv] [--quiet] FILE...
```

Per file it: audits attribute types with the C library (`nc_inq_att`); copies
the file byte-for-byte to `<name>.repair-tmp` in the same directory; in the
copy deletes each NC_CHAR text attribute and re-creates it with
`nc_put_att_string`, re-creating any later attribute of the same variable so
the attribute order is unchanged; verifies the copy against the original
(format, dimensions, variables, chunking, filters, byte order, every attribute's
type, length and raw value, and every data value as raw bytes in slabs of at
most `--slab-mib`, so NaN payloads, signed zeros and fill values are compared
exactly); reads the axis `spacing` attributes of the copy exactly as ORAC does
(`nc_get_att_string`); copies mode and group; optionally hard-links the
original into `--original-dir` (same filesystem, zero copy); and atomically
replaces the original with the verified copy (`rename`).  Any failure leaves
the authoritative file untouched and removes the copy.  A file whose text
attributes are already NC_STRING is reported `compliant` and not modified, so
the utility is idempotent.  Files modified more recently than
`--min-age-seconds` (default 600) or with a fresh working copy are skipped, so
a file still being written by a job is never touched.

NetCDF-4/HDF5 notes: attribute deletion and creation are object-header
operations, so no data chunk moves; NC_STRING values live in the HDF5 global
heap, so a repaired file is a few kilobytes larger; the netCDF-4 format needs
no `nc_redef` for attribute changes (the netCDF4 package handles define mode
regardless).  Because an interrupted HDF5 metadata write could leave a file
unreadable, the utility never edits the authoritative file in place: it edits
a copy and installs it only after verification.

## 5. Regression tests

`tests/test_netcdf_string_attributes.py` checks with the NetCDF C library
(`src/oraclut/io/netcdf_c.py`, a ctypes binding to `libnetcdf`):

- the writer emits NC_STRING for every text attribute (axes, other variables,
  global) and `nc_get_att_string` reads them;
- a plain `setncattr` file reproduces the production failure (`NC_ECHAR`,
  "Attempt to convert between text & numbers");
- the writer's self-check fails the write if text attributes become NC_CHAR;
- the repair converts, preserves attribute order, storage layout, fill values,
  NaN payloads and signed zeros, is idempotent, leaves compliant files
  untouched, honours `--attributes`, handles many affected attributes across
  variables and the global scope, refuses unexpected attribute types (invalid
  UTF-8, user-defined types) without modifying the file, leaves the original
  untouched when conversion or verification fails, detects tampered data and
  metadata, skips recent files and fresh working copies, removes stale working
  copies, retains originals as hard links, reads large variables in bounded
  slabs, and exposes dry-run, check and inventory modes;
- the IDL-written V11 reference is compliant (skipped if the archive is absent).

Run with `python -m pytest tests/test_netcdf_string_attributes.py`.

## 6. Release check for future products

Before a product is released (or copied into an ORAC `sad_dir`), run

```
PYTHONPATH=src python -m oraclut.repair_string_attributes --check <products>
```

which exits 0 only if every text attribute is NC_STRING and the axis
`spacing` attributes read through `nc_get_att_string`.  The writer's self-check
makes a non-compliant product from this repository impossible to write, but
the check also covers products from other writers or older revisions.  The
`validation/netcdf_string_attributes/orac_c_reader/` harness compiles ORAC's
own unmodified `common/nc_get_string_att.c` from a read-only ORAC checkout and
reads the five axis attributes with it, which is the closest available test to
the retrieval itself when no ORAC executable is accessible.
