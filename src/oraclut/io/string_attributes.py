"""Audit and metadata-only repair of text attributes in ORAC LUT NetCDF files.

The NetCDF metadata contract of the ORAC retrieval
-------------------------------------------------
ORAC (ORAC-CC/orac, ``src/read_sad_lut.F90``, ``sad_dimension_read_nc``) reads
the ``spacing`` attribute of the five LUT axis variables ``optical_depth``,
``effective_radius``, ``satellite_zenith``, ``solar_zenith`` and
``relative_azimuth`` through ``common/nc_get_string_att.c``, which calls
``nc_get_att_string``.  That C call accepts only an NC_STRING attribute; on an
NC_CHAR attribute it returns NC_ECHAR and ORAC stops with
``ERROR: ncdf_get_string_att(): NetCDF: Attempt to convert between text & numbers``.
The IDL-written tables (V11, V19, V21) store every text attribute as
NC_STRING; ORAC reads nothing else as text (its ``units`` reads are
debug-only and status-checked).  The contract adopted here is therefore the
IDL representation: every text attribute, global or per variable, is
NC_STRING.  The only NC_CHAR attribute that must stay NC_CHAR is a
``_FillValue`` of a character variable, which is a fill character, not text.

Repair strategy
---------------
The repair never touches the numerical content and never modifies the
authoritative file in place.  For each file it

1. audits attribute types with the NetCDF C library (``oraclut.io.netcdf_c``);
2. copies the file byte-for-byte to ``<name>.repair-tmp`` in the same directory;
3. converts the NC_CHAR text attributes of the copy to NC_STRING with
   ``nc_del_att`` / ``nc_put_att_string`` (through the netCDF4 package),
   re-creating any later attribute of the same variable so attribute order is
   preserved;
4. verifies the copy against the original: dimensions, variables, storage
   layout, every attribute (type, length, raw value) and every data value,
   byte for byte, in bounded-memory slabs;
5. reads the axis ``spacing`` attributes of the copy exactly as ORAC does;
6. optionally hard-links the original into a retention directory, then
   atomically replaces the original with the verified copy.

A failure at any step leaves the authoritative file untouched and removes the
copy.  Running the utility again on a repaired file reports it as compliant
and changes nothing.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import stat
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Sequence

import numpy as np
from netCDF4 import Dataset

from . import netcdf_c as nc

ORAC_AXIS_VARIABLES = ("optical_depth", "effective_radius", "satellite_zenith",
                       "solar_zenith", "relative_azimuth")
ORAC_REQUIRED_STRING_ATTRIBUTE = "spacing"
# NC_CHAR attributes that are not text and must never be converted.
NON_TEXT_CHAR_ATTRIBUTES = frozenset({"_FillValue"})
TEMP_SUFFIX = ".repair-tmp"
DEFAULT_SLAB_BYTES = 64 * 1024 * 1024
DEFAULT_MIN_AGE_SECONDS = 600
_VERSION_PATTERN = re.compile(r"_v(\d+)\.nc$")

Logger = Callable[[str], None]


def _silent(_: str) -> None:
    pass


@dataclass(frozen=True)
class Conversion:
    """One NC_CHAR text attribute to be rewritten as NC_STRING with the same value."""

    variable: str | None
    name: str
    value: str
    char_length: int

    @property
    def location(self) -> str:
        return f"{self.variable or 'GLOBAL'}:{self.name}"

    @property
    def required_by_orac(self) -> bool:
        return self.variable in ORAC_AXIS_VARIABLES and self.name == ORAC_REQUIRED_STRING_ATTRIBUTE


@dataclass
class FileAudit:
    """Attribute inventory of one file and the conversions it needs."""

    path: Path
    file_format: str
    attributes: list[nc.AttributeInfo]
    conversions: list[Conversion]
    problems: list[str]

    @property
    def compliant(self) -> bool:
        return not self.conversions and not self.problems

    @property
    def required_conversions(self) -> list[Conversion]:
        return [c for c in self.conversions if c.required_by_orac]

    @property
    def text_attribute_types(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for attribute in self.attributes:
            if attribute.xtype in (nc.NC_CHAR, nc.NC_STRING) and attribute.name not in NON_TEXT_CHAR_ATTRIBUTES:
                counts[attribute.type_name] = counts.get(attribute.type_name, 0) + 1
        return counts


def audit_file(path: str | os.PathLike, attribute_names: Iterable[str] | None = None) -> FileAudit:
    """Inventory the attributes of ``path`` and list the NC_CHAR text attributes to convert.

    ``attribute_names`` restricts the conversions to the named attributes (for
    example ``{"spacing"}``); by default every NC_CHAR text attribute is
    selected, which is the representation of the IDL-written tables.  Anything
    the repair could not handle exactly is reported in ``problems`` and the
    file is then never modified.
    """

    path = Path(path)
    selected = None if attribute_names is None else set(attribute_names)
    problems: list[str] = []
    conversions: list[Conversion] = []
    with nc.open_file(path) as ncid:
        fmt = nc.inq_format(ncid)
        if nc.inq_ngrps(ncid) != 0:
            problems.append("file contains groups; only flat ORAC LUT files are supported")
    attributes = nc.list_attributes(path)
    if fmt != 3:
        problems.append(f"file format {nc.NC_FORMAT_NAMES.get(fmt, fmt)} cannot hold NC_STRING attributes "
                        "(only NC_FORMAT_NETCDF4 can)")
    with nc.open_file(path) as ncid:
        for attribute in attributes:
            if attribute.xtype == nc.NC_NAT or attribute.xtype >= nc.NC_FIRSTUSERTYPEID:
                problems.append(f"{attribute.location}: unexpected attribute type {attribute.type_name}")
                continue
            if attribute.xtype != nc.NC_CHAR or attribute.name in NON_TEXT_CHAR_ATTRIBUTES:
                continue
            if selected is not None and attribute.name not in selected:
                continue
            varid = nc.NC_GLOBAL if attribute.variable is None else nc.inq_varid(ncid, attribute.variable)
            raw = nc.get_att_text(ncid, varid, attribute.name, attribute.length)
            try:
                value = raw.decode("utf-8")
            except UnicodeDecodeError:
                problems.append(f"{attribute.location}: NC_CHAR value is not valid UTF-8 and cannot "
                                "be represented exactly as NC_STRING")
                continue
            if "\x00" in value:
                problems.append(f"{attribute.location}: NC_CHAR value contains an embedded NUL character")
                continue
            conversions.append(Conversion(attribute.variable, attribute.name, value, attribute.length))
    return FileAudit(path, nc.NC_FORMAT_NAMES.get(fmt, str(fmt)), attributes, conversions, problems)


def check_text_attributes_are_strings(path: str | os.PathLike) -> None:
    """Raise ValueError unless every text attribute of ``path`` is NC_STRING (writer self-check)."""

    audit = audit_file(path)
    if audit.conversions or audit.problems:
        details = [c.location + " is NC_CHAR" for c in audit.conversions] + audit.problems
        raise ValueError(f"{path}: text attributes violate the ORAC NC_STRING contract: " + "; ".join(details))


# ---------------------------------------------------------------------------
# Conversion of the working copy
# ---------------------------------------------------------------------------

def _container(dataset: Dataset, variable: str | None):
    return dataset if variable is None else dataset.variables[variable]


def convert_attributes(path: str | os.PathLike, audit: FileAudit) -> int:
    """Rewrite the audited NC_CHAR text attributes of ``path`` as NC_STRING, in place.

    For every variable (or the global scope) holding a conversion, the
    attributes from the first converted one onward are deleted and re-created
    in their original order, so the attribute order of the file is unchanged.
    Numeric attributes are re-created from their exact in-file values; an
    attribute that could not be re-created exactly raises before anything is
    deleted.  Returns the number of attributes converted.
    """

    if audit.problems:
        raise ValueError("refusing to convert a file with audit problems: " + "; ".join(audit.problems))
    wanted = {(c.variable, c.name): c for c in audit.conversions}
    plans: list[tuple[str | None, list[tuple[nc.AttributeInfo, str]]]] = []
    for scope in [None] + [a.variable for a in audit.attributes if a.variable is not None]:
        if any(p[0] == scope for p in plans):
            continue
        ordered = [a for a in audit.attributes if a.variable == scope]
        starts = [i for i, a in enumerate(ordered) if (scope, a.name) in wanted]
        if not starts:
            continue
        tail: list[tuple[nc.AttributeInfo, str]] = []
        for attribute in ordered[starts[0]:]:
            if (scope, attribute.name) in wanted:
                action = "string"
            elif attribute.xtype == nc.NC_STRING:
                action = "string-array"
            elif attribute.xtype in nc.NUMERIC_TYPES:
                action = "numeric"
            elif attribute.xtype == nc.NC_CHAR and attribute.name not in NON_TEXT_CHAR_ATTRIBUTES:
                action = "char"
            else:
                raise ValueError(f"{attribute.location}: attribute of type {attribute.type_name} follows a "
                                 "converted attribute and cannot be re-created exactly")
            tail.append((attribute, action))
        plans.append((scope, tail))

    converted = 0
    with Dataset(path, "r+") as dataset:
        for scope, tail in plans:
            target = _container(dataset, scope)
            # Capture every value before deleting anything.
            values = []
            for attribute, action in tail:
                if action == "string":
                    values.append(wanted[(scope, attribute.name)].value)
                else:
                    values.append(target.getncattr(attribute.name))
            for attribute, _ in reversed(tail):
                target.delncattr(attribute.name)
            for (attribute, action), value in zip(tail, values):
                if action == "string":
                    target.setncattr_string(attribute.name, value)
                    converted += 1
                elif action == "string-array":
                    target.setncattr_string(attribute.name, value if isinstance(value, list) else [value])
                else:   # numeric (exact numpy dtype) or NC_CHAR text left as NC_CHAR
                    target.setncattr(attribute.name, value)
    return converted


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

@dataclass
class VerificationReport:
    """Outcome of comparing a repaired copy with its original."""

    ok: bool
    differences: list[str] = field(default_factory=list)
    dimensions_compared: int = 0
    variables_compared: int = 0
    attributes_compared: int = 0
    attributes_converted: int = 0
    bytes_compared: int = 0
    slabs_compared: int = 0


def slab_ranges(shape: Sequence[int], itemsize: int, max_bytes: int) -> Iterator[tuple[int, int]]:
    """Yield ``(start, stop)`` ranges along axis 0 so that each slab is at most ``max_bytes``.

    A single row larger than ``max_bytes`` is still yielded as one slab (it is
    the smallest unit along the leading axis); otherwise every slab fits the
    bound and the slabs tile the leading axis exactly.
    """

    if len(shape) == 0:
        yield (0, 1)
        return
    row_bytes = int(np.prod(shape[1:], dtype=np.int64)) * itemsize if len(shape) > 1 else itemsize
    rows = max(1, max_bytes // max(row_bytes, 1))
    for start in range(0, shape[0], rows):
        yield (start, min(start + rows, shape[0]))
    if shape[0] == 0:
        yield (0, 0)


def _attribute_values(path: Path) -> list[tuple[str | None, str, int, int, object]]:
    """Ordered ``(variable, name, xtype, length, raw value)`` for every attribute via the C API."""

    rows = []
    with nc.open_file(path) as ncid:
        for attribute in nc.list_attributes(path):
            varid = nc.NC_GLOBAL if attribute.variable is None else nc.inq_varid(ncid, attribute.variable)
            if attribute.xtype == nc.NC_CHAR:
                value: object = nc.get_att_text(ncid, varid, attribute.name, attribute.length)
            elif attribute.xtype == nc.NC_STRING:
                value = nc.get_att_string(ncid, varid, attribute.name, attribute.length)
            elif attribute.xtype in nc.NUMERIC_TYPES:
                value = nc.get_att_raw(ncid, varid, attribute.name, attribute.xtype, attribute.length)
            else:
                value = ("unsupported", attribute.xtype)
            rows.append((attribute.variable, attribute.name, attribute.xtype, attribute.length, value))
    return rows


def _raw_slab(variable, start: int, stop: int) -> np.ndarray:
    if variable.ndim == 0:
        return np.asarray(variable[...])
    return np.asarray(variable[start:stop, ...])


def verify_preserved(original: str | os.PathLike, repaired: str | os.PathLike,
                     conversions: Sequence[Conversion], *, slab_bytes: int = DEFAULT_SLAB_BYTES,
                     log: Logger = _silent) -> VerificationReport:
    """Check that ``repaired`` equals ``original`` in everything but the listed conversions.

    Compared exactly: file format, dimensions (names, order, lengths, unlimited
    flag), variables (names, order, type, dimensions, shape, chunking, filters,
    byte order), every attribute (name, order, external type, length and raw
    in-file value) and every data value as raw bytes, which also preserves NaN
    payloads and signed zeros.  Data are read in slabs of at most
    ``slab_bytes`` along the leading dimension.
    """

    original, repaired = Path(original), Path(repaired)
    report = VerificationReport(ok=True)
    diffs = report.differences
    expected = {(c.variable, c.name): c for c in conversions}

    with Dataset(original, "r") as a, Dataset(repaired, "r") as b:
        for label, va, vb in (("data_model", a.data_model, b.data_model),
                              ("disk_format", a.disk_format, b.disk_format),
                              ("file_format", a.file_format, b.file_format)):
            if va != vb:
                diffs.append(f"{label}: {va!r} != {vb!r}")
        if list(a.groups) or list(b.groups):
            diffs.append("groups present; verification supports flat files only")
        if list(a.dimensions) != list(b.dimensions):
            diffs.append(f"dimension names/order differ: {list(a.dimensions)} != {list(b.dimensions)}")
        else:
            for name in a.dimensions:
                da, db = a.dimensions[name], b.dimensions[name]
                if len(da) != len(db) or da.isunlimited() != db.isunlimited():
                    diffs.append(f"dimension {name}: length/unlimited differ")
                report.dimensions_compared += 1
        if list(a.variables) != list(b.variables):
            diffs.append(f"variable names/order differ: {list(a.variables)} != {list(b.variables)}")
            report.ok = False
            return report
        for name in a.variables:
            va, vb = a.variables[name], b.variables[name]
            for label, pa, pb in (("dtype", va.dtype, vb.dtype), ("dimensions", va.dimensions, vb.dimensions),
                                  ("shape", va.shape, vb.shape), ("chunking", va.chunking(), vb.chunking()),
                                  ("filters", va.filters(), vb.filters()), ("endian", va.endian(), vb.endian()),
                                  ("datatype", str(va.datatype), str(vb.datatype))):
                if pa != pb:
                    diffs.append(f"variable {name}: {label} differs: {pa!r} != {pb!r}")
            report.variables_compared += 1

        # Attributes: C-API view, exact raw values.
        rows_a, rows_b = _attribute_values(original), _attribute_values(repaired)
        if [(r[0], r[1]) for r in rows_a] != [(r[0], r[1]) for r in rows_b]:
            diffs.append("attribute names or order differ: "
                         f"{[(r[0] or 'GLOBAL', r[1]) for r in rows_a]} != {[(r[0] or 'GLOBAL', r[1]) for r in rows_b]}")
        else:
            for ra, rb in zip(rows_a, rows_b):
                variable, name, xtype_a, len_a, value_a = ra
                _, _, xtype_b, len_b, value_b = rb
                location = f"{variable or 'GLOBAL'}:{name}"
                report.attributes_compared += 1
                conversion = expected.get((variable, name))
                if conversion is not None:
                    if not (xtype_a == nc.NC_CHAR and xtype_b == nc.NC_STRING and len_b == 1
                            and isinstance(value_a, bytes) and value_a.decode("utf-8") == value_b[0]
                            and value_b[0] == conversion.value):
                        diffs.append(f"{location}: conversion not as expected "
                                     f"({nc.type_name(xtype_a)}[{len_a}] {value_a!r} -> "
                                     f"{nc.type_name(xtype_b)}[{len_b}] {value_b!r})")
                    else:
                        report.attributes_converted += 1
                elif (xtype_a, len_a, value_a) != (xtype_b, len_b, value_b):
                    diffs.append(f"{location}: attribute differs ({nc.type_name(xtype_a)}[{len_a}] "
                                 f"{value_a!r} != {nc.type_name(xtype_b)}[{len_b}] {value_b!r})")
        missing = [c.location for c in conversions
                   if (c.variable, c.name) not in {(r[0], r[1]) for r in rows_b}]
        if missing:
            diffs.append("converted attributes missing from the repaired file: " + ", ".join(missing))

        # Data: raw bytes, bounded slabs.
        if not diffs:
            for name in a.variables:
                va, vb = a.variables[name], b.variables[name]
                for v in (va, vb):
                    v.set_auto_maskandscale(False)
                    if v.dtype.kind == "S":
                        v.set_auto_chartostring(False)
                itemsize = va.dtype.itemsize if va.dtype != object else 8
                for start, stop in slab_ranges(va.shape, itemsize, slab_bytes):
                    sa, sb = _raw_slab(va, start, stop), _raw_slab(vb, start, stop)
                    report.slabs_compared += 1
                    if sa.dtype != sb.dtype or sa.shape != sb.shape:
                        diffs.append(f"variable {name}[{start}:{stop}]: slab dtype/shape differ")
                        break
                    if sa.dtype == object:
                        same = bool(np.array_equal(sa, sb))
                    else:
                        same = sa.tobytes() == sb.tobytes()
                        report.bytes_compared += sa.nbytes
                    if not same:
                        diffs.append(f"variable {name}[{start}:{stop}]: data bytes differ")
                        break
                log(f"    verified {name} {va.shape} {va.dtype}")
    report.ok = not diffs
    return report


# ---------------------------------------------------------------------------
# Per-file repair
# ---------------------------------------------------------------------------

@dataclass
class FileResult:
    """Outcome of ``repair_file`` for one path."""

    path: str
    status: str                    # compliant | needs-repair | repaired | skipped | failed
    message: str = ""
    size_bytes: int = 0
    mtime_original: str = ""
    file_format: str = ""
    text_attribute_types: dict = field(default_factory=dict)
    conversions: list[dict] = field(default_factory=list)
    required_conversions: int = 0
    problems: list[str] = field(default_factory=list)
    verification: dict | None = None
    orac_spacing_read: dict = field(default_factory=dict)
    original_retained_as: str | None = None
    elapsed_seconds: float = 0.0

    @property
    def version(self) -> str:
        match = _VERSION_PATTERN.search(Path(self.path).name)
        return match.group(1) if match else ""


def _iso(timestamp: float) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(timestamp))


def _fsync_path(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _copy_file(source: Path, target: Path, log: Logger) -> None:
    started = time.time()
    with open(source, "rb") as src, open(target, "wb") as dst:
        shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
        dst.flush()
        os.fsync(dst.fileno())
    if target.stat().st_size != source.stat().st_size:
        raise OSError(f"copy size mismatch for {source}")
    log(f"    copied {source.stat().st_size} bytes in {time.time() - started:.1f} s")


def _orac_spacing_reads(path: Path) -> dict[str, str]:
    """Read every present axis ``spacing`` attribute exactly as ORAC does."""

    reads: dict[str, str] = {}
    with Dataset(path, "r") as dataset:
        present = [axis for axis in ORAC_AXIS_VARIABLES
                   if axis in dataset.variables and ORAC_REQUIRED_STRING_ATTRIBUTE in dataset.variables[axis].ncattrs()]
    for axis in present:
        reads[axis] = nc.read_string_attribute_like_orac(path, axis, ORAC_REQUIRED_STRING_ATTRIBUTE)
    return reads


def repair_file(path: str | os.PathLike, *, dry_run: bool = False,
                attribute_names: Iterable[str] | None = None, original_dir: str | os.PathLike | None = None,
                min_age_seconds: float = DEFAULT_MIN_AGE_SECONDS, slab_bytes: int = DEFAULT_SLAB_BYTES,
                log: Logger = _silent) -> FileResult:
    """Audit and, unless ``dry_run``, repair one LUT file as described in the module docstring."""

    path = Path(path)
    started = time.time()
    result = FileResult(path=str(path), status="failed")
    try:
        st = path.stat()
    except OSError as exc:
        result.message = f"cannot stat: {exc}"
        return result
    result.size_bytes, result.mtime_original = st.st_size, _iso(st.st_mtime)
    log(f"{path.name}: {st.st_size} bytes, modified {result.mtime_original}")

    try:
        audit = audit_file(path, attribute_names)
    except Exception as exc:   # unreadable or not NetCDF: report, never modify
        result.message = f"audit failed: {exc}"
        result.elapsed_seconds = time.time() - started
        log(f"  FAILED {result.message}")
        return result
    result.file_format = audit.file_format
    result.text_attribute_types = audit.text_attribute_types
    result.conversions = [asdict(c) for c in audit.conversions]
    result.required_conversions = len(audit.required_conversions)
    result.problems = list(audit.problems)
    log(f"  text attributes {audit.text_attribute_types}; {len(audit.conversions)} to convert "
        f"({len(audit.required_conversions)} required by ORAC)")

    if audit.problems:
        result.message = "; ".join(audit.problems)
        result.elapsed_seconds = time.time() - started
        log(f"  FAILED (left unchanged): {result.message}")
        return result
    if not audit.conversions:
        result.status, result.message = "compliant", "all text attributes already NC_STRING"
        try:
            result.orac_spacing_read = _orac_spacing_reads(path)
        except nc.NetCDFError as exc:
            result.status, result.message = "failed", f"ORAC-style spacing read failed: {exc}"
        result.elapsed_seconds = time.time() - started
        log(f"  {result.status}: {result.message}")
        return result
    if dry_run:
        result.status, result.message = "needs-repair", f"{len(audit.conversions)} NC_CHAR text attributes"
        result.elapsed_seconds = time.time() - started
        log(f"  needs repair (dry run, unchanged)")
        return result

    age = time.time() - st.st_mtime
    if age < min_age_seconds:
        result.status = "skipped"
        result.message = f"modified {age:.0f} s ago (< {min_age_seconds:.0f} s); may still be written"
        result.elapsed_seconds = time.time() - started
        log(f"  skipped: {result.message}")
        return result

    temp = path.with_name(path.name + TEMP_SUFFIX)
    retained: Path | None = None
    if original_dir is not None:
        original_dir = Path(original_dir)
        original_dir.mkdir(parents=True, exist_ok=True)
        if original_dir.stat().st_dev != st.st_dev:
            result.status = "skipped"
            result.message = f"retention directory {original_dir} is on another filesystem (hard link impossible)"
            result.elapsed_seconds = time.time() - started
            log(f"  skipped: {result.message}")
            return result
        retained = original_dir / path.name
        if retained.exists():
            result.status = "skipped"
            result.message = f"retention target {retained} already exists"
            result.elapsed_seconds = time.time() - started
            log(f"  skipped: {result.message}")
            return result
    if temp.exists():
        temp_age = time.time() - temp.stat().st_mtime
        if temp_age < min_age_seconds:
            result.status = "skipped"
            result.message = f"working copy {temp.name} is {temp_age:.0f} s old; another repair may be running"
            result.elapsed_seconds = time.time() - started
            log(f"  skipped: {result.message}")
            return result
        log(f"  removing stale working copy {temp.name} ({temp_age:.0f} s old)")
        temp.unlink()

    replaced = False
    try:
        _copy_file(path, temp, log)
        if path.stat().st_mtime != st.st_mtime or path.stat().st_size != st.st_size:
            raise RuntimeError("original changed while it was being copied")
        converted = convert_attributes(temp, audit)
        log(f"    converted {converted} attributes in the working copy")
        verification = verify_preserved(path, temp, audit.conversions, slab_bytes=slab_bytes, log=log)
        result.verification = asdict(verification)
        if not verification.ok:
            raise RuntimeError("verification failed: " + " | ".join(verification.differences[:10]))
        if verification.attributes_converted != len(audit.conversions):
            raise RuntimeError(f"verified {verification.attributes_converted} conversions, "
                               f"expected {len(audit.conversions)}")
        log(f"    verification ok: {verification.variables_compared} variables, "
            f"{verification.attributes_compared} attributes, {verification.bytes_compared} data bytes")
        reads = _orac_spacing_reads(temp)
        for axis, value in reads.items():
            expected_value = next((c.value for c in audit.conversions
                                   if c.variable == axis and c.name == ORAC_REQUIRED_STRING_ATTRIBUTE), None)
            if expected_value is not None and value != expected_value:
                raise RuntimeError(f"ORAC-style read of {axis}:spacing gave {value!r}, expected {expected_value!r}")
        result.orac_spacing_read = reads
        log(f"    ORAC-style nc_get_att_string reads: {reads}")
        shutil.copymode(path, temp)
        try:
            os.chown(temp, -1, st.st_gid)
        except OSError as exc:
            log(f"    warning: could not set group of the working copy: {exc}")
        _fsync_path(temp)
        if path.stat().st_mtime != st.st_mtime or path.stat().st_size != st.st_size:
            raise RuntimeError("original changed during the repair; not replacing it")
        if retained is not None:
            os.link(path, retained)
            result.original_retained_as = str(retained)
        os.replace(temp, path)
        replaced = True
        _fsync_path(path.parent)
        after = audit_file(path, attribute_names)
        if not after.compliant:
            raise RuntimeError("post-replacement audit still reports work: "
                               + "; ".join([c.location for c in after.conversions] + after.problems))
        result.orac_spacing_read = _orac_spacing_reads(path)
        result.status = "repaired"
        result.message = f"{converted} attributes converted to NC_STRING; data verified identical"
        log(f"  repaired: {result.message}")
    except Exception as exc:
        result.status = "failed"
        result.message = ("replaced but post-check failed: " if replaced else "original unchanged: ") + str(exc)
        log(f"  FAILED {result.message}")
    finally:
        if not replaced and temp.exists():
            try:
                temp.unlink()
                log(f"    removed working copy {temp.name}")
            except OSError as exc:   # pragma: no cover - filesystem dependent
                log(f"    warning: could not remove working copy {temp}: {exc}")
    result.elapsed_seconds = time.time() - started
    return result


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

INVENTORY_COLUMNS = ("path", "version", "size_bytes", "mtime_original", "status", "conversions",
                     "required_conversions", "text_attribute_types", "message")


def write_inventory(results: Sequence[FileResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        handle.write("\t".join(INVENTORY_COLUMNS) + "\n")
        for r in results:
            handle.write("\t".join(str(v) for v in (
                r.path, r.version, r.size_bytes, r.mtime_original, r.status, len(r.conversions),
                r.required_conversions, json.dumps(r.text_attribute_types, sort_keys=True),
                r.message.replace("\t", " "))) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m oraclut.repair_string_attributes",
        description="Audit or repair NC_CHAR text attributes of ORAC LUT NetCDF files (metadata only; "
                    "numerical content is never recomputed or rewritten).")
    parser.add_argument("files", nargs="+", type=Path, help="LUT files to audit or repair")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="report what would change; modify nothing")
    mode.add_argument("--check", action="store_true",
                      help="release check: modify nothing, exit 3 if any file is not compliant")
    parser.add_argument("--attributes", default=None,
                        help="comma-separated attribute names to convert (default: every NC_CHAR text attribute)")
    parser.add_argument("--original-dir", type=Path, default=None,
                        help="retain each original as a hard link in this directory (same filesystem)")
    parser.add_argument("--min-age-seconds", type=float, default=DEFAULT_MIN_AGE_SECONDS,
                        help="skip files modified more recently than this (default %(default)s)")
    parser.add_argument("--slab-mib", type=float, default=DEFAULT_SLAB_BYTES / 2**20,
                        help="maximum data slab read at once during verification (default %(default)s MiB)")
    parser.add_argument("--log-dir", type=Path, default=None, help="write a JSON record per file here")
    parser.add_argument("--inventory", type=Path, default=None, help="write a TSV inventory of all files here")
    parser.add_argument("--quiet", action="store_true", help="print one line per file only")
    args = parser.parse_args(argv)

    attribute_names = None if args.attributes is None else [a.strip() for a in args.attributes.split(",") if a.strip()]
    verbose = not args.quiet

    def log(message: str) -> None:
        if verbose or not message.startswith("    "):
            print(message, flush=True)

    print(f"NetCDF C library {nc.library_version()}", flush=True)
    results: list[FileResult] = []
    for file in args.files:
        result = repair_file(file, dry_run=args.dry_run or args.check, attribute_names=attribute_names,
                             original_dir=args.original_dir, min_age_seconds=args.min_age_seconds,
                             slab_bytes=int(args.slab_mib * 2**20), log=log)
        results.append(result)
        if args.log_dir is not None:
            args.log_dir.mkdir(parents=True, exist_ok=True)
            record = asdict(result)
            record["version"] = result.version
            record["recorded"] = _iso(time.time())
            with open(args.log_dir / (Path(file).name + ".repair.json"), "w") as handle:
                json.dump(record, handle, indent=1)
    if args.inventory is not None:
        write_inventory(results, args.inventory)

    counts: dict[str, int] = {}
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1
    print("summary: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), flush=True)
    if any(r.status == "failed" for r in results):
        return 1
    if args.check and any(r.status != "compliant" for r in results):
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
