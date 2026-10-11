"""Minimal ctypes binding to the NetCDF C library for attribute-type inspection.

The netCDF4 Python package reads and writes attributes but does not report the
external NetCDF type (NC_CHAR versus NC_STRING) an attribute is stored with,
and that distinction is exactly what the ORAC retrieval depends on
(ORAC-CC/orac ``common/nc_get_string_att.c`` calls ``nc_get_att_string``, which
fails with NC_ECHAR, "Attempt to convert between text & numbers", on an
NC_CHAR attribute).  This module exposes the handful of C-API calls needed to
audit attribute types with the real library, to read a string attribute the
way ORAC does, and to write attributes of explicit external type in tests.

Only the inspection path is used by the production writer and the repair
utility; the write helpers exist so tests can construct files the netCDF4
package cannot (raw NC_CHAR bytes, user-defined attribute types).
"""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

# External NetCDF types (netcdf.h).
NC_NAT = 0
NC_BYTE = 1
NC_CHAR = 2
NC_SHORT = 3
NC_INT = 4
NC_FLOAT = 5
NC_DOUBLE = 6
NC_UBYTE = 7
NC_USHORT = 8
NC_UINT = 9
NC_INT64 = 10
NC_UINT64 = 11
NC_STRING = 12
NC_FIRSTUSERTYPEID = 32

NC_TYPE_NAMES = {
    NC_NAT: "NC_NAT", NC_BYTE: "NC_BYTE", NC_CHAR: "NC_CHAR", NC_SHORT: "NC_SHORT",
    NC_INT: "NC_INT", NC_FLOAT: "NC_FLOAT", NC_DOUBLE: "NC_DOUBLE", NC_UBYTE: "NC_UBYTE",
    NC_USHORT: "NC_USHORT", NC_UINT: "NC_UINT", NC_INT64: "NC_INT64",
    NC_UINT64: "NC_UINT64", NC_STRING: "NC_STRING",
}
NUMERIC_TYPES = frozenset({NC_BYTE, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE, NC_UBYTE,
                           NC_USHORT, NC_UINT, NC_INT64, NC_UINT64})

NC_GLOBAL = -1
NC_NOWRITE = 0
NC_WRITE = 1
NC_MAX_NAME = 256
NC_NOERR = 0
NC_ECHAR = -56          # "Attempt to convert between text & numbers"
NC_ENOTATT = -43        # "Attribute not found"

NC_FORMAT_NAMES = {1: "NC_FORMAT_CLASSIC", 2: "NC_FORMAT_64BIT_OFFSET", 3: "NC_FORMAT_NETCDF4",
                   4: "NC_FORMAT_NETCDF4_CLASSIC", 5: "NC_FORMAT_64BIT_DATA"}


class NetCDFError(OSError):
    """A NetCDF C-library call returned a non-zero status."""

    def __init__(self, status: int, message: str, context: str = ""):
        self.status = status
        self.message = message
        super().__init__(f"{context}: NetCDF: {message} (status {status})" if context
                         else f"NetCDF: {message} (status {status})")


def type_name(xtype: int) -> str:
    """Human-readable name of an external NetCDF type id."""

    return NC_TYPE_NAMES.get(xtype, f"user-defined type {xtype}")


def _candidate_libraries() -> list[str]:
    candidates = []
    override = os.environ.get("ORACLUT_LIBNETCDF")
    if override:
        candidates.append(override)
    prefix = Path(sys.prefix) / "lib"
    if prefix.is_dir():
        for pattern in ("libnetcdf.so", "libnetcdf.so.*", "libnetcdf*.dylib"):
            candidates.extend(str(path) for path in sorted(prefix.glob(pattern)))
    found = ctypes.util.find_library("netcdf")
    if found:
        candidates.append(found)
    return candidates


_LIBRARY: ctypes.CDLL | None = None


def library() -> ctypes.CDLL:
    """Load (once) and return the NetCDF C library with argument types declared."""

    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    errors = []
    for candidate in _candidate_libraries():
        try:
            lib = ctypes.CDLL(candidate)
        except OSError as exc:  # pragma: no cover - depends on the environment
            errors.append(f"{candidate}: {exc}")
            continue
        _declare(lib)
        _LIBRARY = lib
        return lib
    raise OSError("The NetCDF C library (libnetcdf) could not be loaded: " + "; ".join(errors))


def _declare(lib: ctypes.CDLL) -> None:
    c_int, c_size_t, c_char_p, c_void_p = ctypes.c_int, ctypes.c_size_t, ctypes.c_char_p, ctypes.c_void_p
    p_int, p_size_t = ctypes.POINTER(c_int), ctypes.POINTER(c_size_t)
    lib.nc_strerror.restype = c_char_p
    lib.nc_strerror.argtypes = [c_int]
    lib.nc_inq_libvers.restype = c_char_p
    lib.nc_inq_libvers.argtypes = []
    lib.nc_open.argtypes = [c_char_p, c_int, p_int]
    lib.nc_create.argtypes = [c_char_p, c_int, p_int]
    lib.nc_close.argtypes = [c_int]
    lib.nc_redef.argtypes = [c_int]
    lib.nc_enddef.argtypes = [c_int]
    lib.nc_inq_format.argtypes = [c_int, p_int]
    lib.nc_inq_nvars.argtypes = [c_int, p_int]
    lib.nc_inq_natts.argtypes = [c_int, p_int]
    lib.nc_inq_varname.argtypes = [c_int, c_int, c_char_p]
    lib.nc_inq_varid.argtypes = [c_int, c_char_p, p_int]
    lib.nc_inq_varnatts.argtypes = [c_int, c_int, p_int]
    lib.nc_inq_attname.argtypes = [c_int, c_int, c_int, c_char_p]
    lib.nc_inq_att.argtypes = [c_int, c_int, c_char_p, p_int, p_size_t]
    lib.nc_inq_attlen.argtypes = [c_int, c_int, c_char_p, p_size_t]
    lib.nc_inq_type.argtypes = [c_int, c_int, c_char_p, p_size_t]
    lib.nc_get_att.argtypes = [c_int, c_int, c_char_p, c_void_p]
    lib.nc_get_att_text.argtypes = [c_int, c_int, c_char_p, c_char_p]
    lib.nc_get_att_string.argtypes = [c_int, c_int, c_char_p, ctypes.POINTER(c_char_p)]
    lib.nc_free_string.argtypes = [c_size_t, ctypes.POINTER(c_char_p)]
    lib.nc_put_att.argtypes = [c_int, c_int, c_char_p, c_int, c_size_t, c_void_p]
    lib.nc_put_att_text.argtypes = [c_int, c_int, c_char_p, c_size_t, c_char_p]
    lib.nc_put_att_string.argtypes = [c_int, c_int, c_char_p, c_size_t, ctypes.POINTER(c_char_p)]
    lib.nc_def_vlen.argtypes = [c_int, c_char_p, c_int, p_int]
    lib.nc_inq_grps.argtypes = [c_int, p_int, p_int]


def strerror(status: int) -> str:
    """Return the library's message for a status code."""

    return library().nc_strerror(status).decode()


def check(status: int, context: str = "") -> None:
    """Raise NetCDFError for a non-zero status."""

    if status != NC_NOERR:
        raise NetCDFError(status, strerror(status), context)


def library_version() -> str:
    """Version string of the loaded NetCDF C library."""

    return library().nc_inq_libvers().decode()


@contextmanager
def open_file(path: str | os.PathLike, mode: int = NC_NOWRITE) -> Iterator[int]:
    """Open a NetCDF file with ``nc_open`` and close it on exit."""

    lib = library()
    ncid = ctypes.c_int()
    check(lib.nc_open(os.fsencode(path), mode, ctypes.byref(ncid)), f"nc_open({path})")
    try:
        yield ncid.value
    finally:
        lib.nc_close(ncid.value)


def inq_format(ncid: int) -> int:
    fmt = ctypes.c_int()
    check(library().nc_inq_format(ncid, ctypes.byref(fmt)), "nc_inq_format")
    return fmt.value


def inq_ngrps(ncid: int) -> int:
    count = ctypes.c_int()
    check(library().nc_inq_grps(ncid, ctypes.byref(count), None), "nc_inq_grps")
    return count.value


def inq_nvars(ncid: int) -> int:
    count = ctypes.c_int()
    check(library().nc_inq_nvars(ncid, ctypes.byref(count)), "nc_inq_nvars")
    return count.value


def inq_varname(ncid: int, varid: int) -> str:
    name = ctypes.create_string_buffer(NC_MAX_NAME + 1)
    check(library().nc_inq_varname(ncid, varid, name), "nc_inq_varname")
    return name.value.decode()


def inq_varid(ncid: int, name: str) -> int:
    varid = ctypes.c_int()
    check(library().nc_inq_varid(ncid, name.encode(), ctypes.byref(varid)), f"nc_inq_varid({name})")
    return varid.value


def inq_natts(ncid: int, varid: int) -> int:
    count = ctypes.c_int()
    if varid == NC_GLOBAL:
        check(library().nc_inq_natts(ncid, ctypes.byref(count)), "nc_inq_natts")
    else:
        check(library().nc_inq_varnatts(ncid, varid, ctypes.byref(count)), "nc_inq_varnatts")
    return count.value


def inq_attname(ncid: int, varid: int, index: int) -> str:
    name = ctypes.create_string_buffer(NC_MAX_NAME + 1)
    check(library().nc_inq_attname(ncid, varid, index, name), "nc_inq_attname")
    return name.value.decode()


def inq_att(ncid: int, varid: int, name: str) -> tuple[int, int]:
    """Return ``(xtype, length)`` of an attribute."""

    xtype = ctypes.c_int()
    length = ctypes.c_size_t()
    check(library().nc_inq_att(ncid, varid, name.encode(), ctypes.byref(xtype), ctypes.byref(length)),
          f"nc_inq_att({name})")
    return xtype.value, length.value


def inq_type_size(ncid: int, xtype: int) -> int:
    size = ctypes.c_size_t()
    check(library().nc_inq_type(ncid, xtype, None, ctypes.byref(size)), f"nc_inq_type({xtype})")
    return size.value


def get_att_text(ncid: int, varid: int, name: str, length: int) -> bytes:
    """Raw bytes of an NC_CHAR attribute (``nc_get_att_text``)."""

    buffer = ctypes.create_string_buffer(length + 1)
    check(library().nc_get_att_text(ncid, varid, name.encode(), buffer), f"nc_get_att_text({name})")
    return buffer.raw[:length]


def get_att_string(ncid: int, varid: int, name: str, length: int) -> list[str]:
    """Elements of an NC_STRING attribute (``nc_get_att_string``); fails with NC_ECHAR on NC_CHAR."""

    lib = library()
    strings = (ctypes.c_char_p * max(length, 1))()
    check(lib.nc_get_att_string(ncid, varid, name.encode(), strings), f"nc_get_att_string({name})")
    try:
        return [(strings[i] or b"").decode() for i in range(length)]
    finally:
        lib.nc_free_string(length, strings)


def get_att_raw(ncid: int, varid: int, name: str, xtype: int, length: int) -> bytes:
    """Raw in-file representation of a numeric attribute (``nc_get_att`` with no conversion)."""

    size = inq_type_size(ncid, xtype) * length
    buffer = ctypes.create_string_buffer(max(size, 1))
    check(library().nc_get_att(ncid, varid, name.encode(), buffer), f"nc_get_att({name})")
    return buffer.raw[:size]


def put_att_text(ncid: int, varid: int, name: str, value: bytes) -> None:
    """Write raw bytes as an NC_CHAR attribute (test helper)."""

    check(library().nc_put_att_text(ncid, varid, name.encode(), len(value), value), f"nc_put_att_text({name})")


def put_att_string(ncid: int, varid: int, name: str, values: list[str]) -> None:
    """Write an NC_STRING attribute (``nc_put_att_string``)."""

    array = (ctypes.c_char_p * len(values))(*[v.encode() for v in values])
    check(library().nc_put_att_string(ncid, varid, name.encode(), len(values), array), f"nc_put_att_string({name})")


def def_vlen(ncid: int, name: str, base_type: int) -> int:
    """Define a variable-length user type (test helper for unexpected attribute types)."""

    typeid = ctypes.c_int()
    check(library().nc_def_vlen(ncid, name.encode(), base_type, ctypes.byref(typeid)), "nc_def_vlen")
    return typeid.value


class _VlenT(ctypes.Structure):
    _fields_ = [("len", ctypes.c_size_t), ("p", ctypes.c_void_p)]


def put_att_vlen_int(ncid: int, varid: int, name: str, typeid: int, values: list[int]) -> None:
    """Write one vlen-of-int element as an attribute of user-defined type (test helper)."""

    data = (ctypes.c_int * len(values))(*values)
    element = _VlenT(len(values), ctypes.cast(data, ctypes.c_void_p))
    check(library().nc_put_att(ncid, varid, name.encode(), typeid, 1, ctypes.byref(element)), f"nc_put_att({name})")


def redef(ncid: int) -> None:
    check(library().nc_redef(ncid), "nc_redef")


def enddef(ncid: int) -> None:
    check(library().nc_enddef(ncid), "nc_enddef")


@dataclass(frozen=True)
class AttributeInfo:
    """One attribute as the C library reports it."""

    variable: str | None      # None for a global attribute
    name: str
    xtype: int
    length: int

    @property
    def type_name(self) -> str:
        return type_name(self.xtype)

    @property
    def location(self) -> str:
        return f"{self.variable or 'GLOBAL'}:{self.name}"


def list_attributes(path: str | os.PathLike) -> list[AttributeInfo]:
    """Every attribute in the root group, global first then per variable, in file order."""

    out: list[AttributeInfo] = []
    with open_file(path) as ncid:
        for index in range(inq_natts(ncid, NC_GLOBAL)):
            name = inq_attname(ncid, NC_GLOBAL, index)
            xtype, length = inq_att(ncid, NC_GLOBAL, name)
            out.append(AttributeInfo(None, name, xtype, length))
        for varid in range(inq_nvars(ncid)):
            variable = inq_varname(ncid, varid)
            for index in range(inq_natts(ncid, varid)):
                name = inq_attname(ncid, varid, index)
                xtype, length = inq_att(ncid, varid, name)
                out.append(AttributeInfo(variable, name, xtype, length))
    return out


def read_string_attribute_like_orac(path: str | os.PathLike, variable: str, name: str) -> str:
    """Read a variable attribute exactly as ORAC's ``nc_get_string_att`` does.

    Mirrors ORAC-CC/orac ``common/nc_get_string_att.c`` (nc_open, nc_inq_varid,
    nc_inq_attlen, nc_get_att_string, first element).  Raises NetCDFError with
    status NC_ECHAR for an NC_CHAR attribute, which is the production failure
    "ERROR: ncdf_get_string_att(): NetCDF: Attempt to convert between text & numbers".
    """

    with open_file(path) as ncid:
        varid = inq_varid(ncid, variable)
        length = ctypes.c_size_t()
        check(library().nc_inq_attlen(ncid, varid, name.encode(), ctypes.byref(length)), f"nc_inq_attlen({name})")
        values = get_att_string(ncid, varid, name, length.value)
    return values[0] if values else ""
