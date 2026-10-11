"""Writer for the ORAC LUT level-2 NetCDF representation."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np
from netCDF4 import Dataset

from .string_attributes import check_text_attributes_are_strings


# Global dimension declaration order of the legacy V2 writer (write_v2_lut.pro,
# ncdf_dimdef sequence). The pressure dimension, when present, is declared
# between relative_azimuth and channels; the channel-class dimensions follow the
# string-length dimensions. Dimensions not named here keep their given order
# after these.
V2_DIMENSION_ORDER = (
    "optical_depth", "effective_radius", "satellite_zenith", "solar_zenith",
    "relative_azimuth", "surface_pressure", "channels", "length",
    "st1", "st2", "st3", "st4", "solar_channels", "thermal_channels", "mixed_channels",
)

# Coordinate attributes the legacy writer sets as literal constants. They are
# applied only where the caller supplies no value for that attribute; anything
# the caller does supply wins. ``spacing`` is deliberately absent: it is the
# LUT-definition header word and must come from the grid, never be inferred.
V2_COORDINATE_DEFAULTS: dict[str, dict[str, Any]] = {
    "surface_pressure": {
        "long_name": "surface pressure",
        "units": "hPa",
        "valid_range": np.asarray([900.0, 1100.0], dtype=np.float32),
    },
}


def v2_dimension_order(dimensions: Mapping[str, int]) -> list[str]:
    """Return dimension names in the legacy V2 declaration order."""

    known = [name for name in V2_DIMENSION_ORDER if name in dimensions]
    return known + [name for name in dimensions if name not in V2_DIMENSION_ORDER]


def put_attribute(target: Any, name: str, value: Any) -> None:
    """Write one attribute with an explicit NetCDF external type.

    Text (``str``, or a sequence of ``str``) is written as NC_STRING through
    ``setncattr_string`` (``nc_put_att_string``), which is the type the ORAC
    reader requires for the axis ``spacing`` attributes and the type the
    IDL-written tables use for every text attribute.  The netCDF4 package's
    plain ``setncattr`` writes ``str`` as NC_CHAR, which ORAC's
    ``nc_get_att_string`` rejects with NC_ECHAR; that default produced the
    V22-V25 regression and is deliberately not used for text here.  Non-text
    values are written through ``setncattr`` with their NumPy dtype.  ``bytes``
    is refused so that no text reaches the file with an implicit type.
    """

    if isinstance(value, (bytes, bytearray)):
        raise TypeError(f"attribute {name!r}: pass text as str, not bytes, so its NetCDF type is explicit")
    if isinstance(value, str):
        target.setncattr_string(name, value)
    elif isinstance(value, np.ndarray) and value.dtype.kind in "US":
        target.setncattr_string(name, [str(item) for item in value.ravel().tolist()])
    elif isinstance(value, (list, tuple)) and value and all(isinstance(item, str) for item in value):
        target.setncattr_string(name, list(value))
    else:
        target.setncattr(name, value)


def write_v2_lut(
    path: str | Path,
    *,
    lut_level: int,
    revision: int,
    dimensions: Mapping[str, int],
    variables: Mapping[str, np.ndarray],
    variable_dimensions: Mapping[str, tuple[str, ...]],
    variable_attributes: Mapping[str, Mapping[str, Any]] | None = None,
    global_attributes: Mapping[str, Any] | None = None,
) -> None:
    """Write supplied arrays using the observed ORAC LUT level-2 layout.

    Scientific calculations and array construction remain outside this writer.
    ``lut_level`` and ``revision`` are separate arguments so revision ``21`` is
    not mistaken for the format/schema level.
    """

    if lut_level != 2:
        raise ValueError(f"Only ORAC LUT level 2 is supported, got {lut_level!r}")
    if not isinstance(revision, int) or revision < 0:
        raise ValueError(f"revision must be a non-negative integer, got {revision!r}")
    output = Path(path)
    attrs = variable_attributes or {}
    with Dataset(output, "w", format="NETCDF4") as dataset:
        for name in v2_dimension_order(dimensions):
            size = dimensions[name]
            if not isinstance(size, int) or size < 0:
                raise ValueError(f"Invalid dimension {name!r}: {size!r}")
            dataset.createDimension(name, size)
        for name, values in variables.items():
            if name not in variable_dimensions:
                raise ValueError(f"Missing dimensions for variable {name!r}")
            dims = tuple(variable_dimensions[name])
            if any(dim not in dimensions for dim in dims):
                raise ValueError(f"Variable {name!r} references an unknown dimension: {dims!r}")
            array = np.asarray(values)
            expected = tuple(dimensions[dim] for dim in dims)
            if array.shape != expected:
                raise ValueError(f"Variable {name!r} has shape {array.shape}, expected {expected}")
            fill = attrs.get(name, {}).get("_FillValue", False)
            variable = dataset.createVariable(name, array.dtype, dims, fill_value=fill)
            variable[:] = array
            variable_attrs = {**V2_COORDINATE_DEFAULTS.get(name, {}), **attrs.get(name, {})}
            for attribute, value in variable_attrs.items():
                if attribute != "_FillValue":
                    put_attribute(variable, attribute, value)
        for attribute, value in (global_attributes or {}).items():
            put_attribute(dataset, attribute, value)
    # Self-check with the NetCDF C library: every text attribute must be
    # NC_STRING (the ORAC reader contract), independent of library defaults.
    check_text_attributes_are_strings(output)
