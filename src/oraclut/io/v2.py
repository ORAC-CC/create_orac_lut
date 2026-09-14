"""Writer for the ORAC LUT level-2 NetCDF representation."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np
from netCDF4 import Dataset


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
        for name, size in dimensions.items():
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
            for attribute, value in attrs.get(name, {}).items():
                if attribute != "_FillValue":
                    variable.setncattr(attribute, value)
        for attribute, value in (global_attributes or {}).items():
            dataset.setncattr(attribute, value)
