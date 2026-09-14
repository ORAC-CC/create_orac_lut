"""Reader for the current ORAC V2 NetCDF LUT products.

The reader deliberately keeps NetCDF dimension and variable ordering intact.
It is an inspection layer, not a conversion to a new scientific schema.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError as exc:  # pragma: no cover - depends on the environment
    raise ImportError(
        "Reading ORAC LUTs requires the netCDF4 Python package; "
        "install it through the project Conda environment."
    ) from exc


class LutFormatError(ValueError):
    """Raised when a file is not structurally compatible with V2."""


@dataclass(frozen=True)
class LutFile:
    """An in-memory, ordering-preserving V2 LUT inspection result."""

    path: Path
    file_format: str
    dimensions: Mapping[str, int]
    coordinates: Mapping[str, np.ndarray]
    variables: Mapping[str, np.ndarray]
    variable_dimensions: Mapping[str, tuple[str, ...]]
    variable_dtypes: Mapping[str, str]
    variable_attributes: Mapping[str, Mapping[str, Any]]
    global_attributes: Mapping[str, Any]

    @property
    def variable_names(self) -> tuple[str, ...]:
        """Return variable names in the NetCDF declaration order."""

        return tuple(self.variables)


V2_DIMENSIONS = (
    "optical_depth",
    "effective_radius",
    "satellite_zenith",
    "solar_zenith",
    "relative_azimuth",
    "channels",
)

V2_COORDINATE_VARIABLES = V2_DIMENSIONS[:5] + ("channel_id",)

V2_REQUIRED_VARIABLES = (
    "optical_depth",
    "effective_radius",
    "satellite_zenith",
    "solar_zenith",
    "relative_azimuth",
    "channel_id",
    "central_wavelength",
    "central_wavenumber",
    "solar_channel_flag",
    "mixed_channel_flag",
    "thermal_channel_flag",
    "average_volume_per_particle",
    "extinction_coefficient",
    "extinction_coefficient_ratio",
    "single_scatter_albedo",
    "asymmetry_parameter",
    "T_dv",
    "T_dd",
    "R_dv",
    "R_dd",
    "R_0v",
    "R_0d",
    "T_0d",
    "T_00",
)


def _attributes(obj: Any) -> dict[str, Any]:
    return {name: getattr(obj, name) for name in obj.ncattrs()}


def _validate_structure(
    dimensions: Mapping[str, int],
    variables: Mapping[str, np.ndarray],
    variable_dimensions: Mapping[str, tuple[str, ...]],
) -> None:
    missing_dimensions = [name for name in V2_DIMENSIONS if name not in dimensions]
    if missing_dimensions:
        raise LutFormatError(
            "Missing required V2 dimensions: " + ", ".join(missing_dimensions)
        )
    missing_variables = [name for name in V2_REQUIRED_VARIABLES if name not in variables]
    if missing_variables:
        raise LutFormatError(
            "Missing required V2 variables: " + ", ".join(missing_variables)
        )
    if dimensions.get("thermal_channels", 0) and "E_md" not in variables:
        raise LutFormatError("Thermal-channel LUT is missing required variable 'E_md'")
    for name in V2_DIMENSIONS[:5]:
        if variable_dimensions.get(name) != (name,):
            raise LutFormatError(
                f"Coordinate variable {name!r} must have dimensions ({name!r},), "
                f"got {variable_dimensions.get(name)!r}"
            )
    if variable_dimensions.get("channel_id") != ("channels",):
        raise LutFormatError(
            "Coordinate variable 'channel_id' must have dimensions ('channels',), "
            f"got {variable_dimensions.get('channel_id')!r}"
        )


def read_lut(path: str | Path, *, validate_v21: bool = True) -> LutFile:
    """Read a V2 LUT without transposing or otherwise transforming arrays.

    Parameters
    ----------
    path:
        Existing NetCDF LUT path.
    validate_v21:
        When true (the default), require the dimensions and variables observed
        in current V2 products. Set false only for deliberate exploratory
        inspection of a partial file.
    """

    lut_path = Path(path)
    if not lut_path.is_file():
        raise FileNotFoundError(f"ORAC LUT does not exist: {lut_path}")

    with Dataset(lut_path, mode="r") as dataset:
        dimensions = {name: len(dim) for name, dim in dataset.dimensions.items()}
        variables = {name: np.array(variable[:], copy=True) for name, variable in dataset.variables.items()}
        variable_dimensions = {
            name: tuple(variable.dimensions) for name, variable in dataset.variables.items()
        }
        variable_dtypes = {
            name: str(variable.dtype) for name, variable in dataset.variables.items()
        }
        variable_attributes = {
            name: _attributes(variable) for name, variable in dataset.variables.items()
        }
        coordinates = {name: variables[name] for name in V2_COORDINATE_VARIABLES if name in variables}
        result = LutFile(
            path=lut_path,
            file_format=dataset.file_format,
            dimensions=dimensions,
            coordinates=coordinates,
            variables=variables,
            variable_dimensions=variable_dimensions,
            variable_dtypes=variable_dtypes,
            variable_attributes=variable_attributes,
            global_attributes=_attributes(dataset),
        )

    if validate_v21:
        _validate_structure(result.dimensions, result.variables, result.variable_dimensions)
    return result
