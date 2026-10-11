"""The Python V2 writer must reproduce the legacy surface_pressure representation.

Legacy facts are read from the immutable captured aerosol references rather
than retyped, so the tests follow the legacy product if it is ever recaptured.
The dimension-order and default-attribute tests fail against the pre-fix writer,
which declared dimensions in caller order and wrote no surface_pressure metadata;
the end-to-end product test additionally requires the LUT-definition spacing to
be forwarded through the pipeline metadata table.
"""

from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from oraclut.config import read_lut_grid
from oraclut.io.lut import read_lut
from oraclut.io.v2 import V2_DIMENSION_ORDER, v2_dimension_order, write_v2_lut


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
LEGACY_AEROSOL = (
    ROOT / "validation" / "generated" / "aerosol_legacy_reference"
    / "meteosat10_seviri_m_aerosol_a12_pa79_v21_legacy.nc"
)
LEGACY_AEROSOL_IR2 = (
    ROOT / "validation" / "generated" / "aerosol_visible_ir_comparison" / "legacy"
    / "meteosat-10_seviri_m_aerosol_a12_pa79_v21_legacy.nc"
)
LEGACY_CLOUD = (
    ROOT / "validation" / "generated" / "meteosat10_seviri_liquid_water_stg_cloud_test"
    / "meteosat10_seviri_liquid_water_stg_cloud_test_legacy_reference_v21.nc"
)
PYTHON_AEROSOL = (
    ROOT / "validation" / "generated" / "aerosol_comparison" / "python"
    / "meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc"
)
PYTHON_CLOUD = (
    ROOT / "create_orac_lut" / "luts" / "meteosat-10_seviri_cloud"
    / "meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"
)

needs_legacy = pytest.mark.skipif(not LEGACY_AEROSOL.is_file(), reason="captured legacy aerosol reference missing")


def _write_minimal(path, dimensions, extra_attrs=None):
    variables = {"surface_pressure": np.asarray([950.0, 1013.0, 1050.0], dtype=np.float32)}
    variable_dimensions = {"surface_pressure": ("surface_pressure",)}
    write_v2_lut(
        path, lut_level=2, revision=21, dimensions=dimensions, variables=variables,
        variable_dimensions=variable_dimensions, variable_attributes=extra_attrs,
    )


def test_dimension_order_follows_the_legacy_writer_regardless_of_caller_order():
    caller_order = {
        "optical_depth": 2, "effective_radius": 2, "satellite_zenith": 2, "solar_zenith": 2,
        "relative_azimuth": 2, "channels": 1, "length": 32, "st1": 26, "st2": 11, "st3": 6,
        "st4": 7, "surface_pressure": 3, "solar_channels": 1, "thermal_channels": 1,
    }
    assert v2_dimension_order(caller_order) == [
        "optical_depth", "effective_radius", "satellite_zenith", "solar_zenith",
        "relative_azimuth", "surface_pressure", "channels", "length", "st1", "st2", "st3", "st4",
        "solar_channels", "thermal_channels",
    ]
    # Unknown dimensions are kept, after the known ones, in the order given.
    assert v2_dimension_order({"zzz": 1, "channels": 1, "aaa": 2})[-2:] == ["zzz", "aaa"]


@needs_legacy
def test_written_dimension_order_matches_the_legacy_aerosol_reference(tmp_path):
    with Dataset(LEGACY_AEROSOL) as legacy:
        legacy_order = list(legacy.dimensions)
        legacy_sizes = {name: len(dim) for name, dim in legacy.dimensions.items()}
    # Feed the writer the legacy dimensions deliberately scrambled.
    scrambled = {name: legacy_sizes[name] for name in reversed(legacy_order)}
    path = tmp_path / "order.nc"
    _write_minimal(path, scrambled)
    with Dataset(path) as written:
        assert list(written.dimensions) == legacy_order


@needs_legacy
def test_surface_pressure_defaults_match_the_legacy_constants(tmp_path):
    with Dataset(LEGACY_AEROSOL) as legacy:
        reference = legacy.variables["surface_pressure"]
        expected = {name: reference.getncattr(name) for name in reference.ncattrs()}
        expected_dtype = reference.dtype
        expected_values = reference[:].tolist()
    path = tmp_path / "defaults.nc"
    # The caller supplies only the LUT-definition spacing; everything else is a writer default.
    _write_minimal(path, {"surface_pressure": 3}, {"surface_pressure": {"spacing": expected["spacing"]}})
    with Dataset(path) as written:
        variable = written.variables["surface_pressure"]
        assert variable.dtype == expected_dtype
        assert variable.dimensions == ("surface_pressure",)
        assert variable[:].tolist() == expected_values == [950.0, 1013.0, 1050.0]
        assert set(variable.ncattrs()) == set(expected) == {"long_name", "spacing", "units", "valid_range"}
        assert variable.getncattr("long_name") == expected["long_name"]
        assert variable.getncattr("units") == expected["units"]
        assert variable.getncattr("spacing") == expected["spacing"]
        np.testing.assert_array_equal(variable.getncattr("valid_range"), expected["valid_range"])
        assert variable.getncattr("valid_range").dtype == expected["valid_range"].dtype


def test_caller_attributes_take_precedence_over_defaults(tmp_path):
    path = tmp_path / "override.nc"
    _write_minimal(path, {"surface_pressure": 3}, {"surface_pressure": {"units": "mb", "spacing": "linear"}})
    with Dataset(path) as written:
        variable = written.variables["surface_pressure"]
        assert variable.getncattr("units") == "mb"
        assert variable.getncattr("spacing") == "linear"
        assert variable.getncattr("long_name") == "surface pressure"


def test_reader_keeps_the_pressure_spacing_keyword():
    grid = read_lut_grid(INPUTS / "lut" / "aerosol_test.lut")
    assert grid.surface_pressure_spacing == "uneven_linear"
    assert read_lut_grid(INPUTS / "lut" / "liquid-water-cloud_test.lut").surface_pressure_spacing is None


@needs_legacy
def test_python_aerosol_product_declares_dimensions_like_legacy():
    if not PYTHON_AEROSOL.is_file():
        pytest.skip("regenerated Python aerosol validation product missing")
    with Dataset(LEGACY_AEROSOL) as legacy, Dataset(PYTHON_AEROSOL) as python:
        assert list(python.dimensions) == list(legacy.dimensions)
        legacy_var = legacy.variables["surface_pressure"]
        python_var = python.variables["surface_pressure"]
        assert python_var.dtype == legacy_var.dtype
        assert python_var[:].tolist() == legacy_var[:].tolist()
        for name in ("long_name", "units"):
            assert python_var.getncattr(name) == legacy_var.getncattr(name)
        np.testing.assert_array_equal(python_var.getncattr("valid_range"), legacy_var.getncattr("valid_range"))
        assert "spacing" in python_var.ncattrs(), "spacing must be forwarded from the LUT definition"
        assert python_var.getncattr("spacing") == legacy_var.getncattr("spacing")


@pytest.mark.skipif(not (LEGACY_CLOUD.is_file() and PYTHON_CLOUD.is_file()), reason="cloud references missing")
def test_cloud_products_have_no_surface_pressure_and_keep_legacy_order():
    with Dataset(LEGACY_CLOUD) as legacy, Dataset(PYTHON_CLOUD) as python:
        assert "surface_pressure" not in python.dimensions
        assert "surface_pressure" not in python.variables
        assert list(python.dimensions) == list(legacy.dimensions)


def test_v2_dimension_order_lists_the_pressure_dimension_before_channels():
    order = list(V2_DIMENSION_ORDER)
    assert order.index("surface_pressure") == order.index("relative_azimuth") + 1
    assert order.index("surface_pressure") == order.index("channels") - 1
