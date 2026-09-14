from __future__ import annotations

import numpy as np
import pytest
from netCDF4 import Dataset

from oraclut.io.lut import LutFormatError, read_lut
from oraclut.io.v2 import write_v2_lut


REFERENCE = "/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"


def test_reader_preserves_reference_structure():
    try:
        lut = read_lut(REFERENCE)
    except FileNotFoundError:
        pytest.skip("authorised ORAC reference archive is unavailable")

    assert lut.file_format == "NETCDF4"
    assert tuple(lut.dimensions)[:6] == (
        "optical_depth", "effective_radius", "satellite_zenith",
        "solar_zenith", "relative_azimuth", "channels",
    )
    assert lut.dimensions["channels"] == 11
    np.testing.assert_array_equal(lut.coordinates["effective_radius"], np.arange(1, 40, 2, dtype=np.float32))
    assert lut.variable_dimensions["T_dv"] == (
        "satellite_zenith", "optical_depth", "effective_radius", "channels"
    )
    assert lut.variables["central_wavelength"].dtype == np.dtype("float32")
    assert lut.variable_attributes["central_wavelength"]["units"] == "microns"


def test_reader_reads_synthetic_v21_without_transpose(tmp_path):
    path = tmp_path / "synthetic.nc"
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("optical_depth", 1)
        ds.createDimension("effective_radius", 1)
        ds.createDimension("satellite_zenith", 1)
        ds.createDimension("solar_zenith", 1)
        ds.createDimension("relative_azimuth", 1)
        ds.createDimension("channels", 1)
        ds.createDimension("solar_channels", 1)
        ds.createDimension("thermal_channels", 1)
        ds.createDimension("mixed_channels", 1)
        ds.createDimension("length", 1)
        for name in ("optical_depth", "effective_radius", "satellite_zenith", "solar_zenith", "relative_azimuth"):
            ds.createVariable(name, "f4", (name,))[:] = 0
        ds.createVariable("channel_id", "i2", ("channels",))[:] = 1
        for name in ("central_wavelength", "central_wavenumber", "solar_channel_flag", "mixed_channel_flag", "thermal_channel_flag"):
            ds.createVariable(name, "f4", ("channels",))[:] = 0
        for name in ("average_volume_per_particle",):
            ds.createVariable(name, "f4", ("effective_radius",))[:] = 0
        for name in ("extinction_coefficient", "extinction_coefficient_ratio", "single_scatter_albedo", "asymmetry_parameter"):
            ds.createVariable(name, "f4", ("effective_radius", "channels"))[:] = 0
        for name, dims in {
            "T_dv": ("satellite_zenith", "optical_depth", "effective_radius", "channels"),
            "T_dd": ("optical_depth", "effective_radius", "channels"),
            "R_dv": ("satellite_zenith", "optical_depth", "effective_radius", "channels"),
            "R_dd": ("optical_depth", "effective_radius", "channels"),
            "R_0v": ("relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth", "effective_radius", "solar_channels"),
            "R_0d": ("solar_zenith", "optical_depth", "effective_radius", "solar_channels"),
            "T_0d": ("solar_zenith", "optical_depth", "effective_radius", "solar_channels"),
            "T_00": ("solar_zenith", "optical_depth", "effective_radius", "solar_channels"),
            "E_md": ("satellite_zenith", "optical_depth", "effective_radius", "thermal_channels"),
        }.items():
            ds.createVariable(name, "f4", dims)[:] = np.arange(np.prod([len(ds.dimensions[d]) for d in dims]), dtype="f4").reshape(tuple(len(ds.dimensions[d]) for d in dims))
    lut = read_lut(path)
    assert lut.variable_dimensions["R_0v"] == (
        "relative_azimuth", "satellite_zenith", "solar_zenith", "optical_depth", "effective_radius", "solar_channels"
    )
    np.testing.assert_array_equal(lut.variables["T_dv"], np.zeros((1, 1, 1, 1), dtype="f4"))


def test_reader_rejects_missing_required_variable(tmp_path):
    path = tmp_path / "invalid.nc"
    with Dataset(path, "w") as ds:
        ds.createDimension("optical_depth", 1)
    with pytest.raises(LutFormatError, match="Missing required V2 dimensions"):
        read_lut(path)


def test_v2_writer_preserves_declared_order_and_separates_revision(tmp_path):
    path = tmp_path / "small_v2.nc"
    dimensions = {"optical_depth": 2, "channels": 1}
    variables = {
        "optical_depth": np.array([1.0, 2.0], dtype="f4"),
        "channel_id": np.array([1], dtype="i2"),
        "T_dd": np.array([[3.0], [4.0]], dtype="f4"),
    }
    write_v2_lut(
        path, lut_level=2, revision=21, dimensions=dimensions,
        variables=variables,
        variable_dimensions={
            "optical_depth": ("optical_depth",),
            "channel_id": ("channels",),
            "T_dd": ("optical_depth", "channels"),
        },
    )
    with Dataset(path) as ds:
        assert list(ds.dimensions) == ["optical_depth", "channels"]
        assert ds.variables["T_dd"].dimensions == ("optical_depth", "channels")
        np.testing.assert_array_equal(ds.variables["T_dd"][:], variables["T_dd"])
