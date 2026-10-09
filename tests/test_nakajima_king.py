from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pytest

from validation.nakajima_king.compare_nakajima_king import (
    ComparisonError,
    calculate_displacement,
    check_compatible,
    compare_luts,
    create_figure,
    extract_nakajima_king,
    read_lut_metadata,
    resolve_selection,
)


def _write_lut(
    path: Path,
    *,
    offset: float = 0.0,
    optical_depth=(1.0, 5.0),
    wavelengths=(0.65, 1.6, 2.2),
) -> None:
    with Dataset(path, "w", format="NETCDF4") as dataset:
        dimensions = {
            "optical_depth": len(optical_depth),
            "effective_radius": 3,
            "satellite_zenith": 1,
            "solar_zenith": 1,
            "relative_azimuth": 1,
            "channels": 3,
            "solar_channels": 2,
            "st1": 4,
            "st2": 3,
        }
        for name, size in dimensions.items():
            dataset.createDimension(name, size)
        for name, values in {
            "optical_depth": optical_depth,
            "effective_radius": (5.0, 10.0, 15.0),
            "satellite_zenith": (0.0,),
            "solar_zenith": (30.0,),
            "relative_azimuth": (0.0,),
        }.items():
            dataset.createVariable(name, "f4", (name,))[:] = values
        dataset.createVariable("channel_id", "i2", ("channels",))[:] = (1, 2, 3)
        dataset.createVariable("solar_channel_id", "i2", ("solar_channels",))[:] = (1, 3)
        dataset.createVariable("central_wavelength", "f4", ("channels",))[:] = wavelengths
        dataset.createVariable("platform", "S1", ("st1",))[:] = np.asarray(list("test"), dtype="S1")
        dataset.createVariable("instrument", "S1", ("st2",))[:] = np.asarray(list("abc"), dtype="S1")
        dims = (
            "relative_azimuth", "satellite_zenith", "solar_zenith",
            "optical_depth", "effective_radius", "solar_channels",
        )
        reflectance = np.empty((1, 1, 1, len(optical_depth), 3, 2), dtype=np.float32)
        for tau in range(len(optical_depth)):
            for radius in range(3):
                reflectance[0, 0, 0, tau, radius, 0] = 0.1 * tau + 0.01 * radius + offset
                reflectance[0, 0, 0, tau, radius, 1] = 0.2 * tau + 0.02 * radius + offset
        dataset.createVariable("R_0v", "f4", dims)[:] = reflectance


def test_extracts_channel_grids_and_displacement(tmp_path):
    first = tmp_path / "a.nc"
    second = tmp_path / "b.nc"
    _write_lut(first)
    _write_lut(second, offset=0.01)
    metadata_a = read_lut_metadata(first)
    metadata_b = read_lut_metadata(second)
    check_compatible(metadata_a, metadata_b)
    selection = resolve_selection(
        metadata_a,
        x_channel=1,
        y_channel=3,
        solar_zenith=30.0,
        satellite_zenith=0.0,
        relative_azimuth=0.0,
        surface_pressure=None,
    )
    data_a = extract_nakajima_king(metadata_a, selection)
    data_b = extract_nakajima_king(metadata_b, selection)
    assert data_a.x_reflectance.shape == (2, 3)
    np.testing.assert_allclose(data_a.x_reflectance[1], (0.1, 0.11, 0.12))
    dx, dy, magnitude = calculate_displacement(data_a, data_b)
    np.testing.assert_allclose(dx, 0.01, atol=1.0e-7)
    np.testing.assert_allclose(dy, 0.01, atol=1.0e-7)
    np.testing.assert_allclose(magnitude, np.sqrt(2.0) * 0.01, atol=1.0e-7)


def test_requires_explicit_pair_when_multiple_solar_channel_pairs_exist(tmp_path):
    path = tmp_path / "a.nc"
    _write_lut(path)
    metadata = replace(read_lut_metadata(path), solar_channel_ids=np.asarray([1, 2, 3]))
    with pytest.raises(ComparisonError, match="multiple solar-channel pairs are possible"):
        resolve_selection(
            metadata,
            x_channel=None,
            y_channel=None,
            solar_zenith=None,
            satellite_zenith=None,
            relative_azimuth=None,
            surface_pressure=None,
        )


def test_rejects_different_grids_without_interpolation(tmp_path):
    first = tmp_path / "a.nc"
    second = tmp_path / "b.nc"
    _write_lut(first)
    _write_lut(second, optical_depth=(1.0, 6.0))
    with pytest.raises(ComparisonError, match="coordinate 'optical_depth' differs"):
        check_compatible(read_lut_metadata(first), read_lut_metadata(second))


def test_rejects_materially_different_channel_wavelengths(tmp_path):
    first = tmp_path / "a.nc"
    second = tmp_path / "b.nc"
    _write_lut(first)
    _write_lut(second, wavelengths=(0.65, 1.6, 2.21))
    with pytest.raises(ComparisonError, match="central_wavelength arrays differ"):
        check_compatible(read_lut_metadata(first), read_lut_metadata(second))


def test_writes_four_panel_figure(tmp_path):
    first = tmp_path / "a.nc"
    second = tmp_path / "b.nc"
    output = tmp_path / "comparison.png"
    _write_lut(first)
    _write_lut(second, offset=0.01)
    result = compare_luts(
        first,
        second,
        output=output,
        x_channel=1,
        y_channel=3,
        solar_zenith=30.0,
        satellite_zenith=0.0,
        relative_azimuth=0.0,
    )
    assert result == output
    assert output.stat().st_size > 10_000
    assert plt.get_fignums() == []


def test_contour_lines_have_nakajima_king_orientation_and_labels(tmp_path):
    first = tmp_path / "a.nc"
    second = tmp_path / "b.nc"
    _write_lut(first)
    _write_lut(second, offset=0.01)
    metadata_a = read_lut_metadata(first)
    metadata_b = read_lut_metadata(second)
    selection = resolve_selection(
        metadata_a,
        x_channel=1,
        y_channel=3,
        solar_zenith=30.0,
        satellite_zenith=0.0,
        relative_azimuth=0.0,
        surface_pressure=None,
    )
    data_a = extract_nakajima_king(metadata_a, selection)
    figure = create_figure(data_a, extract_nakajima_king(metadata_b, selection))
    panel_a = figure.axes[0]
    optical_depth_count = data_a.optical_depth.size
    np.testing.assert_allclose(panel_a.lines[0].get_xdata(), data_a.x_reflectance[0, :])
    np.testing.assert_allclose(
        panel_a.lines[optical_depth_count].get_xdata(), data_a.x_reflectance[:, 0]
    )
    labels = [annotation.get_text() for annotation in panel_a.texts]
    assert any(label.startswith("τ=") for label in labels)
    assert any(label.startswith("rₑ=") and label.endswith("µm") for label in labels)
    plt.close(figure)
