"""Regression tests for the reusable V21/V22 LUT comparison utility."""

import importlib.util
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "validation" / "v21_v22_terra_modis" / "compare_luts.py"
SPEC = importlib.util.spec_from_file_location("compare_luts", MODULE_PATH)
COMPARE_LUTS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE_LUTS)


def _write_pair(reference_path, candidate_path):
   with Dataset(reference_path, "w") as dataset:
      dataset.createDimension("x", 2)
      dataset.createDimension("channels", 2)
      dataset.createVariable("x", "f4", ("x",))[:] = [1.0, 2.0]
      dataset.createVariable("channel_id", "i2", ("channels",))[:] = [7, 9]
      values = dataset.createVariable("values", "f4", ("x", "channels"), fill_value=np.float32(-999.0))
      values.units = "1"
      values[:] = [[1.0, 2.0], [3.0, -999.0]]
   with Dataset(candidate_path, "w") as dataset:
      dataset.createDimension("channels", 2)
      dataset.createDimension("x", 2)
      dataset.createVariable("x", "f4", ("x",))[:] = [1.0, 2.0]
      dataset.createVariable("channel_id", "i2", ("channels",))[:] = [7, 9]
      values = dataset.createVariable("values", "f4", ("channels", "x"), fill_value=np.float32(-999.0))
      values.units = "1"
      values[:] = [[1.0, 3.5], [2.0, -999.0]]


def test_compare_aligns_dimensions_and_reports_coordinates(tmp_path):
   reference = tmp_path / "reference.nc"
   candidate = tmp_path / "candidate.nc"
   _write_pair(reference, candidate)
   result, coordinate_rows, calls = COMPARE_LUTS.compare(reference, candidate)
   values = next(row for row in result["comparisons"] if row["variable"] == "values")
   assert values["candidate_transpose_order"] == [1, 0]
   assert values["jointly_valid_elements"] == 3
   assert values["differing_elements"] == 1
   assert values["maximum_absolute_difference"] == 0.5
   assert values["maximum_absolute_location"] == {
      "x": {"index": 1, "coordinate": 2.0},
      "channels": {"index": 0, "coordinate": 7},
   }
   assert values["mask_difference"] == 0
   assert {row["dimension"] for row in coordinate_rows if row["variable"] == "values"} == {"x", "channels"}
   assert calls == []


def test_compare_records_mask_mismatch(tmp_path):
   reference = tmp_path / "reference.nc"
   candidate = tmp_path / "candidate.nc"
   _write_pair(reference, candidate)
   with Dataset(candidate, "a") as dataset:
      dataset.variables["values"][1, 1] = np.float32(4.0)
   result, _, _ = COMPARE_LUTS.compare(reference, candidate)
   values = next(row for row in result["comparisons"] if row["variable"] == "values")
   assert values["reference_masked"] == 1
   assert values["candidate_masked"] == 0
   assert values["mask_difference"] == 1
