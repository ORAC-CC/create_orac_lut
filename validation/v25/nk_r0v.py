"""Nakajima-King R_0v comparison of two compact LUTs (validation only).

Runs validation/nakajima_king/compare_nakajima_king.py unchanged except that
the coordinate and channel vectors are read without netCDF4's valid_range
masking: a compact LUT holding a subset of an instrument's channels has
channel IDs above its channel count, which the legacy writer's
valid_range (0, nchan) would mask.  Arguments as the original tool.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nakajima_king"))
import compare_nakajima_king as nk   # noqa: E402


def _numeric_vector(dataset, name):
   variable = nk._require_variable(dataset, name)
   variable.set_auto_mask(False)
   values = np.asarray(variable[:], dtype=np.float64)
   if values.ndim != 1:
      raise nk.ComparisonError(f"{dataset.filepath()}: variable {name!r} must be one-dimensional")
   return values


nk._numeric_vector = _numeric_vector

if __name__ == "__main__":
   sys.exit(nk.main())
