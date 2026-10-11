"""Numerically compatible subset of the IDL ``TOTAL`` function.

The V21 LUT generator uses scalar ``TOTAL`` reductions and ``/DOUBLE``.  The
dimension and cumulative forms are included because they define the ordering
unambiguously and make the compatibility routine independently testable.  IDL
dimensions are one-based and IDL arrays are flattened with their first index
varying fastest.

``/INTEGER``, ``/NAN`` and ``/PRESERVE_TYPE`` are not used by the V21 LUT
calculation and are deliberately outside this interface.
"""

import numpy as np


def _accumulator_dtype(values, double):
   """Return the dtype used by the supported IDL TOTAL invocation."""

   if values.dtype.kind not in "biufc":
      raise TypeError("IDL TOTAL requires a numerical array")
   if double:
      return np.dtype(np.complex128 if values.dtype.kind == "c" else np.float64)
   if values.dtype == np.dtype(np.complex128):
      return np.dtype(np.complex128)
   if values.dtype.kind == "c":
      return np.dtype(np.complex64)
   if values.dtype == np.dtype(np.float64):
      return np.dtype(np.float64)
   if values.dtype.kind in "iu" and values.dtype.itemsize == 8:
      return np.dtype(np.float64)
   return np.dtype(np.float32)


def idl_total(array, dimension=None, *, cumulative=False, double=False):
   """Return the IDL-compatible total of ``array`` for V21 LUT use.

   Parameters mirror the required IDL forms ``TOTAL(Array [, Dimension],
   [/CUMULATIVE], [/DOUBLE])``.  ``dimension`` is one-based; zero or ``None``
   reduces the complete array in IDL column-major order.  Accumulation is
   sequential through ``numpy.cumsum`` rather than NumPy's pairwise ``sum``.
   """

   values = np.asarray(array)
   if values.size == 0:
      raise ValueError("IDL TOTAL compatibility does not support an empty array")
   dtype = _accumulator_dtype(values, bool(double))

   if dimension is None or dimension == 0:
      flattened = values.ravel(order="F")
      accumulated = np.cumsum(flattened, dtype=dtype)
      if cumulative and values.ndim > 0:
         return accumulated.reshape(values.shape, order="F")
      return accumulated[-1]

   if isinstance(dimension, (bool, np.bool_)) or not isinstance(dimension, (int, np.integer)):
      raise TypeError("IDL TOTAL dimension must be an integer")
   axis = int(dimension) - 1
   if axis < 0 or axis >= values.ndim:
      raise ValueError(f"IDL TOTAL dimension {dimension} is invalid for a {values.ndim}-dimensional array")
   accumulated = np.cumsum(values, axis=axis, dtype=dtype)
   if cumulative:
      return accumulated
   return np.take(accumulated, -1, axis=axis)
