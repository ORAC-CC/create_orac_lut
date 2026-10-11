"""Python counterpart of the IDL built-in INTERPOL(V, X, XOUT) as used for the particle profile.

Both legacy generators call

   scatreltau = INTERPOL(mmstr.rext, mmstr.height, hlayers)

IDL INTERPOL is piecewise linear inside the tabulated range and, because the
segment index is clamped to the end segments, extrapolates linearly from the
two nearest nodes below the first node and above the last one.  Ascending or
descending abscissae are accepted.  The arithmetic is single precision, as in
the legacy run.  This behaviour was established against IDL 8.9 and is the
validated correction recorded in docs/validation_matrix.md (Finding 2).
"""

import numpy as np


def interpol(v, x, xout):
   """Return V interpolated from abscissae X onto XOUT (float32)."""

   x = np.asarray(x, dtype=np.float32)
   v = np.asarray(v, dtype=np.float32)
   xout = np.asarray(xout, dtype=np.float32)
   if x.ndim != 1 or x.size < 2 or v.shape != x.shape:
      raise ValueError("interpol requires at least two (x, v) nodes of equal length")
   order = np.argsort(x, kind="stable")
   x, v = x[order], v[order]
   if np.any(np.diff(x) <= 0):
      raise ValueError("interpol abscissae must be strictly monotonic")
   # IDL: s = VALUE_LOCATE(x, xout) > 0L < (m-2)
   segment = np.clip(np.searchsorted(x, xout, side="right") - 1, 0, x.size - 2)
   x0, x1 = x[segment], x[segment + 1]
   v0, v1 = v[segment], v[segment + 1]
   return (v0 + (xout - x0) * (v1 - v0) / (x1 - x0)).astype(np.float32)


def idl_interpol(v, x, xout):
   """IDL INTERPOL(V, X, XOUT), linear and irregular, in the IDL result type of the operands.

   s = VALUE_LOCATE(X, XOUT) > 0L < (m-2);  p = (XOUT-X[s])*diff/(X[s+1]-X[s]) + V[s],
   diff = V[s+1] - V[s].  The arithmetic type is that of the operands together
   (FLOAT if all are FLOAT, DOUBLE if any is DOUBLE), X may be increasing or
   decreasing (VALUE_LOCATE), and values outside X are extrapolated from the
   end intervals, as in IDL.
   """

   v = np.asarray(v)
   x = np.asarray(x)
   xout = np.asarray(xout)
   dtype = np.result_type(v.dtype, x.dtype, xout.dtype)
   v = v.astype(dtype)
   x = x.astype(dtype)
   xo = xout.astype(dtype)
   m = x.size
   if x.ndim != 1 or m < 2 or v.shape != x.shape:
      raise ValueError("interpol requires at least two (x, v) nodes of equal length")
   if np.all(np.diff(x) > 0):
      s = np.searchsorted(x, xo, side="right") - 1
   elif np.all(np.diff(x) < 0):
      s = np.searchsorted(-x, -xo, side="right") - 1
   else:
      raise ValueError("interpol abscissae must be strictly monotonic")
   s = np.clip(s, 0, m - 2)
   return (xo - x[s]) * (v[s + 1] - v[s]) / (x[s + 1] - x[s]) + v[s]
