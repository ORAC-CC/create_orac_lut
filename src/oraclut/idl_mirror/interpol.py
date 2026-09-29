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
