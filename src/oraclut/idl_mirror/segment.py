"""Python mirror of the legacy IDL ``segment.pro`` routine.

IDL: ``segment, oldx, oldy, x, y, fe, minn=minn``

The IDL routine uses the IDL ``INT_TABULATED`` library function.  That
function first interpolates onto a uniformly spaced grid with a natural cubic
spline, then applies the five-point Newton--Cotes rule in groups of four
intervals.  Keeping that detail is important: replacing it with a trapezoid
integral changes which SRF samples are retained by the segmentation loop.
"""

import numpy as np


def _natural_spline_second_derivatives(x, y):
   """IDL ``SPL_INIT`` with its default natural end conditions."""

   x = np.asarray(x, dtype=np.float32)
   y = np.asarray(y, dtype=np.float32)
   n = x.size
   if n < 2:
      raise ValueError("INT_TABULATED requires at least two points")
   if np.any(np.diff(x) <= 0.0):
      raise ValueError("INT_TABULATED requires strictly ascending x values")

   # Numerical Recipes spline, which is the algorithm used by IDL SPL_INIT.
   y2 = np.zeros(n, dtype=np.float32)
   u = np.zeros(max(n - 1, 1), dtype=np.float32)
   for i in range(1, n - 1):
      sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
      p = sig * y2[i - 1] + np.float32(2.0)
      y2[i] = (sig - np.float32(1.0)) / p
      slope_change = ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                      (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
      u[i] = (np.float32(6.0) * slope_change / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p
   for k in range(n - 2, -1, -1):
      y2[k] = y2[k] * y2[k + 1] + u[k]
   return y2


def _spline_interpolate(x, y, y2, xout):
   """IDL ``SPL_INTERP`` for ascending tabulated points."""

   x = np.asarray(x, dtype=np.float32)
   y = np.asarray(y, dtype=np.float32)
   y2 = np.asarray(y2, dtype=np.float32)
   xout = np.asarray(xout, dtype=np.float32)
   result = np.empty_like(xout, dtype=np.float32)
   for index, value in np.ndenumerate(xout):
      if value <= x[0]:
         klo = 0
         khi = 1
      elif value >= x[-1]:
         klo = x.size - 2
         khi = x.size - 1
      else:
         klo = int(np.searchsorted(x, value, side="right") - 1)
         khi = klo + 1
      h = x[khi] - x[klo]
      a = (x[khi] - value) / h
      b = (value - x[klo]) / h
      result[index] = (a * y[klo] + b * y[khi] +
                       ((a ** 3 - a) * y2[klo] + (b ** 3 - b) * y2[khi]) * h ** 2 / np.float32(6.0))
   return result


def int_tabulated(x, y, sort=False):
   """Mirror IDL ``INT_TABULATED(X, Y)`` for float input arrays.

   ``sort=True`` is the narrowly scoped compatibility path used by the
   EarthCARE QM=2 diagnostic.  IDL accepts the local wavelength reversals in
   that SRF while retaining the original sample order for ``segment``; the
   tabulated integration itself orders each sub-array.  Other callers retain
   the strict ascending-input check by default.
   """

   x = np.asarray(x, dtype=np.float32)
   y = np.asarray(y, dtype=np.float32)
   if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
      raise ValueError("INT_TABULATED requires equal-length one-dimensional arrays")
   if x.size < 2:
      raise ValueError("INT_TABULATED requires at least two points")
   if sort:
      order = np.argsort(x, kind="stable")
      x = x[order]
      y = y[order]
   if np.any(np.diff(x) <= 0.0):
      raise ValueError("INT_TABULATED requires strictly ascending x values")

   # IDL: Xsegments = N_ELEMENTS(X)-1, rounded up to a multiple of four.
   xsegments = x.size - 1
   xsegments += (-xsegments) % 4
   h = (x[-1] - x[0]) / np.float32(xsegments)
   xgrid = h * np.arange(xsegments + 1, dtype=np.float32) + x[0]
   y2 = _natural_spline_second_derivatives(x, y)
   z = _spline_interpolate(x, y, y2, xgrid)

   # IDL: five-point Newton--Cotes / Boole rule, grouped by four intervals.
   end_indices = np.arange(4, z.size, 4, dtype=np.int32)
   value = np.float32(2.0) * h * (
      np.float32(7.0) * (z[end_indices - 4] + z[end_indices]) +
      np.float32(32.0) * (z[end_indices - 3] + z[end_indices - 1]) +
      np.float32(12.0) * z[end_indices - 2]
   ) / np.float32(45.0)
   return np.sum(value, dtype=np.float32)


def segment(oldx, oldy, fe, minn=None, sort=False):
   """Mirror ``create_orac_lut/segment.pro`` and return ``(x, y)``.

   ``sort`` is passed to each IDL-style tabulated integration, not applied to
   the source arrays before segmentation.  This distinction matters for the
   local reversals in the EarthCARE MSI channel-1 SRF.
   """

   oldx = np.asarray(oldx, dtype=np.float32).copy()
   oldy = np.asarray(oldy, dtype=np.float32).copy()
   if oldx.ndim != 1 or oldy.ndim != 1 or oldx.size != oldy.size or oldx.size < 2:
      raise ValueError("segment requires equal-length one-dimensional arrays with at least two points")
   if fe < 0.0:
      raise ValueError("segment requires a non-negative fractional error")

   # IDL: reverse descending input before the segmentation loop and restore it
   # at the end, preserving the original SRF order.
   swap = bool(oldx[0] > oldx[-1])
   if swap:
      oldx = oldx[::-1].copy()
      oldy = oldy[::-1].copy()

   integrate = lambda x, y: int_tabulated(x, y, sort=sort)
   true_total_area = integrate(oldx, oldy)
   true_trapezoidal_area = _trapezoidal_integral(oldx, oldy)
   n = oldx.size - 1
   max_segment_error = np.float32(0.1) * true_total_area

   while True:
      x = [oldx[0]]
      y = [oldy[0]]
      i = 0
      j = 0
      k = 0
      l = 1

      while True:
         j += 1
         area_old = integrate(oldx[i:j + 1], oldy[i:j + 1])
         xt = np.asarray(x + [oldx[j]], dtype=np.float32)
         yt = np.asarray(y + [oldy[j]], dtype=np.float32)
         area_new = integrate(xt[k:l + 1], yt[k:l + 1])
         if abs(area_new - area_old) > max_segment_error:
            # IDL: accept the point immediately before the rejected point.
            x.append(oldx[j - 1])
            y.append(oldy[j - 1])
            k = l
            l += 1
            i = j - 1
            j -= 1
         if j == n:
            break

      x.append(oldx[n])
      y.append(oldy[n])
      x = np.asarray(x, dtype=np.float32)
      y = np.asarray(y, dtype=np.float32)

      trapezoidal_area = _trapezoidal_integral(x, y)
      least_squares = bool(trapezoidal_area == true_trapezoidal_area)
      current_fe = min(
         abs(true_total_area - integrate(x, y)) / true_total_area,
         abs(trapezoidal_area - true_trapezoidal_area) / true_trapezoidal_area,
      )
      max_segment_error *= np.float32(0.9)
      number_ok = True if minn is None else x.size >= int(minn)
      if (current_fe <= fe and number_ok) or least_squares:
         break

   if swap:
      x = x[::-1].copy()
      y = y[::-1].copy()
   return x, y


def _trapezoidal_integral(x, y):
   """IDL: integrate_trapeziodal.pro."""

   x = np.asarray(x, dtype=np.float32)
   y = np.asarray(y, dtype=np.float32)
   return np.sum((x[1:] - x[:-1]) * (np.float32(0.5) * (y[:-1] + y[1:])), dtype=np.float32)
