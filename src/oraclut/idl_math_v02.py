"""IDL-compatible elementary arithmetic used by the V21 LUT generator.

NumPy's float32 and float64 ``cos``, ``exp`` and ``power`` use NumPy's own
(SIMD) implementations, which differ from the C library by an ulp at some
arguments.  IDL 8.9 evaluates these functions with glibc libm (``libidl.so``
imports ``cos``, ``cosf``, ``exp``, ``expf``, ``pow`` and ``powf``); FLOAT
``EXP`` is the DOUBLE ``exp`` rounded to FLOAT (this, not ``expf``, reproduces
the V21 T_dv operators bitwise), and integer powers are evaluated by repeated
squaring.  The routines here reproduce
the IDL arithmetic; each is used where the IDL source has the corresponding
expression (validation/v21_v22_terra_modis/REPORT_single_precision_equivalence.md).

* ``idl_cos_deg(x)``  IDL ``cos(x*!dtor)`` for FLOAT x: FLOAT product with the
  FLOAT system variable !DTOR (0.017453292, bits 0x3C8EFA35), then ``cosf``.
* ``idl_exp(x)``      IDL ``exp(x)``: glibc ``exp``, rounded to FLOAT for FLOAT x.
* ``idl_pow(x, y)``   IDL ``x^y`` for DOUBLE x and a DOUBLE exponent: ``pow``.
* ``idl_ipow(x, n)``  IDL ``x^n`` for an INTEGER exponent n >= 1: repeated
  squaring in the type of x (x^2 = x*x, x^4 = (x*x)*(x*x)).
"""

import ctypes
import math

import numpy as np

IDL_DTOR = np.float32(0.0174532925)      # IDL: !DTOR (FLOAT)

_LIBM = ctypes.CDLL("libm.so.6")         # glibc libm, as called by IDL
_LIBM.cosf.restype = ctypes.c_float
_LIBM.cosf.argtypes = [ctypes.c_float]


def _elementwise(function, values, dtype):
   values = np.asarray(values, dtype=dtype)
   result = np.empty(values.shape, dtype=dtype)
   for index in np.ndindex(values.shape):
      result[index] = function(float(values[index]))
   if result.ndim == 0:
      return dtype(result)
   return result


def idl_cos_deg(angle):
   """Return IDL ``cos(angle*!dtor)`` for FLOAT ``angle`` (a float32 scalar or array)."""

   radians = np.asarray(angle, dtype=np.float32) * IDL_DTOR    # FLOAT product, rounded as in IDL
   return _elementwise(_LIBM.cosf, radians, np.float32)


def idl_exp(x):
   """Return IDL ``exp(x)``: glibc ``exp``, rounded to float32 for float32 ``x`` (not ``expf``)."""

   x = np.asarray(x)
   if x.dtype == np.float32:
      return _elementwise(lambda value: np.float32(math.exp(value)), x, np.float32)
   if x.dtype == np.float64:
      return _elementwise(math.exp, x, np.float64)
   raise TypeError(f"idl_exp supports FLOAT and DOUBLE values, not {x.dtype}")


def idl_pow(x, exponent):
   """Return IDL ``x^exponent`` for float64 ``x`` and a float64 scalar exponent (glibc ``pow``)."""

   exponent = float(exponent)
   return _elementwise(lambda value: math.pow(value, exponent), x, np.float64)


def idl_ipow(x, n):
   """Return IDL ``x^n`` for an integer exponent ``n`` >= 1, by repeated squaring in the type of ``x``."""

   n = int(n)
   if n < 1:
      raise ValueError("idl_ipow supports exponents n >= 1")
   base = np.asarray(x)
   result = None
   while n:
      if n & 1:
         result = base if result is None else result * base
      n >>= 1
      if n:
         base = base * base
   return result
