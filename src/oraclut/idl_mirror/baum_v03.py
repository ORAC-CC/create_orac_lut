"""Reader for the Baum severely-rough ice optical-property tables.

This is the Python counterpart of the legacy ``read_baum_lambda`` routine in
``create_orac_lut/read_baum.pro``.  The supplied tables store fields as
``(effective_diameter, wavelength)`` and P11 as
``(phase_angle, effective_diameter, wavelength)``.  ORAC uses
volume-equivalent/effective radius, so the legacy reader converts diameter to
radius before interpolation.

The arithmetic is that of read_baum.pro in IDL: FLOAT tables, FLOAT fractional
indices from INTERPOL, the FLOAT radius extrapolation, and INTERPOLATE /GRID
evaluated in double precision (validation/v21_v22_terra_modis/
REPORT_single_precision_equivalence.md).
"""

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from .interpol import idl_interpol


def _index(axis, targets):
   """IDL: interpol(indgen(n_elements(axis)), axis, targets), the fractional grid index.

   indgen is INT, so the arithmetic is FLOAT for FLOAT targets and DOUBLE for
   DOUBLE targets (the phase angles), as in IDL.
   """

   return idl_interpol(np.arange(np.asarray(axis).size, dtype=np.int16), axis, targets)


def _grid_interpolate(table, indices):
   """IDL: interpolate(table, i_1, i_2[, i_3], /grid) -- multilinear, edge-clamped.

   IDL converts the location arrays to FLOAT (the DOUBLE phase-angle index
   included), evaluates the interpolation in double precision and returns the
   type of the table (FLOAT); this reproduces the V21 LUT Baum extinction,
   albedo, asymmetry and the radiative-transfer operators of the pagg LUT
   bitwise.
   """

   result = np.asarray(table, dtype=np.float64)
   for axis, index in enumerate(indices):
      n = result.shape[axis]
      position = np.clip(np.asarray(index, dtype=np.float32).astype(np.float64), 0.0, n - 1.0)
      lower = np.minimum(np.floor(position).astype(int), n - 2)
      shape = [1] * result.ndim
      shape[axis] = position.size
      fraction = (position - lower).reshape(shape)
      low = np.take(result, lower, axis=axis)
      high = np.take(result, lower + 1, axis=axis)
      result = low + fraction * (high - low)
   return result.astype(np.float32)


def _extend_radius(field, re_baum, re_last):
   """IDL: temp2[*, n_Re_baum] = (1.-f_Re_extrap)*temp2[*, n_Re_baum-2] + f_Re_extrap*temp2[*, n_Re_baum-1] (FLOAT).

   ``field`` has the radius on axis 1 (the IDL second index).
   """

   n = re_baum.size
   f = (np.float32(re_last) - re_baum[n - 2]) / (re_baum[n - 1] - re_baum[n - 2])
   extra = (np.float32(1.0) - f) * np.take(field, n - 2, axis=1) + f * np.take(field, n - 1, axis=1)
   return np.concatenate([field, np.expand_dims(extra, 1)], axis=1).astype(np.float32)


class BaumTable:
   """In-memory generic Baum table with the legacy interpolation semantics."""

   def __init__(self, filename):
      self.filename = Path(filename)
      if not self.filename.is_file() and not self.filename.is_absolute():
         # ``.mm`` files retain the legacy repository-relative
         # ``input_files/baum/...`` spelling.
         repository_root = Path(__file__).resolve().parents[3]
         for repository_file in (repository_root / self.filename,
                                 repository_root / "create_orac_lut" / self.filename):
            if repository_file.is_file():
               self.filename = repository_file
               break
      if not self.filename.is_file():
         raise FileNotFoundError(f"Baum optical-property file not found: {self.filename}")
      with Dataset(self.filename, "r") as dataset:
         # The tables are FLOAT, as ncdf_varget returns them to IDL; Re_baum = De / 2. (FLOAT)
         self.wavelengths = np.asarray(dataset.variables["wavelengths"][:], dtype=np.float32)
         self.effective_radii = np.asarray(dataset.variables["effective_diameter"][:], dtype=np.float32) / np.float32(2.0)
         self.phase_angles = np.asarray(dataset.variables["phase_angles"][:], dtype=np.float32)
         self.extinction_efficiency = np.asarray(dataset.variables["extinction_efficiency"][:], dtype=np.float32)
         self.single_scattering_albedo = np.asarray(dataset.variables["single_scattering_albedo"][:], dtype=np.float32)
         self.asymmetry_parameter = np.asarray(dataset.variables["asymmetry_parameter"][:], dtype=np.float32)
         self.p11_phase_function = np.asarray(dataset.variables["p11_phase_function"][:], dtype=np.float32)

      expected = (self.effective_radii.size, self.wavelengths.size)
      for name in ("extinction_efficiency", "single_scattering_albedo", "asymmetry_parameter"):
         if getattr(self, name).shape != expected:
            raise ValueError(f"{self.filename}: Baum field {name} has unexpected shape")
      expected_phase = (self.phase_angles.size, *expected)
      if self.p11_phase_function.shape != expected_phase:
         raise ValueError(f"{self.filename}: Baum P11 field has unexpected shape")

   def interpolate(self, wavelengths, effective_radii, phase_angles):
      """Return ``Bext, w, g, P11`` (FLOAT) in ORAC's ``(wavelength, radius)`` and ``(angle, wavelength, radius)`` order.

      IDL read_baum_lambda: the wavelengths and radii are FLOAT, the phase
      angles keep their type (DOUBLE from the generator).  If the largest LUT
      radius exceeds the table, the table is extended to it by linear
      extrapolation in FLOAT before interpolating.
      """

      wavelengths = np.asarray(wavelengths, dtype=np.float32).reshape(-1)
      effective_radii = np.asarray(effective_radii, dtype=np.float32).reshape(-1)
      phase_angles = np.asarray(phase_angles).reshape(-1)
      re_baum = self.effective_radii
      # file order (radius, wavelength) and (angle, radius, wavelength); IDL (wl, Re) and (wl, Re, theta)
      fields = [f.T for f in (self.extinction_efficiency, self.single_scattering_albedo, self.asymmetry_parameter)]
      p11 = np.transpose(self.p11_phase_function, (2, 1, 0))
      if re_baum[-1] < effective_radii[-1]:
         fields = [_extend_radius(f, re_baum, effective_radii[-1]) for f in fields]
         p11 = _extend_radius(p11, re_baum, effective_radii[-1])
         re_baum = np.append(re_baum, effective_radii[-1]).astype(np.float32)
      i_wl = _index(self.wavelengths, wavelengths)
      i_re = _index(re_baum, effective_radii)
      i_theta = _index(self.phase_angles, phase_angles)
      bext, albedo, asymmetry = (_grid_interpolate(f, (i_wl, i_re)) for f in fields)
      phase = np.transpose(_grid_interpolate(p11, (i_wl, i_re, i_theta)), (2, 0, 1))   # IDL: P[i,*,*] = temp[*,*,i]
      return bext, albedo, asymmetry, phase


def read_baum_lambda(filename, wavelengths, effective_radii, phase_angles):
   """Read/interpolate one Baum full-phase-matrix table like IDL ``read_baum_lambda``."""

   return BaumTable(filename).interpolate(wavelengths, effective_radii, phase_angles)
