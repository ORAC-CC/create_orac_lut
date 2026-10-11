"""Reader for the Baum severely-rough ice optical-property tables.

This is the Python counterpart of the legacy ``read_baum_lambda`` routine in
``create_orac_lut/read_baum.pro``.  The supplied tables store fields as
``(effective_diameter, wavelength)`` and P11 as
``(phase_angle, effective_diameter, wavelength)``.  ORAC uses
volume-equivalent/effective radius, so the legacy reader converts diameter to
radius before interpolation.
"""

from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def _linear_indices(axis, targets):
   """Return lower indices and fractions for linear interpolation/extrapolation."""

   axis = np.asarray(axis, dtype=np.float64)
   targets = np.asarray(targets, dtype=np.float64).reshape(-1)
   if axis.ndim != 1 or axis.size < 2 or np.any(np.diff(axis) <= 0):
      raise ValueError("Baum interpolation axes must be strictly increasing")
   lower = np.searchsorted(axis, targets, side="right") - 1
   lower = np.clip(lower, 0, axis.size - 2)
   fraction = (targets - axis[lower]) / (axis[lower + 1] - axis[lower])
   # IDL INTERPOL clamps below the first tabulated point.  The historical
   # reader explicitly extends only the upper effective-radius boundary before
   # calling INTERPOL, so lower-bound extrapolation must not be performed.
   fraction = np.maximum(fraction, 0.0)
   return lower, fraction


def _interpolate_axis(values, axis_values, targets, axis):
   """Linearly interpolate one array axis, including endpoint extrapolation."""

   lower, fraction = _linear_indices(axis_values, targets)
   low = np.take(values, lower, axis=axis)
   high = np.take(values, lower + 1, axis=axis)
   shape = [1] * values.ndim
   shape[axis] = fraction.size
   fraction = fraction.reshape(shape)
   return low + fraction * (high - low)


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
         self.wavelengths = np.asarray(dataset.variables["wavelengths"][:], dtype=np.float64)
         self.effective_radii = np.asarray(dataset.variables["effective_diameter"][:], dtype=np.float64) / 2.0
         self.phase_angles = np.asarray(dataset.variables["phase_angles"][:], dtype=np.float64)
         self.extinction_efficiency = np.asarray(dataset.variables["extinction_efficiency"][:], dtype=np.float64)
         self.single_scattering_albedo = np.asarray(dataset.variables["single_scattering_albedo"][:], dtype=np.float64)
         self.asymmetry_parameter = np.asarray(dataset.variables["asymmetry_parameter"][:], dtype=np.float64)
         self.p11_phase_function = np.asarray(dataset.variables["p11_phase_function"][:], dtype=np.float64)

      expected = (self.effective_radii.size, self.wavelengths.size)
      for name in ("extinction_efficiency", "single_scattering_albedo", "asymmetry_parameter"):
         if getattr(self, name).shape != expected:
            raise ValueError(f"{self.filename}: Baum field {name} has unexpected shape")
      expected_phase = (self.phase_angles.size, *expected)
      if self.p11_phase_function.shape != expected_phase:
         raise ValueError(f"{self.filename}: Baum P11 field has unexpected shape")

   def interpolate(self, wavelengths, effective_radii, phase_angles):
      """Return ``Bext, w, g, P11`` in ORAC's ``(wavelength, radius)`` order."""

      wavelengths = np.asarray(wavelengths, dtype=np.float64).reshape(-1)
      effective_radii = np.asarray(effective_radii, dtype=np.float64).reshape(-1)
      phase_angles = np.asarray(phase_angles, dtype=np.float64).reshape(-1)

      def scalar_field(field):
         by_wavelength = _interpolate_axis(field, self.wavelengths, wavelengths, axis=1)
         by_radius = _interpolate_axis(by_wavelength, self.effective_radii, effective_radii, axis=0)
         return by_radius.T

      bext = scalar_field(self.extinction_efficiency)
      albedo = scalar_field(self.single_scattering_albedo)
      asymmetry = scalar_field(self.asymmetry_parameter)

      by_wavelength = _interpolate_axis(self.p11_phase_function, self.wavelengths, wavelengths, axis=2)
      by_radius = _interpolate_axis(by_wavelength, self.effective_radii, effective_radii, axis=1)
      by_angle = _interpolate_axis(by_radius, self.phase_angles, phase_angles, axis=0)
      phase = np.transpose(by_angle, (0, 2, 1))
      return bext, albedo, asymmetry, phase


def read_baum_lambda(filename, wavelengths, effective_radii, phase_angles):
   """Read/interpolate one Baum full-phase-matrix table like IDL ``read_baum_lambda``."""

   return BaumTable(filename).interpolate(wavelengths, effective_radii, phase_angles)
