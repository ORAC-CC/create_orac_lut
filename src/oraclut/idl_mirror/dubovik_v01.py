"""Readers and numerical mirror of the historical Dubovik T-matrix tables.

These tables are the optical-property backend used by the legacy
``dubovik_lognormal_multiple_eps.pro`` routine.  The files are ASCII and are
deliberately read as records here rather than converted to a new format: this
keeps the external database as the scientific source of truth.

IDL: ``dubovik_read_filenames.pro``, ``dubovik_read_kernel_kext.pro``,
``dubovik_read_kernel_scatt_matrix.pro``, ``dubovik_lognormal.pro`` and
``dubovik_lognormal_multiple_eps.pro``.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


VALID_EPS = np.asarray(
   [0.33490, 0.36690, 0.40190, 0.44030, 0.48230, 0.52830, 0.57870,
    0.63390, 0.69440, 0.76070, 0.83330, 0.91290, 1.00000, 1.09540,
    1.20000, 1.31450, 1.44000, 1.57740, 1.72800, 1.89290, 2.07360,
    2.27150, 2.48832, 2.72580, 2.98600], dtype=np.float64)


@dataclass(frozen=True)
class DubovikFilenames:
   """The seven files associated with each tabulated aspect ratio."""

   kext: tuple[str, ...]
   k11: tuple[str, ...]
   k12: tuple[str, ...]
   k22: tuple[str, ...]
   k33: tuple[str, ...]
   k34: tuple[str, ...]
   k44: tuple[str, ...]


@dataclass(frozen=True)
class DubovikGrid:
   radius: np.ndarray
   wavelength: float


@dataclass(frozen=True)
class DubovikKext:
   ri_real: np.ndarray
   ri_imag: np.ndarray
   radius: np.ndarray
   wvl: float
   epsilon: np.ndarray | float
   kext: np.ndarray
   kabs: np.ndarray
   eps_distribution: np.ndarray | None = None
   key_rd: int | None = None


@dataclass(frozen=True)
class DubovikScatteringMatrix:
   ri_real: np.ndarray
   ri_imag: np.ndarray
   radius: np.ndarray
   theta: np.ndarray
   wvl: float
   epsilon: np.ndarray | float
   k11: np.ndarray
   eps_distribution: np.ndarray | None = None
   key_rd: int | None = None


def _numbers(line):
   """Parse a whitespace-separated numeric record, accepting Fortran D exponents."""

   return [float(word.replace("D", "E").replace("d", "e")) for word in line.split()]


def _first_numbers(line, count, description):
   values = []
   for word in line.split():
      try:
         values.append(float(word.replace("D", "E").replace("d", "e")))
      except ValueError:
         break
      if len(values) == count:
         break
   if len(values) != count:
      raise ValueError(f"Dubovik {description}: expected {count} numeric values")
   return values


class _NumericStream:
   """Consume numeric arrays spanning the wrapped lines used by the tables."""

   def __init__(self, handle):
      self.handle = handle
      self.pending = []

   def take(self, count, description):
      values = []
      while len(values) < count:
         if self.pending:
            needed = count - len(values)
            values.extend(self.pending[:needed])
            self.pending = self.pending[needed:]
            continue
         line = self.handle.readline()
         if line == "":
            raise ValueError(f"Dubovik {description}: unexpected end of file")
         if not line.strip():
            continue
         try:
            self.pending.extend(_numbers(line))
         except ValueError as exc:
            raise ValueError(f"Dubovik {description}: expected numeric data, got {line!r}") from exc
      return np.asarray(values, dtype=np.float64)


def read_dubovik_filenames(file):
   """Read ``ROUTINE/name.dat`` in the order used by the IDL routine."""

   lines = Path(file).read_text().splitlines()
   if not lines:
      raise ValueError(f"Dubovik filename list is empty: {file}")
   try:
      nfiles = int(lines[0].split()[0])
   except (TypeError, ValueError) as exc:
      raise ValueError(f"Dubovik filename count is invalid in {file}") from exc
   # The distributed file contains historical secondary lists after the
   # first 25*7 records.  IDL reads exactly the requested records and leaves
   # that trailing material untouched.
   names = [line.strip() for line in lines[1:] if line.strip()][:7 * nfiles]
   expected = 7 * nfiles
   if len(names) != expected:
      raise ValueError(f"Dubovik filename list has {len(names)} names; expected {expected}")
   records = [[names[7 * i + j] for j in range(7)] for i in range(nfiles)]
   # name.dat stores K11..K44 followed by Kext; the public structure follows
   # the IDL keyword order, which starts with Kext.
   columns = [tuple(record[j] for record in records) for j in (6, 0, 1, 2, 3, 4, 5)]
   return DubovikFilenames(*columns)


def read_dubovik_grid(file):
   """Read ``grid1.dat`` or ``grid1.dat.fix``."""

   values = Path(file).read_text().split()
   if len(values) < 2:
      raise ValueError(f"Dubovik grid is incomplete: {file}")
   count = int(float(values[0]))
   if len(values) < count + 2:
      raise ValueError(f"Dubovik grid has fewer than {count} radii: {file}")
   return DubovikGrid(np.asarray(values[2:count + 2], dtype=np.float64), float(values[1]))


def _read_common_header(handle, grid, scattering):
   first = _first_numbers(handle.readline(), 3, "kernel header")
   _first_numbers(handle.readline(), 1, "interval count")
   if scattering:
      ntheta = int(_first_numbers(handle.readline(), 1, "angle count")[0])
      theta_stream = _NumericStream(handle)
      theta = theta_stream.take(ntheta, "angle grid")
   else:
      theta = None
   return first, theta


def _read_kernel_metadata(handle, fixed):
   if not fixed:
      return None, None, None
   key_rd = int(_first_numbers(handle.readline(), 1, "fixed key_RD")[0])
   _first_numbers(handle.readline(), 2, "fixed header")
   neps = int(_first_numbers(handle.readline(), 1, "fixed epsilon count")[0])
   handle.readline()  # description
   eps_rows = [_first_numbers(handle.readline(), 2, "fixed epsilon row") for _ in range(neps)]
   epsilon = np.asarray([row[0] for row in eps_rows], dtype=np.float64)
   eps_distribution = np.asarray([row[1] for row in eps_rows], dtype=np.float64)
   return key_rd, epsilon, eps_distribution


def _read_ri_header(handle):
   handle.readline()  # real RI description
   handle.readline()  # imaginary RI description
   nreal, nimag = _first_numbers(handle.readline(), 2, "RI grid dimensions")
   return int(nreal), abs(int(nimag))


def read_dubovik_kernel_kext(file, grid1=None):
   """Read one Dubovik extinction/absorption kernel.

   Returned array order is ``[real_RI, imaginary_RI, radius]``, matching the
   IDL structures consumed by ``INTERPOL``.
   """

   file = Path(file)
   fixed = file.name.endswith(".fix")
   if grid1 is None:
      grid1 = file.parent / ("grid1.dat.fix" if fixed else "grid1.dat")
   grid = read_dubovik_grid(grid1)
   with file.open() as handle:
      if fixed:
         key_rd, epsilon, eps_distribution = _read_kernel_metadata(handle, fixed)
         nr = abs(int(_first_numbers(handle.readline(), 1, "radius count")[0]))
      else:
         key_rd = None
         epsilon = _first_numbers(handle.readline(), 3, "kernel header")[2]
         nr = abs(int(_first_numbers(handle.readline(), 1, "radius count")[0]))
         eps_distribution = None
      if nr != grid.radius.size:
         raise ValueError(f"Dubovik radius count mismatch in {file}: {nr} != {grid.radius.size}")
      nreal, nimag = _read_ri_header(handle)
      kext = np.empty((nreal, nimag, nr), dtype=np.float64)
      kabs = np.empty_like(kext)
      ri_real = np.empty(nreal, dtype=np.float32)
      ri_imag = np.empty(nimag, dtype=np.float32)
      stream = _NumericStream(handle)
      wvl = None
      for ir in range(nreal):
         for ik in range(nimag):
            element_line = handle.readline()
            if not element_line:
               raise ValueError(f"Dubovik {file}: missing kernel record")
            _first_numbers(element_line, 1, "element record")
            metadata = _first_numbers(handle.readline(), 3, "wavelength/RI record")
            wvl, ri_real[ir], ri_imag[ik] = metadata
            handle.readline()  # EXTINCTION label
            kext[ir, ik, :] = stream.take(nr, "extinction row")
            handle.readline()  # ABSORPTION label
            kabs[ir, ik, :] = stream.take(nr, "absorption row")
   return DubovikKext(ri_real, ri_imag, grid.radius, float(wvl), epsilon, kext, kabs,
                      eps_distribution, key_rd)


def read_dubovik_kernel_scatt_matrix(file, grid1=None):
   """Read one Dubovik scattering-matrix kernel (normally the K11 file)."""

   file = Path(file)
   fixed = file.name.endswith(".fix")
   if grid1 is None:
      grid1 = file.parent / ("grid1.dat.fix" if fixed else "grid1.dat")
   grid = read_dubovik_grid(grid1)
   with file.open() as handle:
      if fixed:
         key_rd, epsilon, eps_distribution = _read_kernel_metadata(handle, fixed)
         nr = abs(int(_first_numbers(handle.readline(), 1, "radius count")[0]))
         ntheta = int(_first_numbers(handle.readline(), 1, "angle count")[0])
         theta = _NumericStream(handle).take(ntheta, "angle grid")
      else:
         key_rd = None
         epsilon = _first_numbers(handle.readline(), 3, "kernel header")[2]
         nr = abs(int(_first_numbers(handle.readline(), 1, "radius count")[0]))
         ntheta = int(_first_numbers(handle.readline(), 1, "angle count")[0])
         theta = _NumericStream(handle).take(ntheta, "angle grid")
         eps_distribution = None
      if nr != grid.radius.size:
         raise ValueError(f"Dubovik radius count mismatch in {file}: {nr} != {grid.radius.size}")
      # The angle count is repeated in the file header and is retained here
      # rather than inferred from the grid directory.
      nreal, nimag = _read_ri_header(handle)
      k11 = np.empty((nreal, nimag, nr, theta.size), dtype=np.float64)
      ri_real = np.empty(nreal, dtype=np.float32)
      ri_imag = np.empty(nimag, dtype=np.float32)
      stream = _NumericStream(handle)
      wvl = None
      for ir in range(nreal):
         for ik in range(nimag):
            element_line = handle.readline()
            if not element_line:
               raise ValueError(f"Dubovik {file}: missing kernel record")
            _first_numbers(element_line, 1, "element record")
            metadata = _first_numbers(handle.readline(), 3, "wavelength/RI record")
            wvl, ri_real[ir], ri_imag[ik] = metadata
            for radius_index in range(nr):
               k11[ir, ik, radius_index, :] = stream.take(theta.size, "scattering row")
   return DubovikScatteringMatrix(ri_real, ri_imag, grid.radius, theta, float(wvl), epsilon, k11,
                                  eps_distribution, key_rd)


def _bilinear(values, x, y):
   """IDL INTERPOL's bilinear interpolation for an in-range RI coordinate."""

   nx, ny = values.shape[:2]
   ix = int(np.clip(np.floor(x), 0, nx - 2))
   iy = int(np.clip(np.floor(y), 0, ny - 2))
   fx = x - ix
   fy = y - iy
   return ((1.0 - fx) * (1.0 - fy) * values[ix, iy]
           + fx * (1.0 - fy) * values[ix + 1, iy]
           + (1.0 - fx) * fy * values[ix, iy + 1]
           + fx * fy * values[ix + 1, iy + 1])


def _quadratic_interpolate(x, y, xout):
   """Local three-point quadratic interpolation used by IDL /LSQUADRATIC."""

   order = np.argsort(x, kind="stable")
   x = np.asarray(x, dtype=np.float64)[order]
   y = np.asarray(y, dtype=np.float64)[order]
   xout = np.asarray(xout, dtype=np.float64)
   result = np.empty_like(xout)
   for index, value in np.ndenumerate(xout):
      interval = int(np.clip(np.searchsorted(x, value, side="right") - 1, 0, x.size - 2))
      start = int(np.clip(interval - (0 if interval == 0 else 1), 0, x.size - 3))
      nodes = x[start:start + 3]
      coeff = y[start:start + 3]
      result[index] = sum(coeff[j] * np.prod([(value - nodes[k]) / (nodes[j] - nodes[k])
                                                for k in range(3) if k != j]) for j in range(3))
   return result


def dubovik_lognormal(base_path, number_conc, rm, s, wn, ri, eps, dqv=None,
                      renorm_ph=False, silent=False):
   """Mirror ``dubovik_lognormal.pro`` for one aspect ratio."""

   base_path = Path(base_path)
   matches = np.flatnonzero((VALID_EPS - float(eps)) ** 2 < 1.0e-5)
   if matches.size != 1:
      raise ValueError(f"EPS value is not valid: {eps}")
   index = int(matches[0])
   names = read_dubovik_filenames(base_path / "ROUTINE" / "name.dat")
   kernel_dir = base_path / "KERNEL_n22_181"
   extinction = read_dubovik_kernel_kext(kernel_dir / names.kext[index])
   scattering = read_dubovik_kernel_scatt_matrix(kernel_dir / names.k11[index])

   real_ri = float(np.real(ri))
   imag_ri = float(np.imag(ri))
   if imag_ri >= 0.0:
      raise ValueError("Imaginary part of RI must be negative")
   n_grid = (real_ri - extinction.ri_real[0]) / (extinction.ri_real[1] - extinction.ri_real[0])
   k_log = np.log(-imag_ri)
   k_grid = (k_log - np.log(-extinction.ri_imag[0])) / (
      np.log(-extinction.ri_imag[1]) - np.log(-extinction.ri_imag[0]))
   kext_interp = np.asarray([_bilinear(extinction.kext[:, :, i], n_grid, k_grid)
                             for i in range(extinction.radius.size)])
   kabs_interp = np.asarray([_bilinear(extinction.kabs[:, :, i], n_grid, k_grid)
                             for i in range(extinction.radius.size)])
   k11_interp = np.asarray([[_bilinear(scattering.k11[:, :, i, j], n_grid, k_grid)
                             for j in range(scattering.theta.size)]
                            for i in range(scattering.radius.size)])

   rd = extinction.radius
   rm_d = float(rm) * extinction.wvl * float(wn)
   distribution = (float(number_conc) * np.sqrt(8.0 * np.pi) / 3.0 * rd ** 3 / np.log(float(s))
                   * np.exp(-0.5 * (np.log(rd / rm_d) / np.log(float(s))) ** 2))
   scale = (float(wn) * extinction.wvl) ** 2
   text = np.sum(distribution * kext_interp) / scale
   tabs = np.sum(distribution * kabs_interp) / scale
   tsca = text - tabs
   denominator = np.sum(distribution * (kext_interp - kabs_interp))
   phase_table = np.sum(distribution[:, None] * k11_interp, axis=0) / denominator
   theta = np.deg2rad(scattering.theta)
   sine = np.sin(theta)
   cosine = np.cos(theta)
   g = np.sum((sine[:-1] * cosine[:-1] * phase_table[:-1]
               + sine[1:] * cosine[1:] * phase_table[1:]) * np.diff(theta)) * 0.25
   phase = phase_table
   if dqv is not None:
      dqv = np.asarray(dqv, dtype=np.float64)
      phase = np.exp(_quadratic_interpolate(cosine, np.log(phase_table), dqv))
      if dqv.size > theta.size:
         sqv = np.sqrt(1.0 - dqv * dqv)
         dt = np.diff(np.arccos(dqv))
         norm = 1.0 if not renorm_ph else np.sum((sqv[:-1] * phase[:-1] + sqv[1:] * phase[1:]) * dt) / 4.0
         phase = phase / norm
         g = np.sum((sqv[:-1] * dqv[:-1] * phase[:-1] + sqv[1:] * dqv[1:] * phase[1:]) * dt) * 0.25
         g = g / norm
   return text * 1.0e3, tsca * 1.0e3, tsca / text, float(g), phase


def dubovik_lognormal_multiple_eps(base_path, number_conc, rm, s, wn, ri, eps, neps,
                                   dqv=None, renorm_ph=False, silent=False):
   """Mirror ``dubovik_lognormal_multiple_eps.pro``."""

   eps = np.asarray(eps, dtype=np.float64)
   neps = np.asarray(neps, dtype=np.float64)
   if eps.ndim != 1 or neps.shape != eps.shape:
      raise ValueError("EPS and NEPS must have the same number of elements")
   if np.any(neps < 0.0) or np.sum(neps) <= 0.0:
      raise ValueError("NEPS values must be non-negative with a positive total")
   if np.ndim(rm) not in (0, 1) or np.size(rm) not in (1, eps.size):
      raise ValueError("Size of RM and EPS are not consistent")
   if np.ndim(s) not in (0, 1) or np.size(s) not in (1, eps.size):
      raise ValueError("Size of S and EPS are not consistent")
   if np.ndim(ri) not in (0, 1) or np.size(ri) not in (1, eps.size):
      raise ValueError("Size of RI and EPS are not consistent")
   weights = neps / np.sum(neps)
   rms = np.full(eps.size, rm, dtype=np.float64) if np.size(rm) == 1 else np.asarray(rm, dtype=np.float64)
   spreads = np.full(eps.size, s, dtype=np.float64) if np.size(s) == 1 else np.asarray(s, dtype=np.float64)
   ris = np.full(eps.size, ri, dtype=np.complex128) if np.size(ri) == 1 else np.asarray(ri, dtype=np.complex128)
   values = [dubovik_lognormal(base_path, number_conc, rms[i], spreads[i], wn, ris[i], eps[i],
                               dqv=dqv, renorm_ph=renorm_ph, silent=silent)
             for i in range(eps.size) if weights[i] > 0.0]
   positive = weights > 0.0
   bexts = np.asarray([value[0] for value in values])
   bscas = np.asarray([value[1] for value in values])
   gs = np.asarray([value[3] for value in values])
   selected_weights = weights[positive]
   bext = np.sum(bexts * selected_weights)
   bsca = np.sum(bscas * selected_weights)
   g = np.sum(bscas * selected_weights * gs) / bsca
   if dqv is None:
      phase = None
   else:
      phases = np.asarray([value[4] for value in values])
      phase = np.sum(bscas[:, None] * selected_weights[:, None] * phases, axis=0) / bsca
   return bext, bsca, bsca / bext, float(g), phase
