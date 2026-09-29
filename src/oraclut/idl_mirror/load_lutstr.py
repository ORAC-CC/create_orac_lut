"""Python counterpart of load_lutstr.pro: read a LUT grid definition file.

IDL:  Function lut_quadrature, N, spacing, values, Invalid_N_FLAG, max_val=max_val
      pro load_lutstr, file, max_sat_zenith, LUTstr, include_pressure = include_pressure

The file contains pairs of lines: "<N> <spacing>" then the values (two end
points for 'linear'/'logarithmic', all N values for 'uneven_linear'/
'uneven_logarithmic'), for optical depth, effective radius, solar zenith,
satellite zenith and relative azimuth, and, for the aerosol formulation, a
sixth pair for surface pressure.  The satellite-zenith end point is replaced
by the instrument's maximum satellite zenith angle, as in the IDL.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np


def lut_quadrature(n, spacing, values, max_val=None):
   """Expand one grid exactly as the IDL function lut_quadrature does.

   Returns ``(x, invalid_n_flag)``.  Arithmetic is single precision, as the
   IDL's is (float() input, findgen).
   """

   n = int(n)
   values = np.asarray(values, dtype=np.float32)
   findgen = np.arange(n, dtype=np.float32)
   if spacing == "linear":
      if max_val is not None:
         x = values[0] + (np.float32(max_val) - values[0]) * findgen / (n - 1)
      else:
         x = values[0] + (values[1] - values[0]) * findgen / (n - 1)
   elif spacing == "logarithmic":
      if max_val is not None:
         x = np.float32(10) ** (np.log10(values[0]) + np.log10(np.float32(max_val) / values[0]) * findgen / (n - 1))
      else:
         x = np.float32(10) ** (np.log10(values[0]) + np.log10(values[1] / values[0]) * findgen / (n - 1))
   elif spacing in ("uneven_linear", "uneven_logarithmic"):
      x = values
   else:
      raise ValueError("Invalid spacing descriptor in LUT definition. Must be one of: "
                       "'linear', 'uneven_linear', 'logarithmic' or 'uneven_logarithmic'")
   invalid_n_flag = (x.size != n)
   return np.asarray(x, dtype=np.float32), invalid_n_flag


def _grid_lines(file):
   """Return the non-comment lines of the LUT definition (IDL: WHERE Lines NE '' ...).

   Deliberate difference: a trailing "# comment" is removed from every line;
   the IDL only tolerates it on the header lines, where strsplit stops at the
   spacing word.
   """

   lines = []
   for raw in Path(file).read_text().splitlines():
      line = raw.split("#", 1)[0].strip()
      if line != "" and line[0] != "*":
         lines.append(line)
   return lines


def load_lutstr(file, max_sat_zenith, include_pressure=False):
   """Read the LUT grid definition ``file`` and return the structure ``lutstr``."""

   lines = _grid_lines(file)
   if len(lines) < (12 if include_pressure else 10):
      raise ValueError(f"LUT definition {file} is incomplete")

   # First two data lines should be optical depth values
   words = lines[0].split()
   opd_n = int(words[0])
   opd_spacing = words[1].lower()
   opd, invalid_n_flag = lut_quadrature(opd_n, opd_spacing, [float(x) for x in lines[1].split()])
   if invalid_n_flag:
      raise ValueError("Optical Depth N does not match the number of values given")

   # Next two data lines should be effective radius values
   words = lines[2].split()
   efr_n = int(words[0])
   efr_spacing = words[1].lower()
   efr, invalid_n_flag = lut_quadrature(efr_n, efr_spacing, [float(x) for x in lines[3].split()])
   if invalid_n_flag:
      raise ValueError("Effective radius N does not match the number of values given")

   # Next two data lines should be solar zenith angle values
   words = lines[4].split()
   soz_n = int(words[0])
   soz_spacing = words[1].lower()
   soz, invalid_n_flag = lut_quadrature(soz_n, soz_spacing, [float(x) for x in lines[5].split()])
   if invalid_n_flag:
      raise ValueError("Solar zenith angle N does not match the number of values given")

   # Next two data lines should be satellite zenith angle values; the upper end
   # point is the instrument's maximum satellite zenith angle (IDL max_val).
   words = lines[6].split()
   saz_n = int(words[0])
   saz_spacing = words[1].lower()
   saz, invalid_n_flag = lut_quadrature(saz_n, saz_spacing, [float(x) for x in lines[7].split()], max_val=max_sat_zenith)
   if invalid_n_flag:
      raise ValueError("Satellite zenith angle N does not match the number of values given")

   # Next two data lines should be relative azimuth values
   words = lines[8].split()
   raa_n = int(words[0])
   raa_spacing = words[1].lower()
   if raa_spacing != "linear":
      raise ValueError("Relative azimuth spacing must be linear")
   raa, invalid_n_flag = lut_quadrature(raa_n, raa_spacing, [float(x) for x in lines[9].split()])
   if invalid_n_flag:
      raise ValueError("Relative azimuth N does not match the number of values given")

   if include_pressure:
      # if aerosol lut file then pressure also needed
      words = lines[10].split()
      prs_n = int(words[0])
      prs_spacing = words[1].lower()
      prs, invalid_n_flag = lut_quadrature(prs_n, prs_spacing, [float(x) for x in lines[11].split()])
      # The preserved IDL does not test Invalid_N_FLAG for pressure (V2.1 does);
      # the count is checked here like every other grid.
      if invalid_n_flag:
         raise ValueError("Surface pressure N does not match the number of values given")
      lutstr = SimpleNamespace(
         opd_n=opd_n, efr_n=efr_n, soz_n=soz_n, saz_n=saz_n, raa_n=raa_n, prs_n=prs_n,
         opd_spacing=opd_spacing, efr_spacing=efr_spacing, soz_spacing=soz_spacing,
         saz_spacing=saz_spacing, raa_spacing=raa_spacing, prs_spacing=prs_spacing,
         opd=opd, efr=efr, soz=soz, saz=saz, raa=raa, prs=prs,
      )
   else:
      # Build the output structure
      lutstr = SimpleNamespace(
         opd_n=opd_n, efr_n=efr_n, soz_n=soz_n, saz_n=saz_n, raa_n=raa_n,
         opd_spacing=opd_spacing, efr_spacing=efr_spacing, soz_spacing=soz_spacing,
         saz_spacing=saz_spacing, raa_spacing=raa_spacing,
         opd=opd, efr=efr, soz=soz, saz=saz, raa=raa,
      )
   return lutstr
