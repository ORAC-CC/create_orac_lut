"""Python counterpart of load_atmstr.pro: read an atmospheric profile.

IDL:  pro read_presdat, file, presstr     (three columns: height km, pressure hPa, temperature K)
      pro load_atmstr, file, atmospheres, atmstr

atmospheres = 0 selects the old three-column midsatm.dat reader; any other
code selects the RFM/MODTRAN .atm format (*HGT, *PRE, *TEM blocks).  Levels
above 100 km are dropped and the profile is reversed so that level 0 is the
top of the atmosphere and the last level is the surface.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np


def read_presdat(file):
   """Read the old three-column profile (comments start with '#')."""

   rows = []
   for raw in Path(file).read_text().splitlines():
      line = raw.strip()
      if line == "" or line[0] == "#":
         continue
      rows.append([float(x) for x in line.split()[:3]])
   rows = np.asarray(rows, dtype=np.float32)
   presstr = SimpleNamespace(nlevels=rows.shape[0], height=rows[:, 0], pressure=rows[:, 1], temperature=rows[:, 2])
   return presstr


def _read_block(lines, i, nlevels):
   """Read nlevels comma/space separated numbers starting at lines[i]; return (values, next i)."""

   values = []
   while len(values) < nlevels:
      values.extend(float(x) for x in lines[i].replace(",", " ").split() if x.strip() != "")
      i += 1
   return np.asarray(values[:nlevels], dtype=np.float32), i


def load_atmstr(file, atmospheres):
   """Return the structure ``atmstr`` (nlevels, height, pressure, temperature)."""

   if int(atmospheres) == 0:
      return read_presdat(file)

   lines = Path(file).read_text().splitlines()
   i = 0
   # Comment lines start with "!"
   while lines[i].startswith("!"):
      i += 1
   nlevels = int(lines[i].split()[0])
   i += 1
   # IDL: READF label; READF H; READF label; READF P; READF label; READF T
   if not lines[i].lstrip().upper().startswith("*HGT"):
      raise ValueError(f"Atmosphere file {file}: expected *HGT block, found {lines[i]!r}")
   h, i = _read_block(lines, i + 1, nlevels)
   if not lines[i].lstrip().upper().startswith("*PRE"):
      raise ValueError(f"Atmosphere file {file}: expected *PRE block, found {lines[i]!r}")
   p, i = _read_block(lines, i + 1, nlevels)
   if not lines[i].lstrip().upper().startswith("*TEM"):
      raise ValueError(f"Atmosphere file {file}: expected *TEM block, found {lines[i]!r}")
   t, i = _read_block(lines, i + 1, nlevels)

   q = np.flatnonzero(h <= 100)                        # MODTRAN code limited to 100 km
   nlevels = q.size
   atmstr = SimpleNamespace(nlevels=nlevels, height=h[q][::-1].copy(), pressure=p[q][::-1].copy(),
                            temperature=t[q][::-1].copy())
   return atmstr
