"""Python counterpart of load_gasstr.pro: read MODTRAN gas optical-depth profiles.

IDL:  pro read_gasstr, Filename, Onegasstr
      pro load_gasstr, atmospheres, platform, instrument, channelid, gasdir, gasstr

One file per channel: ModtranGasOpd_A<atmosphere>_<platform>_<instrument>_ch<nn>.gas
holding, after two comment lines, four labelled header items (*atmosphere,
*instrument, *channelid, *nlevels in any order) and then nlevels rows of
height (km) and cumulative gas optical depth at that level.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np


def read_gasstr(filename):
   """Read one gas optical-depth file into the structure ``onegasstr``."""

   lines = Path(filename).read_text().splitlines()
   i = 2                                              # skip the two comment lines
   atmosphere = instrument = channelid = None
   nlevels = None
   for _ in range(4):
      label = lines[i][1:].strip().lower()
      value = lines[i + 1].strip()
      if label == "atmosphere":
         atmosphere = value
      elif label == "instrument":
         instrument = value
      elif label == "channelid":
         channelid = value
      elif label == "nlevels":
         nlevels = int(value)
      else:
         raise ValueError(f"Gas file {filename}: unexpected header item {lines[i]!r}")
      i += 2
   if nlevels is None:
      raise ValueError(f"Gas file {filename} does not declare nlevels")
   values = np.asarray([[float(x) for x in lines[i + j].split()[:2]] for j in range(nlevels)], dtype=np.float32)
   onegasstr = SimpleNamespace(atmosphere=atmosphere, instrument=instrument, channelid=channelid,
                               nlevels=nlevels, height=values[:, 0], tau_gas=values[:, 1])
   return onegasstr


def load_gasstr(atmospheres, platform, instrument, channelid, gasdir):
   """Return the list ``gasstr`` with one structure per channel in ``channelid``."""

   gasdir = Path(gasdir)
   gasstr = []
   for channel in channelid:
      gasfile = "ModtranGasOpd_A" + str(atmospheres) + "_" + platform + "_" + instrument + "_ch" + f"{int(channel):02d}" + ".gas"
      if not (gasdir / gasfile).is_file():
         raise FileNotFoundError(f"Gas profile for channel {int(channel)} does not exist: {gasdir / gasfile}")
      gasstr.append(read_gasstr(gasdir / gasfile))
   return gasstr
