"""Python counterpart of load_inststr.pro: read an instrument definition file.

IDL:  pro load_inststr, file, inststr, RequestedChannelID = RequestedChannelID

The IDL routine returns the structure ``inststr``; here the same tags are
returned as attributes of a SimpleNamespace, so ``inststr.channelid`` and
``inststr.solar_channel_flag`` read exactly as in the IDL.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np


def load_inststr(file, requestedchannelid=None):
   """Read the instrument definition file ``file``.

   ``requestedchannelid`` is the IDL RequestedChannelID keyword: when given, only
   those channels are kept in the returned structure.
   """

   file = Path(file)
   instrument_filename = file.name

   # Read in the instrument description; remove leading and trailing spaces;
   # remove blank lines and commented lines (those starting with "*" or "#").
   lines = [line.strip() for line in file.read_text().splitlines()]
   lines = [line for line in lines if line != "" and line[0] not in "*#"]

   # Deconstruct each line into LHS = RHS
   lhs = []
   rhs = []
   for line in lines:
      result = line.split("=")
      lhs.append(result[0])
      rhs.append(result[1])

   # Lines without square brackets are single-valued instrument constants;
   # lines with [n] are channel-specific values for channel n.
   singles = [i for i in range(len(lhs)) if "[" not in lhs[i]]
   multiples = [i for i in range(len(lhs)) if "[" in lhs[i]]

   platform = instrument = instrument_version = None
   view = max_sat_zenith = None
   available_channelid = solar_channelid = thermal_channelid = None
   for i in singles:
      # IDL: strlowcase(strcompress(LHS,/remove_all)) removes all white space
      name = lhs[i].lower().replace(" ", "").replace("\t", "")
      if name == "platform":
         platform = rhs[i].strip().lower()
      elif name == "instrument":
         instrument = rhs[i].strip().lower()
      elif name == "instrumentversion":
         instrument_version = rhs[i].strip().lower()
      elif name == "view":
         view = int(rhs[i])
      elif name == "availablechannels":
         available_channelid = np.asarray([int(x) for x in rhs[i].split()], dtype=np.int16)
      elif name == "solarchannels":
         solar_channelid = np.asarray([int(x) for x in rhs[i].split()], dtype=np.int16)
      elif name == "thermalchannels":
         thermal_channelid = np.asarray([int(x) for x in rhs[i].split()], dtype=np.int16)
      elif name == "maximumsatellitezenith":
         max_sat_zenith = np.float32(rhs[i])
      else:
         raise ValueError("Unknown variable in instrument definition file: " + str(file) + ": " + lhs[i])
   for name, value in (("platform", platform), ("instrument", instrument),
                       ("instrument version", instrument_version), ("view", view),
                       ("available channels", available_channelid),
                       ("maximum satellite zenith", max_sat_zenith)):
      if value is None:
         raise ValueError(f"Instrument definition file {file} does not define '{name}'")
   if solar_channelid is None:
      solar_channelid = np.zeros(0, dtype=np.int16)
   if thermal_channelid is None:
      thermal_channelid = np.zeros(0, dtype=np.int16)

   number_of_nadir_channels = len(available_channelid)
   srf_file = [""] * number_of_nadir_channels
   # TEMPORARY FOR BACK COMPATIBILITY (old SAD file values carried through)
   oldf0 = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldf1 = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldnefr = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldwvn = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldb1 = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldb2 = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldt1 = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldt2 = np.zeros(number_of_nadir_channels, dtype=np.float32)
   oldnebt = np.zeros(number_of_nadir_channels, dtype=np.float32)
   # *******************************
   snr = np.zeros(number_of_nadir_channels, dtype=np.float32)
   rgu = np.zeros(number_of_nadir_channels, dtype=np.float32)
   rou = np.zeros(number_of_nadir_channels, dtype=np.float32)
   rua = np.zeros(number_of_nadir_channels, dtype=np.float32)
   rub = np.zeros(number_of_nadir_channels, dtype=np.float32)
   ruc = np.zeros(number_of_nadir_channels, dtype=np.float32)
   refbt = np.zeros(number_of_nadir_channels, dtype=np.float32)
   nedt = np.zeros(number_of_nadir_channels, dtype=np.float32)
   solar_channel_flag = np.zeros(number_of_nadir_channels, dtype=np.int16)
   thermal_channel_flag = np.zeros(number_of_nadir_channels, dtype=np.int16)

   # Flag the solar and thermal channels among the available channels.
   # Deliberate difference: the IDL tests "Count Eq -1", which never fires (WHERE
   # returns Count = 0 for no match), so an unlisted channel would silently
   # flag the last available channel; here it is an error.
   for channel in solar_channelid:
      j = np.flatnonzero(channel == available_channelid)
      if j.size == 0:
         raise ValueError("Solar channel not included in available channels: " + str(channel))
      solar_channel_flag[j] = 1
   for channel in thermal_channelid:
      j = np.flatnonzero(channel == available_channelid)
      if j.size == 0:
         raise ValueError("Thermal channel not included in available channels: " + str(channel))
      thermal_channel_flag[j] = 1

   # Channel specific values: "name [n] = value"
   for i in multiples:
      lsb = lhs[i].index("[")
      rsb = lhs[i].index("]")
      channelid = int(lhs[i][lsb + 1:rsb])
      j = np.flatnonzero(channelid == available_channelid)
      if j.size == 0:
         raise ValueError("channel not included in available channels: " + str(channelid))
      name = lhs[i][:lsb].lower().replace(" ", "").replace("\t", "")
      value = rhs[i].strip()
      if name == "srf":
         srf_file[j[0]] = value
      # TEMPORARY FOR BACK COMPATIBILITY
      elif name == "oldf0":
         oldf0[j] = np.float32(value)
      elif name == "oldf1":
         oldf1[j] = np.float32(value)
      elif name == "oldnefr":
         oldnefr[j] = np.float32(value)
      elif name == "oldwvn":
         oldwvn[j] = np.float32(value)
      elif name == "oldb1":
         oldb1[j] = np.float32(value)
      elif name == "oldb2":
         oldb2[j] = np.float32(value)
      elif name == "oldt1":
         oldt1[j] = np.float32(value)
      elif name == "oldt2":
         oldt2[j] = np.float32(value)
      elif name == "oldnebt":
         oldnebt[j] = np.float32(value)
      # *******************************
      elif name == "snr":
         snr[j] = np.float32(value)
      elif name == "rgu":
         rgu[j] = np.float32(value)
      elif name == "rou":
         rou[j] = np.float32(value)
      elif name == "rua":
         rua[j] = np.float32(value)
      elif name == "rub":
         rub[j] = np.float32(value)
      elif name == "ruc":
         ruc[j] = np.float32(value)
      elif name == "refbt":
         refbt[j] = np.float32(value)
      elif name == "nedt":
         nedt[j] = np.float32(value)
      else:
         print("Warning: Unknown variable in instrument definition file: " + name)

   # If the RequestedChannelID keyword has been set, keep only those channels.
   # Deliberate difference: the IDL only warns about a requested channel that
   # is not in the file and silently drops it; here it is an error.
   if requestedchannelid is not None and len(requestedchannelid) > 0:
      match = np.zeros(number_of_nadir_channels, dtype=bool)
      for channel in requestedchannelid:
         matchch = np.flatnonzero(int(channel) == available_channelid)
         if matchch.size == 0:
            raise ValueError("channel " + str(channel) + " not found in " + str(file))
         if matchch.size > 1:
            raise ValueError("channel " + str(channel) + " found multiple times in " + str(file))
         match[matchch[0]] = True
      matchi = np.flatnonzero(match)
      number_of_nadir_channels = matchi.size
      available_channelid = available_channelid[matchi]
      solar_channel_flag = solar_channel_flag[matchi]
      thermal_channel_flag = thermal_channel_flag[matchi]
      srf_file = [srf_file[i] for i in matchi]
      # TEMPORARY FOR BACK COMPATIBILITY
      oldf0 = oldf0[matchi]
      oldf1 = oldf1[matchi]
      oldnefr = oldnefr[matchi]
      oldwvn = oldwvn[matchi]
      oldb1 = oldb1[matchi]
      oldb2 = oldb2[matchi]
      oldt1 = oldt1[matchi]
      oldt2 = oldt2[matchi]
      oldnebt = oldnebt[matchi]
      # *******************************
      snr = snr[matchi]
      rgu = rgu[matchi]
      rou = rou[matchi]
      rua = rua[matchi]
      rub = rub[matchi]
      ruc = ruc[matchi]
      refbt = refbt[matchi]
      nedt = nedt[matchi]

   mixed_channel_flag = (solar_channel_flag & thermal_channel_flag).astype(np.int16)

   # If dual view then nadir results are replicated for the forward view
   if view > 0:
      number_of_channels = 2 * number_of_nadir_channels
   else:
      number_of_channels = number_of_nadir_channels

   inststr = SimpleNamespace(
      instrument_filename=instrument_filename,
      platform=platform,
      instrument=instrument,
      instrument_version=instrument_version,
      max_sat_zenith=max_sat_zenith,
      view=view,
      number_of_channels=number_of_channels,
      number_of_nadir_channels=number_of_nadir_channels,
      channelid=available_channelid,
      solar_channel_flag=solar_channel_flag,
      mixed_channel_flag=mixed_channel_flag,
      thermal_channel_flag=thermal_channel_flag,
      srf_file=srf_file,
      # TEMPORARY FOR BACK COMPATIBILITY
      oldf0=oldf0, oldf1=oldf1, oldnefr=oldnefr, oldwvn=oldwvn,
      oldb1=oldb1, oldb2=oldb2, oldt1=oldt1, oldt2=oldt2, oldnebt=oldnebt,
      # *******************************
      snr=snr, rgu=rgu, rou=rou, rua=rua, rub=rub, ruc=ruc, refbt=refbt, nedt=nedt,
   )
   return inststr
