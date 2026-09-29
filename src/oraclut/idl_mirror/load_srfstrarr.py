"""Python counterparts of load_srfstrarr.pro and the routines it calls.

IDL:  pro load_srfstrarr, inststr, solar_spectrum_filename, srfstrarr, QM, nwvl_max
      pro read_srfstr, filename, srfstr                 (input_files/read_srfstr.pro)
      FUNCTION integrate_trapeziodal, x, y              (integrate_trapeziodal.pro)
      function simpson_integral, x, y                   (input_files/simpson_integral.pro)
      pro load_solar_spectrum, solar_spectrum_filename, sssi
      pro bbconstants, filename, b1, b2, t1, t2         (input_files/bbconstants.pro)

For each requested channel the spectral response function (SRF) is read, the
effective centre wavelength/wavenumber, the in-band solar irradiance F0 and the
blackbody fit constants B1, B2, T1, T2 are determined, and the spectral
quadrature for the LUT calculation is set according to QM (srf_quad):
   QM = 0  every SRF point (not ported: unvalidated)
   QM = 1  one point at the centre wavelength with unit weight (the validated path)
   QM = 2  segmented SRF (segment.pro, fe=0.001, minn=12)
"""

from pathlib import Path
from types import SimpleNamespace
import re

import numpy as np

from .segment import segment


def integrate_trapeziodal(x, y):
   """IDL: dx = x(1:*)-x ; ym = 0.5*(y+y(1:*)) ; RETURN, TOTAL(dx*ym)  (single precision)."""

   x = np.asarray(x, dtype=np.float32)
   y = np.asarray(y, dtype=np.float32)
   dx = x[1:] - x[:-1]
   ym = np.float32(0.5) * (y[:-1] + y[1:])
   return np.sum(dx * ym, dtype=np.float32)


def simpson_integral(x, y):
   """IDL simpson_integral.pro as it stands: RETURN, total(x*y).

   The routine computes deltax and Yave but returns total(x*y); this is
   reproduced deliberately because F0 in every validated product depends on it.
   """

   return np.sum(np.asarray(x, dtype=np.float32) * np.asarray(y, dtype=np.float32), dtype=np.float32)


def read_srfstr(filename):
   """Read an RTTOV-format SRF file (4 header lines, then wavenumber, response)."""

   filename = Path(filename)
   lines = filename.read_text().splitlines()
   nrows = int(lines[2].strip())                     # third header line: number of data lines
   data = []
   for line in lines[4:4 + nrows]:
      parts = line.split()
      data.append((float(parts[0]), float(parts[1])))
   data = np.asarray(data, dtype=np.float32)          # IDL: data = fltarr(2, nrows)
   if data.shape[0] != nrows:
      raise ValueError(f"SRF file {filename} declares {nrows} points but holds {data.shape[0]}")
   wvn = data[:, 0]
   srf = data[:, 1]
   wvl = np.float32(1e4) / wvn

   wvl_centre = integrate_trapeziodal(wvl, srf * wvl) / integrate_trapeziodal(wvl, srf)
   wvn_centre = integrate_trapeziodal(wvn, srf * wvn) / integrate_trapeziodal(wvn, srf)

   srfstr = SimpleNamespace(filename=filename.name, nwvl=wvl.size, wvl=wvl, wvn=wvn, srf=srf,
                            wvl_centre=np.float32(wvl_centre), wvn_centre=np.float32(wvn_centre))
   return srfstr


def load_solar_spectrum(solar_spectrum_filename):
   """IDL: loadxy, solar_spectrum_filename, 0, 1, wvl, val  (columns 0 and 1, '*'/'#' comments).

   Deliberate difference: loadxyz.pro reads DOUBLE; the validated Python read
   the table as single precision and that is retained here.
   """

   rows = []
   for raw in Path(solar_spectrum_filename).read_text().splitlines():
      line = raw.strip()
      if line == "" or line[0] in "*#":
         continue
      parts = line.split()
      rows.append((float(parts[0]), float(parts[1])))
   values = np.asarray(rows, dtype=np.float32)
   sssi = SimpleNamespace(name=Path(solar_spectrum_filename).name, wvl=values[:, 0], val=values[:, 1])
   return sssi


_BBCONSTANTS_PATTERN = re.compile(
   r"'([^']+)':\s*begin\s*\n\s*b1=\s*([^&]+)&\s*b2=\s*([^&]+)&\s*t1=\s*([^&]+)&\s*t2=\s*([^\n]+)",
   re.IGNORECASE,
)


def bbconstants(filename, bbconstants_file):
   """Return (b1, b2, t1, t2) for SRF file ``filename`` from bbconstants.pro.

   bbconstants.pro is an IDL CASE statement used as a data table; it is parsed
   here rather than executed.  An SRF file absent from the table is an error,
   as the IDL CASE without ELSE would be.
   """

   text = Path(bbconstants_file).read_text()
   for match in _BBCONSTANTS_PATTERN.finditer(text):
      if match.group(1) == filename:
         return tuple(np.float32(float(match.group(i))) for i in range(2, 6))
   raise ValueError(f"bbconstants: no blackbody constants for SRF file {filename!r} in {bbconstants_file}")


def _qm2_segmentation_input(srfstr, inststr):
   """Return the paired SRF arrays passed to the QM=2 segment routine.

   The EarthCARE MSI channel-1 RTTOV file contains a small number of local
   wavelength reversals after the IDL reader converts its wavenumbers to
   wavelength.  IDL keeps those raw arrays in their original order and its
   tabulated integration orders each sub-array internally.  The Python
   ``segment`` mirror exposes that behavior as ``sort=True``; the raw arrays
   remain untouched for the centre-wavelength, solar-constant, and blackbody
   calculations above.

   Other SRFs retain their original ascending/descending handling, and a
   duplicate wavelength remains an error rather than being silently discarded.
   """

   wvl = np.asarray(srfstr.wvl, dtype=np.float32)
   srf = np.asarray(srfstr.srf, dtype=np.float32)
   differences = np.diff(wvl)
   if np.all(differences > 0.0) or np.all(differences < 0.0):
      return wvl, srf

   is_earthcare_msi = (
      str(inststr.platform).lower() == "earthcare"
      and str(inststr.instrument).lower() == "msi"
   )
   if not is_earthcare_msi:
      # Let segment/int_tabulated raise its normal strict-ordering error for
      # malformed input from every other instrument.
      return wvl, srf

   order = np.argsort(wvl, kind="stable")
   sorted_wvl = wvl[order]
   if np.any(np.diff(sorted_wvl) <= 0.0):
      raise ValueError("QM=2 EarthCARE SRF contains duplicate wavelength values")
   return wvl, srf


def _qm2_uses_idl_tabulated_sort(srfstr, inststr):
   """Whether the narrow EarthCARE QM=2 INT_TABULATED compatibility applies."""

   differences = np.diff(np.asarray(srfstr.wvl, dtype=np.float32))
   return (
      str(inststr.platform).lower() == "earthcare"
      and str(inststr.instrument).lower() == "msi"
      and not (np.all(differences > 0.0) or np.all(differences < 0.0))
   )


def load_srfstrarr(inststr, solar_spectrum_filename, qm, in_path):
   """Build ``srfstrarr`` (one structure per channel) and return (srfstrarr, nwvl_max).

   ``in_path`` locates ``srf/`` and ``bbconstants.pro``; the IDL hard-codes
   'input_files/srf/' relative to the working directory.
   """

   in_path = Path(in_path)
   if qm not in (0, 1, 2):
      raise ValueError("srf_quad (QM) must be 0, 1 or 2")
   if qm == 0:
      raise NotImplementedError("srf_quad=0 (every raw SRF point) is not ported and validated")
   # Work out the maximum number of spectral points that will be used
   nwvl_max = 1
   for i in range(inststr.number_of_nadir_channels):
      srfstr = read_srfstr(in_path / "srf" / inststr.srf_file[i])
      if qm == 0:
         nwvl = srfstr.nwvl
      elif qm == 1:
         nwvl = 1
      elif qm == 2:
         segment_wvl, segment_srf = _qm2_segmentation_input(srfstr, inststr)
         segmented_wvl, _ = segment(
            segment_wvl, segment_srf, np.float32(0.001), minn=12,
            sort=_qm2_uses_idl_tabulated_sort(srfstr, inststr),
         )
         nwvl = segmented_wvl.size
      nwvl_max = max(nwvl_max, nwvl)

   sssi = load_solar_spectrum(solar_spectrum_filename)

   srfstrarr = [
      SimpleNamespace(nwvl=0, wvl_centre=np.float32(0), wvn_centre=np.float32(0), f0=np.float32(0),
                      b1=np.float32(0), b2=np.float32(0), t1=np.float32(0), t2=np.float32(0),
                      wvl=np.zeros(nwvl_max, dtype=np.float32), val=np.zeros(nwvl_max, dtype=np.float32))
      for _ in range(inststr.number_of_nadir_channels)
   ]

   for i in range(inststr.number_of_nadir_channels):
      srfstr = read_srfstr(in_path / "srf" / inststr.srf_file[i])

      srfstrarr[i].wvl_centre = srfstr.wvl_centre
      srfstrarr[i].wvn_centre = srfstr.wvn_centre

      # Calculate the solar constant for the channel (different units for SW and
      # thermal channels).  IDL uses INTERPOL; np.interp agrees within the
      # tabulated solar spectrum, which covers every SRF used.
      if srfstr.wvl_centre > 3:
         # IDL: wvn = 1d4/sssi.wvl
         #      F0 = (simpson_integral(srfstr.wvn, interpol(1d8*sssi.val/wvn^2, wvn, srfstr.wvn)*srfstr.srf)
         #            / simpson_integral(srfstr.wvn, srfstr.srf)) / !pi
         wvn = np.float64(1e4) / sssi.wvl.astype(np.float64)
         irradiance = np.float64(1e8) * sssi.val.astype(np.float64) / wvn ** 2
         order = np.argsort(wvn)                       # np.interp needs ascending abscissae
         solar_on_srf = np.interp(srfstr.wvn.astype(np.float64), wvn[order], irradiance[order]).astype(np.float32)
         srfstrarr[i].f0 = np.float32(
            simpson_integral(srfstr.wvn, solar_on_srf * srfstr.srf)
            / simpson_integral(srfstr.wvn, srfstr.srf) / np.float32(np.pi))
      else:
         # IDL: F0 = (simpson_integral(srfstr.wvl, interpol(sssi.val*10, sssi.wvl, srfstr.wvl)*srfstr.srf)
         #            / simpson_integral(srfstr.wvl, srfstr.srf)) / !pi
         # (the validated Python multiplies by 10 after interpolating; retained)
         solar_on_srf = np.interp(srfstr.wvl.astype(np.float64), sssi.wvl.astype(np.float64),
                                  sssi.val.astype(np.float64)).astype(np.float32)
         srfstrarr[i].f0 = np.float32(
            simpson_integral(srfstr.wvl, solar_on_srf * np.float32(10.0) * srfstr.srf)
            / simpson_integral(srfstr.wvl, srfstr.srf) / np.float32(np.pi))

      b1, b2, t1, t2 = bbconstants(srfstr.filename, in_path / "bbconstants.pro")
      srfstrarr[i].b1 = b1
      srfstrarr[i].b2 = b2
      srfstrarr[i].t1 = t1
      srfstrarr[i].t2 = t2

      if qm == 0:
         srfstrarr[i].nwvl = srfstr.nwvl
         srfstrarr[i].wvl[0:srfstr.nwvl] = srfstr.wvl
         srfstrarr[i].val[0:srfstr.nwvl] = srfstr.srf
      elif qm == 1:
         srfstrarr[i].nwvl = 1
         srfstrarr[i].wvl[0] = srfstrarr[i].wvl_centre
         srfstrarr[i].val[0] = 1
      elif qm == 2:
         # IDL: segment(srfstr.wvl, srfstr.srf, x, y, .001, minn=12),
         # then copy x/y into the SRF structure without reordering them.
         segment_wvl, segment_srf = _qm2_segmentation_input(srfstr, inststr)
         segmented_wvl, segmented_srf = segment(
            segment_wvl, segment_srf, np.float32(0.001), minn=12,
            sort=_qm2_uses_idl_tabulated_sort(srfstr, inststr),
         )
         srfstrarr[i].nwvl = segmented_wvl.size
         srfstrarr[i].wvl[0:segmented_wvl.size] = segmented_wvl
         srfstrarr[i].val[0:segmented_srf.size] = segmented_srf
   return srfstrarr, nwvl_max
