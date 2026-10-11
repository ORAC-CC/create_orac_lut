"""Stage 7: identify the IDL arithmetic of read_baum_lambda (INTERPOL indices + INTERPOLATE /GRID).

For a single Baum component the V21 LUT variables extinction_coefficient,
single_scatter_albedo and asymmetry_parameter are the FLOAT results of

   i_wl = interpol(indgen(n_wl), wl_baum, wl)       (FLOAT arithmetic)
   i_Re = interpol(indgen(n_Re), Re_baum, Re)       (Re_baum = De/2., extended to Re[last])
   x    = interpolate(temp2, i_wl, i_Re, /grid)     (bilinear, edge-clamped)

(read_baum.pro), passed through the double-precision mixing, which returns
them unchanged for one component.  Candidate formulas for the IDL INTERPOLATE
arithmetic are evaluated and compared bitwise with the production V21 pagg LUT
for all 36 Terra MODIS channels and 24 effective radii.

Usage: python baum_interpolate.py   -> results/baum_interpolate.json
"""

import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.idl_mirror import load_inststr, load_lutstr, load_srfstrarr   # noqa: E402

f32 = np.float32
INPUTS = ROOT / "create_orac_lut" / "input_files"
BAUM = INPUTS / "baum" / "AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc"
V21 = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/terra_modis_m_water-ice_a01_pagg_v21.nc")


def interpol_index(x, xout, dtype):
   """IDL interpol(indgen(n), x, xout): (xout-x[s])*1/(x[s+1]-x[s]) + s in the given type."""
   x = np.asarray(x, dtype)
   xout = np.asarray(xout, dtype)
   s = np.clip(np.searchsorted(x, xout, side="right") - 1, 0, x.size - 2)
   return ((xout - x[s]) * dtype(1)) / (x[s + 1] - x[s]) + dtype(1) * s.astype(dtype)


def bilinear(table, ix, iy, form, dtype):
   """INTERPOLATE(table, ix, iy, /GRID) candidates; table[ix_axis, iy_axis] (IDL first index = wavelength)."""
   nx, ny = table.shape
   t = table.astype(dtype)
   ix = np.clip(np.asarray(ix, dtype), 0, nx - 1)
   iy = np.clip(np.asarray(iy, dtype), 0, ny - 1)
   x0 = np.minimum(np.floor(ix).astype(int), nx - 2)
   y0 = np.minimum(np.floor(iy).astype(int), ny - 2)
   fx = (ix - x0.astype(dtype))[:, None]
   fy = (iy - y0.astype(dtype))[None, :]
   a00 = t[x0][:, y0]
   a10 = t[x0 + 1][:, y0]
   a01 = t[x0][:, y0 + 1]
   a11 = t[x0 + 1][:, y0 + 1]
   one = dtype(1)
   if form == "lerp_x_then_y":
      b0 = a00 + fx * (a10 - a00)
      b1 = a01 + fx * (a11 - a01)
      return b0 + fy * (b1 - b0)
   if form == "lerp_y_then_x":
      b0 = a00 + fy * (a01 - a00)
      b1 = a10 + fy * (a11 - a10)
      return b0 + fx * (b1 - b0)
   if form == "weights":
      return (one - fx) * (one - fy) * a00 + fx * (one - fy) * a10 + (one - fx) * fy * a01 + fx * fy * a11
   if form == "weights_x_then_y":
      return (one - fy) * ((one - fx) * a00 + fx * a10) + fy * ((one - fx) * a01 + fx * a11)
   if form == "weights_y_then_x":
      return (one - fx) * ((one - fy) * a00 + fy * a01) + fx * ((one - fy) * a10 + fy * a11)
   raise ValueError(form)


def extended(field, re_baum, re_last):
   """read_baum.pro radius extrapolation: Re_baum += Re[last]; temp2[*, n] = (1.-f)*temp2[*, n-2] + f*temp2[*, n-1]."""
   n = re_baum.size
   if not re_baum[-1] < re_last:
      return field, re_baum
   f = (f32(re_last) - re_baum[n - 2]) / (re_baum[n - 1] - re_baum[n - 2])
   temp2 = np.zeros((field.shape[0], n + 1), f32)
   temp2[:, :n] = field
   temp2[:, n] = (f32(1.0) - f) * temp2[:, n - 2] + f * temp2[:, n - 1]
   return temp2, np.append(re_baum, f32(re_last)).astype(f32)


def main():
   inststr = load_inststr(INPUTS / "inst" / "terra_modis_v1.inst", requestedchannelid=list(range(1, 37)))
   lutstr = load_lutstr(INPUTS / "lut" / "ice-cloud.lut", inststr.max_sat_zenith)
   srfstrarr, _ = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
   wl = np.asarray([s.wvl[0] for s in srfstrarr], f32)
   re = np.asarray(lutstr.efr, f32)
   with Dataset(BAUM) as d:
      wl_baum = np.asarray(d["wavelengths"][:], f32)
      re_baum = np.asarray(d["effective_diameter"][:], f32) / f32(2.0)
      fields = {k: np.asarray(d[k][:], f32).T for k in ("extinction_efficiency", "single_scattering_albedo", "asymmetry_parameter")}
   with Dataset(V21) as d:
      d.set_auto_maskandscale(False)
      reference = {"extinction_efficiency": np.asarray(d["extinction_coefficient"][:], f32).T,
                   "single_scattering_albedo": np.asarray(d["single_scatter_albedo"][:], f32).T,
                   "asymmetry_parameter": np.asarray(d["asymmetry_parameter"][:], f32).T}
   result = {}
   for index_type in (np.float32, np.float64):
      for form in ("lerp_x_then_y", "lerp_y_then_x", "weights", "weights_x_then_y", "weights_y_then_x"):
         for arithmetic in (np.float32, np.float64):
            counts = {}
            for name, field in fields.items():
               temp2, re_ext = extended(field, re_baum, re[-1])
               i_wl = interpol_index(wl_baum, wl, index_type)
               i_re = interpol_index(re_ext, re, index_type)
               value = bilinear(temp2, i_wl, i_re, form, arithmetic).astype(f32)
               counts[name] = int(np.count_nonzero(value.view(np.uint32) != reference[name].view(np.uint32)))
            key = f"index {np.dtype(index_type).name}, {form}, arithmetic {np.dtype(arithmetic).name}"
            result[key] = counts
            print(f"{key:60s} differing of {reference['asymmetry_parameter'].size}: {counts}", flush=True)
   (HERE / "results" / "baum_interpolate.json").write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
   main()
