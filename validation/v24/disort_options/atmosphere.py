"""Representative clear-sky terms for the ORAC forward model (validation only).

In ORAC these come from RTTOV, interpolated to the cloud pressure
(src/interpol_solar.F90, src/interpol_thermal.F90).  They are identical for
every DISORT experiment, so the test only needs representative values.  They
are built from the LUT generator's own inputs, read through the production
loaders:

   * temperature profile: create_orac_lut/input_files/atm/mls.atm (load_atmstr);
   * gas optical depth: the MODTRAN 3.5 mid-latitude summer per-channel
     profiles create_orac_lut/input_files/gas/ModtranGasOpd_A2_<platform>_<instrument>_chNN.gas
     (load_gasstr), cumulative from the top of the atmosphere;
   * Planck function: the band-corrected form of the LUT (B1, B2, T1, T2), as in
     ORAC's T2R.

Solar: Tac, Tbc are nadir transmittances above and below the cloud top
(ORAC raises them to sec(sza), sec(vza) itself: src/fm_solar.F90:749-781).
Thermal, along the view direction (RTTOV convention), non-scattering
Schwarzschild sums over the profile layers (layer temperature = mean of its
levels):
   Tac      transmittance cloud top -> TOA
   Rac_up   atmospheric emission above the cloud reaching the TOA
   Rac_dwn  downwelling atmospheric emission at the cloud top
   Rbc_up   upwelling radiance at the cloud level: surface emission
            eps_s B(Ts) attenuated to the cloud plus the layers in between
   B_c      Planck radiance at the cloud-top temperature
Reflection of downwelling radiation by the surface below the cloud is
neglected (eps_s = 0.98).
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.idl_mirror import load_atmstr, load_gasstr   # noqa: E402
from orac_fm import t2r                                    # noqa: E402

INPUTS = ROOT / "create_orac_lut" / "input_files"
SURFACE_EMISSIVITY = 0.98
CLOUD_TOP_KM = {"liquid": 4.0, "ice": 9.0}


def profile(platform, instrument, channels):
   atm = load_atmstr(INPUTS / "atm" / "mls.atm", 2)
   gas = load_gasstr(2, platform, instrument, channels, INPUTS / "gas")
   tau = {int(g.channelid): np.asarray(g.tau_gas, dtype=np.float64) for g in gas}
   return (np.asarray(atm.height, dtype=np.float64), np.asarray(atm.temperature, dtype=np.float64),
           np.asarray(atm.pressure, dtype=np.float64), tau)


def terms(platform, instrument, channels, planck, cloud, vza):
   """{channel: dict} of solar (Tac, Tbc) and thermal terms for a cloud top and view zenith."""

   height, temp, pressure, tau = profile(platform, instrument, channels)
   ic = int(np.argmin(np.abs(height - CLOUD_TOP_KM[cloud])))
   mu = np.cos(np.deg2rad(vza))
   out = {}
   for ch in channels:
      t = tau[ch]
      entry = {"Tac": float(np.exp(-t[ic])), "Tbc": float(np.exp(-(t[-1] - t[ic]))),
               "cloud_top_km": float(height[ic]), "cloud_top_hPa": float(pressure[ic]), "T_c": float(temp[ic])}
      if ch in planck:
         p = planck[ch]
         b = lambda T: float(t2r(np.asarray(T, dtype=np.float64), p)[0])
         layer_t = 0.5 * (temp[:-1] + temp[1:])
         up = sum(b(layer_t[k]) * (np.exp(-t[k] / mu) - np.exp(-t[k + 1] / mu)) for k in range(ic))
         down = sum(b(layer_t[k]) * (np.exp(-(t[ic] - t[k + 1]) / mu) - np.exp(-(t[ic] - t[k]) / mu)) for k in range(ic))
         below = sum(b(layer_t[k]) * (np.exp(-(t[k] - t[ic]) / mu) - np.exp(-(t[k + 1] - t[ic]) / mu))
                     for k in range(ic, t.size - 1))
         surface = SURFACE_EMISSIVITY * b(temp[-1]) * np.exp(-(t[-1] - t[ic]) / mu)
         entry.update({"Tac_lw": float(np.exp(-t[ic] / mu)), "Rac_up": float(up), "Rac_dwn": float(down),
                       "Rbc_up": float(surface + below), "B_c": b(temp[ic]), "T_s": float(temp[-1])})
      out[ch] = entry
   return out
