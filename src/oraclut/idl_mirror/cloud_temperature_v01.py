"""Cloud temperature profile of the thermal-emission DISORT call of create_orac_cloud_lut (LUT release V25).

Up to V24 the emission call of create_orac_cloud_lut gives every in-cloud
layer one temperature (create_orac_cloud_lut.pro line 561, ``temp = 250.0``)
and divides the emitted radiance by the Planck radiance at that temperature,
so E_md is the emissivity of an isothermal cloud.  LUT release V25 keeps
250 K as the cloud-top temperature and as the normalising temperature
(owner decision 2026-10-15) and lets the temperature rise below the cloud top:

   'cirrostratus'  water-ice LUTs: the supplied cirrostratus profile
                   (references/data/ocalut_cloudprofile_Cirrostratus.dat)
   'wet_adiabat'   liquid-water LUTs: the saturated liquid-water adiabat
   'isothermal'    the V24 calculation (no profile)

This is the method of the historical V25 (commits 7f6b4e2, 13b9df3, 5722280,
990d970; validation/v25/results/REPORT_v25_cloud_temperature_profile.md),
re-implemented for the V24 code with one change decided by the owner: the
cloud top is at 250 K (the V24 temperature), not the historical 240 K.

The supplied profile (origin
/tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf, cloud
type Cirrostratus; P. Watts' vertically inhomogeneous LUT tests for OCA)
tabulates, for 13 total optical depths COT = 0.0625 ... 256, 99 rows of
% extinction f_i, layer optical depth COT f_i / 100, cumulative optical depth,
depth below the cloud top z_i = H(COT) i / 98 and temperature above the cloud
top 8 K/km z_i.  read_cloud_profile verifies all of this on every read: the
extinction shape f_i is common to all COTs and the only COT-dependent quantity
is the cloud depth H(COT) (2 km at COT 0.0625 to 11 km at COT 128 and 256).
Row i is taken as a sample at depth z_i carrying the optical depth
COT f_i / 100, so the cumulative optical-depth fraction at z_i is
F_i = (sum_{j<=i} f_j - f_i / 2) / 100 (midpoint rule; F = 0 at the top and
F = 1 at the base).  This gives a COT-independent depth fraction zeta(F) and

   T(F; tau_055) = T_top + 8 K/km H(tau_055) zeta(F),

with H interpolated linearly in COT between the supplied values, held at
H(0.0625) below COT 0.0625 and not extrapolated above COT 256.  The profile
is mapped through the cumulative optical depth, not through a physical
height: no extinction coefficient enters.

The liquid-water adiabat follows the saturated (pseudo-adiabatic) lapse rate
over liquid water downward from the reference state (T_top, 628 hPa) along
the path length z = x / beta below the cloud top, x the cumulative 0.55 um
optical depth and beta = 20 km^-1 a prescribed extinction coefficient (no
cap, altitude, surface or freezing limit; the pressure is evolved
hydrostatically along the path).  At tau_055 = 256 the path is 12.8 km long;
that deep state is not a physical cloud but contributes little to the
emission (REPORT_v25_temperature_profiles.md).

In the emission call each in-cloud layer is divided into EMISSION_LAYERS
equal-optical-depth sub-layers with the layer's single-scattering albedo and
phase moments; the profile gives the temperature at every sub-layer boundary
(DISORT's Planck source is linear in optical depth within a layer).  The
diffuse and direct-beam calls are not changed.
"""

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

# Profiles known to create_orac_cloud_lut and the substance (.mm) each is defined for
CLOUD_PROFILES = ("isothermal", "cirrostratus", "wet_adiabat")
PROFILE_SUBSTANCE = {"cirrostratus": "water-ice", "wet_adiabat": "liquid-water"}

CIRROSTRATUS_FILE = Path(__file__).resolve().parents[3] / "references" / "data" / "ocalut_cloudprofile_Cirrostratus.dat"
PROFILE_ROWS = 99                        # rows of each supplied profile
LAPSE_RATE = 8.0                         # K/km: temperature above the cloud top = 8 K/km x depth (verified on reading)

BETA_EXT_055_LIQUID = 20.0               # km^-1: path length z = tau_055 / beta of the liquid-water adiabat
P_TOP_LIQUID = 628.0                     # hPa: reference pressure at the cloud top of the liquid-water adiabat
ADIABAT_STEP = 0.01                      # km: Runge-Kutta step of the adiabat

EMISSION_LAYERS = 40                     # sub-layers of each in-cloud layer in the emission call
DISORT_MAX_LAYERS = 46                   # MXCLY of the production DISORT (create_orac_lut/disort2/src/DISORTfunctions.f)

# Thermodynamic constants (SI) of the liquid-water adiabat
G0 = 9.80665                             # standard gravity, m s^-2
R_DRY = 287.04                           # gas constant of dry air, J kg^-1 K^-1
R_VAPOUR = 461.5                         # gas constant of water vapour, J kg^-1 K^-1
C_P_DRY = 1005.7                         # specific heat of dry air at constant pressure, J kg^-1 K^-1
EPSILON = R_DRY / R_VAPOUR               # molecular mass ratio water / dry air (0.622)


def read_cloud_profile(path):
   """Read and verify a supplied cloud-profile file (format of ocalut_cloudprofile_*.dat).

   Returns a structure with the cloud type, origin, SHA-256 digest, the
   total optical depths cot, the common extinction fractions (summing to 1)
   and the cloud depth H (km) of each profile.
   """

   path = Path(path)
   raw = path.read_bytes()
   lines = raw.decode().splitlines()
   head = lines[0].split("#")
   rows = lines[1].split("#")
   if len(head) != 2 or head[1].strip() != "COTs" or len(rows) != 2 or rows[1].strip() != "layers":
      raise ValueError(f"{path}: unexpected header {lines[0]!r} / {lines[1]!r}")
   ncot, nrows = int(head[0]), int(rows[0])
   if nrows != PROFILE_ROWS:
      raise ValueError(f"{path}: {nrows} rows per profile, expected {PROFILE_ROWS}")
   if not lines[2].startswith("Origin:") or not lines[3].startswith("Cloud type:"):
      raise ValueError(f"{path}: the origin and cloud type lines are missing")
   cot = []
   data = []
   i = 4
   while i < len(lines):
      if lines[i].strip() == "":
         i += 1
         continue
      words = lines[i].split()
      if len(words) != 4 or " ".join(words[1:]) != "LUT Cot value" or "% Extinction" not in lines[i + 1]:
         raise ValueError(f"{path}: malformed profile header at line {i + 1}: {lines[i]!r}")
      block = [[float(v) for v in lines[j].split()] for j in range(i + 2, i + 2 + nrows)]
      if len(block) != nrows or any(len(b) != 5 for b in block):
         raise ValueError(f"{path}: the profile at line {i + 1} does not have {nrows} rows of 5 values")
      cot.append(float(words[0]))
      data.append(block)
      i += 2 + nrows
   if len(cot) != ncot:
      raise ValueError(f"{path}: {len(cot)} profiles read, the header says {ncot}")
   cot = np.asarray(cot, dtype=np.float64)
   data = np.asarray(data, dtype=np.float64)                 # (ncot, rows, 5)
   ext, dcot, cum, z, dt = (data[:, :, k] for k in range(5))

   # The structure the interpolation relies on, verified for every profile
   if not np.all(np.isfinite(data)) or np.any(np.diff(cot) <= 0.0):
      raise ValueError(f"{path}: non-finite values or optical depths not increasing")
   if np.any(ext != ext[0]) or np.any(ext < 0.0) or abs(ext[0].sum() - 100.0) > 1e-3:
      raise ValueError(f"{path}: the extinction percentages differ between profiles or do not sum to 100")
   for k in range(cot.size):
      if np.any(dcot[k] < 0.0) or abs(dcot[k].sum() - cot[k]) > 1e-4 * cot[k] or abs(cum[k, -1] - cot[k]) > 1e-4 * cot[k]:
         raise ValueError(f"{path}: the layer optical depths of profile {cot[k]} do not add up to it")
      if np.abs(cum[k] - np.cumsum(dcot[k])).max() > 1e-5 * max(1.0, cot[k]) or np.any(np.diff(cum[k]) <= 0.0):
         raise ValueError(f"{path}: the cumulative optical depth of profile {cot[k]} is not the running sum")
      if z[k, 0] != 0.0 or dt[k, 0] != 0.0 or np.any(np.diff(z[k]) <= 0.0) or np.any(np.diff(dt[k]) <= 0.0):
         raise ValueError(f"{path}: depth and temperature of profile {cot[k]} do not rise from 0 at the cloud top")
      if np.abs(dt[k] - LAPSE_RATE * z[k]).max() > 1e-4:
         raise ValueError(f"{path}: profile {cot[k]} does not follow dT = {LAPSE_RATE} K/km x depth")
      if np.abs(z[k] / z[k, -1] - np.arange(nrows) / (nrows - 1.0)).max() > 1e-5:
         raise ValueError(f"{path}: the rows of profile {cot[k]} are not equally spaced in depth")
   if np.any(np.diff(z[:, -1]) < 0.0):
      raise ValueError(f"{path}: the cloud depth decreases with optical depth")
   return SimpleNamespace(path=path, sha256=hashlib.sha256(raw).hexdigest(),
                          origin=lines[2].split(":", 1)[1].strip(), cloud_type=lines[3].split(":", 1)[1].strip(),
                          cot=cot, extinction=ext[0] / 100.0, depth=z[:, -1].copy())


def profile_depth_fraction(profstr, fraction):
   """zeta(F): depth fraction z / H at cumulative optical-depth fraction F (midpoint rule)."""

   midpoints = np.cumsum(profstr.extinction) - 0.5 * profstr.extinction
   zeta = np.arange(PROFILE_ROWS) / (PROFILE_ROWS - 1.0)
   return np.interp(fraction, np.r_[0.0, midpoints, 1.0], np.r_[0.0, zeta, 1.0])


def profile_depth(profstr, tau_055):
   """H(tau_055) (km), linear in COT between the supplied profiles; held below the smallest, never extrapolated above."""

   if tau_055 > profstr.cot[-1] * (1.0 + 1e-9):
      raise ValueError(f"the optical depth {tau_055} exceeds the largest supplied profile ({profstr.cot[-1]})")
   return float(np.interp(tau_055, profstr.cot, profstr.depth))


def saturation_vapour_pressure_liquid(t):
   """Saturation vapour pressure over liquid water (Pa); Murphy and Koop (2005, QJRMS 131, 1539) Eq. 10."""

   return np.exp(54.842763 - 6763.22 / t - 4.210 * np.log(t) + 0.000367 * t
                 + np.tanh(0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 * np.log(t) + 0.014025 * t))


def latent_heat_vaporisation(t):
   """Latent heat of vaporisation (J kg^-1); Rogers and Yau (1989)."""

   return 2.501e6 - 2370.0 * (t - 273.15)


def wet_adiabatic_lapse_rate(t, p):
   """Saturated (pseudo-adiabatic) lapse rate over liquid water, -dT/dz (K m^-1); AMS Glossary."""

   e_s = saturation_vapour_pressure_liquid(t)
   if not np.isfinite(e_s) or e_s >= p:
      raise ValueError(f"the saturation vapour pressure {e_s:.4g} Pa reaches the pressure {p:.4g} Pa at {t:.1f} K")
   lv = latent_heat_vaporisation(t)
   r_s = EPSILON * e_s / (p - e_s)
   return G0 * (1.0 + lv * r_s / (R_DRY * t)) / (C_P_DRY + lv ** 2 * r_s * EPSILON / (R_DRY * t ** 2))


def wet_adiabat(t_top, p_top, z_max):
   """(z km, T K, p hPa) of the saturated liquid-water adiabat from (t_top, p_top) down to z_max.

   Fourth-order Runge-Kutta in steps of ADIABAT_STEP of
      dT/dz = Gamma_w(T, p),   dp/dz = p g / (R_d T)   (hydrostatic, dry air).
   z is the path length below the cloud top, not an altitude.
   """

   nsteps = max(1, int(np.ceil(z_max / ADIABAT_STEP - 1e-9)))
   z = np.linspace(0.0, max(z_max, ADIABAT_STEP), nsteps + 1)
   t = np.empty(nsteps + 1)
   p = np.empty(nsteps + 1)
   t[0] = t_top
   p[0] = 100.0 * p_top

   def rates(tk, pa):
      return wet_adiabatic_lapse_rate(tk, pa), pa * G0 / (R_DRY * tk)

   for i in range(nsteps):
      h = 1000.0 * (z[i + 1] - z[i])
      k1t, k1p = rates(t[i], p[i])
      k2t, k2p = rates(t[i] + 0.5 * h * k1t, p[i] + 0.5 * h * k1p)
      k3t, k3p = rates(t[i] + 0.5 * h * k2t, p[i] + 0.5 * h * k2p)
      k4t, k4p = rates(t[i] + h * k3t, p[i] + h * k3p)
      t[i + 1] = t[i] + h * (k1t + 2.0 * k2t + 2.0 * k3t + k4t) / 6.0
      p[i + 1] = p[i] + h * (k1p + 2.0 * k2p + 2.0 * k3p + k4p) / 6.0
      if not (np.isfinite(t[i + 1]) and np.isfinite(p[i + 1])):
         raise ValueError(f"the wet adiabat is not finite {z[i + 1]:.2f} km below the cloud top")
   return z, t, p / 100.0


def setup_cloud_profile(name, max_opd, t_top):
   """The cloud profile of one LUT (a structure, like setup_disort), for optical depths up to max_opd.

   'cirrostratus' reads the supplied profile file; 'wet_adiabat' tabulates
   the adiabat from (t_top, P_TOP_LIQUID) to z = max_opd / BETA_EXT_055_LIQUID.
   """

   if name == "cirrostratus":
      profstr = read_cloud_profile(CIRROSTRATUS_FILE)
      profile_depth(profstr, float(max_opd))
      return SimpleNamespace(name=name, t_top=float(t_top), profstr=profstr)
   if name == "wet_adiabat":
      z, t, p = wet_adiabat(float(t_top), P_TOP_LIQUID, float(max_opd) / BETA_EXT_055_LIQUID)
      return SimpleNamespace(name=name, t_top=float(t_top), z=z, t=t, p=p)
   raise ValueError(f"unknown cloud profile {name!r}; known: {CLOUD_PROFILES}")


def cloud_temperature(cloudprof, tau_055, fraction):
   """Temperature (K) at cumulative optical-depth fraction(s) F of a cloud of 0.55 um optical depth tau_055."""

   fraction = np.clip(np.asarray(fraction, dtype=np.float64), 0.0, 1.0)
   if cloudprof.name == "cirrostratus":
      profstr = cloudprof.profstr
      return cloudprof.t_top + LAPSE_RATE * profile_depth(profstr, tau_055) * profile_depth_fraction(profstr, fraction)
   z = fraction * tau_055 / BETA_EXT_055_LIQUID
   if np.any(z > cloudprof.z[-1] * (1.0 + 1e-9)):
      raise ValueError(f"the wet adiabat is tabulated to {cloudprof.z[-1]:.3f} km; {z.max():.3f} km requested")
   return np.interp(z, cloudprof.z, cloudprof.t)


def describe_cloud_profile(cloudprof):
   """One line describing the profile (printed and recorded in the provenance)."""

   if cloudprof.name == "cirrostratus":
      profstr = cloudprof.profstr
      return (f"cirrostratus: supplied profile {profstr.path.name} (sha256 {profstr.sha256[:12]}, {profstr.cloud_type}, "
              f"{profstr.cot.size} optical depths {profstr.cot[0]:g}-{profstr.cot[-1]:g}); T_top {cloudprof.t_top:g} K, "
              f"T = T_top + {LAPSE_RATE:g} K/km x depth, cloud depth {profstr.depth[0]:g}-{profstr.depth[-1]:g} km "
              f"linear in optical depth; {EMISSION_LAYERS} emission layers per cloud layer; E_md / B(T_top)")
   return (f"wet_adiabat: saturated liquid-water adiabat from T_top {cloudprof.t_top:g} K, {P_TOP_LIQUID:g} hPa along "
           f"z = tau_055 / {BETA_EXT_055_LIQUID:g} km^-1 (to {cloudprof.z[-1]:g} km: {cloudprof.t[-1]:.2f} K, "
           f"{cloudprof.p[-1]:.1f} hPa); {EMISSION_LAYERS} emission layers per cloud layer; E_md / B(T_top)")


def emission_layers(emtau, emssa, empmo, scatreltau, tau_055, cloudprof, temp):
   """Sub-layered inputs (emtau, emssa, empmo, emtemp) of the emission call for the in-cloud layers.

   Each in-cloud layer (top down) becomes EMISSION_LAYERS sub-layers of equal
   optical depth with the layer's single-scattering albedo and moments.  The
   boundary temperatures follow the cumulative particle optical-depth
   fraction (scatreltau) of the 0.55 um optical depth tau_055 through the
   cloud profile, so every channel sees the same cloud.  ``temp`` is the
   emission temperature of the calling code, which must be the cloud-top
   temperature of the profile.
   """

   f32 = np.float32
   n = EMISSION_LAYERS
   if f32(temp) != f32(cloudprof.t_top):
      raise ValueError(f"the emission temperature {temp} differs from the cloud-top temperature {cloudprof.t_top} of the profile")
   if emtau.size * n > DISORT_MAX_LAYERS:
      raise ValueError(f"{emtau.size} in-cloud layers x {n} sub-layers exceed the {DISORT_MAX_LAYERS} layers of DISORT (MXCLY)")
   weights = np.asarray(scatreltau, dtype=np.float64)
   if np.any(weights <= 0.0):
      raise ValueError("every in-cloud layer must hold particles")
   subtau = np.repeat(emtau / f32(n), n).astype(f32)
   subssa = np.repeat(emssa, n).astype(f32)
   subpmo = np.asfortranarray(np.repeat(empmo, n, axis=1).astype(f32))
   fraction = np.concatenate([[0.0], np.cumsum(np.repeat(weights / weights.sum() / n, n))])
   fraction[-1] = 1.0
   subtemp = cloud_temperature(cloudprof, float(tau_055), fraction).astype(f32)
   return subtau, subssa, subpmo, subtemp
