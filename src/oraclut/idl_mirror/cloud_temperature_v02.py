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

Emission layers (emission_layers).  The production cloud is one DISORT layer
of optical depth D (in the channel), single-scattering albedo and phase
moments; the profile gives B(F) = B_nu[T(F)], the Planck radiance of the
channel's DISORT window at cumulative optical-depth fraction F (the same
temperature profile for every channel).  DISORT's thermal source is linear in
optical depth within a layer, with slope XR1 = (B_bottom - B_top) / dtau; in
single precision a layer whose source changes much over a small optical
depth loses accuracy (round-off of about 2e-8 (dB / B) / dtau^2 relative;
validation/v25/temperature_profiles/conditioning.py): the 40 equal
sub-layers of the historical method gave E_md = -1.5 at tau = 1e-10.  The
emission call is therefore built as follows.

* Linear layering (sublayer_bounds): n = min(EMISSION_LAYERS, floor(D / h))
  sub-layers, equal while D <= n h, otherwise a geometric progression in
  optical depth from a top sub-layer h thick, so that the emitting top of a
  thick cloud is resolved; h = TOP_SUBLAYER_TAU, doubled until every
  sub-layer k satisfies (|dB_k| / B_mean) / dtau_k^2 <= PLANCK_SLOPE_LIMIT
  (the sub-layers are resolvable by DISORT), and used if that leaves at least
  STAIRCASE_CELLS sub-layers.  The source is the L2
  projection of B(F) onto continuous piecewise-linear functions
  (planck_projection), which preserves the integral of B; DISORT receives the
  temperatures whose window Planck radiance is the projected value at each
  boundary.
* Isothermal staircase, when no linear layering is resolvable (an optically
  thin cloud whose temperature still spans the profile):
  K = min(STAIRCASE_CELLS, max(1, floor(D / STAIRCASE_CELL_TAU))) equal
  sub-layers, each isothermal at the temperature of its mean Planck
  radiance, joined by sub-layers of zero optical depth that carry the
  temperature steps (DISORT gives a layer without optical depth no source
  slope, XR1 = 0).  An isothermal layer has no source slope, so the call is
  free of the cancellation above; cells thinner than STAIRCASE_CELL_TAU only
  add DISORT's ordinary thin-layer round-off, so below it the cloud is one
  layer with the optical-depth mean B_mean.  The staircase preserves the
  integral of B over every cell; for an optically thin cloud the emergent
  radiance (1 - w) int_0^D B dt / mu + O(D^2) depends on the source only
  through B_mean = int_0^1 B(F) dF, so E_md -> (1 - w) D B_mean / (mu B_top)
  -> 0 as D -> 0.
"""

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ..radiative_transfer.legacy_disort import plkavg

# Profiles known to create_orac_cloud_lut and the substance (.mm) each is defined for
CLOUD_PROFILES = ("isothermal", "cirrostratus", "wet_adiabat")
PROFILE_SUBSTANCE = {"cirrostratus": "water-ice", "wet_adiabat": "liquid-water"}

CIRROSTRATUS_FILE = Path(__file__).resolve().parents[3] / "references" / "data" / "ocalut_cloudprofile_Cirrostratus.dat"
PROFILE_ROWS = 99                        # rows of each supplied profile
LAPSE_RATE = 8.0                         # K/km: temperature above the cloud top = 8 K/km x depth (verified on reading)

BETA_EXT_055_LIQUID = 20.0               # km^-1: path length z = tau_055 / beta of the liquid-water adiabat
P_TOP_LIQUID = 628.0                     # hPa: reference pressure at the cloud top of the liquid-water adiabat
ADIABAT_STEP = 0.01                      # km: Runge-Kutta step of the adiabat

EMISSION_LAYERS = 46                     # largest number of linear sub-layers: MXCLY of the production DISORT
DISORT_MAX_LAYERS = 46                   # MXCLY of the production DISORT (create_orac_lut/disort2/src/DISORTfunctions.f)
TOP_SUBLAYER_TAU = 0.01                  # first (top) sub-layer of the geometric layering (channel optical depth)
PLANCK_SLOPE_LIMIT = 1000.0              # largest (|dB_k| / B_mean) / dtau_k^2 of a linear sub-layer
STAIRCASE_CELLS = 23                     # largest number of isothermal cells of the staircase (23 cells + 22 joints <= 46 layers)
STAIRCASE_CELL_TAU = 1e-3                # thinnest isothermal cell (channel optical depth)
QUADRATURE_PANELS = 64                   # equal panels of [0, 1] added to the quadrature of planck_projection
GAUSS_X, GAUSS_W = np.polynomial.legendre.leggauss(3)

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
              f"linear in optical depth; emission layers: up to {EMISSION_LAYERS} geometric (top {TOP_SUBLAYER_TAU:g}) with projected "
              f"Planck source, or {STAIRCASE_CELLS} isothermal cells; E_md / B(T_top)")
   return (f"wet_adiabat: saturated liquid-water adiabat from T_top {cloudprof.t_top:g} K, {P_TOP_LIQUID:g} hPa along "
           f"z = tau_055 / {BETA_EXT_055_LIQUID:g} km^-1 (to {cloudprof.z[-1]:g} km: {cloudprof.t[-1]:.2f} K, "
           f"{cloudprof.p[-1]:.1f} hPa); emission layers: up to {EMISSION_LAYERS} geometric (top {TOP_SUBLAYER_TAU:g}) with projected "
           f"Planck source, or {STAIRCASE_CELLS} isothermal cells; E_md / B(T_top)")


def profile_breakpoints(cloudprof):
   """Fractions F at which the profile temperature has kinks (the row midpoints of the supplied profile)."""

   if cloudprof.name == "cirrostratus":
      extinction = cloudprof.profstr.extinction
      return np.cumsum(extinction) - 0.5 * extinction
   return np.zeros(0)


def inverse_planck(radiance, wnlo, wnhi, t):
   """Temperature whose DISORT window Planck radiance (PLKAVG) is ``radiance``; Newton iteration from ``t``.

   PLKAVG takes a float32 temperature and resolves the radiance to about
   1e-5, the precision with which DISORT itself evaluates its TEMPER; the
   iterate closest to ``radiance`` is returned once within 1e-5 (3e-5 at
   worst).
   """

   best, misfit = t, np.inf
   for _ in range(20):
      value = plkavg(wnlo, wnhi, t)
      if abs(value - radiance) < misfit:
         best, misfit = t, abs(value - radiance)
      if misfit <= 1e-5 * radiance:
         return best
      slope = (plkavg(wnlo, wnhi, t + 0.01) - plkavg(wnlo, wnhi, t - 0.01)) / 0.02
      t -= (value - radiance) / slope
   if misfit <= 3e-5 * radiance:
      return best
   raise ValueError(f"no temperature with window Planck radiance {radiance}")


def sublayer_bounds(dtau, top):
   """Boundaries (cumulative optical-depth fractions) of the linear emission sub-layers for a top sub-layer ``top`` thick.

   n = min(EMISSION_LAYERS, floor(dtau / top)) sub-layers: equal ones while
   dtau <= n top, else a geometric progression from a top sub-layer ``top``
   thick (ratio by bisection).
   """

   n = min(EMISSION_LAYERS, int(dtau / top))
   if n <= 1:
      return np.asarray([0.0, 1.0])
   if dtau <= n * top:
      return np.arange(n + 1) / n
   low, high = 1.0, 2.0
   while top * (high ** n - 1.0) / (high - 1.0) < dtau:
      high *= 2.0
   for _ in range(200):
      r = 0.5 * (low + high)
      low, high = (r, high) if top * (r ** n - 1.0) / (r - 1.0) < dtau else (low, r)
   r = 0.5 * (low + high)
   bounds = top * (r ** np.arange(n + 1) - 1.0) / (r - 1.0) / dtau
   bounds[-1] = 1.0
   return bounds


def planck_projection(cloudprof, tau_055, bounds, wnlo, wnhi):
   """L2 projection of the window Planck radiance B(F) of the profile onto piecewise-linear functions on ``bounds``.

   Returns the radiances at the boundaries (planck_means gives the
   projection onto cell constants).  The integrals of B against the hat functions are
   evaluated with 3-point
   Gauss-Legendre quadrature on panels bounded by the sub-layer boundaries,
   the profile's kinks and QUADRATURE_PANELS equal steps.
   """

   edges = np.unique(np.concatenate([bounds, profile_breakpoints(cloudprof), np.linspace(0.0, 1.0, QUADRATURE_PANELS + 1)]))
   width = np.diff(edges)
   f = (edges[:-1, None] + 0.5 * width[:, None] * (1.0 + GAUSS_X[None, :])).ravel()
   w = (0.5 * width[:, None] * GAUSS_W[None, :]).ravel()
   b = np.asarray([plkavg(wnlo, wnhi, float(t)) for t in cloud_temperature(cloudprof, tau_055, f)])
   n = bounds.size - 1
   h = np.diff(bounds)
   cell = np.clip(np.searchsorted(bounds, f, side="right") - 1, 0, n - 1)
   x = (f - bounds[cell]) / h[cell]                          # position within the sub-layer, 0 .. 1
   load = np.zeros(n + 1)
   np.add.at(load, cell, w * b * (1.0 - x))
   np.add.at(load, cell + 1, w * b * x)
   mass = np.zeros((n + 1, n + 1))
   for k in range(n):
      mass[k:k + 2, k:k + 2] += h[k] * np.asarray([[2.0, 1.0], [1.0, 2.0]]) / 6.0
   return np.linalg.solve(mass, load)


def planck_means(cloudprof, tau_055, bounds, wnlo, wnhi):
   """Mean window Planck radiance of the profile in every cell between ``bounds`` (the projection onto cell constants)."""

   edges = np.unique(np.concatenate([bounds, profile_breakpoints(cloudprof), np.linspace(0.0, 1.0, QUADRATURE_PANELS + 1)]))
   width = np.diff(edges)
   f = (edges[:-1, None] + 0.5 * width[:, None] * (1.0 + GAUSS_X[None, :])).ravel()
   w = (0.5 * width[:, None] * GAUSS_W[None, :]).ravel()
   b = np.asarray([plkavg(wnlo, wnhi, float(t)) for t in cloud_temperature(cloudprof, tau_055, f)])
   cell = np.clip(np.searchsorted(bounds, f, side="right") - 1, 0, bounds.size - 2)
   return np.bincount(cell, weights=w * b, minlength=bounds.size - 1) / np.diff(bounds)


def emission_layers(emtau, emssa, empmo, tau_055, cloudprof, temp, wnlo, wnhi):
   """Inputs (emtau, emssa, empmo, emtemp) of the emission call for the in-cloud layer (module docstring).

   ``emtau``, ``emssa``, ``empmo[:, 0]`` describe the one in-cloud layer in
   the channel; ``tau_055`` is the node's 0.55 um optical depth, which fixes
   the cloud's temperature profile; ``wnlo`` and ``wnhi`` are the channel's
   DISORT window; ``temp`` is the emission temperature of the calling code,
   which must be the cloud-top temperature of the profile.  The sub-layers
   keep the layer's single-scattering albedo and moments; their optical
   depths add up to the layer's.
   """

   f32 = np.float32
   if f32(temp) != f32(cloudprof.t_top):
      raise ValueError(f"the emission temperature {temp} differs from the cloud-top temperature {cloudprof.t_top} of the profile")
   if emtau.size != 1:
      raise ValueError(f"the emission layering is defined for one in-cloud layer, not {emtau.size}")
   dtau = float(emtau[0])
   tau_055 = float(tau_055)
   mean = float(planck_means(cloudprof, tau_055, np.asarray([0.0, 1.0]), wnlo, wnhi)[0])
   # linear layering: the finest top sub-layer whose sub-layers DISORT resolves
   top = TOP_SUBLAYER_TAU
   subtemp = None
   while top <= dtau:
      bounds = sublayer_bounds(dtau, top)
      if bounds.size - 1 < STAIRCASE_CELLS:
         break                                                 # coarser than the staircase: use the staircase
      radiance = planck_projection(cloudprof, tau_055, bounds, wnlo, wnhi)
      thickness = dtau * np.diff(bounds)
      if np.max(np.abs(np.diff(radiance)) / mean / thickness ** 2) <= PLANCK_SLOPE_LIMIT:
         if np.any(radiance <= 0.0):
            raise ValueError("the projected Planck radiance is not positive")
         profile = cloud_temperature(cloudprof, tau_055, bounds)
         subtemp = np.asarray([inverse_planck(b, wnlo, wnhi, float(t)) for b, t in zip(radiance, profile)], dtype=f32)
         break
      top *= 2.0
   if subtemp is None:
      # isothermal staircase: equal cells at their mean Planck radiance, joined by layers without optical depth
      ncell = min(STAIRCASE_CELLS, max(1, int(dtau / STAIRCASE_CELL_TAU)))
      bounds = np.arange(ncell + 1) / ncell
      cells = cloud_temperature(cloudprof, tau_055, 0.5 * (bounds[:-1] + bounds[1:]))
      level = np.asarray([inverse_planck(b, wnlo, wnhi, float(t))
                          for b, t in zip(planck_means(cloudprof, tau_055, bounds, wnlo, wnhi), cells)])
      thickness = np.zeros(2 * ncell - 1)
      thickness[0::2] = dtau * np.diff(bounds)
      subtemp = np.concatenate([[level[0]], np.repeat(level, 2)[1:]]).astype(f32)
   n = thickness.size
   if n > DISORT_MAX_LAYERS:
      raise ValueError(f"{n} emission sub-layers exceed the {DISORT_MAX_LAYERS} layers of DISORT (MXCLY)")
   subtau = thickness.astype(f32)
   last = int(np.flatnonzero(subtau > 0.0)[-1]) if np.any(subtau > 0.0) else n - 1
   subtau[last] = f32(dtau - float(np.sum(np.delete(subtau, last), dtype=np.float64)))
   if np.any(subtau < 0.0):
      raise ValueError("an emission sub-layer has a negative optical depth")
   return subtau, np.repeat(emssa, n).astype(f32), np.asfortranarray(np.repeat(empmo, n, axis=1).astype(f32)), subtemp
