"""Cloud temperature profile for the V25 thermal-emission calculation (create_orac_cloud_lut).

Up to V24 the emission DISORT call gives every in-cloud layer one temperature
(create_orac_cloud_lut.pro line 561, ``temp = 250.0``), so the LUT emissivity
E_md = UU / B(temp) is that of an isothermal cloud and does not depend on the
value of temp.  V25 keeps a fixed cloud-top reference temperature,
T_TOP_K = 240 K, and lets the temperature increase downward through the cloud
along a saturated adiabat.  E_md stays normalised by B(T_TOP_K), so it carries
the differential emission of the warmer lower cloud.  The LUT grids are
unchanged: cloud temperature is not a LUT dimension and no cloud-base
temperature is retrieved.

The V25 model (its constants are representative, deliberate choices, not
microphysical results):

* T_top = T_TOP_K at the cloud top of the production atmosphere, the upper
  boundary of the first layer that holds particles (the .mm profile; for the
  production profiles and mls.atm this is 4 km, 628 hPa).
* Nominal geometric depth H_nominal = tau_055 / beta_ext_055 with the
  representative volume extinction coefficient beta_ext_055 of the substance
  (LIQUID_WATER: 20 km^-1; WATER_ICE: 1 km^-1), capped at H_max (2.5 km; 6 km):
  H = min(H_nominal, H_max).
* The whole cloud optical depth is mapped onto 0 <= z <= H: a cumulative
  0.55-um optical depth x below the cloud top (0 <= x <= tau_055) lies at
  z(x) = H x / tau_055, and z = 0 everywhere when tau_055 = 0.  A capped cloud
  therefore spans 0..H continuously, never an isothermal remainder.  tau_055 is
  the LUT optical-depth coordinate (lutstr.opd): the particle optical depth of
  every channel is tau_055 * BextRat, BextRat = Bext / Bext(0.55 um), so one
  physical cloud structure serves every spectral channel.
* T(z) is the saturated adiabat of the substance, integrated downward from
  (T_top, p_top) with the pressure of the production atmosphere at height
  z_top - z (log-linear in height; continued below the lowest level with the
  lowest layer's scale height, which only an ice cloud deeper than the cloud
  top height reaches).

beta_ext_055 is prescribed.  It is not derived from the per-particle extinction
cross-section or the normalised size distribution of
generate_scattering_properties: the LUT fixes no number concentration, so the
microphysics does not define a physical thickness.

Thermodynamics.  The repository had no saturated-adiabat routine, so this
module holds a compact one:

* saturation vapour pressure over liquid water (valid for supercooled water)
  and over ice, and the latent heat of sublimation: Murphy and Koop (2005,
  Q. J. R. Meteorol. Soc. 131, 1539-1565), Eqs. 10, 7 and 5;
* latent heat of vaporisation L_v = 2.501e6 - 2370 (T - 273.15) J kg^-1
  (Rogers and Yau 1989);
* pseudo-adiabatic (saturated) lapse rate
      Gamma_s = g (1 + L r_s / (R_d T)) / (c_pd + L^2 r_s eps / (R_d T^2)),
      r_s = eps e_s / (p - e_s)
  (AMS Glossary of Meteorology, "moist-adiabatic lapse rate"), with e_sw and
  L_v for liquid-water clouds and e_si and L_s for ice clouds.  Mixed-phase
  physics is not represented.

The emission DISORT call subdivides each in-cloud layer into EMISSION_SUBLAYERS
equal-optical-depth sub-layers (same single-scattering albedo and phase
moments) and gives DISORT the adiabat temperature at every sub-layer boundary;
DISORT takes the Planck source as linear in optical depth within a layer and
warns above a 10 K step (CHEKIN).  The production cloud is one atmosphere
layer, so without the subdivision the whole cloud would be one linear ramp.
"""

from dataclasses import dataclass

import numpy as np

# ---------------------------------------------------------------------------
# V25 model constants
# ---------------------------------------------------------------------------
T_TOP_K = 240.0                          # cloud-top reference temperature, as the existing calculation

# substance (.mm "substance") -> (phase, beta_ext_055 km^-1, H_max km)
CLOUD_MODELS = {
   "liquid-water": ("liquid", 20.0, 2.5),
   "water-ice": ("ice", 1.0, 6.0),
}

# Sub-layers of every in-cloud layer in the emission DISORT call.  The largest
# temperature span is the capped ice cloud (6 km, 42 K), so 12 sub-layers keep
# every DISORT temperature step below 5 K (4.6 K in the first 0.5 km of ice;
# below 2 K for liquid water), well inside DISORT's 10 K guidance.
EMISSION_SUBLAYERS = 12

# DISORT 2.0 as compiled for production accepts at most MXCLY = 46 computational
# layers (create_orac_lut/disort2/src/DISORTfunctions.f), so the sub-layered
# emission call is limited to that many layers in total.
DISORT_MAX_LAYERS = 46

# Integration step of the adiabat (km)
PROFILE_STEP_KM = 0.01

# ---------------------------------------------------------------------------
# Thermodynamic constants (SI)
# ---------------------------------------------------------------------------
G0 = 9.80665                             # standard gravity, m s^-2
R_DRY = 287.04                           # gas constant of dry air, J kg^-1 K^-1
R_VAPOUR = 461.5                         # gas constant of water vapour, J kg^-1 K^-1
C_P_DRY = 1005.7                         # specific heat of dry air at constant pressure, J kg^-1 K^-1
EPSILON = R_DRY / R_VAPOUR               # molecular mass ratio water / dry air (0.622)
M_WATER = 0.018015                       # molar mass of water, kg mol^-1
T_FREEZE = 273.15


def saturation_vapour_pressure_liquid(temperature_k):
   """Saturation vapour pressure over liquid water (Pa); Murphy and Koop (2005) Eq. 10, 123-332 K."""

   t = np.asarray(temperature_k, dtype=np.float64)
   return np.exp(54.842763 - 6763.22 / t - 4.210 * np.log(t) + 0.000367 * t
                 + np.tanh(0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 * np.log(t) + 0.014025 * t))


def saturation_vapour_pressure_ice(temperature_k):
   """Saturation vapour pressure over ice (Pa); Murphy and Koop (2005) Eq. 7, above 110 K."""

   t = np.asarray(temperature_k, dtype=np.float64)
   return np.exp(9.550426 - 5723.265 / t + 3.53068 * np.log(t) - 0.00728332 * t)


def latent_heat_vaporisation(temperature_k):
   """Latent heat of vaporisation (J kg^-1); Rogers and Yau (1989)."""

   return 2.501e6 - 2370.0 * (np.asarray(temperature_k, dtype=np.float64) - T_FREEZE)


def latent_heat_sublimation(temperature_k):
   """Latent heat of sublimation (J kg^-1); Murphy and Koop (2005) Eq. 5 (J mol^-1) / M_WATER."""

   t = np.asarray(temperature_k, dtype=np.float64)
   return (46782.5 + 35.8925 * t - 0.07414 * t ** 2 + 541.5 * np.exp(-(t / 123.75) ** 2)) / M_WATER


def saturated_lapse_rate(temperature_k, pressure_pa, phase):
   """Saturated (pseudo-)adiabatic lapse rate -dT/dz (K m^-1) over liquid water or ice."""

   t = float(temperature_k)
   p = float(pressure_pa)
   if phase == "liquid":
      e_s = float(saturation_vapour_pressure_liquid(t))
      latent = float(latent_heat_vaporisation(t))
   elif phase == "ice":
      e_s = float(saturation_vapour_pressure_ice(t))
      latent = float(latent_heat_sublimation(t))
   else:
      raise ValueError(f"phase must be 'liquid' or 'ice', got {phase!r}")
   if e_s >= p:
      raise ValueError(f"saturation vapour pressure {e_s:.1f} Pa reaches the pressure {p:.1f} Pa at {t:.1f} K")
   r_s = EPSILON * e_s / (p - e_s)
   return G0 * (1.0 + latent * r_s / (R_DRY * t)) / (C_P_DRY + latent ** 2 * r_s * EPSILON / (R_DRY * t ** 2))


def pressure_interpolator(height_km, pressure_hpa):
   """Return p(z_km) in hPa: log-linear in height between the profile levels.

   Below the lowest level the lowest layer's log-pressure gradient is continued
   (the production atmospheres end at the surface, 0 km).
   """

   h = np.asarray(height_km, dtype=np.float64)
   lp = np.log(np.asarray(pressure_hpa, dtype=np.float64))
   order = np.argsort(h)
   h, lp = h[order], lp[order]
   if h.size < 2:
      raise ValueError("the atmosphere needs at least two levels")
   slope_bottom = (lp[1] - lp[0]) / (h[1] - h[0])

   def pressure_hpa_at(z_km):
      z = float(z_km)
      if z < h[0]:
         return float(np.exp(lp[0] + slope_bottom * (z - h[0])))
      return float(np.exp(np.interp(z, h, lp)))

   return pressure_hpa_at


def adiabatic_profile(t_top_k, cloud_top_km, height_km, pressure_hpa, phase, h_max_km, step_km=PROFILE_STEP_KM):
   """Integrate the saturated adiabat downward from the cloud top.

   Returns ``(depth_km, temperature_k)`` on 0, step, ..., h_max (the last node
   is exactly h_max), by fourth-order Runge-Kutta in depth with the pressure of
   the atmosphere at height cloud_top_km - depth.
   """

   if h_max_km <= 0.0:
      raise ValueError("h_max_km must be positive")
   pressure_at = pressure_interpolator(height_km, pressure_hpa)
   nsteps = int(np.ceil(h_max_km / step_km - 1e-9))
   depth = np.linspace(0.0, h_max_km, nsteps + 1)
   temperature = np.empty(nsteps + 1, dtype=np.float64)
   temperature[0] = float(t_top_k)

   def rate(z_km, t):                                   # dT/dz with z the depth below the cloud top (m)
      return saturated_lapse_rate(t, 100.0 * pressure_at(cloud_top_km - z_km), phase)

   for i in range(nsteps):
      z0, t0 = depth[i], temperature[i]
      dz_km = depth[i + 1] - z0
      dz_m = 1000.0 * dz_km
      k1 = rate(z0, t0)
      k2 = rate(z0 + 0.5 * dz_km, t0 + 0.5 * dz_m * k1)
      k3 = rate(z0 + 0.5 * dz_km, t0 + 0.5 * dz_m * k2)
      k4 = rate(z0 + dz_km, t0 + dz_m * k3)
      temperature[i + 1] = t0 + dz_m * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
   return depth, temperature


@dataclass(frozen=True)
class CloudTemperatureModel:
   """The V25 adiabatic cloud of one LUT: its constants and the tabulated adiabat T(depth)."""

   substance: str
   phase: str
   t_top_k: float
   beta_ext_055_per_km: float
   h_max_km: float
   cloud_top_km: float
   cloud_top_hpa: float
   depth_km: np.ndarray
   temperature_k: np.ndarray

   def geometric_depth_km(self, tau_055):
      """H = min(tau_055 / beta_ext_055, H_max); 0 for tau_055 = 0."""

      tau = float(tau_055)
      if tau < 0.0:
         raise ValueError(f"the optical depth must not be negative, got {tau}")
      return min(tau / self.beta_ext_055_per_km, self.h_max_km)

   def depth_of_optical_depth_km(self, tau_055, x):
      """z(x) = H x / tau_055 for cumulative 0.55-um optical depth x below the cloud top (0 when tau_055 = 0)."""

      tau = float(tau_055)
      x = np.asarray(x, dtype=np.float64)
      if tau <= 0.0:
         return np.zeros_like(x)
      if np.any(x < 0.0) or np.any(x > tau * (1.0 + 1e-6)):
         raise ValueError("the cumulative optical depth must lie between 0 and tau_055")
      return self.geometric_depth_km(tau) * np.clip(x / tau, 0.0, 1.0)

   def temperature_at_depth_k(self, depth_km):
      """T(z) from the tabulated adiabat (linear interpolation; z is clipped to 0..H_max)."""

      z = np.clip(np.asarray(depth_km, dtype=np.float64), 0.0, self.h_max_km)
      return np.interp(z, self.depth_km, self.temperature_k)

   def boundary_temperatures_k(self, tau_055, fractions):
      """Temperatures at layer boundaries given as fractions 0..1 of the cloud optical depth."""

      f = np.asarray(fractions, dtype=np.float64)
      if np.any(f < 0.0) or np.any(f > 1.0 + 1e-9):
         raise ValueError("fractions must lie between 0 and 1")
      return self.temperature_at_depth_k(self.geometric_depth_km(tau_055) * np.clip(f, 0.0, 1.0))

   def cloud_base_temperature_k(self, tau_055):
      return float(self.temperature_at_depth_k(self.geometric_depth_km(tau_055)))

   def describe(self):
      return (f"adiabatic (V25): T_top {self.t_top_k:.1f} K at {self.cloud_top_km:.2f} km / "
              f"{self.cloud_top_hpa:.1f} hPa; {self.phase}-saturated adiabat; beta_ext(0.55 um) "
              f"{self.beta_ext_055_per_km:g} /km; H_max {self.h_max_km:g} km "
              f"(T {self.cloud_base_temperature_k(np.inf):.1f} K at H_max); {EMISSION_SUBLAYERS} emission sub-layers per cloud layer")

   def global_attributes(self):
      """NetCDF global attributes recording the V25 cloud temperature treatment."""

      return {
         "cloud_temperature_profile": "adiabatic",
         "cloud_top_temperature_K": np.float32(self.t_top_k),
         "cloud_adiabat": f"saturated adiabat over {'liquid water' if self.phase == 'liquid' else 'ice'}",
         "cloud_extinction_coefficient_055um_per_km": np.float32(self.beta_ext_055_per_km),
         "cloud_maximum_geometric_depth_km": np.float32(self.h_max_km),
         "cloud_top_height_km": np.float32(self.cloud_top_km),
         "cloud_top_pressure_hPa": np.float32(self.cloud_top_hpa),
         "cloud_emission_sublayers": np.int32(EMISSION_SUBLAYERS),
         "cloud_temperature_profile_note": (
            "Thermal emission (E_md) of a cloud whose temperature rises from the cloud top "
            "(cloud_top_temperature_K) along the saturated adiabat; geometric depth "
            "H = min(tau_055 / beta, H_max) with the whole optical depth mapped linearly onto 0..H; "
            "E_md is normalised by the Planck radiance at cloud_top_temperature_K and may exceed 1."),
      }


def cloud_temperature_model(substance, atmstr, scatreltau, t_top_k=T_TOP_K):
   """Build the V25 model of a LUT from its .mm substance, atmosphere and particle layer profile.

   The cloud top is the upper level of the first atmosphere layer (top down)
   whose relative particle optical depth is positive.
   """

   key = str(substance).lower()
   if key not in CLOUD_MODELS:
      raise ValueError(f"the V25 adiabatic cloud temperature profile is defined for substances "
                       f"{sorted(CLOUD_MODELS)}, not {substance!r}")
   phase, beta, h_max = CLOUD_MODELS[key]
   layers = np.flatnonzero(np.asarray(scatreltau) > 0.0)
   if layers.size == 0:
      raise ValueError("the particle profile has no layer with particles")
   cloud_top_km = float(atmstr.height[layers[0]])
   cloud_top_hpa = pressure_interpolator(atmstr.height, atmstr.pressure)(cloud_top_km)
   depth, temperature = adiabatic_profile(t_top_k, cloud_top_km, atmstr.height, atmstr.pressure, phase, h_max)
   return CloudTemperatureModel(substance=key, phase=phase, t_top_k=float(t_top_k), beta_ext_055_per_km=beta,
                                h_max_km=h_max, cloud_top_km=cloud_top_km, cloud_top_hpa=cloud_top_hpa,
                                depth_km=depth, temperature_k=temperature)


def emission_layers(dtau, ssalb, pmom, scatreltau, tau_055, model, nsub=EMISSION_SUBLAYERS):
   """Sub-layered emission-call inputs (dtau, ssalb, pmom, temper) for the in-cloud layers.

   ``dtau``, ``ssalb``, ``pmom[:, layer]`` and ``scatreltau`` are the in-cloud
   layers, top down.  Every layer becomes ``nsub`` sub-layers of equal optical
   depth with its own single-scattering albedo and moments; the boundary
   temperatures follow the cumulative particle optical-depth fraction (the same
   at every wavelength) through the model's depth mapping and adiabat.
   """

   f32 = np.float32
   dtau = np.asarray(dtau, dtype=f32)
   ssalb = np.asarray(ssalb, dtype=f32)
   pmom = np.asarray(pmom, dtype=f32)
   weights = np.asarray(scatreltau, dtype=np.float64)
   nsub = int(nsub)
   if nsub < 1:
      raise ValueError("nsub must be at least 1")
   if dtau.ndim != 1 or ssalb.shape != dtau.shape or pmom.shape[1] != dtau.size or weights.shape != dtau.shape:
      raise ValueError("dtau, ssalb, pmom[:, layer] and scatreltau must describe the same layers")
   if dtau.size == 0:
      raise ValueError("no in-cloud layer")
   if np.any(weights <= 0.0):
      raise ValueError("every in-cloud layer must have a positive particle optical depth")
   if dtau.size * nsub > DISORT_MAX_LAYERS:
      raise ValueError(f"{dtau.size} in-cloud layers x {nsub} sub-layers exceed the {DISORT_MAX_LAYERS} layers "
                       "the production DISORT accepts (MXCLY)")

   emtau = np.repeat(dtau / f32(nsub), nsub).astype(f32)
   emssa = np.repeat(ssalb, nsub).astype(f32)
   empmo = np.asfortranarray(np.repeat(pmom, nsub, axis=1).astype(f32))
   # Cumulative particle optical-depth fraction at every sub-layer boundary, 0 .. 1 exactly
   fraction = np.concatenate([[0.0], np.cumsum(np.repeat(weights / weights.sum() / nsub, nsub))])
   fraction[-1] = 1.0
   temper = model.boundary_temperatures_k(tau_055, fraction).astype(f32)
   return emtau, emssa, empmo, temper
