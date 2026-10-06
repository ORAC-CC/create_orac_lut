"""Cloud temperature profile for the V25 thermal-emission calculation (create_orac_cloud_lut).

Up to V24 the emission DISORT call gives every in-cloud layer one temperature
(create_orac_cloud_lut.pro line 561, ``temp = 250.0``), so the LUT emissivity
E_md = UU / B(temp) is that of an isothermal cloud and does not depend on the
value of temp.  V25 keeps a fixed cloud-top reference temperature,
T_TOP_K = 240 K, and lets the temperature increase downward through the cloud
along the saturated adiabat of the LUT's substance.  E_md stays normalised by
B(T_TOP_K), so it carries the differential emission of the warmer lower
cloud.  The LUT grids are unchanged: cloud temperature is not a LUT dimension
and no cloud-base temperature is retrieved.

The V25 model (its constants are representative, deliberate choices, not
microphysical results):

* Reference cloud-top state T_top = T_TOP_K, p_top = P_TOP_HPA, the same for
  liquid water and ice.
* Thermodynamic path-length coordinate: a cumulative 0.55-um optical depth x
  below the cloud top lies at the distance

      z(x) = x / beta_ext_055,   i.e.  dz = d(tau_055) / beta_ext_055,

  along the reference adiabat, with the representative volume extinction
  coefficient beta_ext_055 of the substance (LIQUID_WATER: 20 km^-1;
  WATER_ICE: 1 km^-1).  z is not limited, is not an absolute altitude and is
  not a LUT variable; the LUT's cloud position in the atmosphere does not
  enter.  tau_055 is the LUT optical-depth coordinate (lutstr.opd): the
  particle optical depth of every channel is tau_055 * BextRat,
  BextRat = Bext / Bext(0.55 um), so one thermodynamic structure serves every
  spectral channel of a tau_055 node.
* T(z) is the saturated adiabat of the substance, integrated downward from
  the reference state.  The lapse rate needs the pressure, which is evolved
  hydrostatically along the same path, dp/dz = p g / (R_d T) (dry air), from
  p_top; no atmospheric profile is consulted.
* The LUT is phase-pure: a liquid-water LUT follows the liquid-water adiabat
  throughout and an ice LUT the ice adiabat throughout, whatever temperature
  is reached.  There is no melting, no freezing-point ceiling and no phase
  switch; extreme reference states at large optical depth are a limitation of
  this simplified parameterisation and are reported, not limited.

beta_ext_055 is prescribed.  It is not derived from the per-particle
extinction cross-section or the normalised size distribution of
generate_scattering_properties: the LUT fixes no number concentration, so the
microphysics does not define a physical thickness.

Thermodynamics.  The repository had no saturated-adiabat routine, so this
module holds a compact one:

* saturation vapour pressure over liquid water (valid for supercooled water)
  and over ice, and the latent heat of sublimation: Murphy and Koop (2005,
  Q. J. R. Meteorol. Soc. 131, 1539-1565), Eqs. 10, 7 and 5;
* latent heat of vaporisation L_v = 2.501e6 - 2370 (T - T_FREEZE) J kg^-1
  (Rogers and Yau 1989);
* pseudo-adiabatic (saturated) lapse rate
      Gamma_s = g (1 + L r_s / (R_d T)) / (c_pd + L^2 r_s eps / (R_d T^2)),
      r_s = eps e_s / (p - e_s)
  (AMS Glossary of Meteorology, "moist-adiabatic lapse rate"), with e_sw and
  L_v for liquid-water clouds and e_si and L_s for ice clouds.

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
# Reference pressure of the cloud-top state (hPa).  The pressure enters the
# saturated lapse rate only through the saturation mixing ratio; this value is
# the mid-latitude-summer pressure at 4 km, the state from which the first V25
# formulation started, kept so that the reference adiabat is unchanged there.
P_TOP_HPA = 628.0

# substance (.mm "substance") -> (phase, beta_ext_055 km^-1)
CLOUD_MODELS = {
   "liquid-water": ("liquid", 20.0),
   "water-ice": ("ice", 1.0),
}

# Sub-layers of every in-cloud layer in the emission DISORT call.
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
T_FREEZE = 273.15                        # reference temperature of the latent-heat fit only


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
   if not np.isfinite(e_s) or e_s >= p:
      raise ValueError(f"saturation vapour pressure {e_s:.4g} Pa reaches the pressure {p:.4g} Pa at {t:.1f} K: "
                       "the saturated adiabat is no longer defined")
   r_s = EPSILON * e_s / (p - e_s)
   return G0 * (1.0 + latent * r_s / (R_DRY * t)) / (C_P_DRY + latent ** 2 * r_s * EPSILON / (R_DRY * t ** 2))


def adiabatic_profile(t_top_k, p_top_hpa, phase, z_max_km, step_km=PROFILE_STEP_KM):
   """Integrate the saturated adiabat downward from the reference cloud-top state.

   Returns ``(depth_km, temperature_k, pressure_hpa)`` on 0, step, ..., z_max
   (the last node is exactly z_max), by fourth-order Runge-Kutta in the path
   length z of the coupled system
      dT/dz = Gamma_s(T, p),   dp/dz = p g / (R_d T)   (hydrostatic, dry air).
   No atmospheric profile enters: z is the distance below the reference cloud
   top along the adiabat, not an altitude.
   """

   if z_max_km < 0.0:
      raise ValueError("z_max_km must not be negative")
   nsteps = max(1, int(np.ceil(z_max_km / step_km - 1e-9)))
   depth = np.linspace(0.0, max(z_max_km, step_km), nsteps + 1)
   temperature = np.empty(nsteps + 1, dtype=np.float64)
   pressure = np.empty(nsteps + 1, dtype=np.float64)
   temperature[0] = float(t_top_k)
   pressure[0] = 100.0 * float(p_top_hpa)

   def rates(t, p):                                     # (dT/dz, dp/dz) per metre
      return saturated_lapse_rate(t, p, phase), p * G0 / (R_DRY * t)

   for i in range(nsteps):
      t0, p0 = temperature[i], pressure[i]
      h = 1000.0 * (depth[i + 1] - depth[i])
      k1t, k1p = rates(t0, p0)
      k2t, k2p = rates(t0 + 0.5 * h * k1t, p0 + 0.5 * h * k1p)
      k3t, k3p = rates(t0 + 0.5 * h * k2t, p0 + 0.5 * h * k2p)
      k4t, k4p = rates(t0 + h * k3t, p0 + h * k3p)
      temperature[i + 1] = t0 + h * (k1t + 2.0 * k2t + 2.0 * k3t + k4t) / 6.0
      pressure[i + 1] = p0 + h * (k1p + 2.0 * k2p + 2.0 * k3p + k4p) / 6.0
      if not (np.isfinite(temperature[i + 1]) and np.isfinite(pressure[i + 1])):
         raise ValueError(f"the {phase} adiabat is not finite at {depth[i + 1]:.2f} km below the cloud top")
   return depth, temperature, pressure / 100.0


@dataclass(frozen=True)
class CloudTemperatureModel:
   """The V25 adiabatic cloud of one LUT: its constants and the tabulated reference adiabat."""

   substance: str
   phase: str
   t_top_k: float
   p_top_hpa: float
   beta_ext_055_per_km: float
   depth_km: np.ndarray          # path length below the cloud top, 0 .. the deepest node of the LUT
   temperature_k: np.ndarray
   pressure_hpa: np.ndarray

   def depth_of_optical_depth_km(self, x):
      """z(x) = x / beta_ext_055 for cumulative 0.55-um optical depth x below the cloud top."""

      x = np.asarray(x, dtype=np.float64)
      if np.any(x < 0.0):
         raise ValueError("the cumulative optical depth must not be negative")
      return x / self.beta_ext_055_per_km

   def temperature_at_depth_k(self, depth_km):
      """T(z) from the tabulated adiabat (linear interpolation; z must lie within the table)."""

      z = np.asarray(depth_km, dtype=np.float64)
      if np.any(z < 0.0) or np.any(z > self.depth_km[-1] * (1.0 + 1e-9)):
         raise ValueError(f"the reference adiabat is tabulated to {self.depth_km[-1]:.3f} km; requested "
                          f"{float(np.max(z)):.3f} km")
      return np.interp(z, self.depth_km, self.temperature_k)

   def pressure_at_depth_hpa(self, depth_km):
      z = np.asarray(depth_km, dtype=np.float64)
      return np.interp(z, self.depth_km, self.pressure_hpa)

   def boundary_temperatures_k(self, tau_055, fractions):
      """Temperatures at layer boundaries given as fractions 0..1 of the node's optical depth tau_055."""

      tau = float(tau_055)
      if tau < 0.0:
         raise ValueError(f"the optical depth must not be negative, got {tau}")
      f = np.asarray(fractions, dtype=np.float64)
      if np.any(f < 0.0) or np.any(f > 1.0 + 1e-9):
         raise ValueError("fractions must lie between 0 and 1")
      return self.temperature_at_depth_k(self.depth_of_optical_depth_km(np.clip(f, 0.0, 1.0) * tau))

   def cloud_base_temperature_k(self, tau_055):
      return float(self.temperature_at_depth_k(self.depth_of_optical_depth_km(tau_055)))

   def describe(self):
      z = float(self.depth_km[-1])
      return (f"adiabatic (V25): T_top {self.t_top_k:.1f} K, p_top {self.p_top_hpa:.1f} hPa (reference state); "
              f"{self.phase}-saturated adiabat throughout; z = tau_055 / {self.beta_ext_055_per_km:g} km "
              f"(tabulated to {z:g} km: {float(self.temperature_k[-1]):.1f} K, {float(self.pressure_hpa[-1]):.4g} hPa); "
              f"{EMISSION_SUBLAYERS} emission sub-layers per cloud layer")

   def global_attributes(self):
      """NetCDF global attributes recording the V25 cloud temperature treatment."""

      return {
         "cloud_temperature_profile": "adiabatic",
         "cloud_top_temperature_K": np.float32(self.t_top_k),
         "cloud_top_reference_pressure_hPa": np.float32(self.p_top_hpa),
         "cloud_adiabat": f"saturated adiabat over {'liquid water' if self.phase == 'liquid' else 'ice'}",
         "cloud_extinction_coefficient_055um_per_km": np.float32(self.beta_ext_055_per_km),
         "cloud_emission_sublayers": np.int32(EMISSION_SUBLAYERS),
         "cloud_temperature_profile_note": (
            "Thermal emission (E_md) of a cloud whose temperature rises below the cloud top "
            "(cloud_top_temperature_K, cloud_top_reference_pressure_hPa) along the phase-pure saturated "
            "adiabat, with the path length below the cloud top z = tau_055 / cloud_extinction_coefficient_055um_per_km "
            "(not limited; not an altitude; pressure evolved hydrostatically along the adiabat); "
            "E_md is normalised by the Planck radiance at cloud_top_temperature_K and may exceed 1."),
      }


def cloud_temperature_model(substance, max_tau_055, t_top_k=T_TOP_K, p_top_hpa=P_TOP_HPA):
   """Build the V25 model of a LUT from its .mm substance and the largest optical depth of its grid.

   The reference adiabat is tabulated to z = max_tau_055 / beta_ext_055, the
   deepest point any emission call of the LUT evaluates.
   """

   key = str(substance).lower()
   if key not in CLOUD_MODELS:
      raise ValueError(f"the V25 adiabatic cloud temperature profile is defined for substances "
                       f"{sorted(CLOUD_MODELS)}, not {substance!r}")
   phase, beta = CLOUD_MODELS[key]
   if float(max_tau_055) < 0.0:
      raise ValueError("max_tau_055 must not be negative")
   depth, temperature, pressure = adiabatic_profile(t_top_k, p_top_hpa, phase, float(max_tau_055) / beta)
   return CloudTemperatureModel(substance=key, phase=phase, t_top_k=float(t_top_k), p_top_hpa=float(p_top_hpa),
                                beta_ext_055_per_km=beta, depth_km=depth, temperature_k=temperature,
                                pressure_hpa=pressure)


def emission_layers(dtau, ssalb, pmom, scatreltau, tau_055, model, nsub=EMISSION_SUBLAYERS):
   """Sub-layered emission-call inputs (dtau, ssalb, pmom, temper) for the in-cloud layers.

   ``dtau``, ``ssalb``, ``pmom[:, layer]`` and ``scatreltau`` are the in-cloud
   layers, top down.  Every layer becomes ``nsub`` sub-layers of equal optical
   depth with its own single-scattering albedo and moments; the boundary
   temperatures follow the cumulative particle optical-depth fraction (the same
   at every wavelength) of the node's tau_055 through the model's path-length
   coordinate and adiabat.
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
