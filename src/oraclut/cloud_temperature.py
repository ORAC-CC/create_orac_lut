"""Vertically inhomogeneous cloud profile for the V25 thermal-emission calculation (create_orac_cloud_lut).

Up to V24 the emission DISORT call gives every in-cloud layer one temperature
(create_orac_cloud_lut.pro line 561, ``temp = 250.0``), so the LUT emissivity
E_md = UU / B(temp) is that of an isothermal cloud and does not depend on the
value of temp.  V25 keeps a fixed cloud-top reference temperature,
T_TOP_K = 240 K, and gives the cloud the vertical structure of a supplied
cloud profile: the distribution of the cloud optical depth with depth below
the cloud top and the temperature departure from the cloud top at each
depth.  E_md stays normalised by B(T_TOP_K), so it carries the differential
emission of the warmer lower cloud.  The LUT grids are unchanged: cloud
temperature is not a LUT dimension and no cloud-base temperature is retrieved.

The profile data (references/data/ocalut_cloudprofile_Cirrostratus.dat;
origin /tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf,
cloud type Cirrostratus; the vertical profiles of P. Watts' vertically
inhomogeneous LUT tests for OCA / EUMETSAT) tabulate, for 13 total cloud
optical depths COT = 0.0625 ... 256 (doubling), 99 rows of

   % extinction  f_i      layer optical depth  dCOT_i = COT f_i / 100
   cumulative optical depth  C_i = sum_{j <= i} dCOT_j
   depth below the cloud top  z_i = H(COT) i / 98   (i = 0 ... 98)
   temperature above the cloud top  dT_i = 8 K/km x z_i.

Established numerically from the file (validation/v25): the extinction shape
f_i is the same for all 13 COTs and sums to 100; the rows are equally spaced
in depth from the cloud top (z = 0) to the cloud base (z = H); dT = 8 z holds
in every row; the only COT-dependent quantity is the cloud depth H(COT), which
rises from 2.0 km (COT 0.0625) to 11.0 km (COT 128 and 256).  The supplied
model is therefore

   dtau_i(COT) = COT f_i,   z_i(COT) = H(COT) i / 98,   dT_i = 8 K/km z_i,

with f_i and H(COT) read from the file.  The cloud does not become 256 km
deep at COT 256; it stays 11 km deep while the optical depth per unit depth
grows.

Interpretation of the 99 rows (P. Watts' DISORT example was not available):
row i is a sample at depth z_i carrying the optical depth dCOT_i, so the
optical depth accumulated at z_i is taken as C_i - dCOT_i / 2 (the midpoint
rule; the cloud top has C = 0 at z = 0 and the base C = COT at z = H).  This
defines a COT-independent relation between the cumulative optical-depth
fraction F and the depth fraction zeta = z / H, and hence
dT(F; COT) = 8 K/km H(COT) zeta(F).

Interpolation to the ORAC optical-depth grid: H(COT) is interpolated
linearly in COT between the 13 native values (the file is exactly linear in
COT from 0.0625 to 1 and nearly so to 8; linear-in-log(COT) interpolation
differs by at most 0.23 km, 1.9 K at the base, at the intermediate Grid B
nodes), held at H(0.0625) = 2 km below COT 0.0625 (no extrapolation; E_md
tends to zero with the optical depth), and not extrapolated above COT 256
(the Grid B maximum).  f_i needs no interpolation.  Every native profile is
reproduced exactly at its own COT.

DISORT representation: the production kernel accepts MXCLY = 46 computational
layers, so the 99-row profile is remapped conservatively in cumulative
optical-depth space onto EMISSION_LAYERS equal-optical-depth sub-layers of
each in-cloud atmosphere layer, with the boundary temperatures
T_TOP_K + dT(F_k; COT) at the sub-layer boundaries F_k = k / N.  The total
optical depth and the cumulative distribution at every boundary are exact;
the thermal source is linear in optical depth within each sub-layer (DISORT's
own assumption).  The single-scattering albedo and phase moments of the LUT
particle model are the same in every sub-layer (the profile prescribes
extinction and temperature structure, not microphysics), so the diffuse and
direct-beam DISORT calls, which see one homogeneous cloud layer, are
unchanged: in plane-parallel radiative transfer the geometric distribution of
extinction within a homogeneous layer does not enter the reflection and
transmission operators; the profile matters only through the pairing of
optical depth with temperature in the emission call.
"""

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# V25 constants
# ---------------------------------------------------------------------------
T_TOP_K = 240.0                          # cloud-top reference temperature (the file gives departures from it)
LAPSE_K_PER_KM = 8.0                     # dT = 8 K/km z in the supplied profiles (verified, not assumed)
PROFILE_FILES = {
   "cirrostratus": Path(__file__).resolve().parents[2] / "references" / "data" / "ocalut_cloudprofile_Cirrostratus.dat",
}
# Sub-layers of every in-cloud atmosphere layer in the emission DISORT call
# (chosen from the layer-count convergence of validation/v25).
EMISSION_LAYERS = 40
# DISORT 2.0 as compiled for production accepts at most MXCLY = 46 computational
# layers (create_orac_lut/disort2/src/DISORTfunctions.f).
DISORT_MAX_LAYERS = 46
NATIVE_ROWS = 99


@dataclass(frozen=True)
class CloudProfileSet:
   """The validated native profiles of one supplied cloud-profile file."""

   cloud_type: str
   origin: str
   path: str
   md5: str
   cot: np.ndarray                 # (13,) total cloud optical depths
   extinction_percent: np.ndarray  # (99,) the common extinction shape, sums to 100
   dcot: np.ndarray                # (13, 99) layer optical depths as tabulated
   cumulative: np.ndarray          # (13, 99) cumulative optical depth as tabulated
   depth_km: np.ndarray            # (13, 99) depth below the cloud top
   delta_t_k: np.ndarray           # (13, 99) temperature above the cloud top

   @property
   def depth_scale_km(self):
      """H(COT): the cloud depth of each native profile."""
      return self.depth_km[:, -1]

   @property
   def optical_fraction_midpoints(self):
      """Cumulative optical-depth fraction at each row (midpoint rule), common to all COTs."""
      f = self.extinction_percent / 100.0
      return np.cumsum(f) - 0.5 * f

   def depth_fraction(self, optical_fraction):
      """zeta(F): depth fraction z/H at cumulative optical-depth fraction F, 0 <= F <= 1."""
      F = np.asarray(optical_fraction, dtype=np.float64)
      if np.any(F < 0.0) or np.any(F > 1.0 + 1e-9):
         raise ValueError("the cumulative optical-depth fraction must lie between 0 and 1")
      zeta = np.arange(NATIVE_ROWS) / (NATIVE_ROWS - 1.0)
      return np.interp(np.clip(F, 0.0, 1.0), np.r_[0.0, self.optical_fraction_midpoints, 1.0], np.r_[0.0, zeta, 1.0])

   def depth_at(self, cot):
      """H(COT), interpolated linearly in COT; held at the smallest native value below it."""
      cot = float(cot)
      if cot < 0.0:
         raise ValueError(f"the optical depth must not be negative, got {cot}")
      if cot > self.cot[-1] * (1.0 + 1e-9):
         raise ValueError(f"the optical depth {cot} exceeds the largest supplied profile ({self.cot[-1]}); "
                          "the profiles are not extrapolated")
      return float(np.interp(cot, self.cot, self.depth_scale_km))

   def delta_t_at(self, cot, optical_fraction):
      """Temperature above the cloud top at cumulative optical-depth fraction F of a cloud of optical depth COT."""
      return LAPSE_K_PER_KM * self.depth_at(cot) * self.depth_fraction(optical_fraction)

   def layer_boundaries(self, cot, nlayers):
      """(fractions, dT) at the N + 1 boundaries of N equal-optical-depth layers."""
      n = int(nlayers)
      if n < 1:
         raise ValueError("nlayers must be at least 1")
      fractions = np.arange(n + 1) / n
      return fractions, self.delta_t_at(cot, fractions)

   def native_layer_optical_depths(self, cot):
      """dCOT_i = COT f_i / 100 for an arbitrary COT (the tabulated values at a native COT, to their rounding)."""
      return float(cot) * self.extinction_percent / 100.0


def load_cloud_profiles(path):
   """Parse and validate a supplied cloud-profile file (format of ocalut_cloudprofile_*.dat)."""

   path = Path(path)
   text = path.read_text()
   lines = text.splitlines()
   md5 = hashlib.md5(text.encode()).hexdigest()
   head = re.match(r"\s*(\d+)\s+# COTs", lines[0])
   rows_line = re.match(r"\s*(\d+)\s+# layers", lines[1])
   if not head or not rows_line:
      raise ValueError(f"{path}: unexpected header {lines[0]!r} / {lines[1]!r}")
   ncot, nrows = int(head.group(1)), int(rows_line.group(1))
   if nrows != NATIVE_ROWS:
      raise ValueError(f"{path}: {nrows} rows per profile, expected {NATIVE_ROWS}")
   origin = lines[2].split(":", 1)[1].strip() if lines[2].startswith("Origin:") else ""
   cloud_type = lines[3].split(":", 1)[1].strip() if lines[3].startswith("Cloud type:") else ""
   cots, blocks = [], []
   i = 4
   while i < len(lines):
      if lines[i].strip() == "":
         i += 1
         continue
      m = re.match(r"\s*([0-9.Ee+-]+)\s+LUT Cot value", lines[i])
      if not m or i + 1 >= len(lines) or "% Extinction" not in lines[i + 1]:
         raise ValueError(f"{path}: malformed profile header at line {i + 1}: {lines[i]!r}")
      if i + 2 + nrows > len(lines):
         raise ValueError(f"{path}: the profile at line {i + 1} is truncated")
      rows = []
      for j in range(i + 2, i + 2 + nrows):
         values = lines[j].split()
         if len(values) != 5:
            raise ValueError(f"{path}: line {j + 1} has {len(values)} values, expected 5")
         rows.append([float(v) for v in values])
      cots.append(float(m.group(1)))
      blocks.append(np.asarray(rows, dtype=np.float64))
      i += 2 + nrows
   if len(cots) != ncot:
      raise ValueError(f"{path}: {len(cots)} profiles read, header says {ncot}")
   cot = np.asarray(cots)
   data = np.stack(blocks)                                   # (ncot, 99, 5)
   ext, dcot, cum, z, dt = (data[:, :, k] for k in range(5))
   # -- validation of the supplied data, every profile --
   if not np.all(np.isfinite(data)):
      raise ValueError(f"{path}: non-finite values")
   if np.any(np.diff(cot) <= 0.0):
      raise ValueError(f"{path}: the COTs are not increasing")
   if np.abs(ext - ext[0]).max() > 0.0:
      raise ValueError(f"{path}: the extinction shape differs between profiles")
   if abs(ext[0].sum() - 100.0) > 1e-3 or np.any(ext < 0.0):
      raise ValueError(f"{path}: the extinction percentages do not sum to 100")
   for k in range(len(cot)):
      if np.any(dcot[k] < 0.0) or abs(dcot[k].sum() - cot[k]) > 1e-4 * cot[k]:
         raise ValueError(f"{path}: layer dCOT of COT {cot[k]} do not sum to the total")
      if np.abs(cum[k] - np.cumsum(dcot[k])).max() > 1e-5 * max(1.0, cot[k]) or np.any(np.diff(cum[k]) <= 0.0):
         raise ValueError(f"{path}: cumulative COT of profile {cot[k]} is not the running sum or not monotonic")
      if abs(cum[k, -1] - cot[k]) > 1e-4 * cot[k]:
         raise ValueError(f"{path}: cumulative COT of profile {cot[k]} does not end at the total")
      if z[k, 0] != 0.0 or np.any(np.diff(z[k]) <= 0.0) or np.any(np.diff(dt[k]) <= 0.0) or dt[k, 0] != 0.0:
         raise ValueError(f"{path}: depth or temperature of profile {cot[k]} is not monotonic from 0 at the top")
      if np.abs(dt[k] - LAPSE_K_PER_KM * z[k]).max() > 1e-4:
         raise ValueError(f"{path}: profile {cot[k]} does not follow dT = {LAPSE_K_PER_KM} K/km z")
      if np.abs(z[k] / z[k, -1] - np.arange(nrows) / (nrows - 1.0)).max() > 1e-5:
         raise ValueError(f"{path}: the rows of profile {cot[k]} are not equally spaced in depth")
   if np.any(np.diff(z[:, -1]) < 0.0):
      raise ValueError(f"{path}: the cloud depth does not increase with COT")
   return CloudProfileSet(cloud_type=cloud_type, origin=origin, path=str(path), md5=md5, cot=cot,
                          extinction_percent=ext[0].copy(), dcot=dcot, cumulative=cum, depth_km=z, delta_t_k=dt)


@dataclass(frozen=True)
class CloudVerticalProfile:
   """The V25 vertical cloud model of one LUT: a supplied profile set at the reference cloud-top temperature."""

   name: str
   profiles: CloudProfileSet
   t_top_k: float
   max_tau_055: float

   def boundary_temperatures_k(self, tau_055, fractions):
      """T_top + dT at layer boundaries given as cumulative optical-depth fractions 0..1."""
      tau = float(tau_055)
      if tau < 0.0:
         raise ValueError(f"the optical depth must not be negative, got {tau}")
      return self.t_top_k + self.profiles.delta_t_at(tau, fractions)

   def cloud_base_temperature_k(self, tau_055):
      return float(self.boundary_temperatures_k(tau_055, [1.0])[0])

   def describe(self):
      p = self.profiles
      return (f"{self.name} (V25): supplied vertically inhomogeneous profile ({p.cloud_type}, {Path(p.path).name}, "
              f"md5 {p.md5[:8]}; {p.cot.size} native COTs {p.cot[0]:g}-{p.cot[-1]:g}, {NATIVE_ROWS} rows); "
              f"T_top {self.t_top_k:.1f} K, dT = {LAPSE_K_PER_KM:g} K/km x depth, cloud depth {p.depth_scale_km[0]:g}-"
              f"{p.depth_scale_km[-1]:g} km (linear in COT); {EMISSION_LAYERS} emission layers per cloud layer")

   def global_attributes(self):
      """NetCDF global attributes recording the V25 cloud profile treatment."""
      p = self.profiles
      return {
         "cloud_vertical_profile": self.name,
         "cloud_vertical_profile_file": Path(p.path).name,
         "cloud_vertical_profile_md5": p.md5,
         "cloud_vertical_profile_origin": p.origin,
         "cloud_vertical_profile_type": p.cloud_type,
         "cloud_top_temperature_K": np.float32(self.t_top_k),
         "cloud_emission_layers": np.int32(EMISSION_LAYERS),
         "cloud_vertical_profile_note": (
            "Thermal emission (E_md) of a cloud with the supplied vertical distribution of optical depth and "
            "temperature departure from the cloud top (dT = 8 K/km x depth; depth H(COT) interpolated linearly "
            "in COT between the supplied profiles, 2-11 km), at the reference cloud-top temperature "
            "cloud_top_temperature_K; remapped onto cloud_emission_layers equal-optical-depth layers; "
            "E_md is normalised by the Planck radiance at cloud_top_temperature_K and may exceed 1; "
            "reflection and transmission operators are those of the homogeneous cloud (V24)."),
      }


def cloud_vertical_profile_model(name, max_tau_055, t_top_k=T_TOP_K, path=None):
   """Build the V25 model of a LUT from the named supplied profile and the largest optical depth of its grid."""

   key = str(name).lower()
   if key not in PROFILE_FILES:
      raise ValueError(f"unknown cloud vertical profile {name!r}; available: {sorted(PROFILE_FILES)}")
   profiles = load_cloud_profiles(PROFILE_FILES[key] if path is None else path)
   if float(max_tau_055) > profiles.cot[-1] * (1.0 + 1e-9):
      raise ValueError(f"the LUT optical depth {max_tau_055} exceeds the largest supplied profile ({profiles.cot[-1]})")
   return CloudVerticalProfile(name=key, profiles=profiles, t_top_k=float(t_top_k), max_tau_055=float(max_tau_055))


def emission_layers(dtau, ssalb, pmom, scatreltau, tau_055, model, nlayers=EMISSION_LAYERS):
   """Sub-layered emission-call inputs (dtau, ssalb, pmom, temper) for the in-cloud layers.

   ``dtau``, ``ssalb``, ``pmom[:, layer]`` and ``scatreltau`` are the in-cloud
   atmosphere layers, top down.  Every layer becomes ``nlayers`` sub-layers of
   equal optical depth with its own single-scattering albedo and moments; the
   boundary temperatures follow the cumulative particle optical-depth fraction
   of the node's tau_055 (the same at every wavelength) through the supplied
   profile.
   """

   f32 = np.float32
   dtau = np.asarray(dtau, dtype=f32)
   ssalb = np.asarray(ssalb, dtype=f32)
   pmom = np.asarray(pmom, dtype=f32)
   weights = np.asarray(scatreltau, dtype=np.float64)
   n = int(nlayers)
   if n < 1:
      raise ValueError("nlayers must be at least 1")
   if dtau.ndim != 1 or ssalb.shape != dtau.shape or pmom.shape[1] != dtau.size or weights.shape != dtau.shape:
      raise ValueError("dtau, ssalb, pmom[:, layer] and scatreltau must describe the same layers")
   if dtau.size == 0:
      raise ValueError("no in-cloud layer")
   if np.any(weights <= 0.0):
      raise ValueError("every in-cloud layer must have a positive particle optical depth")
   if dtau.size * n > DISORT_MAX_LAYERS:
      raise ValueError(f"{dtau.size} in-cloud layers x {n} sub-layers exceed the {DISORT_MAX_LAYERS} layers "
                       "the production DISORT accepts (MXCLY)")

   emtau = np.repeat(dtau / f32(n), n).astype(f32)
   emssa = np.repeat(ssalb, n).astype(f32)
   empmo = np.asfortranarray(np.repeat(pmom, n, axis=1).astype(f32))
   # Cumulative particle optical-depth fraction at every sub-layer boundary, 0 .. 1 exactly
   fraction = np.concatenate([[0.0], np.cumsum(np.repeat(weights / weights.sum() / n, n))])
   fraction[-1] = 1.0
   temper = model.boundary_temperatures_k(tau_055, fraction).astype(f32)
   return emtau, emssa, empmo, temper
