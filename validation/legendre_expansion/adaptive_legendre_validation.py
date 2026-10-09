"""Coefficient-level validation of the adaptive Legendre expansion (Phase 3).

Validation-only code.  For representative liquid-water (production
3.5 x effective-radius integration), ice-sphere (0.001-100 um) and aerosol
(log-normal classes, create_range mixing) cases it runs the production
per-wavelength routine generate_scattering_properties._adaptive_mie_wavelength

   production  at the rule-based Gauss-Legendre order Nq
   reference   at a higher order (2 Nq, or 1.5 Nq where the Mie kernel's
               10000-angle limit requires), giving both the rule's L at that
               order (stability) and every exact coefficient up to each phase
               function's degree bound D_r (the deliberately long reference)

and compares, in double precision:
   * the production series of L terms with the directly calculated averaged
     phase function at ~1400 dense independent angles (captured from the
     reference run: 0-1 deg every 0.005, 1-30 every 0.05, 30-180 every 0.25),
     in 0-1, 1-30 and 30-180 deg;
   * low-order (l <= 60) and high-order (61 <= l < L) chi_l = omega_l/(2l+1)
     against the long reference, and the largest reference coefficient
     dropped beyond L;
   * omega_1/3 against the Mie asymmetry parameter;
   * L at the higher order.

   PYTHONPATH=src python validation/legendre_expansion/adaptive_legendre_validation.py
"""

import csv
import importlib
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"

gsp = importlib.import_module("oraclut.idl_mirror.generate_scattering_properties")
le = importlib.import_module("oraclut.idl_mirror.legendre_expansion")
from oraclut.idl_mirror.load_mmdat import load_mmdat   # noqa: E402

DENSE_THETA = np.unique(np.concatenate((np.arange(0.0, 1.0, 0.005), np.arange(1.0, 30.0, 0.05), np.arange(30.0, 180.0001, 0.25))))
# MODIS 0.412 / 0.646 / 2.114 / 3.785 / 11.03 um, SLSTR 0.555 um (float32 centre wavelengths of the .inst files)
LIQUID = ("liquid-water_stg.mm", [5.0, 10.0, 20.0, 40.0], [0.4118525, 0.5548, 0.6462791, 2.1142, 3.785, 11.0262])
ICE = ("water-ice_sph.mm", [41.0, 93.0], [0.4118525, 0.6462791])
AEROSOL = ("aerosol_a76.mm", [0.1, 1.0, 3.3598182201, 10.0], [0.55, 0.6382, 10.7963])


def class_inputs(mmfile, radii):
   mmstr = load_mmdat(INPUTS / "microphysics" / mmfile, INPUTS)
   radii = np.asarray(radii, dtype=np.float64)
   if mmstr.distname[0] == "log_normal":
      lut_mrat, lut_rm = gsp.create_range(mmstr.mrat, mmstr.rm, mmstr.s, radii)
   else:
      lut_mrat = np.ones((1, radii.size))
      lut_rm = radii[None, :].copy()
   return mmstr, lut_mrat, lut_rm


def run_case(mmfile, radii, wavelength):
   mmstr, lut_mrat, lut_rm = class_inputs(mmfile, radii)
   ri = np.array([gsp._interpol_complex(c.cm, c.wl, [wavelength])[0] for c in mmstr.comp])
   production = gsp._adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, ri, wavelength, f"{mmfile} {wavelength}")
   nq0 = production[4].shape[0]
   factor = 2.0 if 2 * nq0 + DENSE_THETA.size <= le.MIE_MAX_ANGLES else 1.5
   if int(factor * nq0) + DENSE_THETA.size > le.MIE_MAX_ANGLES:
      raise SystemExit(f"{mmfile} {wavelength}: no reference order fits the Mie kernel limit")
   captured = []
   original_order, original_length, original_theta = gsp.initial_quadrature_order, gsp.expansion_length, gsp.CHECK_THETA

   def capture(omega, degree, check_mu, check_phase):
      captured.append((int(degree), np.array(check_phase)))
      return original_length(omega, degree, check_mu, check_phase)

   gsp.initial_quadrature_order = lambda degree: int(factor * original_order(degree))
   gsp.expansion_length = capture
   gsp.CHECK_THETA = DENSE_THETA
   try:
      reference = gsp._adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, ri, wavelength, f"{mmfile} {wavelength} reference")
   finally:
      gsp.initial_quadrature_order, gsp.expansion_length, gsp.CHECK_THETA = original_order, original_length, original_theta
   if len(captured) != len(radii):
      raise SystemExit(f"{mmfile} {wavelength}: the reference run repeated a radius ({len(captured)} passes for "
                       f"{len(radii)} radii); its dense phase functions cannot be matched to radii")
   mu = np.cos(np.deg2rad(DENSE_THETA))
   regions = {"0_1": DENSE_THETA < 1.0, "1_30": (DENSE_THETA >= 1.0) & (DENSE_THETA < 30.0), "30_180": DENSE_THETA >= 30.0}
   rows = []
   for r, re in enumerate(radii):
      degree, direct = captured[r]
      omega = production[5][:, r]
      omega_ref = reference[5][:degree + 1, r]
      length = int(production[6][r])
      series = le.reconstruct_phase_function(omega, length, mu)
      rel = np.abs(series / direct - 1.0)
      chi = omega[:length] / (2.0 * np.arange(length) + 1.0)
      chi_ref = omega_ref / (2.0 * np.arange(degree + 1) + 1.0)
      mratbextw = lut_mrat[:, r] * production[0][:, r] * production[1][:, r]
      g_class = float(np.sum(mratbextw * production[2][:, r]) / np.sum(mratbextw))
      row = {"microphysics": mmfile, "wavelength_um": wavelength, "effective_radius_um": re, "degree_bound": degree,
             "nq": nq0, "L": length, "nq_reference": reference[4].shape[0], "L_at_reference_nq": int(reference[6][r]),
             "noise": le.rounding_noise(omega, degree),
             "tail_max_abs_omega": float(np.max(np.abs(omega[length:degree + 1]))) if length <= degree else 0.0,
             "recon_max_rel": float(rel.max()), "recon_max_abs_over_p0": float(np.max(np.abs(series - direct)) / direct[0]),
             "max_dchi_0_60": float(np.max(np.abs(chi[:min(61, length)] - chi_ref[:min(61, length)]))),
             "max_dchi_61_L": float(np.max(np.abs(chi[61:] - chi_ref[61:length]))) if length > 61 else 0.0,
             "max_ref_chi_dropped": float(np.max(np.abs(chi_ref[length:]))) if length <= degree else 0.0,
             "omega1_over_3_minus_g": float(omega[1] / 3.0 - g_class)}
      for name, mask in regions.items():
         row[f"recon_max_rel_{name}"] = float(rel[mask].max())
      rows.append(row)
   return rows


def main():
   rows = []
   for mmfile, radii, wavelengths in (AEROSOL, LIQUID, ICE):
      for wavelength in wavelengths:
         started = time.time()
         rows += run_case(mmfile, radii, wavelength)
         print(f"{mmfile} {wavelength}: L = {[row['L'] for row in rows[-len(radii):]]}, {time.time() - started:.0f} s", flush=True)
   out = HERE / "results" / "adaptive_legendre_validation.csv"
   out.parent.mkdir(parents=True, exist_ok=True)
   with open(out, "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)
   print("written", out, "; all finite:", all(np.isfinite(v) for row in rows for v in row.values() if isinstance(v, float)))


if __name__ == "__main__":
   main()
