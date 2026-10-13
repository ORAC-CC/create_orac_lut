"""Convergence of the V24 radius integration and adaptive Legendre expansion (working-tree code, read-only).

   python validation/v24/scattering_refinements/convergence.py <case> <out.json>

<case> is "<model>:<wavelength um>:<effective radius um>", e.g.
liquid-water_240:0.646:21 (MODIS SRF centre wavelengths are used by the caller).
Every quantity is calculated by the production routines
(generate_scattering_properties.radius_upper_factor / refined_xres,
create_bwgp.create_bwgp, legendre_expansion) and compared with:

  radius integration   the production grid against the same limits with two
                       further nested refinement levels (refined_xres / 4),
                       and the legacy V23 grid (0.001-100 um, step 0.4)
                       against that reference;
  integration limit    3.5 r_e against 5 r_e on the production refinement
                       (liquid water only);
  Legendre expansion   the production order Nq against 2 Nq: King's length L,
                       the moments chi_l for l < L, the coefficients dropped
                       beyond L, the reconstruction error, |omega_0 - 1|,
                       |omega_1 / 3 - g| and max |chi_l|.

Bulk quantities are compared as relative extinction and absolute w and g, the
phase function at the CHECK_THETA angles as relative differences in 0-5,
5-150 and 150-180 degrees.
"""

import json
import sys
import time

import numpy as np

import codebase

codebase.select("v24")
mirror = codebase.module("oraclut.idl_mirror")
gsp = codebase.module("oraclut.idl_mirror.generate_scattering_properties")
cb = codebase.module("oraclut.idl_mirror.create_bwgp")
le = codebase.module("oraclut.idl_mirror.legendre_expansion")

model, wavelength, radius = sys.argv[1].split(":")
wavelength, radius = float(wavelength), float(radius)
inputs = codebase.ROOT / "create_orac_lut" / "input_files"
mm = mirror.load_mmdat(inputs / "microphysics" / f"{model}.mm", inputs)
c = 0
ri = gsp._interpol_complex(mm.comp[c].cm, mm.comp[c].wl, [wavelength])
factor, refine = gsp.radius_upper_factor(mm, c), gsp.refined_xres(mm, c)
check_mu = np.cos(np.deg2rad(le.CHECK_THETA))
wavenumber = float(np.float32(1.0) / np.float32(wavelength))


def bwgp(dqv, factor, refine):
   start = time.perf_counter()
   bext, w, g, phi, vavg = cb.create_bwgp(mm.distname[c], radius, mm.s[c], ri, [wavelength], dqv, scode="mie",
                                          radius_upper_factor=factor, refined_xres=refine, phase_dtype=np.float64)
   return float(bext[0]), float(w[0]), float(g[0]), phi[:, 0], time.perf_counter() - start


def nodes(factor, refine):
   params, npts = cb.mie_integration_limits(mm.distname[c], radius, mm.s[c], wavenumber, factor, refine)
   rl, ru, truncated = cb.size_integration_limits(mm.distname[c], params, wavenumber)
   return (npts if npts is not None else cb.legacy_npts(rl, ru, wavenumber, cb.MIE_XRES)), ru


def compare(a, b):
   """b against the reference a: bulk and phase-function differences."""
   theta = le.CHECK_THETA
   rel = np.abs(b[3] / a[3] - 1.0)
   return {"extinction_rel": abs(b[0] / a[0] - 1.0), "w_abs": abs(b[1] - a[1]), "g_abs": abs(b[2] - a[2]),
           "phase_rel_0_5": float(rel[theta < 5].max()), "phase_rel_5_150": float(rel[(theta >= 5) & (theta < 150)].max()),
           "phase_rel_150_180": float(rel[theta >= 150].max())}


out = {"model": model, "wavelength": wavelength, "effective_radius": radius, "radius_upper_factor": factor,
       "refined_xres": refine}
# ---- radius integration (bulk and the phase function at the check angles)
production = bwgp(check_mu, factor, refine)
fine = bwgp(check_mu, factor, None if refine is None else refine / 4.0)
legacy = bwgp(check_mu, None, None)
out["nodes"] = {"production": nodes(factor, refine), "reference": nodes(factor, None if refine is None else refine / 4.0),
                "v23_legacy": nodes(None, None)}
out["production_vs_reference"] = compare(fine, production)
out["v23_grid_vs_reference"] = compare(fine, legacy)
out["production_bulk"] = {"extinction": production[0], "w": production[1], "g": production[2], "seconds": production[4]}
if factor is not None:
   longer = bwgp(check_mu, 5.0, refine)
   out["limit_3.5_vs_5.0"] = compare(longer, production)

# ---- adaptive Legendre expansion: production order and twice that order
lut_mrat = np.ones((1, 1))
lut_rm = np.asarray([[radius]])
start = time.perf_counter()
result = gsp._adaptive_mie_wavelength(mm, lut_mrat, lut_rm, ri, wavelength, f"{model} {wavelength} um {radius} um")
out["adaptive_seconds"] = time.perf_counter() - start
bext_c, w_c, g_c, vavg_c, phs, omega, lmom, diagnostics = result
length, degree, order = int(lmom[0]), int(diagnostics["degree"][0]), int(diagnostics["order"][0])
x, weights = le.gauss_legendre(2 * order)
phase2 = bwgp(np.concatenate((-x, check_mu)), factor, refine)[3]
omega2 = le.legendre_coefficients(-x, weights, phase2[:2 * order])
length2, noise2, error2 = le.expand(omega2, degree, check_mu, phase2[2 * order:])
chi = omega[:length, 0] / (2.0 * np.arange(length) + 1.0)
chi2 = omega2[:length] / (2.0 * np.arange(length) + 1.0)
out["legendre"] = {
   "degree_bound": degree, "order": order, "L": length, "noise": float(diagnostics["noise"][0]),
   "reconstruction_error": float(diagnostics["error"][0]),
   "L_at_2Nq": length2, "noise_at_2Nq": noise2, "reconstruction_error_at_2Nq": error2,
   "max_chi_difference_2Nq": float(np.max(np.abs(chi - chi2))),
   "max_omega_dropped_beyond_L": float(np.max(np.abs(omega[length:degree + 1, 0]))) if degree + 1 > length else 0.0,
   "omega0_minus_1": float(omega[0, 0] - 1.0), "omega1_over_3_minus_g": float(omega[1, 0] / 3.0 - g_c[0, 0]),
   "max_abs_chi": float(np.max(np.abs(chi))),
}
with open(sys.argv[2], "w") as handle:
   json.dump(out, handle, indent=1)
print(json.dumps({k: out[k] for k in ("model", "wavelength", "effective_radius", "nodes", "legendre")}))
