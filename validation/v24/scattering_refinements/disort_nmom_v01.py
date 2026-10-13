"""DISORT with more than 1000 Legendre moments: acceptance, use, padding and convergence (read-only).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v24/scattering_refinements/disort_nmom.py <model> <wavelength um> <r_e um> <out.json>

The phase function is expanded by the production V24 routine
(generate_scattering_properties._adaptive_mie_wavelength, working tree), which
also returns every coefficient up to the order Nq, so the coefficients that
King's criterion drops beyond L are available.  A single homogeneous cloud
layer (optical depths 1, 8, 64; single-scattering albedo of the class; no
Rayleigh, no gas) is then calculated by the production setup_disort /
call_disort with the original single-precision Intel DISORT kernel, NSTR = 60,
the Grid B view geometry (21 satellite zeniths, both hemispheres, 37 relative
azimuths 0-180), the direct beam (fbeam = 100) at solar zeniths 30 and 60
degrees and the diffuse source (fisot = 100), for these numbers of moments
NMOM (array rows, highest order NMOM - 1):

   L            the expansion (production when L > 1000); max(L, 61) when L < 61,
                because DISORT refuses NMOM < NSTR (CHEKIN) for a scattering layer
   max(1000, L) the production rule
   L + 512      zero padded
   max(L, 61)   the historical dbcc42c rule (when L < 1000)
   N < L        truncations (diagnostic only: production never truncates)
   L < N <= D+1 the coefficients below King's 1e-9 threshold added back

and reports bitwise equalities, the largest differences of the reflectance
R = pi UU / fbeam and of the fluxes, and the time per call.
"""

import json
import sys
import time

import numpy as np

import codebase

codebase.select("v24")
mirror = codebase.module("oraclut.idl_mirror")
gsp = codebase.module("oraclut.idl_mirror.generate_scattering_properties")
legacy_disort = codebase.module("oraclut.radiative_transfer.legacy_disort")
if legacy_disort.disort_precision() != "single" or legacy_disort.single_kernel() != "intel_dlm":
   raise SystemExit("set ORAC_LUT_DISORT_PRECISION=single and ORAC_LUT_DISORT_SP_KERNEL=intel_dlm (the V23/V24 kernel)")

model, wavelength, radius, out_file = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
inputs = codebase.ROOT / "create_orac_lut" / "input_files"
mm = mirror.load_mmdat(inputs / "microphysics" / f"{model}.mm", inputs)
ri = gsp._interpol_complex(mm.comp[0].cm, mm.comp[0].wl, [wavelength])
bext_c, w_c, g_c, vavg_c, phs, omega, lmom, diagnostics = gsp._adaptive_mie_wavelength(
   mm, np.ones((1, 1)), np.asarray([[radius]]), ri, wavelength, f"{model} {wavelength} um {radius} um")
length, degree = int(lmom[0]), int(diagnostics["degree"][0])
omega = omega[:, 0]
chi_all = (omega / (2.0 * np.arange(omega.size) + 1.0)).astype(np.float32)        # every computed moment
ssa = np.float32(w_c[0, 0])

f32 = np.float32
saz = np.arange(21, dtype=f32) * f32(3.75)
raa = np.arange(37, dtype=f32) * f32(5.0)
umu = np.zeros(42, f32)
umu[:21] = -np.cos(np.deg2rad(saz)).astype(f32)
umu[21:] = -umu[20::-1]
bd = np.abs(umu) == 1.0
umu[bd] = f32(0.99999) * umu[bd] / np.abs(umu[bd])
up = slice(41, 20, -1)


def run(nmom_rows, tau, source):
   """DISORT for one layer with the first nmom_rows moments of the expansion (zero beyond L or the array)."""
   pm = np.zeros(nmom_rows, f32)
   keep = min(nmom_rows, length) if source != "tail" else min(nmom_rows, chi_all.size)
   pm[:keep] = chi_all[:keep]
   pmom = np.asfortranarray(pm[:, None])
   disort_vars = mirror.setup_disort(60, 1, 21, 37, nmom_rows)
   dtau = np.asarray([tau], f32)
   ssalb = np.asarray([ssa], f32)
   utau = np.asarray([0.0, tau], f32)
   results = {}
   start = time.perf_counter()
   for name, fbeam, umu0, fisot in (("diffuse", 0.0, np.cos(np.deg2rad(50.0)), 100.0),
                                    ("sza30", 100.0, np.cos(np.deg2rad(30.0)), 0.0),
                                    ("sza60", 100.0, np.cos(np.deg2rad(60.0)), 0.0)):
      rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = mirror.call_disort(
         disort_vars, dtau.copy(), ssalb.copy(), pmom.copy(), utau.copy(), umu, raa, fbeam, f32(umu0), fisot)
      results[name] = {"R": np.asarray(uu[up, 0, :], np.float64) * np.pi / max(fbeam, fisot),
                       "flux": np.concatenate((rfldir, rfldn, flup)).astype(np.float64)}
   return results, (time.perf_counter() - start) / 3.0


def bitwise(a, b):
   return all(np.array_equal(a[k]["R"], b[k]["R"]) and np.array_equal(a[k]["flux"], b[k]["flux"]) for k in a)


def max_diff(a, b, key):
   return float(max(np.max(np.abs(a[k][key] - b[k][key])) for k in a))


report = {"model": model, "wavelength": wavelength, "effective_radius": radius, "L": length, "degree_bound": degree,
          "order": int(diagnostics["order"][0]), "ssa": float(ssa), "nstr": 60, "cases": {}}
exact = max(length, 61)                 # DISORT needs NMOM >= NSTR (60) when there is scattering
report["exact_rows"] = exact
for tau in (1.0, 8.0, 64.0):
   base, t_base = run(exact, tau, "expansion")
   case = {"seconds_per_call_at_L": t_base}
   production, t_prod = run(max(1000, length), tau, "expansion")
   case["production_max(1000,L)_bitwise_equal_L"] = bitwise(base, production)
   case["seconds_per_call_production"] = t_prod
   padded, _ = run(exact + 512, tau, "expansion")
   case["zero_padded_L+512_bitwise_equal_L"] = bitwise(base, padded)
   if length < 1000:
      historical, _ = run(max(length, 61), tau, "expansion")
      case["historical_max(L,61)_bitwise_equal_production"] = bitwise(historical, production)
   one_thousand, t_1000 = run(1000, tau, "expansion")
   case["seconds_per_call_1000"] = t_1000
   truncations = {}
   for n in (200, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000):
      if 61 <= n < length:
         result, _ = run(n, tau, "expansion")
         truncations[n] = {"max_abs_R_difference": max_diff(result, base, "R"),
                           "max_abs_R_difference_sza30": float(np.max(np.abs(result["sza30"]["R"] - base["sza30"]["R"]))),
                           "fluxes_bitwise_equal": all(np.array_equal(result[k]["flux"], base[k]["flux"]) for k in base)}
   case["truncated_to_N"] = truncations
   tails = {}
   for n in sorted({min(length + 100, degree + 1), min(length + 400, degree + 1), degree + 1}):
      n = max(n, 61)
      if n > length and n != exact or n > exact:
         result, _ = run(n, tau, "tail")
         tails[n] = {"max_abs_R_difference": max_diff(result, base, "R"),
                     "max_abs_chi_added": float(np.max(np.abs(chi_all[length:n])))}
   case["below_king_threshold_added_to_N"] = tails
   case["max_R"] = float(max(np.max(np.abs(base[k]["R"])) for k in base))
   report["cases"][str(tau)] = case
with open(out_file, "w") as handle:
   json.dump(report, handle, indent=1)
print(json.dumps(report, indent=1)[:3000])
