"""Checks of the mirrored ORAC forward model (validation only).

   python check_fm.py      (uses the modis_ice_agg baseline LUT of run_variant.py)

1. bicubic(): node values exact; gradients at nodes equal Int_LUT_Common's
   corner finite differences; continuity across a cell boundary.
2. solar(): analytic Jacobians (ORAC derivative_wrt_crp_parameter_brdf_0v,
   equation form 3) against centred finite differences of the full
   interpolate-then-forward-model chain.
3. Black surface: Ref = Tac_0v R_0v.
4. Lambertian surface: equation form 3 equals McGarragh et al. (2018)
   Eqs. (39)-(40) once the paper's mu0-normalised reflectance is converted
   (ORAC's Ref includes cos(sza): Ref = R_0v + mu0 [surface terms]).
5. Surface reciprocity of the pre-processor BRDF: rho_0d(theta) = rho_dv(theta).
6. thermal_radiance(): Jacobians against finite differences.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import orac_fm as fm                       # noqa: E402
from propagate import lut, tables          # noqa: E402

ok = True


def check(flag, message):
   global ok
   ok &= bool(flag)
   print(("PASS  " if flag else "FAIL  ") + message)


t = tables(lut("modis_ice_agg", "baseline"))
xt, xr = np.log10(t["axes"]["optical_depth"]), t["axes"]["effective_radius"]
G = t["R_dd"][:, :, 0]

# 1. bicubic
vals = [fm.bicubic(G, xt, xr, xt[i], xr[j])[0] for i in range(xt.size - 1) for j in range(xr.size - 1)]
check(np.allclose(vals, [G[i, j] for i in range(xt.size - 1) for j in range(xr.size - 1)], rtol=0, atol=1e-12),
      "bicubic reproduces node values")
i, j = 2, 2
f, dt, dr = fm.bicubic(G, xt, xr, xt[i], xr[j])
check(abs(dt - (G[i + 1, j] - G[i - 1, j]) / (xt[i + 1] - xt[i - 1])) < 1e-12 and
      abs(dr - (G[i, j + 1] - G[i, j - 1]) / (xr[j + 1] - xr[j - 1])) < 1e-12,
      "bicubic gradient at an interior node = ORAC's central corner difference")
eps = 1e-7
a = fm.bicubic(G, xt, xr, xt[3] - eps, 0.5 * (xr[1] + xr[2]))[0]
b = fm.bicubic(G, xt, xr, xt[3] + eps, 0.5 * (xr[1] + xr[2]))[0]
check(abs(a - b) < 1e-6, "bicubic continuous across a tau cell boundary")

# 2./3. solar Jacobians and black surface
s = 0
sza = t["axes"]["solar_zenith"][:, None, None]
vza = t["axes"]["satellite_zenith"][None, :, None]
raa = t["axes"]["relative_azimuth"][None, None, :]
rho = {"0v": np.full((4, 4, 7), 0.12), "0d": np.full((4, 4, 7), 0.10), "dv": np.full((4, 4, 7), 0.11), "dd": np.full((4, 4, 7), 0.09)}


def ref_at(lt, re, rho):
   crp, dcrp = {}, {}
   for name, table in (("R_0v", t["R_0v"][..., s]), ("T_00", t["T_00"][..., s][:, :, :, None, None]),
                       ("T_vv", t["T_00"][..., s][:, :, None, :, None]), ("T_0d", t["T_0d"][..., s][:, :, :, None, None]),
                       ("T_dv", t["T_dv"][..., 0][:, :, None, :, None]), ("R_dd", t["R_dd"][..., 0][:, :, None, None, None])):
      f, d1, d2 = fm.bicubic(table, xt, xr, lt, re)
      crp[name], dcrp[name] = f, (d1, d2)
   return fm.solar(crp, dcrp, rho, 0.95, 0.9, sza, vza), crp


lt0, re0 = 0.5 * (xt[2] + xt[3]), 0.5 * (xr[1] + xr[2])
(ref, grads), crp = ref_at(lt0, re0, rho)
h = 1e-5
fd_t = (ref_at(lt0 + h, re0, rho)[0][0] - ref_at(lt0 - h, re0, rho)[0][0]) / (2 * h)
fd_r = (ref_at(lt0, re0 + h, rho)[0][0] - ref_at(lt0, re0 - h, rho)[0][0]) / (2 * h)
check(np.allclose(grads[0], fd_t, rtol=1e-5, atol=1e-8) and np.allclose(grads[1], fd_r, rtol=1e-5, atol=1e-8),
      f"solar analytic Jacobians = finite differences (max |diff| {max(np.abs(grads[0]-fd_t).max(), np.abs(grads[1]-fd_r).max()):.1e})")
zero = {k: np.zeros((4, 4, 7)) for k in rho}
(ref0, _), crp0 = ref_at(lt0, re0, zero)
tac0v = 0.95 ** (1 / np.cos(np.deg2rad(sza))) * 0.95 ** (1 / np.cos(np.deg2rad(vza)))
check(np.allclose(ref0, tac0v * crp0["R_0v"], rtol=1e-12, atol=1e-15), "black surface: Ref = Tac_0v R_0v")

# 4. Lambertian: equation form 3 vs McGarragh et al. (2018) Eqs. (39)-(40)
r_l = 0.2
lam = {k: np.full((4, 4, 7), r_l) for k in rho}
(ref_l, _), c = ref_at(lt0, re0, lam)
mu0 = np.cos(np.deg2rad(sza))
tbc, tac = 0.9, 0.95
tbc0, tbcv, tbcd = tbc ** (1 / mu0), tbc ** (1 / np.cos(np.deg2rad(vza))), tbc
num = (tbc0 * c["T_00"] * r_l + tbcd * c["T_0d"] * r_l)
paper = (c["T_00"] * r_l * c["T_vv"] * tbc0 * tbcv + tbcd * c["T_0d"] * r_l * c["T_vv"] * tbcv
         + num / (1 - r_l * c["R_dd"] * tbcd**2) * (tbcd * c["T_dv"] + c["R_dd"] * tbcd**2 * r_l * c["T_vv"] * tbcv))
check(np.allclose(ref_l, tac0v * (c["R_0v"] + mu0 * paper), rtol=1e-10),
      "Lambertian: equation form 3 = paper Eqs. (39)-(40) with Ref = R_0v + mu0 x surface terms")

# 5. surface reciprocity
f = np.array([0.3, 0.15, 0.02])
th = np.array([0.0, 30.0, 60.0, 75.0])
r1 = fm.rtls_rho(f, th, 20.0, 40.0)["0d"]
r2 = fm.rtls_rho(f, 20.0, th, 40.0)["dv"]
check(np.allclose(r1, r2, rtol=1e-6), f"pre-processor BRDF: rho_0d(theta) = rho_dv(theta) {np.round(r1, 4)}")

# 6. thermal Jacobians
k = t["thermal"].index(31)
ich = t["all"].index(31)
atm = {"Rbc_up": 80.0, "B_c": 60.0, "Rac_dwn": 10.0, "Tac": 0.85, "Rac_up": 12.0}


def rad_at(lt, re):
   crp, dcrp = {}, {}
   for name, table in (("T_dv", t["T_dv"][..., ich]), ("R_dv", t["R_dv"][..., ich]), ("E_md", t["E_md"][..., k])):
      f_, d1, d2 = fm.bicubic(table, xt, xr, lt, re)
      crp[name], dcrp[name] = f_, (d1, d2)
   return fm.thermal_radiance(crp, dcrp, atm)


r, g = rad_at(lt0, re0)
fd = (rad_at(lt0 + h, re0)[0] - rad_at(lt0 - h, re0)[0]) / (2 * h)
check(np.allclose(g[0], fd, rtol=1e-5, atol=1e-8), "thermal analytic Jacobian = finite difference")
print("ALL PASSED" if ok else "FAILURES")
sys.exit(0 if ok else 1)
