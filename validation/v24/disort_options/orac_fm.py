"""The production ORAC fast forward model, mirrored for DISORT-option tests (validation only).

Mirrors ORAC-CC/orac commit eb9a1233c1136ef95fa0c37f650fe519e9134ed4 (master,
2026-08-03; read-only clone outside this repository).  Paths below are
relative to that repository; every formula cites its routine and lines.

Configuration mirrored: one cloud layer (Ctrl%Approach /= AppCld2L), overcast
(X(IFr) = 1), the cloud default i_equation_form = 3 and use_full_brdf = true
(src/read_driver.F90:481, 620, 1250-1271), NetCDF V2 LUTs, so
get_T_dv_from_T_0d = false (read_driver.F90:1292-1297; src/read_sad.F90:62-80),
bicubic (tau, re) interpolation, the default LUTIntSelm (read_driver.F90:624)
in the conda build that defines INCLUDE_NR (config/lib.conda.inc:37-39).

LUT -> ORAC operators (src/read_sad_lut.F90:979-1103; src/fm_routines.F90:32-41;
src/set_crp_solar.F90:128-173; src/set_crp_thermal.F90:102-110):
   R_0v -> Rbd (IR_0v)  interpolated in (tau, satzen, solzen, relazi, re)
   R_dd -> Rfd (IR_dd)  (tau, re)
   T_00 -> Tb  (IT_00)  (tau, solzen, re);  also IT_vv = the same table at satzen
   T_0d -> Tfbd (IT_0d) (tau, solzen, re)
   T_dv -> Td  (IT_dv)  (tau, satzen, re)
   R_dv -> Rd  (IR_dv)  (tau, satzen, re)          thermal only
   E_md -> Em  (IEm)    (tau, satzen, re)          thermal only
   (R_0d and T_dd are used only by the two-layer forms.)

Measurement: solar Ref = pi L / F0 including cos(sza) (the LUT's R_0v
normalisation); thermal and mixed channels: brightness temperature from R2T
(src/r2t.F90); mixed: BT of Rad + f0 Ref, f0 = the LUT's F0 (= F/pi)
(src/fm.F90 FM, mixed-channel block; src/get_spixel.F90:276).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Interpolation: src/int_lut_routines.F90 Int_LUT_Common (bicubic branch,
# lines 163-213) with GZero stencils (src/set_gzero.F90:95-174).  The grid
# coordinates are log10(tau) (src/read_sad_lut.F90:870-878, our LUTs carry
# spacing = "uneven_logarithmic") and re.  Numerical Recipes bcuint is the
# bicubic Hermite patch through the corner values, first derivatives and
# cross derivative; it is written here in Hermite form (identical polynomial).
# ---------------------------------------------------------------------------


def _locate(x, value):
   """Index i0 with x[i0] <= value < x[i0+1], clamped to [0, n-2] (ORAC locate with lims)."""

   i = int(np.searchsorted(x, value, side="right") - 1)
   return min(max(i, 0), x.size - 2)


def _stencil(i0, n):
   """(m1, 0, 1, p1) of Set_GZero lines 117-139 (one-sided at the edges)."""

   i1 = i0 + 1
   if i0 == 0:
      return i0, i0, i1, i1 + 1
   if i1 == n - 1:
      return i0 - 1, i0, i1, i1
   return i0 - 1, i0, i1, i1 + 1


def _hermite(t):
   return np.array([2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t, -2 * t**3 + 3 * t**2, t**3 - t**2])


def _dhermite(t):
   return np.array([6 * t**2 - 6 * t, 3 * t**2 - 4 * t + 1, -6 * t**2 + 6 * t, 3 * t**2 - 2 * t])


def bicubic(G, xt, xr, log_tau, re):
   """F, dF/dlog10(tau), dF/dre from G[tau, re, ...] (trailing axes vectorised).

   G is the angle-interpolated table on the full (tau, re) grid; xt = log10 tau
   nodes, xr = re nodes.  Corner derivatives exactly as Int_LUT_Common
   lines 165-200: central differences over (m1, 1) at the lower corner and
   (0, p1) at the upper corner of each axis.
   """

   i0 = _locate(xt, log_tau)
   j0 = _locate(xr, re)
   tm, t0, t1, tp = _stencil(i0, xt.size)
   rm, r0, r1, rp = _stencil(j0, xr.size)
   d1 = xt[t1] - xt[t0]
   d2 = xr[r1] - xr[r0]
   t = (log_tau - xt[t0]) / d1
   u = (re - xr[r0]) / d2
   y = {(0, 0): G[t0, r0], (1, 0): G[t1, r0], (1, 1): G[t1, r1], (0, 1): G[t0, r1]}
   dyt = {(0, 0): (G[t1, r0] - G[tm, r0]) / (xt[t1] - xt[tm]), (1, 0): (G[tp, r0] - G[t0, r0]) / (xt[tp] - xt[t0]),
          (1, 1): (G[tp, r1] - G[t0, r1]) / (xt[tp] - xt[t0]), (0, 1): (G[t1, r1] - G[tm, r1]) / (xt[t1] - xt[tm])}
   dyr = {(0, 0): (G[t0, r1] - G[t0, rm]) / (xr[r1] - xr[rm]), (1, 0): (G[t1, r1] - G[t1, rm]) / (xr[r1] - xr[rm]),
          (1, 1): (G[t1, rp] - G[t1, r0]) / (xr[rp] - xr[r0]), (0, 1): (G[t0, rp] - G[t0, r0]) / (xr[rp] - xr[r0])}
   dyy = {(0, 0): (G[t1, r1] - G[t1, rm] - G[tm, r1] + G[tm, rm]) / ((xt[t1] - xt[tm]) * (xr[r1] - xr[rm])),
          (1, 0): (G[tp, r1] - G[tp, rm] - G[t0, r1] + G[t0, rm]) / ((xt[tp] - xt[t0]) * (xr[r1] - xr[rm])),
          (1, 1): (G[tp, rp] - G[tp, r0] - G[t0, rp] + G[t0, r0]) / ((xt[tp] - xt[t0]) * (xr[rp] - xr[r0])),
          (0, 1): (G[t1, rp] - G[t1, r0] - G[tm, rp] + G[tm, r0]) / ((xt[t1] - xt[tm]) * (xr[rp] - xr[r0]))}
   ht, hu, dht, dhu = _hermite(t), _hermite(u), _dhermite(t), _dhermite(u)
   f = dft = dfr = 0.0
   for (a, b) in y:
      # Hermite basis indices: value at node a uses h[0] (a = 0) or h[2] (a = 1); slope h[1] or h[3]
      va, sa = (0, 1) if a == 0 else (2, 3)
      vb, sb = (0, 1) if b == 0 else (2, 3)
      terms = ((y[a, b], va, vb, 1.0), (dyt[a, b] * d1, sa, vb, 1.0), (dyr[a, b] * d2, va, sb, 1.0),
               (dyy[a, b] * d1 * d2, sa, sb, 1.0))
      for value, it, iu, _ in terms:
         f = f + value * ht[it] * hu[iu]
         dft = dft + value * dht[it] * hu[iu] / d1
         dfr = dfr + value * ht[it] * dhu[iu] / d2
   return f, dft, dfr


# ---------------------------------------------------------------------------
# Solar: src/fm_solar.F90 FM_Solar, i_layer = 0, use_full_brdf, i_equation_form = 3
#   a = 1 - rho_dd R_dd Tbc_dd                                          (883)
#   b = T_00 rho_0d Tbc_0 + T_0d rho_dd Tbc_d                           (885-886)
#   c = T_dv Tbc_d + R_dd rho_dv T_vv Tbc_dd Tbc_v                      (924-925)
#   e = R_0v + (T_00 rho_0v T_vv Tbc_0v + T_0d rho_dv T_vv Tbc_dv
#               + b c / a) / sec(sza)                                   (926-929)
#   Ref = Fr Tac_0v e + (1 - Fr) Ref_clear,  Fr = 1                     (940-945)
# Transmittances (lines 749-804): Tac_0 = Tac^sec(sza), Tac_v = Tac^sec(vza),
# Tac_d = Tac^dif_trans_fac with dif_trans_fac = 1 (line 643); likewise Tbc.
# Derivatives: derivative_wrt_crp_parameter_brdf_0v, i_equation_form = 3
# (lines 269-286), e_l with the interpolated d_CRP.
# ---------------------------------------------------------------------------


def solar(crp, d_crp, rho, tac, tbc, sza, vza):
   """Ref and dRef/dx (x = log10 tau, re) for one channel.  crp/d_crp: dicts of operators."""

   sec0, secv = 1.0 / np.cos(np.deg2rad(sza)), 1.0 / np.cos(np.deg2rad(vza))
   tac_0, tac_v = tac**sec0, tac**secv
   tbc_0, tbc_v, tbc_d = tbc**sec0, tbc**secv, tbc
   tac_0v = tac_0 * tac_v
   tbc_0v, tbc_dv, tbc_dd = tbc_0 * tbc_v, tbc_d * tbc_v, tbc_d * tbc_d
   r0v, rdd, t00, tvv, t0d, tdv = (crp[k] for k in ("R_0v", "R_dd", "T_00", "T_vv", "T_0d", "T_dv"))
   a = 1.0 - rho["dd"] * rdd * tbc_dd
   b = t00 * rho["0d"] * tbc_0 + t0d * rho["dd"] * tbc_d
   c = tdv * tbc_d + rdd * rho["dv"] * tvv * tbc_dd * tbc_v
   e = r0v + (t00 * rho["0v"] * tvv * tbc_0v + t0d * rho["dv"] * tvv * tbc_dv + b * c / a) / sec0
   ref = tac_0v * e
   grads = []
   for k in range(2):
      d = {name: d_crp[name][k] for name in d_crp}
      e_l = (d["R_0v"]
             + ((d["T_00"] * rho["0v"] * tvv + t00 * rho["0v"] * d["T_vv"]) * tbc_0v
                + (d["T_0d"] * rho["dv"] * tvv + t0d * rho["dv"] * d["T_vv"]) * tbc_dv
                + ((d["T_00"] * rho["0d"] * tbc_0 + d["T_0d"] * rho["dd"] * tbc_d) * c
                   + b * (d["T_dv"] * tbc_d + d["R_dd"] * rho["dv"] * tvv * tbc_dd * tbc_v
                          + rdd * rho["dv"] * d["T_vv"] * tbc_dd * tbc_v)) / a
                + b * c * rho["dd"] * d["R_dd"] * tbc_dd / (a * a)) / sec0)
      grads.append(tac_0v * e_l)
   return ref, grads


# ---------------------------------------------------------------------------
# Thermal: src/fm_thermal.F90 FM_Thermal, single layer (lines 253-275, 416-422)
#   R = (Rbc_up T_dv + B(Tc) E_md + Rac_dwn R_dv) Tac + Rac_up,  Fr = 1
#   BT = R2T(R);  d_BT = dT/dR d_R
# R2T / T2R: src/r2t.F90:71-82, src/t2r.F90 (band-corrected Planck with the
# LUT's B1, B2, T1, T2).
# ---------------------------------------------------------------------------


def t2r(T, planck):
   B1, B2, T1, T2 = planck
   t_eff = T * T2 + T1
   c = np.exp(B2 / t_eff)
   return B1 / (c - 1.0), B1 * (B2 / t_eff) * c * T2 / (t_eff * (c - 1.0) ** 2)


def r2t(R, planck):
   B1, B2, T1, T2 = planck
   R = np.maximum(R, 1e-6)
   c = np.log(B1 / R + 1.0)
   T = (B2 / c - T1) / T2
   return T, B1 * B2 / (T2 * c * c * R * (R + B1))


def thermal_radiance(crp, d_crp, atm):
   """Overcast TOA radiance and dR/dx (x = log10 tau, re) for one channel."""

   r = (atm["Rbc_up"] * crp["T_dv"] + atm["B_c"] * crp["E_md"] + atm["Rac_dwn"] * crp["R_dv"]) * atm["Tac"] + atm["Rac_up"]
   grads = [atm["Tac"] * (atm["Rbc_up"] * d_crp["T_dv"][k] + atm["B_c"] * d_crp["E_md"][k] + atm["Rac_dwn"] * d_crp["R_dv"][k])
            for k in range(2)]
   return r, grads


# ---------------------------------------------------------------------------
# Measurement uncertainty: src/get_measurements.F90:174-404 (SySelm = SelmAux,
# Homog = Coreg = true for cloud: src/read_driver.F90:521-523) for NetCDF LUTs:
#   solar   sigma_L = L / snr  (MODIS)  or  sqrt(rua^2 L^2 + rub^2 L + ruc^2)
#           (SLSTR/AATSR, ru2 = ru^2: src/read_sad_lut.F90:1223-1238);
#           + 1 % (Homog) and 2 % (Coreg) of the signal, in quadrature
#   thermal sigma_BT = NEBT dR/dT(T0) / dR/dT(Tm), T0 = refbt; + 0.5 K, 0.15 K;
#   mixed (day) additionally (1 % Rad / dR/dT) and (2 % Rad / dR/dT), Rad = T2R(Tm)
#   (USE_OLD_MIXED_UNCERTAINTY = false, THERMAL_ALWAYS_HAS_HOMOG_COREG = true:
#   get_measurements.F90:106-107).
# ---------------------------------------------------------------------------


def sigma_solar(ref, noise_rel, production=True):
   s2 = (noise_rel * ref) ** 2
   if production:
      s2 = s2 + (0.01 * ref) ** 2 + (0.02 * ref) ** 2
   return np.sqrt(s2)


def sigma_bt(bt, nebt, refbt, planck, mixed=False, production=True):
   _, d0 = t2r(refbt, planck)
   rad, dm = t2r(bt, planck)
   s2 = (nebt * d0 / dm) ** 2
   if production:
      s2 = s2 + 0.5**2 + 0.15**2
      if mixed:
         s2 = s2 + (0.01 * rad / dm) ** 2 + (0.02 * rad / dm) ** 2
   return np.sqrt(s2)


# ---------------------------------------------------------------------------
# Surface: pre_processing/ross_thick_li_sparse_r.F90 (land BRDF of ORAC's
# pre-processor): Ross-thick + Li-sparse-R kernels with p_ross = (0,),
# p_li_r = (b/r = 1, h/b = 2) (lines 557-558); rho_0v direct (614-636);
# rho_0d, rho_dv by 4 x 4 Gauss-Legendre quadrature over (theta, phi) (641-735);
# rho_dd by 4 x 4 x 4 quadrature (742-787).  Relative azimuth 0 = hot spot.
# ---------------------------------------------------------------------------


def _rtls(theta1, theta2, phi, f):
   cos1, cos2 = np.cos(theta1), np.cos(theta2)
   sin1, sin2 = np.sin(theta1), np.sin(theta2)
   cos_phi, sin_phi = np.cos(phi), np.sin(phi)
   # Ross thick (lines 346-403)
   cos_ksi = np.clip(cos1 * cos2 + sin1 * sin2 * cos_phi, -1.0, 1.0)
   ksi = np.arccos(cos_ksi)
   k_ross = ((np.pi / 2 - ksi) * cos_ksi + np.sin(ksi)) / (cos1 + cos2) - np.pi / 4
   # Li sparse R (lines 186-321), b/r = 1, h/b = 2
   tan_i, tan_r = np.tan(theta1), np.tan(theta2)
   th_i, th_r = np.arctan(tan_i), np.arctan(tan_r)
   a = np.cos(th_i) * np.cos(th_r)
   b = np.sin(th_i) * np.sin(th_r)
   c = tan_i**2 + tan_r**2
   d = tan_i * tan_r * 2.0
   e = tan_i**2 * tan_r**2
   r = 1.0 / np.cos(th_i) + 1.0 / np.cos(th_r)
   g = 2.0 / r
   cos_ksi_p = a + b * cos_phi
   p_ = (1.0 + cos_ksi_p) / (np.cos(th_i) * np.cos(th_r))
   dd = np.sqrt(np.maximum(c - d * cos_phi, 0.0))
   h = np.sqrt(dd * dd + e * sin_phi**2)
   cos_t = g * h
   q = np.where(cos_t > 1.0, 1.0,
                1.0 - (np.arccos(np.minimum(cos_t, 1.0)) - np.sqrt(np.maximum(1.0 - cos_t**2, 0.0)) * cos_t) / np.pi)
   k_li = 0.5 * p_ - q * r
   return f[0] + f[1] * k_ross + f[2] * k_li


def rtls_rho(f, sza, vza, raa):
   """rho_0v, rho_0d, rho_dv, rho_dd of the ORAC pre-processor for kernel weights f."""

   xt, wt = np.polynomial.legendre.leggauss(4)
   qt = (xt + 1.0) * np.pi / 4.0
   wt = wt * np.pi / 4.0
   xp, wp = np.polynomial.legendre.leggauss(4)
   qp = (xp + 1.0) * np.pi
   wp = wp * np.pi
   cs = np.cos(qt) * np.sin(qt) * wt
   s0, sv, ra = np.deg2rad(sza), np.deg2rad(vza), np.deg2rad(raa)
   rho_0v = _rtls(s0, sv, ra, f)
   rho_0d = sum(cs[k] * sum(_rtls(s0, qt[k], qp[l], f) * wp[l] for l in range(4)) for k in range(4)) / np.pi
   rho_dv = sum(cs[k] * sum(_rtls(qt[k], sv, qp[l], f) * wp[l] for l in range(4)) for k in range(4)) / np.pi
   rho_dd = 2.0 * sum(cs[k] * sum(cs[l] / np.pi * sum(_rtls(qt[k], qt[l], qp[m], f) * wp[m] for m in range(4))
                                  for l in range(4)) for k in range(4))
   return {"0v": rho_0v, "0d": rho_0d, "dv": rho_dv, "dd": rho_dd}
