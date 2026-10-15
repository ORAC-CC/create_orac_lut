"""Planck source of the V25 emission call against analytic radiative transfer (read-only for the code).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v25/temperature_profiles/planck_check.py   -> results/planck_check.json

For a non-scattering cloud (single-scattering albedo 0) on a black, cold lower
boundary (BTEMP = 0, as in the DISORT adapter), the radiance leaving the top
in direction mu is the integral of the Planck source over optical depth.  For
the production profiles (cirrostratus, wet adiabat; cloud top 250 K) at
several 0.55 um optical depths (used as the channel optical depth) and the
Terra MODIS channel 20, 31 and 32 windows (centre wavenumber +-0.5 %, as
create_orac_cloud_lut):

* disort: the production DISORT call on the EMISSION_LAYERS sub-layers made by
  emission_layers, against the exact integral of a source linear in optical
  depth between the same boundary Planck radiances (the DISORT assumption):
  checks that DISORT receives and uses the temperature of every boundary, top
  down;
* reversed: the same integral with the boundary order reversed (a test that
  the comparison can tell the order);
* layers: the piecewise-linear integral for N = 10, 20, 40, 46 equal-optical-
  depth sub-layers against N = 4000 (the discretisation error of the emission
  representation without scattering), as a fraction of B(250 K); graded40:
  the same error for 40 sub-layers with boundaries at fractions (k / 40)^3,
  thin at the cloud top (an alternative not used by V25, for comparison);
* E_md: the normalised top radiance UU / B(250 K) at nadir and its isothermal
  value (1 - exp(-tau / mu));
* thin: the same DISORT/exact comparison for optically thin cirrostratus
  clouds (tau_055 = 1e-10 ... 0.3), whose sub-layers are optically very thin
  while their boundary temperatures differ (the cloud depth is held at 2 km).
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.idl_mirror import call_disort, cloud_temperature as ct, setup_disort   # noqa: E402
from oraclut.radiative_transfer.legacy_disort import disort_precision, plkavg, single_kernel   # noqa: E402

if disort_precision() != "single" or single_kernel() != "intel_dlm":
   raise SystemExit("set ORAC_LUT_DISORT_PRECISION=single and ORAC_LUT_DISORT_SP_KERNEL=intel_dlm")
f32 = np.float32
CHANNELS = {20: 3.788211, 31: 11.01842, 32: 12.03196}           # Terra MODIS centre wavelengths (um), as the products
SAZ = np.asarray([0.0, 30.0, 60.0, 75.0], f32)


def top_radiance(tau, planck, mu):
   """Exact top radiance of a non-scattering stack, Planck source linear in optical depth in each layer.

   Layer k contributes exp(-t_k / mu) [B_k (1 - e^-x) + (B_k+1 - B_k) (1 - e^-x (1 + x)) / x],
   x = dtau / mu, evaluated without cancellation for small x.
   """

   radiance, top = 0.0, 0.0
   for k, dtau in enumerate(tau):
      x = dtau / mu
      absorbed = -np.expm1(-x)
      ramp = x / 2.0 - x * x / 3.0 + x ** 3 / 8.0 if x < 1e-3 else (absorbed - x * np.exp(-x)) / x
      radiance += np.exp(-top / mu) * (planck[k] * absorbed + (planck[k + 1] - planck[k]) * ramp)
      top += dtau
   return radiance


def piecewise(cloudprof, tau_055, tau, wnlo, wnhi, mu, n, graded=False):
   fraction = (np.arange(n + 1) / n) ** (3 if graded else 1)
   planck = np.asarray([plkavg(wnlo, wnhi, float(t)) for t in ct.cloud_temperature(cloudprof, tau_055, fraction)])
   return top_radiance(tau * np.diff(fraction), planck, mu)


umu = np.concatenate([-np.cos(np.radians(SAZ)), np.cos(np.radians(SAZ))[::-1]]).astype(f32)
umu[np.abs(umu) == 1.0] *= f32(0.99999)
up = umu[SAZ.size:]
n = ct.EMISSION_LAYERS
disort_vars = setup_disort(60, n, SAZ.size, 1, 61)
pmom = np.zeros((61, 1), f32)
pmom[0] = 1.0
results = []
for name, taus in (("cirrostratus", (1.0, 4.0, 16.0, 64.0, 256.0)), ("wet_adiabat", (1.0, 16.0, 64.0, 256.0))):
   cloudprof = ct.setup_cloud_profile(name, 256.0, 250.0)
   for tau_055 in taus:
      for channel, wvl in CHANNELS.items():
         wn = f32(1e4) / f32(wvl)
         wnlo, wnhi = float(f32(0.995) * wn), float(f32(1.005) * wn)
         bbe = plkavg(wnlo, wnhi, 250.0)
         # the thermal optical depth of the channel is taken as tau_055 (a check of the call, not of the optics)
         subtau, subssa, subpmo, subtemp = ct.emission_layers(np.asarray([tau_055], f32), np.asarray([0.0], f32),
                                                              np.asfortranarray(pmom), np.asarray([1.0], f32),
                                                              tau_055, cloudprof, 250.0)
         uu = call_disort(disort_vars, subtau, subssa, subpmo, np.asarray([0.0, subtau.sum()], f32), umu,
                          np.asarray([0.0], f32), 0.0, f32(1.0), 0.0, plank=True, wnlo=wnlo, wnhi=wnhi, temp=subtemp,
                          nlayer=n)[5][SAZ.size:, 0, 0]
         planck = np.asarray([plkavg(wnlo, wnhi, float(t)) for t in subtemp])
         exact = np.asarray([top_radiance(subtau.astype(np.float64), planck, float(mu)) for mu in up])
         flipped = np.asarray([top_radiance(subtau.astype(np.float64), planck[::-1], float(mu)) for mu in up])
         reference = np.asarray([piecewise(cloudprof, tau_055, float(subtau.sum()), wnlo, wnhi, float(mu), 4000) for mu in up])
         layers = {str(m): float(np.max(np.abs([piecewise(cloudprof, tau_055, float(subtau.sum()), wnlo, wnhi, float(mu), m)
                                                for mu in up] - reference)) / bbe) for m in (10, 20, 40, 46)}
         layers["graded40"] = float(np.max(np.abs([piecewise(cloudprof, tau_055, float(subtau.sum()), wnlo, wnhi, float(mu), 40, True)
                                                   for mu in up] - reference)) / bbe)
         nadir = int(np.argmax(up))
         results.append({"profile": name, "tau_055": tau_055, "channel": channel,
                         "t_top": float(subtemp[0]), "t_base": float(subtemp[-1]),
                         "disort_vs_exact_max_rel": float(np.max(np.abs(uu / exact - 1.0))),
                         "reversed_vs_exact_min_rel": float(np.min(np.abs(flipped / exact - 1.0))),
                         "layers_error_over_b250": layers,
                         "e_md_nadir": float(uu[nadir] / bbe),
                         "e_md_isothermal_nadir": float(1.0 - np.exp(-float(subtau.sum()) / float(up[nadir])))})
         r = results[-1]
         print(f"{name:12s} tau {tau_055:6g} ch {channel}: T {r['t_top']:.1f}-{r['t_base']:.1f} K; DISORT/exact-1 "
               f"{r['disort_vs_exact_max_rel']:.1e}; reversed {r['reversed_vs_exact_min_rel']:.2f}; E_md nadir "
               f"{r['e_md_nadir']:.4f} (isothermal {r['e_md_isothermal_nadir']:.4f}); layers 20/40/46 "
               f"{layers['20']:.1e}/{layers['40']:.1e}/{layers['46']:.1e}; graded 40 {layers['graded40']:.1e}")
# Optically thin clouds: the cirrostratus cloud keeps its 2 km depth (16 K) below COT 0.0625, so each of
# the EMISSION_LAYERS sub-layers is optically very thin while its boundary temperatures differ.
thin = []
cloudprof = ct.setup_cloud_profile("cirrostratus", 256.0, 250.0)
for tau_055 in (1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1):
   for channel in (20, 31):
      wn = f32(1e4) / f32(CHANNELS[channel])
      wnlo, wnhi = float(f32(0.995) * wn), float(f32(1.005) * wn)
      subtau, subssa, subpmo, subtemp = ct.emission_layers(np.asarray([tau_055], f32), np.asarray([0.0], f32),
                                                           np.asfortranarray(pmom), np.asarray([1.0], f32),
                                                           tau_055, cloudprof, 250.0)
      uu = call_disort(disort_vars, subtau, subssa, subpmo, np.asarray([0.0, subtau.sum()], f32), umu,
                       np.asarray([0.0], f32), 0.0, f32(1.0), 0.0, plank=True, wnlo=wnlo, wnhi=wnhi, temp=subtemp,
                       nlayer=n)[5][SAZ.size:, 0, 0]
      planck = np.asarray([plkavg(wnlo, wnhi, float(t)) for t in subtemp])
      exact = np.asarray([top_radiance(subtau.astype(np.float64), planck, float(mu)) for mu in up])
      thin.append({"tau_055": tau_055, "channel": channel, "sublayer_tau": float(subtau[0]),
                   "disort_vs_exact_max_rel": float(np.max(np.abs(uu / exact - 1.0))),
                   "e_md_nadir_disort": float(uu[int(np.argmax(up))] / plkavg(wnlo, wnhi, 250.0)),
                   "e_md_nadir_exact": float(exact[int(np.argmax(up))] / plkavg(wnlo, wnhi, 250.0))})
      r = thin[-1]
      print(f"thin cirrostratus tau {tau_055:7.0e} ch {channel}: sub-layer tau {r['sublayer_tau']:.1e}; DISORT/exact-1 "
            f"{r['disort_vs_exact_max_rel']:.1e}; E_md nadir DISORT {r['e_md_nadir_disort']:+.3e}, exact {r['e_md_nadir_exact']:.3e}")
(HERE / "results" / "planck_check.json").write_text(json.dumps({"profiles": results, "thin": thin}, indent=1) + "\n")
