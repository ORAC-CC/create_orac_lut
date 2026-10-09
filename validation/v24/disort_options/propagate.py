"""Propagate each DISORT experiment's complete operator set through the ORAC forward model (validation only).

For every case and experiment (run_variant.py products, reused unchanged)
the complete, internally consistent set of LUT operators of that experiment
is put through the production ORAC fast forward model (orac_fm.py, mirroring
ORAC-CC/orac eb9a1233), for:

   * surfaces (solar):  black (diagnostic), dark Lambertian 0.05, vegetated land
     (ORAC pre-processor Ross-thick/Li-sparse-R BRDF, representative MODIS-like
     kernel weights), bright snow-like Lambertian;
   * every compact-grid geometry (solar zenith x view zenith x relative azimuth
     nodes; all are Grid B nodes, so no angular interpolation is involved);
   * cloud states at every (tau, re) node and at every cell mid-point in
     (log10 tau, re), with ORAC's bicubic interpolation and its gradients;
   * clear-sky terms from atmosphere.py (cloud top 4 km liquid, 9 km ice).

Measurements: solar reflectance; thermal and mixed (3.7 um) brightness
temperature, as ORAC.  Compared with the numerical reference (double
precision, NSTR = 236) and with the production baseline:
   r = |dy| / sigma, sigma = the production ORAC Sy (noise + homogeneity +
   co-registration) and, separately, instrument noise only;
   Jacobians dy/dlog10(tau), dy/dre: |dJ| / sigma (per unit state) and relative.
Operator cancellation: for each measurement the total dy is compared with
the sum of the magnitudes of the linear contributions of the individual
operators (sum_k |dy/dop_k d op_k|): ratio < 1 means cancellation.

   python propagate.py           -> validation/v24/results/disort_options/fm_*.csv, fm_arrays_*.npz
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
import orac_fm as fm                          # noqa: E402
from atmosphere import terms                  # noqa: E402
from run_variant import CASES, EXPERIMENTS, TMP   # noqa: E402
from oraclut.io.lut import read_lut           # noqa: E402

RESULTS = ROOT / "validation" / "v24" / "results" / "disort_options"
CHANNELS_OF = {case: set(spec[3]) for case, spec in CASES.items()}
CURRENT_CASE = [None]
REFERENCE = "double_nstr236"
PLATFORM = {"modis_liquid": ("terra", "modis", "liquid"), "modis_ice_agg": ("terra", "modis", "ice"),
            "slstr_liquid": ("sentinel-3a", "slstr", "liquid"), "slstr_ice_sph": ("sentinel-3a", "slstr", "ice")}
PLATFORM.update({case + "_bs": spec for case, spec in list(PLATFORM.items())})
# Representative surfaces.  Vegetated land: Ross-thick/Li-sparse-R weights (f_iso, f_vol, f_geo)
# by wavelength band; ORAC's pre-processor uses 1 - emissivity at 3.7 um over land.
VEGETATION = [(0.50, (0.03, 0.01, 0.004)), (0.70, (0.04, 0.02, 0.006)), (1.00, (0.30, 0.15, 0.02)),
              (1.80, (0.22, 0.08, 0.03)), (2.50, (0.12, 0.04, 0.02))]
SNOW = [(0.50, 0.90), (0.70, 0.88), (1.00, 0.80), (1.80, 0.10), (2.50, 0.05)]
LONGWAVE_RHO = 0.03


def lut(case, experiment):
   path = TMP / "runs" / case / experiment
   files = sorted(path.glob("*.nc"))
   return read_lut(files[0]) if files else None


def tables(l):
   """Operator tables as arrays (tau, re, ...) per channel, from a LUT."""

   v = lambda n: np.asarray(l.variables[n], dtype=np.float64)
   out = {"solar": [int(c) for c in v("solar_channel_id")], "thermal": [int(c) for c in v("thermal_channel_id")],
          "all": [int(c) for c in v("channel_id")]}
   out["R_0v"] = np.transpose(v("R_0v"), (3, 4, 2, 1, 0, 5))        # tau, re, sza, vza, raa, solar
   out["T_00"] = np.transpose(v("T_00"), (1, 2, 0, 3))              # tau, re, sza, solar
   out["T_0d"] = np.transpose(v("T_0d"), (1, 2, 0, 3))
   out["T_dv"] = np.transpose(v("T_dv"), (1, 2, 0, 3))              # tau, re, vza, channel
   out["R_dv"] = np.transpose(v("R_dv"), (1, 2, 0, 3))
   out["R_dd"] = v("R_dd")                                          # tau, re, channel
   out["E_md"] = np.transpose(v("E_md"), (1, 2, 0, 3))              # tau, re, vza, thermal
   out["axes"] = {n: v(n) for n in ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")}
   return out


def channel_constants(l, instrument):
   v = lambda n: np.asarray(l.variables[n], dtype=np.float64)
   solar = [int(c) for c in v("solar_channel_id")]
   thermal = [int(c) for c in v("thermal_channel_id")]
   allc = [int(c) for c in v("channel_id")]
   wn = v("central_wavenumber")
   const = {}
   for i, ch in enumerate(solar):
      noise = float(v("rua")[i]) if instrument == "slstr" else (1.0 / float(v("snr")[i]) if float(v("snr")[i]) > 0 else None)
      const.setdefault(ch, {}).update({"F0": float(v("F0")[i]), "noise_rel": noise,
                                       "noise_source": "rua" if instrument == "slstr" else "1/snr"})
   for i, ch in enumerate(thermal):
      const.setdefault(ch, {}).update({"planck": tuple(float(v(k)[i]) for k in ("B1", "B2", "T1", "T2")),
                                       "nedt": float(v("nedt")[i]), "refbt": float(v("refbt")[i])})
   for ch in allc:
      const[ch]["wavelength_um"] = 1e4 / float(wn[allc.index(ch)])
   return const


def surface(name, wavelength, sza, vza, raa):
   shape = np.broadcast(sza, vza, raa).shape
   if name == "black":
      return {k: np.zeros(shape) for k in ("0v", "0d", "dv", "dd")}
   if name == "dark":
      return {k: np.full(shape, 0.05) for k in ("0v", "0d", "dv", "dd")}
   if wavelength > 3.0:
      return {k: np.full(shape, LONGWAVE_RHO if name == "vegetation" else 0.02) for k in ("0v", "0d", "dv", "dd")}
   if name == "snow":
      value = next(r for w, r in SNOW if wavelength <= w) if wavelength <= SNOW[-1][0] else SNOW[-1][1]
      return {k: np.full(shape, value) for k in ("0v", "0d", "dv", "dd")}
   f = next(k for w, k in VEGETATION if wavelength <= w) if wavelength <= VEGETATION[-1][0] else VEGETATION[-1][1]
   rho = fm.rtls_rho(np.asarray(f), sza, vza, raa)
   return {k: np.broadcast_to(np.asarray(r, dtype=np.float64), shape) for k, r in rho.items()}


def states(axes):
   lt = np.log10(axes["optical_depth"])
   re = axes["effective_radius"]
   nodes = [(a, b, "node") for a in lt for b in re]
   mids = [(0.5 * (lt[i] + lt[i + 1]), 0.5 * (re[j] + re[j + 1]), "midpoint") for i in range(lt.size - 1) for j in range(re.size - 1)]
   return nodes + mids


def evaluate(t, const, atm, surf_name):
   """Measurements y[state, channel-key, sza, vza, raa] and Jacobians; and per-operator sensitivities."""

   axes = t["axes"]
   xt, xr = np.log10(axes["optical_depth"]), axes["effective_radius"]
   sza = axes["solar_zenith"][:, None, None]
   vza = axes["satellite_zenith"][None, :, None]
   raa = axes["relative_azimuth"][None, None, :]
   out = {}
   # nadir channels only: the dual-view (SLSTR) oblique-view channels carry the
   # nadir operators copied unchanged by the generator (_replicate_dual_view)
   for ch in [c for c in t["all"] if c in CHANNELS_OF[CURRENT_CASE[0]]]:
      c = const[ch]
      is_solar = ch in t["solar"]
      is_thermal = ch in t["thermal"]
      kind = "mixed" if is_solar and is_thermal else ("solar" if is_solar else "thermal")
      rho = surface(surf_name, c["wavelength_um"], sza, vza, raa) if is_solar else None
      ys, js = [], []
      for lt, re, _ in states(axes):
         rad = 0.0
         drad = [0.0, 0.0]
         if is_solar:
            s = t["solar"].index(ch)
            crp, dcrp = {}, {}
            for name, table, shape in (("R_0v", t["R_0v"][..., s], None), ("T_00", t["T_00"][..., s][:, :, :, None, None], None),
                                       ("T_vv", t["T_00"][..., s][:, :, None, :, None], None),
                                       ("T_0d", t["T_0d"][..., s][:, :, :, None, None], None),
                                       ("T_dv", t["T_dv"][..., t["all"].index(ch)][:, :, None, :, None], None),
                                       ("R_dd", t["R_dd"][..., t["all"].index(ch)][:, :, None, None, None], None)):
               f, dt, dr = fm.bicubic(table, xt, xr, lt, re)
               crp[name], dcrp[name] = f, (dt, dr)
            ref, dref = fm.solar(crp, dcrp, rho, atm[ch]["Tac"], atm[ch]["Tbc"], sza, vza)
            ref = np.broadcast_to(ref, (sza.size, vza.size, raa.size))
            dref = [np.broadcast_to(d, ref.shape) for d in dref]
         if is_thermal:
            k = t["thermal"].index(ch)
            crp, dcrp = {}, {}
            for name, table in (("T_dv", t["T_dv"][..., t["all"].index(ch)]), ("R_dv", t["R_dv"][..., t["all"].index(ch)]),
                                ("E_md", t["E_md"][..., k])):
               f, dt, dr = fm.bicubic(table, xt, xr, lt, re)
               crp[name], dcrp[name] = f, (dt, dr)
            a = atm[ch]
            r, dr_ = fm.thermal_radiance(crp, dcrp, {"Rbc_up": a["Rbc_up"], "B_c": a["B_c"], "Rac_dwn": a["Rac_dwn"],
                                                    "Tac": a["Tac_lw"], "Rac_up": a["Rac_up"]})
            rad = np.broadcast_to(r[None, :, None], (sza.size, vza.size, raa.size))
            drad = [np.broadcast_to(d[None, :, None], rad.shape) for d in dr_]
         if kind == "solar":
            y, j = ref, dref
         else:
            total = rad + (c["F0"] * ref if kind == "mixed" else 0.0)
            dtot = [drad[i] + (c["F0"] * dref[i] if kind == "mixed" else 0.0) for i in range(2)]
            bt, dtdr = fm.r2t(total, c["planck"])
            y, j = bt, [dtdr * d for d in dtot]
         ys.append(np.array(y, dtype=np.float64))
         js.append(np.array(j, dtype=np.float64))
      out[ch] = {"kind": kind, "y": np.array(ys), "J": np.array(js)}       # y[state, sza, vza, raa], J[state, 2, ...]
   return out


def sigmas(res, const, production):
   s = {}
   for ch, r in res.items():
      c = const[ch]
      if r["kind"] == "solar":
         # ORAC sets Sy = 1e6 where the measurement is <= 0 (src/get_measurements.F90:270-271)
         s[ch] = np.where(r["y"] > 0.0, fm.sigma_solar(r["y"], c["noise_rel"], production), 1e3) if c["noise_rel"] else None
      else:
         s[ch] = fm.sigma_bt(r["y"], c["nedt"], c["refbt"], c["planck"], mixed=r["kind"] == "mixed", production=production)
   return s


def main():
   global RESULTS, REFERENCE
   # optional: python propagate.py [--reference NAME] [--cases a,b] [--out DIR]  (debugging / subsets)
   args = sys.argv[1:]
   only = None
   while args:
      key, value = args[0], args[1]
      args = args[2:]
      if key == "--reference":
         REFERENCE = value
      elif key == "--cases":
         only = value.split(",")
      elif key == "--out":
         RESULTS = Path(value)
   RESULTS.mkdir(parents=True, exist_ok=True)
   rows, jac_rows = [], []
   for case, (instfile, mmfile, grid, channels, nmom) in CASES.items():
      if (only and case not in only) or (not only and (case.endswith("_bs") or case.startswith("defect"))):
         continue
      base = lut(case, "baseline")
      ref = lut(case, REFERENCE)
      if base is None or ref is None:
         print(case, "skipped (baseline or reference missing)")
         continue
      platform, instrument, cloud = PLATFORM[case]
      CURRENT_CASE[0] = case
      const = channel_constants(base, instrument)
      planck = {ch: const[ch]["planck"] for ch in const if "planck" in const[ch]}
      axes = tables(base)["axes"]
      st = states(axes)
      state_kind = np.array([k for _, _, k in st])
      state_tau = 10 ** np.array([a for a, _, _ in st])
      state_re = np.array([b for _, b, _ in st])
      sza = axes["solar_zenith"][:, None, None]
      vza = axes["satellite_zenith"][None, :, None]
      raa = axes["relative_azimuth"][None, None, :]
      cos_t = -np.cos(np.deg2rad(sza)) * np.cos(np.deg2rad(vza)) - np.sin(np.deg2rad(sza)) * np.sin(np.deg2rad(vza)) * np.cos(np.deg2rad(raa))
      theta = np.rad2deg(np.arccos(np.clip(cos_t, -1, 1)))
      # thermal terms depend on the view zenith: evaluate per VZA node and assemble
      atm_by_vza = [terms(platform, instrument, channels, planck, cloud, float(v)) for v in axes["satellite_zenith"]]
      for surf in ("black", "dark", "vegetation", "snow"):
         results = {}
         for experiment in EXPERIMENTS:
            l = lut(case, experiment)
            if l is None:
               continue
            t = tables(l)
            per_vza = [evaluate(t, const, atm, surf) for atm in atm_by_vza]
            # take the VZA column iv from the evaluation made with that VZA's clear-sky terms
            res = {}
            for ch in per_vza[0]:
               y = np.stack([per_vza[iv][ch]["y"][:, :, iv, :] for iv in range(len(per_vza))], axis=2)
               j = np.stack([per_vza[iv][ch]["J"][:, :, :, iv, :] for iv in range(len(per_vza))], axis=3)
               res[ch] = {"kind": per_vza[0][ch]["kind"], "y": y, "J": j}
            results[experiment] = res
         reference = results[REFERENCE]
         # every experiment's measurements and Jacobians (float32) for the detailed analyses (scratch)
         full = TMP / "fm_full"
         full.mkdir(parents=True, exist_ok=True)
         np.savez_compressed(full / f"{case}_{surf}.npz",
                             **{f"{e}__{ch}__y": results[e][ch]["y"].astype(np.float32) for e in results for ch in results[e]},
                             **{f"{e}__{ch}__J": results[e][ch]["J"].astype(np.float32) for e in results for ch in results[e]},
                             state_tau=state_tau, state_re=state_re, state_kind=state_kind, theta=theta,
                             sza=axes["solar_zenith"], vza=axes["satellite_zenith"], raa=axes["relative_azimuth"],
                             kinds=np.array([f"{ch}:{results[REFERENCE][ch]['kind']}" for ch in results[REFERENCE]]))
         sig_p = sigmas(reference, const, True)
         sig_n = sigmas(reference, const, False)
         np.savez_compressed(RESULTS / f"fm_arrays_{case}_{surf}.npz",
                             **{f"{e}__{ch}__y": results[e][ch]["y"] for e in results for ch in results[e]
                                if e in ("baseline", REFERENCE, "nstr032", "nstr100", "nstr200", "double", "nocorint")},
                             **{f"{e}__{ch}__J": results[e][ch]["J"] for e in results for ch in results[e]
                                if e in ("baseline", REFERENCE, "nstr032", "nstr100", "nstr200", "double", "nocorint")},
                             state_tau=state_tau, state_re=state_re, state_kind=state_kind,
                             sza=axes["solar_zenith"], vza=axes["satellite_zenith"], raa=axes["relative_azimuth"])
         for experiment, res in results.items():
            for against, other in (("reference", reference), ("baseline", results.get("baseline"))):
               if other is None or (against == "reference" and experiment == REFERENCE) or \
                     (against == "baseline" and experiment == "baseline"):
                  continue
               for ch, r in res.items():
                  o = other[ch]
                  dy = np.abs(r["y"] - o["y"])
                  for label, sig in (("production_sigma", sig_p[ch]), ("noise_sigma", sig_n[ch])):
                     if sig is None:
                        continue
                     rr = dy / sig
                     solar_like = r["kind"] == "solar"
                     # mixed (3.7 um) channels carry the solar term and its backscatter peak too
                     glory = np.broadcast_to(theta >= 170.0, rr.shape) if r["kind"] != "thermal" else np.zeros(rr.shape, bool)
                     k = np.unravel_index(int(np.nanargmax(rr)), rr.shape)
                     rows.append({"case": case, "experiment": experiment, "against": against, "surface": surf, "channel": ch,
                                  "kind": r["kind"], "sigma": label, "max_abs_dy": float(np.nanmax(dy)), "n_nonfinite_excluded": int((~np.isfinite(rr)).sum()),
                                  "max_r": float(np.nanmax(rr)), "rms_r": float(np.sqrt(np.nanmean(rr**2))),
                                  "p99_r": float(np.nanpercentile(rr, 99)),
                                  "max_r_off_glory": float(np.nanmax(rr[~glory])) if (~glory).any() else None,
                                  "max_r_midpoints": float(np.nanmax(rr[state_kind == "midpoint"])),
                                  "max_r_signal_ge_0.01": float(np.nanmax(rr[o["y"] >= 0.01])) if solar_like and (o["y"] >= 0.01).any()
                                  else float(np.nanmax(rr)),
                                  "n_reference_nonpositive": int((o["y"] <= 0).sum()) if solar_like else 0,
                                  "worst_tau": float(state_tau[k[0]]), "worst_re": float(state_re[k[0]]),
                                  "worst_state": state_kind[k[0]], "worst_sza": float(axes["solar_zenith"][k[1]]),
                                  "worst_vza": float(axes["satellite_zenith"][k[2]]), "worst_raa": float(axes["relative_azimuth"][k[3]]),
                                  "worst_theta": float(theta[k[1], k[2], k[3]])})
                     if label == "production_sigma":
                        for ix, xname in ((0, "dlog10tau"), (1, "dre")):
                           dj = np.abs(r["J"][:, ix] - o["J"][:, ix])
                           jref = np.abs(o["J"][:, ix])
                           floor = 0.05 * np.nanmax(jref, axis=0, keepdims=True)
                           rel = dj / np.maximum(jref, np.maximum(floor, 1e-12))
                           kj = np.unravel_index(int(np.nanargmax(dj / sig)), dj.shape)
                           jac_rows.append({"case": case, "experiment": experiment, "against": against, "surface": surf,
                                            "channel": ch, "kind": r["kind"], "derivative": xname,
                                            "max_dJ_over_sigma": float(np.nanmax(dj / sig)),
                                            "rms_dJ_over_sigma": float(np.sqrt(np.nanmean((dj / sig) ** 2))),
                                            "max_rel_dJ": float(np.nanmax(rel)), "p99_rel_dJ": float(np.nanpercentile(rel, 99)),
                                            "worst_tau": float(state_tau[kj[0]]), "worst_re": float(state_re[kj[0]]),
                                            "worst_sza": float(axes["solar_zenith"][kj[1]]),
                                            "worst_vza": float(axes["satellite_zenith"][kj[2]]),
                                            "worst_raa": float(axes["relative_azimuth"][kj[3]])})
         print(case, surf, "done", flush=True)
   # Non-finite reference values (a double-precision NSTR 236 breakdown, see the report) are
   # excluded by the nan-aware statistics above and counted in n_nonfinite_excluded.
   for name, table in (("fm_measurement.csv", rows), ("fm_jacobians.csv", jac_rows)):
      with open(RESULTS / name, "w", newline="") as handle:
         writer = csv.DictWriter(handle, fieldnames=list(table[0]))
         writer.writeheader()
         writer.writerows(table)
   print("written", RESULTS)


if __name__ == "__main__":
   main()
