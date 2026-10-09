"""Analyse the compact-grid V24 / V25 runs of the profile-based formulation (validation only).

Inputs: validation/tmp/v25/compact/<case>/<mode>/ (run_compact.py) and the
finished V24 products.  Outputs, in validation/v25/results/:

   compact_v24_reproduction.csv   isothermal run against the V24 product at the matching
                                  Grid B nodes (every variable; bitwise or not)
   compact_operator_changes.csv   cirrostratus against isothermal: every LUT variable
   compact_e_md.csv               E_md isothermal / cirrostratus / constant-beta by channel,
                                  tau, re (view zenith 0; worst view)
   compact_layers.csv             E_md with 12, 24, 36, 40 emission layers against 46
   compact_ttop.csv               normalised E_md at reference cloud-top temperatures 220 and 260 K
                                  against 240 K (the diagnostic of the report's section 8)
   compact_toa.csv                top-of-atmosphere brightness-temperature difference of the complete
                                  ORAC thermal / mixed forward model, cirrostratus and constant-beta
                                  against V24, by channel, tau, re and view zenith
   compact_timing.csv             DISORT call counts and CPU times by mode
   profiles.csv, profile_grid_b.csv   the supplied profiles and their interpolation to Grid B
   figures/profile_supplied.png, e_md_toa.png, layers_ttop.png, nk_bt.png

   python validation/v25/compare_compact.py
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "v24" / "disort_options"))

from oraclut.io.lut import read_lut                              # noqa: E402
from oraclut import cloud_temperature as ct                      # noqa: E402
import orac_fm as fm                                             # noqa: E402
import atmosphere as atmo                                        # noqa: E402
import propagate as pg                                           # noqa: E402
import run_variant as rv                                         # noqa: E402

TMP = ROOT / "validation" / "tmp" / "v25" / "compact"
RESULTS = HERE / "results"
FIGURES = RESULTS / "figures"
LUT_DIR = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
CASES = {
   "modis_ice_agg": ("terra", "modis", "ice", "terra_modis_m_water-ice_a01_pagg_v24.nc", rv.MODIS_CHANNELS),
   "modis_ice_sph": ("terra", "modis", "ice", "terra_modis_m_water-ice_a01_psph_v24.nc", rv.MODIS_CHANNELS),
   "slstr_ice_sph": ("sentinel-3a", "slstr", "ice", "sentinel-3a_slstr_m_water-ice_a01_psph_v24.nc", rv.SLSTR_CHANNELS),
   "slstr_ice_agg": ("sentinel-3a", "slstr", "ice", "sentinel-3a_slstr_m_water-ice_a01_pagg_v24.nc", rv.SLSTR_CHANNELS),
   "modis_liquid": ("terra", "modis", "liquid", "terra_modis_m_liquid-water_a01_pstg_v24.nc", rv.MODIS_CHANNELS),
   "slstr_liquid": ("sentinel-3a", "slstr", "liquid", "sentinel-3a_slstr_m_liquid-water_a01_pstg_v24.nc", rv.SLSTR_CHANNELS),
}
COORDS = ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth")
PG = {"modis_ice_agg": "modis_ice_agg", "modis_ice_sph": "modis_ice_agg", "slstr_ice_sph": "slstr_ice_sph", "slstr_ice_agg": "slstr_ice_sph",
      "modis_liquid": "modis_liquid", "slstr_liquid": "slstr_liquid"}


def lut(case, mode):
   if mode == "v25":                                     # the V25 run of the case's backend
      for name in ("v25", "cirrostratus", "wet_adiabat"):
         files = sorted((TMP / case / name).glob("*.nc"))
         if files:
            return read_lut(files[0])
      return None
   files = sorted((TMP / case / mode).glob("*.nc"))
   return read_lut(files[0]) if files else None


def node_index(small, full):
   index = {}
   for name in COORDS:
      s = np.asarray(small.variables[name], dtype=np.float64)
      f = np.asarray(full.variables[name], dtype=np.float64)
      index[name] = np.array([int(np.argmin(np.abs(f - v))) for v in s])
      assert np.allclose(f[index[name]], s, rtol=1e-6, atol=1e-6), name
   for dim, var in (("channels", "channel_id"), ("solar_channels", "solar_channel_id"),
                    ("thermal_channels", "thermal_channel_id"), ("mixed_channels", "mixed_channel_id")):
      if var in small.variables:
         a = [int(c) for c in np.asarray(small.variables[var])]
         b = [int(c) for c in np.asarray(full.variables[var])]
         index[dim] = np.array([b.index(c) for c in a])
   return index


def compare_variables(a, b, index=None):
   rows = []
   for name in a.variable_names:
      x = np.asarray(a.variables[name])
      if x.dtype.kind != "f" or name in COORDS:
         continue
      dims = a.variable_dimensions[name]
      y = np.asarray(b.variables[name])
      if index is not None:
         if not all(d in index for d in dims):
            continue
         y = y[np.ix_(*[index[d] for d in dims])]
      d = np.abs(x.astype(np.float64) - y.astype(np.float64))
      rows.append({"variable": name, "max_abs_diff": float(np.nanmax(d)) if d.size else 0.0,
                   "bitwise": bool(np.array_equal(x, y, equal_nan=True)), "values": int(x.size)})
   return rows


def write(name, rows):
   if not rows:
      return
   with open(RESULTS / name, "w", newline="") as handle:
      w = csv.DictWriter(handle, fieldnames=list(rows[0]))
      w.writeheader()
      w.writerows(rows)


def thermal_arrays(l):
   v = lambda n: np.asarray(l.variables[n], dtype=np.float64)
   return {"thermal": [int(c) for c in v("thermal_channel_id")], "all": [int(c) for c in v("channel_id")],
           "solar": [int(c) for c in v("solar_channel_id")],
           "E_md": np.transpose(v("E_md"), (1, 2, 0, 3)), "T_dv": np.transpose(v("T_dv"), (1, 2, 0, 3)),
           "R_dv": np.transpose(v("R_dv"), (1, 2, 0, 3)), "R_0v": np.transpose(v("R_0v"), (3, 4, 2, 1, 0, 5)),
           "axes": {n: v(n) for n in COORDS}}


def profile_figures():
   model = ct.cloud_vertical_profile_model("cirrostratus", 256.0)
   p = model.profiles
   rows = []
   for k, cot in enumerate(p.cot):
      for i in range(99):
         rows.append({"cot": float(cot), "row": i, "extinction_percent": float(p.extinction_percent[i]), "dcot": float(p.dcot[k, i]),
                      "cumulative": float(p.cumulative[k, i]), "depth_km": float(p.depth_km[k, i]), "delta_t_K": float(p.delta_t_k[k, i])})
   write("profiles.csv", rows)
   grid = np.array([1e-10, 0.0078125, 0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 2.82842712474619, 4, 5.65685424949238, 8,
                    11.3137084989848, 16, 22.6274169979695, 32, 45.254833995939, 64, 90.509667991878, 128, 181.019335983756, 256])
   write("profile_grid_b.csv", [{"tau_055": float(t), "depth_km": p.depth_at(t), "delta_t_base_K": model.cloud_base_temperature_k(t) - 240.0,
                                 "depth_km_log_interpolation": float(np.interp(np.log(max(t, p.cot[0])), np.log(p.cot), p.depth_scale_km)),
                                 "native": bool(t in p.cot)} for t in grid])

   fig, axes = plt.subplots(2, 3, figsize=(17, 9))
   zeta = np.arange(99) / 98.0
   axes[0, 0].plot(p.extinction_percent, zeta, "k-")
   axes[0, 0].invert_yaxis()
   axes[0, 0].set_xlabel("% extinction per row (sums to 100; common to all 13 COTs)")
   axes[0, 0].set_ylabel("normalised depth below the cloud top z / H")
   axes[0, 0].set_title("supplied extinction profile", fontsize=10)
   axes[0, 1].plot(np.r_[0.0, p.optical_fraction_midpoints, 1.0], np.r_[0.0, zeta, 1.0], "k-")
   axes[0, 1].set_xlabel("cumulative optical-depth fraction F")
   axes[0, 1].set_ylabel("z / H")
   axes[0, 1].invert_yaxis()
   axes[0, 1].set_title("cumulative optical depth against normalised depth", fontsize=10)
   axes[0, 2].semilogx(p.cot, p.depth_scale_km, "ko", label="supplied H(COT)")
   fine = np.logspace(np.log10(0.0078125), np.log10(256), 300)
   axes[0, 2].semilogx(fine, [p.depth_at(t) for t in fine], "C0-", lw=1, label="linear in COT (used)")
   axes[0, 2].semilogx(fine, np.interp(np.log(np.maximum(fine, p.cot[0])), np.log(p.cot), p.depth_scale_km), "C1--", lw=1, label="linear in log COT")
   axes[0, 2].set_xlabel("total COT (0.55 µm)")
   axes[0, 2].set_ylabel("cloud depth H (km)")
   axes[0, 2].set_title("physical depth against total COT", fontsize=10)
   axes[0, 2].legend(fontsize=8)
   axes[1, 0].semilogx(p.cot, p.delta_t_k[:, -1], "ko")
   axes[1, 0].semilogx(fine, [model.cloud_base_temperature_k(t) - 240.0 for t in fine], "C0-", lw=1)
   axes[1, 0].set_xlabel("total COT (0.55 µm)")
   axes[1, 0].set_ylabel("ΔT at the cloud base (K)")
   axes[1, 0].set_title("base temperature departure against total COT", fontsize=10)
   for k in (0, 4, 8, 12):
      axes[1, 1].plot(p.delta_t_k[k], p.depth_km[k], label=f"COT {p.cot[k]:g}")
   axes[1, 1].invert_yaxis()
   axes[1, 1].set_xlabel("ΔT above the cloud top (K)")
   axes[1, 1].set_ylabel("depth below the cloud top (km)")
   axes[1, 1].set_title("ΔT = 8 K/km × depth", fontsize=10)
   axes[1, 1].legend(fontsize=8)
   for tau, colour in ((0.25, "C0"), (4.0, "C1"), (16.0, "C2"), (256.0, "C3")):
      f, dt = p.layer_boundaries(tau, ct.EMISSION_LAYERS)
      axes[1, 2].step(dt, f, where="post", color=colour, label=f"τ {tau:g}: {ct.EMISSION_LAYERS} layers")
      F = np.linspace(0, 1, 400)
      axes[1, 2].plot(p.delta_t_at(tau, F), F, color=colour, ls=":", lw=1)
   axes[1, 2].invert_yaxis()
   axes[1, 2].set_xlabel("ΔT at the layer boundaries (K)")
   axes[1, 2].set_ylabel("cumulative optical-depth fraction")
   axes[1, 2].set_title("remapped DISORT emission layers (steps) and the supplied relation (dotted)", fontsize=9)
   axes[1, 2].legend(fontsize=8)
   for ax in axes.ravel():
      ax.grid(True, alpha=0.3)
   fig.tight_layout()
   fig.savefig(FIGURES / "profile_supplied.png", dpi=140)
   return model


def main():
   RESULTS.mkdir(parents=True, exist_ok=True)
   FIGURES.mkdir(parents=True, exist_ok=True)
   model = profile_figures()
   liquid = ct.cloud_vertical_profile_model("wet_adiabat", 256.0)
   models = {"ice": model, "liquid": liquid}

   def depth_of(case, tau):
      m = models[CASES[case][2]]
      return m.profiles.depth_at(tau) if CASES[case][2] == "ice" else float(m.depth_of_optical_depth_km(tau))

   def base_of(case, tau):
      return models[CASES[case][2]].cloud_base_temperature_k(tau)
   repro, changes, emd_rows, layer_rows, ttop_rows, toa_rows, timing = [], [], [], [], [], [], []
   bt_curves = {}
   for case, (platform, instrument, cloud, product, nadir) in CASES.items():
      iso, new = lut(case, "isothermal"), lut(case, "v25")
      if iso is None or new is None:
         print(case, "missing runs")
         continue
      beta = lut(case, "constantbeta")
      full = read_lut(LUT_DIR / product)
      for r in compare_variables(iso, full, node_index(iso, full)):
         repro.append({"case": case, "product": product, **r})
      for r in compare_variables(new, iso):
         changes.append({"case": case, **r})
      ref46 = lut(case, "layers46")
      for n in (12, 24, 36, 40):
         l = lut(case, "v25" if n == ct.EMISSION_LAYERS else f"layers{n}")
         if l is None or ref46 is None:
            continue
         e, er = np.asarray(l.variables["E_md"], float), np.asarray(ref46.variables["E_md"], float)
         ch = [int(c) for c in np.asarray(l.variables["thermal_channel_id"])]
         tau = np.asarray(l.variables["optical_depth"])
         for k, c in enumerate(ch):
            if c not in nadir:
               continue
            d = np.abs(e[..., k] - er[..., k])
            i = np.unravel_index(np.argmax(d), d.shape)
            layer_rows.append({"case": case, "channel": c, "layers": n, "max_abs_dE_md_vs_46": float(d.max()),
                               "tau_of_max": float(tau[i[1]]), "E_md_46_there": float(er[..., k][i]), "rms": float(np.sqrt((d ** 2).mean()))})
      for t_top in (220, 260):
         l = lut(case, f"ttop{t_top}")
         if l is None:
            continue
         e, e240 = np.asarray(l.variables["E_md"], float), np.asarray(new.variables["E_md"], float)
         ch = [int(c) for c in np.asarray(l.variables["thermal_channel_id"])]
         tau = np.asarray(l.variables["optical_depth"])
         for k, c in enumerate(ch):
            if c not in nadir:
               continue
            d = e[..., k] - e240[..., k]
            i = np.unravel_index(np.argmax(np.abs(d)), d.shape)
            ttop_rows.append({"case": case, "channel": c, "t_top_K": t_top, "max_abs_dE_md_vs_240": float(np.abs(d).max()),
                              "tau_of_max": float(tau[i[1]]), "E_md_240_there": float(e240[..., k][i]), "E_md_there": float(e[..., k][i]),
                              "max_rel_dE_md": float(np.max(np.abs(d) / np.maximum(e240[..., k], 1e-3)))})
      for mode in ("isothermal", "cirrostratus", "v25", "constantbeta", "layers12", "layers24", "layers36", "layers46", "ttop220", "ttop260"):
         rec = TMP / case / mode / "record.json"
         if rec.exists():
            timing.append({"case": case, **json.loads(rec.read_text())})

      a, b = thermal_arrays(iso), thermal_arrays(new)
      c3 = thermal_arrays(beta) if beta is not None else None
      ax = a["axes"]
      const = pg.channel_constants(iso, instrument)
      planck = {ch: const[ch]["planck"] for ch in const if "planck" in const[ch]}
      for k, ch in enumerate(a["thermal"]):
         if ch not in nadir:
            continue
         kind = "mixed" if ch in a["solar"] else "thermal"
         for i, tau in enumerate(ax["optical_depth"]):
            for j, re in enumerate(ax["effective_radius"]):
               e0, e1 = a["E_md"][i, j, :, k], b["E_md"][i, j, :, k]
               row = {"case": case, "channel": ch, "kind": kind, "wavelength_um": const[ch]["wavelength_um"], "tau": float(tau),
                      "re": float(re), "depth_km": depth_of(case, tau), "dT_base_K": base_of(case, tau) - 240.0,
                      "E_md_v24_vza0": float(e0[0]), "E_md_v25_vza0": float(e1[0]), "dE_md_vza0": float(e1[0] - e0[0]),
                      "ratio_vza0": float(e1[0] / e0[0]) if e0[0] > 0 else None, "max_abs_dE_md_views": float(np.max(np.abs(e1 - e0))),
                      "max_E_md_v25": float(e1.max()), "E_md_v25_vza75": float(e1[-1])}
               if c3 is not None:
                  row["E_md_constantbeta_vza0"] = float(c3["E_md"][i, j, 0, k])
               emd_rows.append(row)
      pg.CURRENT_CASE[0] = PG[case]
      for iv, vza in enumerate(ax["satellite_zenith"]):
         terms = atmo.terms(platform, instrument, [c for c in a["thermal"] if c in nadir], planck, cloud, float(vza))
         for k, ch in enumerate(a["thermal"]):
            if ch not in nadir:
               continue
            kind = "mixed" if ch in a["solar"] else "thermal"
            ci = a["all"].index(ch)
            t = dict(terms[ch])
            t["B_c"] = float(fm.t2r(np.asarray(ct.T_TOP_K), planck[ch])[0])
            for i, tau in enumerate(ax["optical_depth"]):
               for j, re in enumerate(ax["effective_radius"]):
                  bts = []
                  for arr in ((a, b, c3) if c3 is not None else (a, b)):
                     crp = {"E_md": arr["E_md"][i, j, iv, k], "T_dv": arr["T_dv"][i, j, iv, ci], "R_dv": arr["R_dv"][i, j, iv, ci]}
                     rad, _ = fm.thermal_radiance(crp, {n: (0.0, 0.0) for n in crp}, t)
                     if kind == "mixed":
                        si = int(np.argmin(np.abs(ax["solar_zenith"] - 30.0)))
                        pi_ = int(np.argmin(np.abs(ax["relative_azimuth"] - 90.0)))
                        sci = arr["solar"].index(ch)
                        ref = arr["R_0v"][i, j, si, iv, pi_, sci] * t["Tac"] ** (1.0 / np.cos(np.deg2rad(30.0)) + 1.0 / max(np.cos(np.deg2rad(vza)), 1e-3))
                        rad = rad + const[ch]["F0"] * ref
                     bts.append(float(fm.r2t(np.asarray([rad]), planck[ch])[0][0]))
                  sigma = float(fm.sigma_bt(np.asarray([bts[0]]), const[ch]["nedt"], const[ch]["refbt"], planck[ch], mixed=kind == "mixed", production=True)[0])
                  row = {"case": case, "channel": ch, "kind": kind, "tau": float(tau), "re": float(re), "vza": float(vza),
                         "BT_v24_K": bts[0], "BT_v25_K": bts[1], "dBT_K": bts[1] - bts[0], "sigma_production_K": sigma,
                         "r": abs(bts[1] - bts[0]) / sigma}
                  if len(bts) == 3:
                     row["BT_constantbeta_K"] = bts[2]
                     row["dBT_constantbeta_K"] = bts[2] - bts[0]
                  toa_rows.append(row)
                  if iv == 0:
                     bt_curves.setdefault((case, ch), {})[(float(tau), float(re))] = bts[:2]
   write("compact_v24_reproduction.csv", repro)
   write("compact_operator_changes.csv", changes)
   write("compact_e_md.csv", emd_rows)
   write("compact_layers.csv", layer_rows)
   write("compact_ttop.csv", ttop_rows)
   write("compact_toa.csv", toa_rows)
   if timing:
      write("compact_timing.csv", [{"case": r["case"], "mode": r["mode"], "emission_layers": r["emission_layers"], "total_cpu_s": r["total_cpu_s"],
                                    **{f"calls_{k}": v for k, v in r["disort_calls"].items()},
                                    **{f"cpu_{k}_s": v for k, v in r["disort_cpu_s"].items()}} for r in timing])

   cases = [c for c in CASES if any(r["case"] == c for r in emd_rows)]
   if not cases:
      print("no complete case")
      return
   fig, axes = plt.subplots(2, len(cases), figsize=(4.5 * len(cases), 8), squeeze=False)
   for n, case in enumerate(cases):
      rows = [r for r in emd_rows if r["case"] == case]
      for ch in sorted({r["channel"] for r in rows}):
         sel = sorted((r["tau"], r["re"], r["dE_md_vza0"]) for r in rows if r["channel"] == ch)
         res = sorted({s[1] for s in sel})
         re_mid = res[len(res) // 2]
         pts = [s for s in sel if s[1] == re_mid]
         wl = next(r["wavelength_um"] for r in rows if r["channel"] == ch)
         axes[0, n].semilogx([p[0] for p in pts], [p[2] for p in pts], "o-", ms=3, label=f"ch {ch} ({wl:.2f} µm), r_e {re_mid:g}")
      axes[0, n].set_title(f"{case}: E_md(V25 profile) − E_md(V24), VZA 0", fontsize=9)
      axes[0, n].set_xlabel("τ (0.55 µm)")
      axes[0, n].grid(True, alpha=0.3)
      axes[0, n].legend(fontsize=6)
      rows = [r for r in toa_rows if r["case"] == case and r["vza"] == 0.0]
      for ch in sorted({r["channel"] for r in rows}):
         sel = sorted((r["tau"], r["re"], r["dBT_K"], r.get("dBT_constantbeta_K")) for r in rows if r["channel"] == ch)
         res = sorted({s[1] for s in sel})
         re_mid = res[len(res) // 2]
         pts = [s for s in sel if s[1] == re_mid]
         line, = axes[1, n].semilogx([p[0] for p in pts], [p[2] for p in pts], "o-", ms=3, label=f"ch {ch}, r_e {re_mid:g} (profile)")
         if pts[0][3] is not None:
            axes[1, n].semilogx([p[0] for p in pts], [p[3] for p in pts], ":", color=line.get_color(), lw=1, label=f"ch {ch} (constant β, superseded)")
      axes[1, n].set_title(f"{case}: TOA BT − BT(V24), VZA 0", fontsize=9)
      axes[1, n].set_xlabel("τ (0.55 µm)")
      axes[1, n].set_ylabel("ΔBT (K)")
      axes[1, n].grid(True, alpha=0.3)
      axes[1, n].legend(fontsize=6)
   fig.tight_layout()
   fig.savefig(FIGURES / "e_md_toa.png", dpi=140)

   fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
   for n, case in enumerate(cases):
      for ch in sorted({r["channel"] for r in layer_rows if r["case"] == case}):
         pts = sorted((r["layers"], r["max_abs_dE_md_vs_46"]) for r in layer_rows if r["case"] == case and r["channel"] == ch)
         if pts:
            axes[0].semilogy([p[0] for p in pts], [max(p[1], 1e-7) for p in pts], "o-", ms=3, color=f"C{n}", alpha=0.8, label=f"{case} ch {ch}")
   axes[0].set_xlabel("emission layers")
   axes[0].set_ylabel("max |E_md − E_md(46 layers)|")
   axes[0].set_title("layer-count convergence", fontsize=10)
   axes[0].legend(fontsize=5, ncol=2)
   for n, case in enumerate(cases):
      chans = sorted({r["channel"] for r in ttop_rows if r["case"] == case})
      for ch in chans:
         for t_top, marker in ((220, "v"), (260, "^")):
            pts = [r for r in ttop_rows if r["case"] == case and r["channel"] == ch and r["t_top_K"] == t_top]
            if pts:
               axes[1].plot([ch], [pts[0]["max_abs_dE_md_vs_240"]], marker, color=f"C{n}", label=f"{case} T_top {t_top} K" if ch == chans[0] else None)
   axes[1].set_xlabel("channel")
   axes[1].set_ylabel("max |E_md(T_top) − E_md(240 K)|")
   axes[1].set_title("reference cloud-top temperature diagnostic", fontsize=10)
   axes[1].legend(fontsize=6)
   for ax in axes:
      ax.grid(True, alpha=0.3)
   fig.tight_layout()
   fig.savefig(FIGURES / "layers_ttop.png", dpi=140)

   pairs = {"modis": [(31, 32), (20, 31)], "slstr": [(8, 9), (7, 8)]}
   fig, axes = plt.subplots(2, len(cases), figsize=(4.5 * len(cases), 9), squeeze=False)
   for n, case in enumerate(cases):
      for p_, (cx, cy) in enumerate(pairs[CASES[case][1]]):
         ax_ = axes[p_, n]
         gx, gy = bt_curves.get((case, cx)), bt_curves.get((case, cy))
         if not gx or not gy:
            continue
         taus = sorted({k[0] for k in gx})
         res = sorted({k[1] for k in gx})
         for v, (ls, lab) in enumerate((("-", "V24"), ("--", "V25 profile"))):
            for tau in taus:
               ax_.plot([gx[(tau, re)][v] for re in res], [gy[(tau, re)][v] for re in res], ls, color="C0", lw=1)
            for re in res:
               ax_.plot([gx[(tau, re)][v] for tau in taus], [gy[(tau, re)][v] for tau in taus], ls, color="C3", lw=1, label=lab if re == res[0] else None)
         for tau in taus:
            ax_.annotate(f"τ={tau:g}", (gx[(tau, res[-1])][0], gy[(tau, res[-1])][0]), fontsize=6, color="C0")
         ax_.set_xlabel(f"channel {cx} BT (K)")
         ax_.set_ylabel(f"channel {cy} BT (K)")
         ax_.set_title(f"{case}: VZA 0", fontsize=9)
         ax_.grid(True, alpha=0.3)
         ax_.legend(fontsize=7)
   fig.tight_layout()
   fig.savefig(FIGURES / "nk_bt.png", dpi=140)
   print("written", RESULTS)


if __name__ == "__main__":
   main()
