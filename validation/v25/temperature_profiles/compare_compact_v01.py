"""Compare the compact V25 variants with the sealed V24 release and with historical V25 (read-only).

   python validation/v25/temperature_profiles/compare_compact.py   -> results/compact_<model>.json, results/summary.md

For each model (water-ice_agg, liquid-water_240, liquid-water_old) and
variant (compact_lut.py) against the sealed V24 product (work/<model>_v24):

* identity: the variables equal to V24 bit for bit and those that differ;
  the isothermal variant must equal V24 in every variable, the others in
  every variable but E_md;
* isothermal limit (constant): |E_md - E_md(V24)|;
* profile: E_md - E_md(V24) by thermal channel (20, 31, 32), its sign (a
  cloud warmer below its top emits more: >= 0), the physical bound
  E_md B(250 K) <= B(T_base(tau)), the emitted brightness temperature at
  nadir, the number of interior extrema of E_md along tau (smoothness);
* mirror (cooling with depth): E_md - E_md(V24) <= 0;
* layers20 / layers46: |E_md - E_md(profile)|;
* ttop240: the historical 240 K cloud top: |E_md - E_md(profile)|, and
  against historical V25: the full Terra MODIS V25 products of 2026-10-06
  (revision 990d970, shared ORAC LUT directory) at the compact nodes, and for
  water-ice_agg the historical compact results
  (validation/v25/results/compact_e_md.csv: historical V24 and V25 E_md at
  nadir for r_e 5, 33, 93 um).
A difference of at most NOISE (float32 rounding of E_md ~ 1) counts as zero.
Statistics are given for tau >= THICK (the profile statistics) and, for the
optically thin nodes tau < THICK, separately (thin_cloud: the V25 / V24 ratio,
whose optically thin limit is the optical-depth mean of B(T) / B(250 K)).
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.idl_mirror import cloud_temperature as ct   # noqa: E402
from oraclut.radiative_transfer.legacy_disort import plkavg   # noqa: E402

WORK = HERE / "work"
SHARED = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
PRODUCT = {"water-ice_agg": "terra_modis_m_water-ice_a01_pagg", "liquid-water_240": "terra_modis_m_liquid-water_a01_p240",
           "liquid-water_old": "terra_modis_m_liquid-water_a01_pold"}
PROFILE = {"water-ice_agg": "cirrostratus", "liquid-water_240": "wet_adiabat", "liquid-water_old": "wet_adiabat"}
VARIANTS = ("isothermal", "profile", "constant", "mirror", "ttop240", "layers20", "layers46")
NOISE = 1e-6
THICK = 0.0625
f32 = np.float32


def read(path):
   with Dataset(path) as d:
      d.set_auto_maskandscale(False)
      data = {name: np.asarray(v[:]) for name, v in d.variables.items()}
      data["E_md valid_range"] = np.asarray(d["E_md"].getncattr("valid_range"))
   return data


def window(wavelength):
   wn = f32(1e4) / f32(wavelength)
   return float(f32(0.995) * wn), float(f32(1.005) * wn)


def brightness_temperature(radiance, wnlo, wnhi):
   lo, hi = 100.0, 500.0
   for _ in range(60):
      mid = 0.5 * (lo + hi)
      lo, hi = (mid, hi) if plkavg(wnlo, wnhi, mid) < radiance else (lo, mid)
   return 0.5 * (lo + hi)


def interior_extrema(curve):
   slope = np.diff(curve)
   slope = slope[np.abs(slope) > NOISE]
   return int(np.count_nonzero(np.diff(np.sign(slope)) != 0))


def stats(delta):
   i = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
   return {"max_abs": float(np.abs(delta).max()), "min": float(delta.min()), "max": float(delta.max()), "at": [int(k) for k in i]}


summary = []
for model, stem in PRODUCT.items():
   v24_path = WORK / f"{model}_v24" / f"{stem}_v24.nc"
   if not v24_path.is_file():
      print(f"{model}: no V24 product yet")
      continue
   v24 = read(v24_path)
   tau, re, vza = v24["optical_depth"].astype(np.float64), v24["effective_radius"], v24["satellite_zenith"]
   thermal = np.flatnonzero(v24["thermal_channel_flag"])
   channels = [int(c) for c in v24["channel_id"][thermal]]
   wavelength = [float(w) for w in v24["central_wavelength"][thermal]]
   record = {"model": model, "thermal_channels": channels, "variants": {}}
   e24 = v24["E_md"].astype(np.float64)                       # (vza, tau, re, thermal channel)
   cloudprof = ct.setup_cloud_profile(PROFILE[model], 256.0, 250.0)
   e = {"v24": e24}
   for variant in VARIANTS:
      path = WORK / f"{model}_{variant}" / f"{stem}_v25.nc"
      if not path.is_file():
         continue
      v25 = read(path)
      names = [n for n in v24 if n != "E_md valid_range"]
      differing = sorted(n for n in names if not np.array_equal(v24[n], v25[n]))
      e[variant] = v25["E_md"].astype(np.float64)
      timing = json.loads((path.parent / "timing.json").read_text())[-1]
      record["variants"][variant] = {"differing_variables": differing, "E_md_valid_range": v25["E_md valid_range"].tolist(),
                                     "seconds": timing["seconds"], "emission_calls": timing["emission_calls"],
                                     "emission_seconds": timing["emission_seconds"]}
   v = record["variants"]
   ok_identity = ("isothermal" not in v or v["isothermal"]["differing_variables"] == []) and \
                 all(v[k]["differing_variables"] in ([], ["E_md"]) for k in v if k != "isothermal")
   record["identity_ok"] = ok_identity
   if "constant" in e:
      d = e["constant"] - e24
      big = np.abs(e24) > 1e-3
      record["isothermal_limit"] = {**stats(d), "max_rel_where_E_md_above_1e-3": float(np.max(np.abs(d[big]) / np.abs(e24[big])))}
   thick = tau >= THICK
   if "profile" in e:
      ratio = np.where(e24 != 0.0, e["profile"] / np.where(e24 != 0.0, e24, 1.0), np.nan)
      fraction = (np.arange(4000) + 0.5) / 4000
      record["thin_cloud"] = {}
      for k, (channel, wvl) in enumerate(zip(channels, wavelength)):
         wnlo, wnhi = window(wvl)
         rows = []
         for t in np.flatnonzero(~thick):
            limit = np.mean([plkavg(wnlo, wnhi, float(x)) for x in ct.cloud_temperature(cloudprof, float(tau[t]), fraction[::40])]) / plkavg(wnlo, wnhi, 250.0)
            rows.append({"tau": float(tau[t]), "E_md_v24_range": [float(e24[:, t, :, k].min()), float(e24[:, t, :, k].max())],
                         "E_md_v25_range": [float(e["profile"][:, t, :, k].min()), float(e["profile"][:, t, :, k].max())],
                         "ratio_range": [float(np.nanmin(ratio[:, t, :, k])), float(np.nanmax(ratio[:, t, :, k]))] if t > 0 else None,
                         "thin_limit_ratio": float(limit)})
         record["thin_cloud"][str(channel)] = rows
      d = (e["profile"] - e24)[:, thick]
      per_channel = {}
      for k, (channel, wvl) in enumerate(zip(channels, wavelength)):
         wnlo, wnhi = window(wvl)
         b250 = plkavg(wnlo, wnhi, 250.0)
         t_base = np.asarray([ct.cloud_temperature(cloudprof, float(t), [1.0])[0] for t in tau])
         bound = np.asarray([plkavg(wnlo, wnhi, float(t)) for t in t_base]) / b250
         margin = (bound[None, :, None] - e["profile"][..., k])[:, thick]
         bt = {f"r_e {float(r):g}": [round(brightness_temperature(float(x) * b250, wnlo, wnhi), 2) for x in e["profile"][0, :, j, k]]
               for j, r in enumerate(re)}
         bt24 = {f"r_e {float(r):g}": [round(brightness_temperature(float(x) * b250, wnlo, wnhi), 2) for x in e24[0, :, j, k]]
                 for j, r in enumerate(re)}
         extrema = [interior_extrema(e["profile"][i, thick, j, k]) for i in range(vza.size) for j in range(re.size)]
         extrema24 = [interior_extrema(e24[i, thick, j, k]) for i in range(vza.size) for j in range(re.size)]
         per_channel[str(channel)] = {"wavelength_um": wvl, "dE": stats(d[..., k]), "E_md_v25_max": float(e["profile"][:, thick, :, k].max()),
                                      "finite": bool(np.all(np.isfinite(e["profile"][..., k]))),
                                      "negative_dE_beyond_noise": int(np.count_nonzero(d[..., k] < -NOISE)),
                                      "bound_margin_min": float(margin.min()),
                                      "nadir_brightness_temperature_K": bt, "nadir_brightness_temperature_v24_K": bt24,
                                      "interior_extrema_along_tau_max": max(extrema), "interior_extrema_v24_max": max(extrema24),
                                      "dE_nadir_along_tau_mid_re": [float(x) for x in d[0, :, re.size // 2, k]]}
      record["profile"] = per_channel
   for variant in ("mirror",):
      if variant in e:
         d = (e[variant] - e24)[:, thick]
         record[variant] = {**stats(d), "positive_beyond_noise": int(np.count_nonzero(d > NOISE))}
   for variant in ("layers20", "layers46", "ttop240"):
      if variant in e and "profile" in e:
         d = (e[variant] - e["profile"])[:, thick]
         record[f"{variant}_minus_profile"] = {str(c): stats(d[..., k]) for k, c in enumerate(channels)}
   # historical V25 (990d970, 240 K cloud top): the full Terra products at the compact nodes
   shared = SHARED / f"{stem}_v25.nc"
   if "ttop240" in e and shared.is_file():
      with Dataset(shared) as h:
         h.set_auto_maskandscale(False)
         index = [[int(np.argmin(np.abs(np.asarray(h[name][:], np.float64) - x))) for x in grid]
                  for name, grid in (("satellite_zenith", vza), ("optical_depth", tau), ("effective_radius", re))]
         hch = list(np.asarray(h["thermal_channel_id"][:]))
         assert all(np.allclose(np.asarray(h[name][:])[i], grid) for name, grid, i in
                    (("satellite_zenith", vza, index[0]), ("optical_depth", tau, index[1]), ("effective_radius", re, index[2])))
         eh = np.asarray(h["E_md"][:], np.float64)[np.ix_(index[0], index[1], index[2], [hch.index(c) for c in channels])]
      record["historical_v25_products"] = {
         "file": str(shared),
         "ttop240_minus_historical": {str(c): stats((e["ttop240"] - eh)[:, thick, :, k]) for k, c in enumerate(channels)},
         "ttop240_minus_historical_thin": {str(c): stats((e["ttop240"] - eh)[:, ~thick, :, k]) for k, c in enumerate(channels)},
         "profile_minus_historical": {str(c): stats((e["profile"] - eh)[:, thick, :, k]) for k, c in enumerate(channels)},
         "v24_minus_historical_v25": {str(c): stats((e24 - eh)[:, thick, :, k]) for k, c in enumerate(channels)},
         "historical_E_md_at_tau_min": {str(c): [float(eh[:, 0, :, k].min()), float(eh[:, 0, :, k].max())] for k, c in enumerate(channels)},
         "historical_minus_v24_nadir_mid_re": {str(c): [float(x) for x in (eh - e24)[0, :, re.size // 2, k]] for k, c in enumerate(channels)},
         "ttop240_minus_v24_nadir_mid_re": {str(c): [float(x) for x in (e["ttop240"] - e24)[0, :, re.size // 2, k]] for k, c in enumerate(channels)},
      }
   # historical compact results (water-ice_agg: Terra MODIS Baum aggregate)
   if model == "water-ice_agg" and "ttop240" in e:
      rows = [r for r in csv.DictReader((ROOT / "validation/v25/results/compact_e_md.csv").open()) if r["case"] == "modis_ice_agg"]
      table = []
      for r in rows:
         t, rr, c = float(r["tau"]), float(r["re"]), int(r["channel"])
         if not (np.any(np.isclose(tau, t)) and np.any(np.isclose(re, rr)) and c in channels):
            continue
         i, j, k = int(np.argmin(np.abs(tau - t))), int(np.argmin(np.abs(re - rr))), channels.index(c)
         table.append({"channel": c, "tau": t, "re": rr,
                       "historical_v24": float(r["E_md_v24_vza0"]), "historical_v25": float(r["E_md_v25_vza0"]),
                       "v24": e24[0, i, j, k], "ttop240": e["ttop240"][0, i, j, k], "profile": e["profile"][0, i, j, k]})
      dh = np.asarray([x["historical_v25"] - x["historical_v24"] for x in table])
      dn = np.asarray([x["ttop240"] - x["v24"] for x in table])
      record["historical_compact"] = {
         "nodes": len(table),
         "profile_effect_ttop240_minus_historical_max_abs": float(np.abs(dn - dh).max()),
         "profile_effect_relative_max": float(np.max(np.abs(dn - dh) / np.maximum(np.abs(dh), 1e-3))),
         "isothermal_v24_minus_historical_v24_max_abs": float(np.max(np.abs([x["v24"] - x["historical_v24"] for x in table]))),
         "ttop240_minus_historical_v25_max_abs": float(np.max(np.abs([x["ttop240"] - x["historical_v25"] for x in table]))),
         "rows": table}
   (HERE / "results" / f"compact_{model}.json").write_text(json.dumps(record, indent=1, default=float) + "\n")
   summary.append(record)
   print(f"== {model}: identity {'OK' if ok_identity else 'FAILED'}; variants {sorted(v)}")
   for k, x in v.items():
      print(f"   {k:10s} differing {x['differing_variables']}; {x['seconds']:.0f} s, {x['emission_calls']} emission calls "
            f"{x['emission_seconds']:.1f} s")
   if "isothermal_limit" in record:
      print(f"   isothermal limit: max |dE| {record['isothermal_limit']['max_abs']:.2e} (rel {record['isothermal_limit']['max_rel_where_E_md_above_1e-3']:.2e})")
   for c, rows in record.get("thin_cloud", {}).items():
      print(f"   thin ch {c}: " + "; ".join(f"tau {r['tau']:.3g} V25 {r['E_md_v25_range'][0]:+.3g}..{r['E_md_v25_range'][1]:+.3g}"
                                          + (f" ratio {r['ratio_range'][0]:.3f}..{r['ratio_range'][1]:.3f} (limit {r['thin_limit_ratio']:.3f})" if r['ratio_range'] else "")
                                          for r in rows))
   for c, x in record.get("profile", {}).items():
      print(f"   ch {c}: dE {x['dE']['min']:+.4f}..{x['dE']['max']:+.4f}, V25 max {x['E_md_v25_max']:.3f}, dE<0: "
            f"{x['negative_dE_beyond_noise']}, bound margin {x['bound_margin_min']:+.3f}, extrema {x['interior_extrema_along_tau_max']} "
            f"(V24 {x['interior_extrema_v24_max']})")
   if "mirror" in record:
      print(f"   mirror: dE {record['mirror']['min']:+.4f}..{record['mirror']['max']:+.2e}, positive {record['mirror']['positive_beyond_noise']}")
   for variant in ("layers20", "layers46", "ttop240"):
      if f"{variant}_minus_profile" in record:
         print(f"   {variant} - profile: " + ", ".join(f"ch {c} {s['max_abs']:.2e}" for c, s in record[f"{variant}_minus_profile"].items()))
   if "historical_v25_products" in record:
      h = record["historical_v25_products"]
      print("   historical V25 products (tau >= 0.0625): ttop240 - hist " + ", ".join(f"ch {c} {s['max_abs']:.3e}" for c, s in h["ttop240_minus_historical"].items())
            + "; thin: " + ", ".join(f"ch {c} {s['max_abs']:.3e}" for c, s in h["ttop240_minus_historical_thin"].items())
            + "; historical E_md at tau 1e-10: " + ", ".join(f"ch {c} {v[0]:+.3f}..{v[1]:+.3f}" for c, v in h["historical_E_md_at_tau_min"].items())
            + "; V24 - hist V25 " + ", ".join(f"ch {c} {s['max_abs']:.3e}" for c, s in h["v24_minus_historical_v25"].items()))
   if "historical_compact" in record:
      h = record["historical_compact"]
      print(f"   historical compact ({h['nodes']} nodes): profile effect ttop240 vs historical {h['profile_effect_ttop240_minus_historical_max_abs']:.2e} "
            f"(rel {h['profile_effect_relative_max']:.2e}); V24 vs historical V24 {h['isothermal_v24_minus_historical_v24_max_abs']:.2e}; "
            f"ttop240 vs historical V25 {h['ttop240_minus_historical_v25_max_abs']:.2e}")
