"""Markdown tables of the V24 compact validation results (results/*.json -> results/summary.md).

   python validation/v24/scattering_refinements/summarize.py
"""

import glob
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE / "results"
load = lambda pattern: [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(R / pattern)))]
lines = []

lines += ["### Radius-integration and Legendre convergence (convergence.py)", "",
          "| Case (model, um, r_e) | nodes V24 / reference / V23 | V24 vs 4x finer: ext, g, phase 5-150, 150-180 | V23 grid vs 4x finer: ext, g, phase 150-180 | 3.5 vs 5 r_e: ext | L (2 Nq) | max dchi at 2 Nq | reconstruction | noise |",
          "|---|---|---|---|---|---|---|---|---|"]
for d in load("convergence_*.json"):
   p, v, t, L, n = d["production_vs_reference"], d["v23_grid_vs_reference"], d.get("limit_3.5_vs_5.0"), d["legendre"], d["nodes"]
   lines.append(f"| {d['model'].replace('liquid-water_', 'liquid ').replace('water-ice_sph', 'ice sph')} {d['wavelength']} {d['effective_radius']:g} "
                f"| {n['production'][0]} / {n['reference'][0]} / {n['v23_legacy'][0]} "
                f"| {p['extinction_rel']:.1e}, {p['g_abs']:.1e}, {p['phase_rel_5_150']:.1e}, {p['phase_rel_150_180']:.1e} "
                f"| {v['extinction_rel']:.1e}, {v['g_abs']:.1e}, {v['phase_rel_150_180']:.1e} "
                f"| {t['extinction_rel']:.1e} | {L['L']} ({L['L_at_2Nq']}) | {L['max_chi_difference_2Nq']:.0e} | {L['reconstruction_error']:.0e} | {L['noise']:.0e} |"
                if t else
                f"| {d['model'].replace('liquid-water_', 'liquid ').replace('water-ice_sph', 'ice sph')} {d['wavelength']} {d['effective_radius']:g} "
                f"| {n['production'][0]} / {n['reference'][0]} / {n['v23_legacy'][0]} "
                f"| {p['extinction_rel']:.1e}, {p['g_abs']:.1e}, {p['phase_rel_5_150']:.1e}, {p['phase_rel_150_180']:.1e} "
                f"| {v['extinction_rel']:.1e}, {v['g_abs']:.1e}, {v['phase_rel_150_180']:.1e} "
                f"| - | {L['L']} ({L['L_at_2Nq']}) | {L['max_chi_difference_2Nq']:.0e} | {L['reconstruction_error']:.0e} | {L['noise']:.0e} |")

lines += ["", "### DISORT with NMOM = max(1000, L) (disort_nmom.py; original Intel kernel, NSTR 60)", "",
          "| Case (um, r_e) | L | max(1000,L) = L bitwise | L+512 padded bitwise | max(L,61) bitwise (L<1000) | fluxes bitwise under truncation | max dR at N = 1000 / 1500 / 2000 / 2500 | s/call at 1000 / production |",
          "|---|---|---|---|---|---|---|---|"]
for d in load("disort_nmom_*.json"):
   cases = d["cases"].values()
   tr = lambda n: max((c["truncated_to_N"].get(str(n), {}).get("max_abs_R_difference", float("nan")) for c in cases), default=float("nan"))
   lines.append(f"| {d['wavelength']} {d['effective_radius']:g} | {d['L']} | {all(c['production_max(1000,L)_bitwise_equal_L'] for c in cases)} "
                f"| {all(c['zero_padded_L+512_bitwise_equal_L'] for c in cases)} "
                f"| {all(c.get('historical_max(L,61)_bitwise_equal_production', True) for c in cases) if d['L'] < 1000 else 'n/a'} "
                f"| {all(v['fluxes_bitwise_equal'] for c in cases for v in c['truncated_to_N'].values())} "
                f"| {' / '.join('-' if tr(n) != tr(n) else f'{tr(n):.1e}' for n in (1000, 1500, 2000, 2500))} "
                f"| {max(c['seconds_per_call_1000'] for c in cases):.3f} / {max(c['seconds_per_call_production'] for c in cases):.3f} |")

lines += ["", "### Class scattering properties, V24 - V23 (compare_optics.py)", "",
          "| Model | max abs rel. extinction | max abs rel. extinction ratio | max abs dw | max abs dg | L range | max abs dchi (l < 1000) | consistency |",
          "|---|---|---|---|---|---|---|---|"]
for d in load("optics_*.json"):
   rows, c = d["rows"], d["consistency"]
   m = lambda key: max(abs(r[key]) for r in rows)
   lines.append(f"| {d['model']} | {m('extinction_rel'):.1e} | {m('extinction_ratio_rel'):.1e} | {m('w_abs'):.1e} | {m('g_abs'):.1e} "
                f"| {min(r['L_v24'] for r in rows)}-{max(r['L_v24'] for r in rows)} | {m('max_chi_change_l_lt_1000'):.1e} "
                f"| w in [0,1] {c['w_in_0_1']}, chi_0 = 1 {c['chi0_is_1']}, max abs chi {c['max_abs_chi']:.3g}, abs(chi_1 - g) <= {c['max_abs_chi1_minus_g']:.0e} |")

lines += ["", "### Compact LUTs, V24 - V23 (compare_compact.py; identical Grid B nodes and settings)", ""]
for d in load("compact_*.json"):
   lines += [f"**{d['model']}**: variables {d['same_variables']}, dimensions {d['same_dimensions']}, coordinates {d['same_coordinates']}, "
             f"V24 finite {d['v24_finite']}, bitwise identical variables {len(d['bitwise_identical_variables'])}, "
             f"isolated R_0v call outliers {len(d['r0v_call_outliers'])}", "",
             "| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (1, 5, 11, 21, 39 / ice 5, 41, 93) |", "|---|---|---|---|---|---|---|"]
   for group in ("optics", "operators"):
      for name, e in d[group].items():
         lines.append(f"| {name} | {e['median']:.1e} | {e['p99']:.1e} | {e['max']:.1e} | {e['max_at']} | {e['max_rel']:.1e} | "
                      f"{', '.join(f'{v:.0e}' for v in e['max_by_radius'].values())} |")
   lines.append("")

lines += ["### Mie-step cost, 24 Grid B liquid radii (timing.py)", "",
          "| um | V23 | + 3.5 r_e | + 3.5 r_e + refinement | + 3.5 r_e + adaptive Legendre | V24 (all three) | V24 / V23 |", "|---|---|---|---|---|---|---|"]
for d in sorted(load("timing_*.json"), key=lambda d: d["wavelength"]):
   t = d["totals_seconds"]
   lines.append(f"| {d['wavelength']} | {t['v23']:.0f} s | {t['m1']:.0f} s | {t['m1_m3']:.0f} s | {t['m1_m2']:.0f} s | {t['v24']:.0f} s | {t['v24'] / t['v23']:.1f} |")
(R / "summary.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
