"""Accuracy of the V25 emission layering against independent references (read-only for the code).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v25/temperature_profiles/layering.py   -> results/layering.json

Cases: the emission-call inputs recorded by emission_inputs.py (work/inputs/:
in-cloud optical depth, single-scattering albedo and phase moments of Terra
MODIS channels 20, 31, 32 at Grid B optical depths, ice agg and liquid water).
E_md = UU / B(250 K) at view zeniths 0, 60, 75 deg.

* exact (any optical depth): with the single-scattering albedo set to 0 (the
  same optical depths, profiles and windows), the exact radiance of a source
  linear in optical depth on a 4000-node mesh refined at the cloud top (first
  spacing 1e-4 of the channel optical depth or less): an independent
  measure of how well the layering represents the profile's Planck source.
* convergence (thick clouds, with scattering): DISORT cannot be given more
  than 46 layers, and in single precision it cannot reproduce a steep source
  (hat functions of a superposition reference fail: reference_check), so the
  scheme is compared with itself for 12, 23, 34 and 46 layers and other
  thresholds; E_inf = E_46 + (E_46 - E_34) / ((46/34)^2 - 1) (second-order
  Richardson estimate).
* thin: for tau_055 <= 0.25, (a) the case with single-scattering albedo 0
  against the exact non-scattering radiance of the profile (any optical
  depth), (b) the real case against the first-order thin-cloud radiance
  (1 - w) D B_mean / mu.
* schemes: historical (40 equal sub-layers, temperatures sampled at the
  boundaries; the candidate of 2026-10-15) and the V25 layering
  (cloud_temperature.emission_layers) for other TOP_SUBLAYER_TAU,
  PLANCK_SLOPE_LIMIT and largest layer numbers.
* sweep: E_md for 160 channel optical depths D between 1e-5 and 5 (ice agg
  r_e 5 and 33 um, liquid water r_e 11 um; channels 20 and 31; D and tau_055
  scaled together), the error without scattering, and the change of E_md
  across every change of representation (staircase cells, staircase to
  linear, number of linear sub-layers, top sub-layer), located by bisection,
  for PLANCK_SLOPE_LIMIT 300 and 1000.
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
SAZ = np.asarray([0.0, 60.0, 75.0], f32)
UMU = np.concatenate([-np.cos(np.radians(SAZ)), np.cos(np.radians(SAZ))[::-1]]).astype(f32)
UMU[np.abs(UMU) == 1.0] *= f32(0.99999)
UP = UMU[SAZ.size:][::-1].astype(np.float64)                  # view zeniths 0, 60, 75 deg
REFERENCE_NODES = 600
REFERENCE_TOP = 2e-3
DEFAULTS = (ct.EMISSION_LAYERS, ct.PLANCK_SLOPE_LIMIT, ct.TOP_SUBLAYER_TAU)


def emit(case, dtau, temper, ssa=None):
   """E_md (view zeniths 0, 60, 75) of the emission call on sub-layers dtau with boundary temperatures temper."""

   n = dtau.size
   nmom = int(case["nmom"])
   pmom = np.asfortranarray(np.repeat(case["empmo"][:nmom, None], n, axis=1).astype(f32))
   w = case["emssa"] if ssa is None else ssa
   uu = call_disort(setup_disort(60, n, SAZ.size, 1, nmom), np.asarray(dtau, f32), np.full(n, w, f32), pmom,
                    np.asarray([0.0, np.sum(np.asarray(dtau, f32), dtype=f32)], f32), UMU, np.asarray([0.0], f32), 0.0,
                    f32(1.0), 0.0, plank=True, wnlo=case["wnlo"], wnhi=case["wnhi"], temp=np.asarray(temper, f32), nlayer=n)[5]
   return uu[SAZ.size:, 0, 0][::-1].astype(np.float64) / case["b250"]


def planck(case, temperature):
   return np.asarray([plkavg(case["wnlo"], case["wnhi"], float(t)) for t in np.atleast_1d(temperature)])


def profile_planck(case, fraction):
   return planck(case, ct.cloud_temperature(case["cloudprof"], case["tau_055"], fraction))


def exact_nonscattering(tau, b, mu):
   """Exact top radiance of a non-scattering stack, source linear in optical depth in each layer (stable for small x)."""

   radiance, top = 0.0, 0.0
   for k, dtau in enumerate(tau):
      x = dtau / mu
      absorbed = -np.expm1(-x)
      ramp = x / 2.0 - x * x / 3.0 + x ** 3 / 8.0 if x < 1e-3 else (absorbed - x * np.exp(-x)) / x
      radiance += np.exp(-top / mu) * (b[k] * absorbed + (b[k + 1] - b[k]) * ramp)
      top += dtau
   return radiance


def fine_mesh(dtau, nodes=4000, top=1e-4):
   """A mesh of nodes + 1 boundaries on [0, dtau], geometric from a first spacing top (uniform if that is coarser)."""

   if dtau <= nodes * top:
      return np.linspace(0.0, dtau, nodes + 1)
   low, high = 1.0, 1.0 + 1.0 / nodes
   while top * np.expm1(nodes * np.log(high)) / (high - 1.0) < dtau:
      high = 1.0 + 2.0 * (high - 1.0)
   for _ in range(200):
      r = 0.5 * (low + high)
      low, high = (r, high) if top * np.expm1(nodes * np.log(r)) / (r - 1.0) < dtau else (low, r)
   t = top * np.expm1(np.arange(nodes + 1) * np.log(r)) / (r - 1.0)
   t[-1] = dtau
   return t


def exact(case):
   """E_md of the profile without scattering (exact; source linear between 4001 nodes refined at the top)."""

   t = fine_mesh(case["dtau"], top=min(1e-4, case["dtau"] / 4000))
   b = profile_planck(case, t / case["dtau"])
   return np.asarray([exact_nonscattering(np.diff(t), b, mu) for mu in UP]) / case["b250"]


def hat_responses(case, t):
   """Superposition reference that fails in single precision (kept as the record of reference_check)."""

   peak = 300.0
   bpeak = plkavg(case["wnlo"], case["wnhi"], peak)
   responses = np.zeros((t.size, UP.size))
   for j in range(t.size):
      edges = sorted({0.0, t[max(j - 1, 0)], t[j], t[min(j + 1, t.size - 1)], t[-1]})
      temper = np.asarray([peak if e == t[j] else 0.0 for e in edges])
      responses[j] = emit(case, np.diff(edges), temper) * case["b250"] / bpeak
   return responses


def historical(case, ssa=None):
   bounds = np.arange(41) / 40
   temper = ct.cloud_temperature(case["cloudprof"], case["tau_055"], bounds)
   return emit(case, np.full(40, case["dtau"] / 40), temper, ssa)


def v25(case, layers=None, limit=None, top=None, ssa=None):
   ct.EMISSION_LAYERS, ct.PLANCK_SLOPE_LIMIT, ct.TOP_SUBLAYER_TAU = (layers or DEFAULTS[0], limit or DEFAULTS[1], top or DEFAULTS[2])
   try:
      subtau, subssa, subpmo, subtemp = ct.emission_layers(np.asarray([case["dtau"]], f32), np.asarray([case["emssa"]], f32),
                                                           np.asfortranarray(case["empmo"][:int(case["nmom"]), None]),
                                                           case["tau_055"], case["cloudprof"], 250.0, case["wnlo"], case["wnhi"])
   finally:
      ct.EMISSION_LAYERS, ct.PLANCK_SLOPE_LIMIT, ct.TOP_SUBLAYER_TAU = DEFAULTS
   return emit(case, subtau, subtemp, ssa), subtau.size


def representation(case, limit):
   """(kind, number of layers, top sub-layer h): the choice made by emission_layers (its selection rule, repeated)."""

   d, tau = float(np.float32(case["dtau"])), case["tau_055"]       # emission_layers sees the float32 optical depth
   mean = ct.planck_means(case["cloudprof"], tau, np.asarray([0.0, 1.0]), case["wnlo"], case["wnhi"])[0]
   top = ct.TOP_SUBLAYER_TAU
   while top <= d:
      bounds = ct.sublayer_bounds(d, top)
      if bounds.size - 1 < ct.STAIRCASE_CELLS:
         break
      radiance = ct.planck_projection(case["cloudprof"], tau, bounds, case["wnlo"], case["wnhi"])
      if np.max(np.abs(np.diff(radiance)) / mean / (d * np.diff(bounds)) ** 2) <= limit:
         return ("linear", int(bounds.size - 1), top)
      top *= 2.0
   cells = min(ct.STAIRCASE_CELLS, max(1, int(d / ct.STAIRCASE_CELL_TAU)))
   return ("staircase", 2 * cells - 1, 0.0)


def cases():
   profiles = {}
   for path in sorted((HERE / "work" / "inputs").glob("*.npz")):
      data = np.load(path)
      name = str(data["profile"])
      profiles.setdefault(name, ct.setup_cloud_profile(name, 256.0, 250.0))
      for k in range(data["tau_055"].size):
         case = {"model": str(data["model"]), "re": float(data["effective_radius"]), "channel": int(data["channel"][k]),
                 "tau_055": float(data["tau_055"][k]), "wnlo": float(data["wnlo"][k]), "wnhi": float(data["wnhi"][k]),
                 "dtau": float(data["emtau"][k]), "emssa": float(data["emssa"][k]), "nmom": int(data["nmom"][k]),
                 "empmo": data["empmo"][k], "cloudprof": profiles[name]}
         case["b250"] = plkavg(case["wnlo"], case["wnhi"], 250.0)
         yield case


def label(case):
   return f"{case['model']} r_e {case['re']:g} ch {case['channel']} tau {case['tau_055']:g}"


out = {"reference_check": [], "thick": [], "thin": [], "sweep": {}}
all_cases = list(cases())
# 1. a superposition reference fails in single precision: a source exactly linear in optical depth, built from hats
for case in [c for c in all_cases if c["re"] in (33.0, 11.0) and c["channel"] == 20 and c["tau_055"] in (1.0, 16.0)]:
   t = fine_mesh(case["dtau"], nodes=200, top=2e-3)
   linear = case["b250"] + (planck(case, 330.0)[0] - case["b250"]) * t / case["dtau"]
   hats = (hat_responses(case, t) * linear[:, None]).sum(axis=0) / case["b250"]
   whole = emit(case, np.asarray([case["dtau"]]), [250.0, 330.0])
   out["reference_check"].append({"case": label(case), "hat_superposition_minus_undivided_max_abs": float(np.max(np.abs(hats - whole))),
                                  "undivided_E_md": whole.tolist()})
   print(f"superposition check {label(case)}: |hats - undivided| {out['reference_check'][-1]['hat_superposition_minus_undivided_max_abs']:.1e} "
         f"(E_md {whole[0]:.3f})", flush=True)
# 2. thick clouds: exact without scattering; convergence with scattering
SCHEMES = {"v25": {}, "top0.003": {"top": 0.003}, "top0.03": {"top": 0.03}, "limit300": {"limit": 300.0},
           "layers34": {"layers": 34}, "layers23": {"layers": 23}, "layers28": {"layers": 28}}
for case in [c for c in all_cases if c["tau_055"] in (1.0, 4.0, 16.0, 64.0, 128.0, 256.0)]:
   e0 = exact(case)
   row = {"case": label(case), "model": case["model"], "channel": case["channel"], "tau_055": case["tau_055"], "re": case["re"],
          "ssa": case["emssa"], "dtau": case["dtau"], "exact_nonscattering_E_md": e0.tolist(),
          "nonscattering_error": {"historical": float(np.max(np.abs(historical(case, 0.0) - e0)))}}
   for name, options in SCHEMES.items():
      row["nonscattering_error"][name] = float(np.max(np.abs(v25(case, ssa=0.0, **options)[0] - e0)))
   values = {name: v25(case, **options)[0] for name, options in SCHEMES.items()}
   values["historical"] = historical(case)
   e_inf = values["v25"] + (values["v25"] - values["layers34"]) / ((46.0 / 34.0) ** 2 - 1.0)
   row["E_md_v25"] = values["v25"].tolist()
   row["richardson_E_md"] = e_inf.tolist()
   row["scattering_minus_richardson"] = {name: float(np.max(np.abs(v - e_inf))) for name, v in values.items()}
   out["thick"].append(row)
   ns, sc = row["nonscattering_error"], row["scattering_minus_richardson"]
   print(f"{label(case)} (w {case['emssa']:.3f}, D {case['dtau']:.3g}): no scattering: historical {ns['historical']:.1e}, V25 {ns['v25']:.1e} "
         f"(top 0.003 {ns['top0.003']:.1e}, 0.03 {ns['top0.03']:.1e}); scattering vs Richardson: historical {sc['historical']:.1e}, "
         f"V25 {sc['v25']:.1e}, 34/28/23 {sc['layers34']:.1e}/{sc['layers28']:.1e}/{sc['layers23']:.1e}, top 0.003/0.03 "
         f"{sc['top0.003']:.1e}/{sc['top0.03']:.1e}, limit 300 {sc['limit300']:.1e}", flush=True)
# 3. thin clouds
for case in [c for c in all_cases if c["tau_055"] <= 0.25]:
   e, n = v25(case)
   ehist = historical(case)
   bmean = ct.planck_means(case["cloudprof"], case["tau_055"], np.asarray([0.0, 1.0]), case["wnlo"], case["wnhi"])[0]
   first_order = (1.0 - case["emssa"]) * case["dtau"] * bmean / UP / case["b250"]
   isothermal = emit(case, np.asarray([case["dtau"]]), [250.0, 250.0])
   # non-scattering: V25 layering against the exact radiance of the profile
   ct_bounds = np.linspace(0.0, 1.0, 4001)
   exact0 = np.asarray([exact_nonscattering(np.full(4000, case["dtau"] / 4000), profile_planck(case, ct_bounds), mu) for mu in UP]) / case["b250"]
   subtau, _, _, subtemp = ct.emission_layers(np.asarray([case["dtau"]], f32), np.asarray([0.0], f32),
                                             np.asfortranarray(case["empmo"][:int(case["nmom"]), None]), case["tau_055"],
                                             case["cloudprof"], 250.0, case["wnlo"], case["wnhi"])
   e0 = emit(case, subtau, subtemp, ssa=0.0)
   row = {"case": label(case), "model": case["model"], "channel": case["channel"], "tau_055": case["tau_055"], "re": case["re"],
          "dtau": case["dtau"], "ssa": case["emssa"], "layers": n, "E_md": e.tolist(), "E_md_historical": ehist.tolist(),
          "E_md_isothermal": isothermal.tolist(), "first_order": first_order.tolist(),
          "ratio_to_first_order": (e / first_order).tolist() if case["dtau"] > 0 else None,
          "nonscattering_rel_error": float(np.max(np.abs(e0 / exact0 - 1.0)))}
   out["thin"].append(row)
   print(f"{label(case)}: n {n}; E_md {e[0]:+.4e} (historical {ehist[0]:+.4e}, isothermal {isothermal[0]:.4e}); "
         f"/first order {e[0] / first_order[0]:.4f}..{e[-1] / first_order[-1]:.4f}; non-scattering rel. error {row['nonscattering_rel_error']:.1e}", flush=True)
# 4. sweep in optical depth: E_md along D, the change of E_md at every change of representation, and the error
#    without scattering, for the channel 20 and 31 inputs at tau_055 = 1 (D and tau_055 scaled together)
def sweep(base, limit):
   ratio = base["dtau"] / base["tau_055"]
   at = lambda d: dict(base, dtau=float(d), tau_055=float(d / ratio))
   grid = np.geomspace(1e-5, 5.0, 160)
   points, signatures = [], [representation(at(d), limit) for d in grid]
   for d, sig in zip(grid, signatures):
      e0 = v25(at(d), limit=limit, ssa=0.0)[0]
      points.append({"dtau": float(d), "representation": sig, "E_md": v25(at(d), limit=limit)[0].tolist(),
                     "nonscattering_error": float(np.max(np.abs(e0 - exact(at(d)))))})
   jumps = []
   for k in range(grid.size - 1):
      if signatures[k] == signatures[k + 1]:
         continue
      low, high = np.log(grid[k]), np.log(grid[k + 1])
      for _ in range(40):
         mid = 0.5 * (low + high)
         low, high = (mid, high) if representation(at(np.exp(mid)), limit) == signatures[k] else (low, mid)
      below, above = at(np.exp(low) * (1.0 - 1e-6)), at(np.exp(high) * (1.0 + 1e-6))   # distinct in float32
      eb, ea = v25(below, limit=limit)[0], v25(above, limit=limit)[0]
      jumps.append({"dtau": float(np.exp(high)), "from": representation(below, limit), "to": representation(above, limit),
                    "change": float(np.max(np.abs(ea - eb))), "relative": float(np.max(np.abs(ea / eb - 1.0)))})
   return {"case": label(base), "limit": limit, "points": points, "jumps": jumps}


for model, re_ in (("water-ice_agg", 5.0), ("water-ice_agg", 33.0), ("liquid-water_240", 11.0)):
   for channel in (20, 31):
      base = next(c for c in all_cases if c["model"] == model and c["re"] == re_ and c["channel"] == channel and c["tau_055"] == 1.0)
      for limit in (300.0, 1000.0):
         result = sweep(base, limit)
         out["sweep"][f"{label(base)} limit {limit:g}"] = result
         worst = max(result["jumps"], key=lambda j: j["relative"]) if result["jumps"] else None
         error = max(p["nonscattering_error"] for p in result["points"])
         print(f"sweep {label(base)} limit {limit:g}: {len(result['jumps'])} changes of representation; largest change "
               + (f"{worst['change']:.1e} ({worst['relative']:.1e} relative) at D {worst['dtau']:.3g} {worst['from']} -> {worst['to']}" if worst else "-")
               + f"; largest error without scattering {error:.1e}", flush=True)
(HERE / "results" / "layering.json").write_text(json.dumps(out, indent=1) + "\n")
