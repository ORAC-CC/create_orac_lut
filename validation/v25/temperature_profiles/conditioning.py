"""Conditioning of the single-precision DISORT emission call with optically thin sub-layers (read-only).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v25/temperature_profiles/conditioning.py   -> results/conditioning.json

A homogeneous slab of optical depth D whose Planck radiance is exactly linear
in optical depth, B(t) = B(250 K) + (B(T_base) - B(250 K)) t / D, is
represented exactly by DISORT's linear-in-tau source on any subdivision.  So
the emergent radiance must not depend on how the slab is divided; any change
is round-off in the production kernel.  The slab is split into k sub-layers
of optical depth h at the top and the rest in equal layers (46 in all), and
compared with the undivided slab, for h = 1e-8 ... 1e-1, single-scattering
albedos 0 ... 0.999, a Henyey-Greenstein phase function (g 0.85), the
channel 20 and 31 windows, a 16 K and an 80 K temperature span, and view
zeniths 0-75 deg.  The result fixes the smallest sub-layer optical depth
that the V25 emission layering may use.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.idl_mirror import call_disort, setup_disort   # noqa: E402
from oraclut.radiative_transfer.legacy_disort import disort_precision, plkavg, single_kernel   # noqa: E402

if disort_precision() != "single" or single_kernel() != "intel_dlm":
   raise SystemExit("set ORAC_LUT_DISORT_PRECISION=single and ORAC_LUT_DISORT_SP_KERNEL=intel_dlm")
f32 = np.float32
SAZ = np.asarray([0.0, 30.0, 60.0, 75.0], f32)
umu = np.concatenate([-np.cos(np.radians(SAZ)), np.cos(np.radians(SAZ))[::-1]]).astype(f32)
umu[np.abs(umu) == 1.0] *= f32(0.99999)
NMOM = 61


def inverse_planck(b, wnlo, wnhi, t):
   for _ in range(30):
      db = (plkavg(wnlo, wnhi, t + 0.01) - plkavg(wnlo, wnhi, t - 0.01)) / 0.02
      step = (plkavg(wnlo, wnhi, t) - b) / db
      t -= step
      if abs(step) < 1e-7:
         break
   return t


def emitted(dtau, ssa, temper, wnlo, wnhi):
   n = dtau.size
   pmom = np.asfortranarray(np.repeat((0.85 ** np.arange(NMOM, dtype=np.float64)).astype(f32)[:, None], n, axis=1))
   uu = call_disort(setup_disort(60, n, SAZ.size, 1, NMOM), dtau.astype(f32), np.full(n, ssa, f32), pmom,
                    np.asarray([0.0, dtau.sum(dtype=f32)], f32), umu, np.asarray([0.0], f32), 0.0, f32(1.0), 0.0,
                    plank=True, wnlo=wnlo, wnhi=wnhi, temp=temper.astype(f32), nlayer=n)[5]
   return uu[SAZ.size:, 0, 0].astype(np.float64)


results = []
for channel, wvl in ((20, 3.788211), (31, 11.01842)):
   wn = f32(1e4) / f32(wvl)
   wnlo, wnhi = float(f32(0.995) * wn), float(f32(1.005) * wn)
   b250 = plkavg(wnlo, wnhi, 250.0)
   for span in (16.0, 80.0):
      b_base = plkavg(wnlo, wnhi, 250.0 + span)
      for total in (0.1, 1.0, 16.0):
         for ssa in (0.0, 0.5, 0.9, 0.99, 0.999):
            def temperatures(bounds):
               return np.asarray([inverse_planck(b250 + (b_base - b250) * t / total, wnlo, wnhi, 250.0 + span * t / total)
                                  for t in bounds])
            whole = emitted(np.asarray([total]), ssa, temperatures(np.asarray([0.0, total])), wnlo, wnhi)
            for h in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 2.5e-3, 1e-2, 3e-2):
               k = 10
               if k * h >= total:
                  continue
               rest = np.full(46 - k, (total - k * h) / (46 - k))
               dtau = np.concatenate([np.full(k, h), rest])
               bounds = np.concatenate([[0.0], np.cumsum(dtau)])
               bounds[-1] = total
               split = emitted(dtau, ssa, temperatures(bounds), wnlo, wnhi)
               results.append({"channel": channel, "span_K": span, "slab_tau": total, "ssa": ssa, "sublayer_tau": h,
                                "max_rel_change": float(np.max(np.abs(split / whole - 1.0))),
                                "max_change_over_b250": float(np.max(np.abs(split - whole)) / b250)})
summary = {}
for r in results:
   key = f"{r['sublayer_tau']:g}"
   summary[key] = max(summary.get(key, 0.0), r["max_change_over_b250"])
for key, value in summary.items():
   print(f"sub-layer tau {key:>8s}: largest change of the emergent radiance {value:.1e} B(250 K)")
for h in ("0.0001", "0.001", "0.0025", "0.01"):
   worst = max((r for r in results if f"{r['sublayer_tau']:g}" == f"{float(h):g}"), key=lambda r: r["max_change_over_b250"])
   print(f"   worst at {h}: {worst}")
(HERE / "results" / "conditioning.json").write_text(json.dumps({"by_sublayer_tau": summary, "cases": results}, indent=1) + "\n")
