"""Estimate the wall time of the full V24 Terra MODIS liquid-water products (read-only).

   python validation/v24/scattering_refinements/runtime_estimate.py -> results/runtime_estimate.json

Anchored on measured times (atmnode009-010, one CPU):
* the V23 production jobs (validation/slurm/v23_*_49542[234].out start/end);
* the Mie step of V23 and V24 for all 24 Grid B radii at 0.469, 0.646, 1.64,
  3.79 and 11.03 um (timing.py), interpolated linearly in log(time) against
  log(wavelength) over the 36 Terra MODIS channel centres and 0.55 um;
* the DISORT time per call against the number of moments (disort_nmom.py:
  1000 and L for each case), fitted as a + b NMOM.
The V24 DISORT time is the V23 DISORT time scaled, channel by channel and
radius by radius, by the per-call cost at max(1000, L), with L bounded by the
degree bound 2 NStop(x_max) + 1 of the 3.5 r_e integration (an upper bound:
infrared expansions stop well below it).
"""

import datetime
import glob
import json
import re
from pathlib import Path

import numpy as np

import codebase

codebase.select("v24")
mirror = codebase.module("oraclut.idl_mirror")
le = codebase.module("oraclut.idl_mirror.legendre_expansion")
cb = codebase.module("oraclut.idl_mirror.create_bwgp")
HERE = Path(__file__).resolve().parent
ROOT = codebase.ROOT

# V23 production job wall times (hours)
v23_hours = {}
for log in sorted(glob.glob(str(ROOT / "validation" / "slurm" / "v23_*_49542[234].out"))):
   text = Path(log).read_text()
   start = re.search(r"Start time:\s+(.+)", text).group(1).strip()
   end = re.search(r"End time:\s+(.+)", text).group(1).strip()
   parse = lambda s: datetime.datetime.strptime(s.replace("BST ", ""), "%a %d %b %H:%M:%S %Y")
   v23_hours[Path(log).name.split("_")[1]] = (parse(end) - parse(start)).total_seconds() / 3600.0

# Mie step per wavelength (24 radii)
timing = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(HERE / "results" / "timing_liquid-water_240_*.json")))]
lam = np.log([t["wavelength"] for t in timing])
order = np.argsort(lam)
log_v23 = np.log([t["totals_seconds"]["v23"] for t in timing])[order]
log_v24 = np.log([t["totals_seconds"]["v24"] for t in timing])[order]
lam = lam[order]
inputs = ROOT / "create_orac_lut" / "input_files"
inst = mirror.load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=list(range(1, 37)))
srf, nwvl = mirror.load_srfstrarr(inst, inputs / "sun" / "Gueymard2018.sssi", 1, inputs)
wavelengths = np.asarray([0.55] + [s.wvl_centre for s in srf], dtype=np.float64)
mie_v23 = float(np.sum(np.exp(np.interp(np.log(wavelengths), lam, log_v23))))
mie_v24 = float(np.sum(np.exp(np.interp(np.log(wavelengths), lam, log_v24))))

# DISORT per-call time against NMOM
points = []
for f in glob.glob(str(HERE / "results" / "disort_nmom_*.json")):
   d = json.loads(Path(f).read_text())
   for case in d["cases"].values():
      points.append((1000, case["seconds_per_call_1000"]))
      points.append((max(1000, d["L"]), case["seconds_per_call_production"]))
n, t = np.asarray(points, dtype=np.float64).T
slope, intercept = np.polyfit(n, t, 1)
radii = mirror.load_lutstr(inputs / "lut" / "liquid-water-cloud.lut", inst.max_sat_zenith).efr.astype(np.float64)
factors = []
for s in srf:
   wavenumber = float(np.float32(1.0) / np.float32(s.wvl_centre))
   for radius in radii:
      upper = cb.lattice_upper_radius(radius, 3.5, wavenumber)[0]
      nmom = max(1000, le.degree_bound(2.0 * np.pi * upper * wavenumber) + 1)
      factors.append((intercept + slope * nmom) / (intercept + slope * 1000.0))
disort_factor = float(np.mean(factors))

estimate = {}
for product in ("pold", "p240"):
   v23 = v23_hours[product]
   disort_v23 = v23 - mie_v23 / 3600.0
   estimate[product] = {"v23_hours": v23, "v23_mie_hours": mie_v23 / 3600.0, "v24_mie_hours": mie_v24 / 3600.0,
                        "disort_cost_factor": disort_factor, "v24_hours": disort_v23 * disort_factor + mie_v24 / 3600.0}
estimate["pagg"] = {"v23_hours": v23_hours.get("pagg"), "v24_hours": v23_hours.get("pagg"),
                    "note": "Baum class: calculated as in V23 (compact products bitwise identical)"}
estimate["disort_seconds_per_call_fit"] = {"intercept": float(intercept), "per_moment": float(slope)}
(HERE / "results" / "runtime_estimate.json").write_text(json.dumps(estimate, indent=1) + "\n")
print(json.dumps(estimate, indent=1))
