"""Compact Terra MODIS cloud LUT for the V24/V25 cloud-temperature comparison (read-only for the code).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v25/temperature_profiles/compact_lut.py <model> <variant> <out_dir> <channel> ...

<variant>:
  v24          the sealed V24 release (releases/lut_v24.release, verified
               files loaded by the release launcher's import finder)
  isothermal   the V25 working tree, cloud_vertical_profile not set
  profile      the V25 working tree with the production profile of the
               substance (liquid water: wet_adiabat; water ice: cirrostratus)
  constant     as profile, with every sub-layer boundary temperature set to
               the cloud top (250 K): the isothermal limit of the
               sub-layered emission call
  mirror       as profile, reflected about the cloud top, T = 2 T_top - T(F):
               a cloud cooling with depth
  ttop240      as profile with the historical V25 cloud top 240 K (profile
               and Planck normalisation B(240 K))
  layersN      as profile with N emission sub-layers (N <= 46)

The variants other than v24 replace create_orac_luts.emission_layers (and,
for ttop240, plkavg) in the driver's namespace for this process only; the
production code is not changed.  Settings: those of the standard Terra
MODIS run files (atmosphere 2, no gas, Rayleigh on, SRF quadrature 1, NSTR 60,
nmom 1000) on compact_liquid_grid_b.lut or compact_ice_grid_b.lut.  The wall
time and the time spent in the emission DISORT calls (plank=True) are
appended to <out_dir>/timing.json.
"""

import hashlib
import importlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
V25_FILES = ("create_orac_luts.py", "src/oraclut/idl_mirror/cloud_temperature.py", "src/oraclut/idl_mirror/write_v2_lut.py")
PRODUCTION_PROFILE = {"liquid-water": "wet_adiabat", "water-ice": "cirrostratus"}

model, variant, out_dir, channels = sys.argv[1], sys.argv[2], Path(sys.argv[3]), [int(c) for c in sys.argv[4:]]
substance = model.split("_")[0]
if variant == "v24":
   spec = importlib.util.spec_from_file_location("oraclut_release_launcher_v04", ROOT / "create_orac_luts_release_v04.py")
   launcher = importlib.util.module_from_spec(spec)
   spec.loader.exec_module(launcher)
   manifest = launcher.read_manifest(ROOT / "releases" / "lut_v24.release")
   launcher.verify_manifest_files(manifest, ("driver", "package", "module", "config", "adapter", "kernel"))
   sys.meta_path.insert(0, launcher.ReleaseFinder(manifest))
   driver = importlib.import_module("create_orac_luts")
   launcher.substitute_config_files()
   print(f"V24 sealed release: {driver.__file__}", flush=True)
else:
   sys.path.insert(0, str(ROOT / "src"))
   sys.path.insert(0, str(ROOT))
   driver = importlib.import_module("create_orac_luts")
   for name in V25_FILES:
      print(f"V25 source {name}: sha256 {hashlib.sha256((ROOT / name).read_bytes()).hexdigest()}", flush=True)
legacy_disort = importlib.import_module("oraclut.radiative_transfer.legacy_disort")
if legacy_disort.disort_precision() != "single" or legacy_disort.single_kernel() != "intel_dlm":
   raise SystemExit("set ORAC_LUT_DISORT_PRECISION=single and ORAC_LUT_DISORT_SP_KERNEL=intel_dlm (the V24 kernel)")

options = {}
if variant not in ("v24", "isothermal"):
   options["cloud_vertical_profile"] = PRODUCTION_PROFILE[substance]
   cloud_temperature = importlib.import_module("oraclut.idl_mirror.cloud_temperature")
   emission_layers = driver.emission_layers
   if variant == "constant":
      def constant_layers(emtau, emssa, empmo, scatreltau, tau_055, cloudprof, temp):
         subtau, subssa, subpmo, subtemp = emission_layers(emtau, emssa, empmo, scatreltau, tau_055, cloudprof, temp)
         return subtau, subssa, subpmo, np.full_like(subtemp, np.float32(temp))
      driver.emission_layers = constant_layers
   elif variant == "mirror":
      def mirror_layers(emtau, emssa, empmo, scatreltau, tau_055, cloudprof, temp):
         subtau, subssa, subpmo, subtemp = emission_layers(emtau, emssa, empmo, scatreltau, tau_055, cloudprof, temp)
         return subtau, subssa, subpmo, (2.0 * float(temp) - subtemp.astype(np.float64)).astype(np.float32)
      driver.emission_layers = mirror_layers
   elif variant == "ttop240":
      profile240 = cloud_temperature.setup_cloud_profile(PRODUCTION_PROFILE[substance], 256.0, 240.0)
      plkavg = driver.plkavg

      def layers240(emtau, emssa, empmo, scatreltau, tau_055, cloudprof, temp):
         if temp != 250.0:
            raise RuntimeError(f"unexpected emission temperature {temp}")
         return emission_layers(emtau, emssa, empmo, scatreltau, tau_055, profile240, 240.0)

      def plkavg240(wnlo, wnhi, temp):
         if temp != 250.0:
            raise RuntimeError(f"unexpected normalising temperature {temp}")
         return plkavg(wnlo, wnhi, 240.0)
      driver.emission_layers = layers240
      driver.plkavg = plkavg240
   elif variant.startswith("layers"):
      cloud_temperature.EMISSION_LAYERS = int(variant[len("layers"):])
   elif variant != "profile":
      raise SystemExit(f"unknown variant {variant!r}")

# Time the emission calls (plank=True) separately
call_disort = driver.call_disort
emission = {"calls": 0, "seconds": 0.0}


def timed_call_disort(*args, **kwargs):
   if not kwargs.get("plank"):
      return call_disort(*args, **kwargs)
   start = time.perf_counter()
   result = call_disort(*args, **kwargs)
   emission["seconds"] += time.perf_counter() - start
   emission["calls"] += 1
   return result


driver.call_disort = timed_call_disort
grid = HERE / ("compact_ice_grid_b.lut" if substance == "water-ice" else "compact_liquid_grid_b.lut")
out_dir.mkdir(parents=True, exist_ok=True)
start = time.perf_counter()
status = driver.create_orac_cloud_lut(ROOT / "create_orac_lut" / "input_files", "terra_modis_v1.inst", f"{model}.mm",
                                      str(grid), out_dir, 2, channelid=channels, gas=0, no_rayleigh=0, srf_quad=1,
                                      version=24 if variant == "v24" else 25, nstreams=60, nmom=1000, **options)
elapsed = time.perf_counter() - start
if status != 0:
   raise SystemExit(f"create_orac_cloud_lut returned {status}")
timing = out_dir / "timing.json"
records = json.loads(timing.read_text()) if timing.is_file() else []
records.append({"variant": variant, "model": model, "channels": channels, "seconds": elapsed,
                "emission_calls": emission["calls"], "emission_seconds": emission["seconds"]})
timing.write_text(json.dumps(records, indent=1) + "\n")
print(f"{variant} {model}: {elapsed:.0f} s ({emission['calls']} emission calls, {emission['seconds']:.1f} s) -> {out_dir}")
