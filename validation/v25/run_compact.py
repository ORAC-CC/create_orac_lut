"""Compact-grid V24 / V25 comparison runs of the production generator (validation only).

Cases on the compact Grid B subsets validation/v25/compact_v25_ice.lut and
compact_v25_liquid.lut (tau 0.0625, 0.25, 1, 4, 16, 64, 256; r_e 5-93 um ice,
3-35 um liquid; SZA/VZA 0-75; 7 azimuths),
all channels of the DISORT-options cases, production kernel, NSTR 60, from a
scattering cache (the cache depends on the channels and radii only):

   modis_ice_agg   Terra MODIS, Baum aggregate     (cache of validation/v24/disort_options)
   modis_ice_sph   Terra MODIS, ice spheres        (own cache)
   slstr_ice_sph   Sentinel-3A SLSTR, ice spheres  (cache of validation/v24/disort_options)
   slstr_ice_agg   Sentinel-3A SLSTR, Baum aggregate (own cache)
   modis_liquid    Terra MODIS, liquid water stg       (cache of validation/v24/disort_options)
   slstr_liquid    Sentinel-3A SLSTR, liquid water stg (cache of validation/v24/disort_options)

Modes:
   isothermal      the V24 calculation (cloud_vertical_profile absent)
   v25             the V25 calculation of the case's backend (cirrostratus for ice,
                   wet_adiabat for liquid water; production EMISSION_LAYERS); the
                   backend names are accepted as synonyms
   layersN         cirrostratus with N emission layers (layer-count convergence; N <= 46)
   ttopT           cirrostratus at reference cloud-top temperature T K (diagnostic of section 8)
   constantbeta    the superseded z = tau/beta saturated-adiabat formulation of revision 13b9df3
                   (validation/tmp/v25/old_constant_beta_module.py), 12 sub-layers

   python validation/v25/run_compact.py CASE MODE
Outputs: validation/tmp/v25/compact/<case>/<mode>/*.nc with record.json.
"""

import importlib.util
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "v24" / "disort_options"))

import run_variant as rv                                  # noqa: E402  (MIRROR, make_mirror, channel lists)

TMP = ROOT / "validation" / "tmp" / "v25"
OUT = TMP / "compact"
GRIDS = {"ice": "compact_v25_ice.lut", "liquid": "compact_v25_liquid.lut"}
# case: (instrument file, microphysics, channels, nmom, scattering cache directory, grid family, V25 backend)
CASES = {
   "modis_ice_agg": ("terra_modis_v1.inst", "water-ice_agg.mm", rv.MODIS_CHANNELS, 1000, rv.TMP / "scat" / "modis_ice_agg", "ice", "cirrostratus"),
   "modis_ice_sph": ("terra_modis_v1.inst", "water-ice_sph.mm", rv.MODIS_CHANNELS, None, TMP / "scat" / "modis_ice_sph", "ice", "cirrostratus"),
   "slstr_ice_sph": ("sentinel-3a_slstr_v1.inst", "water-ice_sph.mm", rv.SLSTR_CHANNELS, None, rv.TMP / "scat" / "slstr_ice_sph", "ice", "cirrostratus"),
   "slstr_ice_agg": ("sentinel-3a_slstr_v1.inst", "water-ice_agg.mm", rv.SLSTR_CHANNELS, 1000, TMP / "scat" / "slstr_ice_agg", "ice", "cirrostratus"),
   "modis_liquid": ("terra_modis_v1.inst", "liquid-water_stg.mm", rv.MODIS_CHANNELS, None, rv.TMP / "scat" / "modis_liquid", "liquid", "wet_adiabat"),
   "slstr_liquid": ("sentinel-3a_slstr_v1.inst", "liquid-water_stg.mm", rv.SLSTR_CHANNELS, None, rv.TMP / "scat" / "slstr_liquid", "liquid", "wet_adiabat"),
}


def scattering(case):
   import create_orac_luts
   instfile, mmfile, channels, nmom, scat, family, _ = CASES[case]
   if not (scat / "scatfile.npz").exists():
      scat.mkdir(parents=True, exist_ok=True)
      create_orac_luts.create_orac_cloud_lut(rv.MIRROR, instfile, mmfile, GRIDS[family], scat, 2, channelid=channels, srf_quad=1,
                                             scat_only=1, version=None, nstreams=60, nmom=nmom, work_path=scat)
   return scat / "scatfile.npz"


def main():
   case, mode = sys.argv[1], sys.argv[2]
   rv.make_mirror()
   for grid in GRIDS.values():
      (rv.MIRROR / "lut" / grid).write_text((HERE / grid).read_text())
   import create_orac_luts
   from oraclut import cloud_temperature
   instfile, mmfile, channels, nmom, _, family, backend = CASES[case]
   GRID = GRIDS[family]
   cache = scattering(case)
   out = OUT / case / mode
   if out.exists():
      shutil.rmtree(out)
   out.mkdir(parents=True)
   shutil.copy(cache, out / "scatfile.npz")
   options = {}
   layers = 1
   if mode in ("cirrostratus", "wet_adiabat", "v25"):
      options["cloud_vertical_profile"] = backend
      layers = cloud_temperature.EMISSION_LAYERS
   elif mode.startswith("layers"):
      layers = int(mode[len("layers"):])
      options["cloud_vertical_profile"] = backend
      cloud_temperature.EMISSION_LAYERS = layers
      create_orac_luts.emission_layers.__defaults__ = (layers,)
   elif mode.startswith("ttop"):
      t_top = float(mode[len("ttop"):])
      options["cloud_vertical_profile"] = backend
      layers = cloud_temperature.EMISSION_LAYERS
      cloud_temperature.T_TOP_K = t_top
      create_orac_luts.cloud_vertical_profile_model.__defaults__ = (t_top, None)
   elif mode == "constantbeta":
      spec = importlib.util.spec_from_file_location("old_cloud_temperature", TMP / "old_constant_beta_module.py")
      old = importlib.util.module_from_spec(spec)
      spec.loader.exec_module(old)
      options["cloud_vertical_profile"] = backend
      layers = old.EMISSION_SUBLAYERS
      # the superseded formulation: its model and sub-layering replace the profile model
      create_orac_luts.cloud_vertical_profile_model = lambda name, max_tau, **_: old.cloud_temperature_model(
         "water-ice" if family == "ice" else "liquid-water", max_tau)
      create_orac_luts.emission_layers = lambda dtau, ssalb, pmom, srt, tau, model: old.emission_layers(dtau, ssalb, pmom, srt, tau, model)
   elif mode != "isothermal":
      raise SystemExit("MODE must be isothermal, v25 (cirrostratus / wet_adiabat), layersN, ttopT or constantbeta")

   timing = {"diffuse": [0, 0.0], "direct": [0, 0.0], "emission": [0, 0.0]}
   production_call = create_orac_luts.call_disort

   def timed(disort_vars, dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kwargs):
      kind = "emission" if kwargs.get("plank") else ("direct" if fbeam > 0 else "diffuse")
      started = time.process_time()
      result = production_call(disort_vars, dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, **kwargs)
      timing[kind][0] += 1
      timing[kind][1] += time.process_time() - started
      return result

   create_orac_luts.call_disort = timed
   started_cpu, started_wall = time.process_time(), time.time()
   status = create_orac_luts.create_orac_cloud_lut(rv.MIRROR, instfile, mmfile, GRID, out, 2, channelid=channels,
                                                   srf_quad=1, reuse_scat=1, version=None, nstreams=60, nmom=nmom,
                                                   work_path=out / "work", **options)
   record = {"case": case, "mode": mode, "status": status, "emission_layers": layers,
             "total_cpu_s": time.process_time() - started_cpu, "total_wall_s": time.time() - started_wall,
             "disort_calls": {k: v[0] for k, v in timing.items()}, "disort_cpu_s": {k: v[1] for k, v in timing.items()}}
   (out / "record.json").write_text(json.dumps(record, indent=2) + "\n")
   print(json.dumps(record), flush=True)


if __name__ == "__main__":
   main()
