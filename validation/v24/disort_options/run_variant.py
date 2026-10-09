"""Run one DISORT-option experiment through the production LUT path (validation only).

The production cloud generator create_orac_luts.create_orac_cloud_lut is run
on a compact grid of Grid B nodes (compact_disort_liquid.lut /
compact_disort_ice.lut) with:

   * the scattering properties computed once by the production code
     (scat_only) and reused by every experiment (reuse_scat), so that only
     the radiative transfer differs;
   * nstreams passed through the generator's own argument;
   * every other DISORT change applied by replacing the DISORT kernel that
     oraclut.idl_mirror.call_disort calls (the name ``disort`` in that module)
     with an experimental build (build_variants.py: MXCMU 256, settable
     ACCUR / DELTAM / CORINT, optional double precision), or by giving the
     generator a scattering cache whose moments are truncated (NMOM test).

Experiment "baseline" uses the production kernel unchanged.  Every DISORT
call is timed (process CPU time) by wrapping create_orac_luts.call_disort.

   env -i PATH=/usr/bin:/bin HOME=/tmp OMP_NUM_THREADS=1 python run_variant.py CASE EXPERIMENT [REPEAT]

Products: validation/tmp/disort_options/runs/CASE/EXPERIMENT[_rREPEAT]/.
"""

import ctypes
import importlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
TMP = ROOT / "validation" / "tmp" / "disort_options"
INPUTS = ROOT / "create_orac_lut" / "input_files"
MIRROR = TMP / "inputs"

MODIS_CHANNELS = [8, 3, 1, 2, 6, 7, 20, 31, 32]       # 0.41 0.47 0.65 0.86 1.63 2.13 3.79 11.0 12.0 um
SLSTR_CHANNELS = [1, 3, 5, 6, 7, 8, 9]                # 0.55 0.87 1.61 2.25 3.74 10.85 12.0 um
CASES = {
   # name: (instrument file, microphysics, compact grid, channels, nmom (run-file value))
   "modis_liquid": ("terra_modis_v1.inst", "liquid-water_stg.mm", "compact_disort_liquid.lut", MODIS_CHANNELS, None),
   "modis_ice_agg": ("terra_modis_v1.inst", "water-ice_agg.mm", "compact_disort_ice.lut", MODIS_CHANNELS, 1000),
   "slstr_liquid": ("sentinel-3a_slstr_v1.inst", "liquid-water_stg.mm", "compact_disort_liquid.lut", SLSTR_CHANNELS, None),
   "slstr_ice_sph": ("sentinel-3a_slstr_v1.inst", "water-ice_sph.mm", "compact_disort_ice.lut", SLSTR_CHANNELS, None),
}

# Reproductions of corrupted direct-beam calls found in V24 products (severe_negatives.py)
CASES["defect_slstr240"] = ("sentinel-3a_slstr_v1.inst", "liquid-water_240.mm", "defect_slstr240.lut", [2], None)
CASES["defect_modis_agg"] = ("aqua_modis_v1.inst", "water-ice_agg.mm", "defect_modis_agg.lut", [24], 1000)
CASES["defect_modis_ghm5"] = ("terra_modis_v1.inst", "water-ice_ghm.mm", "defect_modis_ghm5.lut", [5], 1000)
CASES["defect_modis_ghm17"] = ("terra_modis_v1.inst", "water-ice_ghm.mm", "defect_modis_ghm17.lut", [17], 1000)
# Smooth negative R_0v of the Baum GHM model at large radius in the 4.5 um band (not an isolated call)
CASES["ringing_modis_ghm24"] = ("terra_modis_v1.inst", "water-ice_ghm.mm", "ringing_modis_ghm24.lut", [24], 1000)

# Near-backscatter geometry sets (Grid B nodes: SZA, VZA 30 and 33.75 deg, RAA 0-180 every 5 deg),
# sharing the scattering cache of the matching case
for _case in list(CASES):
   _inst, _mm, _grid, _ch, _nmom = CASES[_case]
   CASES[_case + "_bs"] = (_inst, _mm, _grid.replace(".lut", "_bs.lut"), _ch, _nmom)
SCAT_OF = {case: case[:-3] if case.endswith("_bs") else case for case in CASES}

# name: (nstreams, kernel, accur, deltam, corint, moments)
#   kernel "production" = the production library; "single"/"double" = build_variants.py
#   moments "production" = the generator's own (max(L, NSTR + 1)); "nstr+1" = truncated to NSTR + 1
# Stream counts: only values that DISORT accepts at every Grid B solar zenith
# (SETDIS refuses |mu0 - mu_i| / mu0 < 1e-4 for a computational node mu_i; every
# NSTR = 2 mod 4 puts a node at mu = 0.5 (60 deg), and NSTR >= 240 at 0 deg).
EXPERIMENTS = {"baseline": (60, "production", 1e-8, 1, 1, "production"),
               "kernel_check": (60, "single", 1e-8, 1, 1, "production")}
for n in (8, 16, 24, 32, 40, 48, 56, 68, 72, 80, 88, 96, 100):
   EXPERIMENTS[f"nstr{n:03d}"] = (n, "production", 1e-8, 1, 1, "production")
for n in (112, 136, 152, 168, 200, 236):
   EXPERIMENTS[f"nstr{n:03d}"] = (n, "single", 1e-8, 1, 1, "production")
for accur in (1e-2, 1e-3, 1e-4, 1e-6, 0.0):
   EXPERIMENTS[f"accur{accur:g}"] = (60, "single", accur, 1, 1, "production")
EXPERIMENTS.update({
   "nocorint": (60, "single", 1e-8, 1, 0, "production"),
   "nodeltam_nocorint": (60, "single", 1e-8, 0, 0, "production"),
   "moments_nstr1": (60, "single", 1e-8, 1, 1, "nstr+1"),
   "double": (60, "double", 1e-8, 1, 1, "production"),
   # targeted interactions
   "nocorint_nstr136": (136, "single", 1e-8, 1, 0, "production"),
   "nocorint_nstr236": (236, "single", 1e-8, 1, 0, "production"),
   "moments_nstr1_nstr136": (136, "single", 1e-8, 1, 1, "nstr+1"),
   "double_nstr136": (136, "double", 1e-8, 1, 1, "production"),
   "double_nstr236": (236, "double", 1e-8, 1, 1, "production"),     # reference
   "nstr032_accur1e-3": (32, "single", 1e-3, 1, 1, "production"),
   "nstr040_accur1e-3": (40, "single", 1e-3, 1, 1, "production"),
   # timing calibration: the MXCMU = 256 build at NSTR 100 (per-call array-zeroing overhead)
   "nstr100_k256": (100, "single", 1e-8, 1, 1, "production"),
})


def make_mirror():
   (MIRROR / "lut").mkdir(parents=True, exist_ok=True)
   for item in INPUTS.iterdir():
      link = MIRROR / item.name
      if item.name != "lut" and not link.is_symlink():
         try:
            link.symlink_to(item)
         except FileExistsError:                     # created by a concurrent run
            pass
   for grid in ("compact_disort_liquid.lut", "compact_disort_ice.lut", "compact_disort_liquid_bs.lut", "compact_disort_ice_bs.lut",
                "defect_slstr240.lut", "defect_modis_agg.lut", "defect_modis_ghm5.lut",
                "defect_modis_ghm17.lut", "ringing_modis_ghm24.lut"):
      (MIRROR / "lut" / grid).write_text((HERE / grid).read_text())


def experimental_kernel(kernel, accur, deltam, corint):
   """A drop-in replacement for oraclut.radiative_transfer.legacy_disort.disort."""

   library = ctypes.CDLL(str(TMP / "build" / kernel / f"libdisort_{kernel}.so"))
   real = np.float64 if kernel == "double" else np.float32
   c_real = ctypes.c_double if kernel == "double" else ctypes.c_float
   vec = np.ctypeslib.ndpointer(real, flags="C_CONTIGUOUS")
   fort = np.ctypeslib.ndpointer(real, flags="F_CONTIGUOUS")
   i32 = ctypes.c_int32
   library.oraclut_disort.argtypes = [i32, vec, vec, i32, fort, vec, c_real, c_real, i32, i32, vec, i32, i32, vec, i32, vec,
                                      c_real, c_real, c_real, vec, vec, vec, vec, vec, fort, vec, vec]
   library.oraclut_disort.restype = ctypes.c_int
   library.oraclut_set_options.argtypes = [ctypes.c_double, i32, i32]
   library.oraclut_set_options(accur, deltam, corint)

   def disort(dtauc, single_scatter_albedo, phase_moments, utau, umu, phi, fbeam, umu0, fisot, *, nstreams=60, plank=False,
              wavenumber_low=0.0, wavenumber_high=0.0, temperature=None):
      # same argument handling as the production binding, in the kernel's precision
      dtauc = np.ascontiguousarray(np.asarray(dtauc, dtype=np.float32).ravel(), dtype=real)
      ssalb = np.ascontiguousarray(np.asarray(single_scatter_albedo, dtype=np.float32).ravel(), dtype=real)
      nlyr = dtauc.size
      phase = np.asfortranarray(np.asarray(phase_moments, dtype=np.float32), dtype=real)
      utau = np.ascontiguousarray(np.asarray(utau, dtype=np.float32).ravel(), dtype=real)
      umu = np.ascontiguousarray(np.asarray(umu, dtype=np.float32).ravel(), dtype=real)
      phi = np.ascontiguousarray(np.asarray(phi, dtype=np.float32).ravel(), dtype=real)
      temper = np.ascontiguousarray(np.zeros(nlyr, dtype=np.float32) if temperature is None
                                    else np.asarray(temperature, dtype=np.float32).ravel(), dtype=real)
      ntau, numu, nphi = utau.size, umu.size, phi.size
      out = {k: np.zeros(ntau, dtype=real) for k in ("rfldir", "rfldn", "flup", "dfdt", "uavg")}
      uu = np.zeros((numu, ntau, nphi), dtype=real, order="F")
      albmed = np.zeros(numu, dtype=real)
      trnmed = np.zeros(numu, dtype=real)
      status = library.oraclut_disort(nlyr, dtauc, ssalb, phase.shape[0] - 1, phase, temper,
                                      float(np.float32(wavenumber_low)), float(np.float32(wavenumber_high)), int(plank), ntau,
                                      utau, int(nstreams), numu, umu, nphi, phi, float(np.float32(fbeam)),
                                      float(np.float32(umu0)), float(np.float32(fisot)), out["rfldir"], out["rfldn"],
                                      out["flup"], out["dfdt"], out["uavg"], uu, albmed, trnmed)
      if status:
         raise RuntimeError(f"experimental DISORT returned error code {status}")
      result = {k: v.astype(np.float32) for k, v in out.items()}
      result.update(uu=np.asfortranarray(uu.astype(np.float32)), albmed=albmed.astype(np.float32),
                    trnmed=trnmed.astype(np.float32))
      return result

   return disort


def scattering(case):
   instfile, mmfile, grid, channels, nmom = CASES[SCAT_OF[case]]
   scat = TMP / "scat" / SCAT_OF[case]
   if not (scat / "scatfile.npz").exists():
      import create_orac_luts
      scat.mkdir(parents=True, exist_ok=True)
      create_orac_luts.create_orac_cloud_lut(MIRROR, instfile, mmfile, grid, scat, 2, channelid=channels, srf_quad=1,
                                             scat_only=1, version=None, nstreams=60, nmom=nmom, work_path=scat)
   return scat / "scatfile.npz"


def run(case, experiment, repeat=None):
   import create_orac_luts
   instfile, mmfile, grid, channels, nmom = CASES[case]
   nstr, kernel, accur, deltam, corint, moments = EXPERIMENTS[experiment]
   cache = scattering(case)
   out = TMP / "runs" / case / (experiment if repeat is None else f"{experiment}_r{repeat}")
   if out.exists():
      shutil.rmtree(out)
   out.mkdir(parents=True)
   if moments == "production":
      shutil.copy(cache, out / "scatfile.npz")
   else:
      saved = dict(np.load(cache))
      lmom = saved["lmom"]
      saved["lmom"] = np.minimum(lmom, nstr + 1)
      amom = saved["amom"].copy()
      amom[nstr + 1:] = 0.0
      saved["amom"] = amom
      np.savez(out / "scatfile.npz", **saved)
   if kernel != "production":
      # Timing repeats with NSTR <= 100 use the MXCMU = 100 build of the same
      # kernel (identical results; DISORT zeroes MXCMU-sized arrays on every
      # call, so the MXCMU = 256 build carries a fixed per-call overhead).
      # (the double-precision MXCMU = 100 build corrupts the heap; double precision is timed with
      # the MXCMU = 256 build, an upper bound on its cost)
      build = kernel + "100" if repeat is not None and nstr <= 100 and kernel == "single" and \
         not experiment.endswith("_k256") else kernel
      importlib.import_module("oraclut.idl_mirror.call_disort").disort = experimental_kernel(build, accur, deltam, corint)
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
   status = create_orac_luts.create_orac_cloud_lut(MIRROR, instfile, mmfile, grid, out, 2, channelid=channels, srf_quad=1,
                                                   reuse_scat=1, version=None, nstreams=nstr, nmom=nmom,
                                                   work_path=out / "work")
   record = {"case": case, "experiment": experiment, "repeat": repeat, "nstreams": nstr, "kernel": kernel,
             "kernel_build": "production" if kernel == "production" else build, "accur": accur,
             "deltam": deltam, "corint": corint, "moments": moments, "status": status,
             "total_cpu_s": time.process_time() - started_cpu, "total_wall_s": time.time() - started_wall,
             "disort_calls": {k: v[0] for k, v in timing.items()}, "disort_cpu_s": {k: v[1] for k, v in timing.items()},
             "disort_cpu_total_s": sum(v[1] for v in timing.values())}
   (out / "record.json").write_text(json.dumps(record, indent=2) + "\n")
   print(json.dumps(record), flush=True)


def main():
   make_mirror()
   case, experiment = sys.argv[1], sys.argv[2]
   if experiment == "scat":
      scattering(case)
      return
   run(case, experiment, int(sys.argv[3]) if len(sys.argv) > 3 else None)


if __name__ == "__main__":
   main()
