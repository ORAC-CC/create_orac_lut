"""Capture the selected direct-beam call from the corrected V22 production path."""

from __future__ import annotations

import importlib
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
INPUTS = CASE / "run" / "input_files"
RESULTS = CASE / "results"

sys.path.insert(0, str(ROOT))
import create_orac_luts_v02 as production


class CaptureComplete(RuntimeError):
   """Stop the compact production calculation after the selected solver call."""


def main() -> int:
   call_module = importlib.import_module("oraclut.idl_mirror.call_disort")
   original_disort = call_module.disort

   def capture_disort(dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot,
                      *, nstreams, plank, wavenumber_low, wavenumber_high, temperature):
      result = original_disort(
         dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot,
         nstreams=nstreams, plank=plank, wavenumber_low=wavenumber_low,
         wavenumber_high=wavenumber_high, temperature=temperature,
      )
      target_umu0 = np.float32(np.cos(np.deg2rad(np.float32(80.0))))
      if np.float32(fbeam) != np.float32(100.0) or np.float32(umu0) != target_umu0:
         return result

      dtauc = np.ascontiguousarray(np.asarray(dtauc, dtype=np.float32).ravel())
      ssalb = np.ascontiguousarray(np.asarray(ssalb, dtype=np.float32).ravel())
      pmom = np.asfortranarray(np.asarray(pmom, dtype=np.float32))
      utau = np.ascontiguousarray(np.asarray(utau, dtype=np.float32).ravel())
      umu = np.ascontiguousarray(np.asarray(umu, dtype=np.float32).ravel())
      phi = np.ascontiguousarray(np.asarray(phi, dtype=np.float32).ravel())
      temper = np.ascontiguousarray(np.asarray(temperature, dtype=np.float32).ravel())
      nlyr = np.int32(dtauc.size)
      nmom = np.int32(pmom.shape[0] - 1)
      ntau = np.int32(utau.size)
      numu = np.int32(umu.size)
      nphi = np.int32(phi.size)

      RESULTS.mkdir(parents=True, exist_ok=True)
      np.savez(
         RESULTS / "corrected_v22_disort_boundary.npz",
         nlyr=nlyr, dtauc=dtauc, ssalb=ssalb, nmom=nmom, pmom=pmom,
         temper=temper, wvnmlo=np.float32(wavenumber_low),
         wvnmhi=np.float32(wavenumber_high), usrtau=np.int32(1), ntau=ntau,
         utau=utau, nstr=np.int32(nstreams), usrang=np.int32(1), numu=numu,
         umu=umu, nphi=nphi, phi=phi, ibcnd=np.int32(0),
         fbeam=np.float32(fbeam), umu0=np.float32(umu0), phi0=np.float32(0.0),
         fisot=np.float32(fisot), lamber=np.int32(1), albedo=np.float32(0.0),
         btemp=np.float32(0.0), ttemp=np.float32(0.0), temis=np.float32(0.0),
         plank=np.int32(int(plank)), onlyfl=np.int32(0), accur=np.float32(1.0e-8),
         prnt=np.zeros(5, dtype=np.int32), maxcly=nlyr, maxulv=ntau,
         maxumu=numu, maxphi=nphi, maxmom=nmom,
         rfldir=result["rfldir"], rfldn=result["rfldn"], flup=result["flup"],
         dfdt=result["dfdt"], uavg=result["uavg"], uu=result["uu"],
         albmed=result["albmed"], trnmed=result["trnmed"],
         pmom_fortran_contiguous=np.bool_(pmom.flags.f_contiguous),
         uu_fortran_contiguous=np.bool_(result["uu"].flags.f_contiguous),
      )
      raise CaptureComplete

   call_module.disort = capture_disort
   try:
      production.create_orac_cloud_lut(
         INPUTS, "terra_modis_v1.inst", "liquid-water_240.mm",
         "compact_p240_node.lut", CASE / "run" / "python_output", "2",
         channelid=[6], gas=0, no_rayleigh=0, srf_quad=1, reuse_scat=0,
         scat_only=0, tmatrix_path=None, version=22, driver=None,
         nstreams=60, nmom=1000, work_path=CASE / "run" / "python_work",
      )
   except CaptureComplete:
      print("ORACLUT corrected V22 DISORT boundary capture complete")
      return 0
   raise RuntimeError("selected corrected V22 DISORT call was not captured")


if __name__ == "__main__":
   raise SystemExit(main())
