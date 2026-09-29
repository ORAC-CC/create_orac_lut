"""Python counterpart of call_disort.pro.

IDL:  pro call_disort, DTAUC, SSALB, PMOM, UTAU, UMU, PHI, FBEAM, UMU0, FISOT,
                       RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED,
                       plank=plank, wnlo=wnlo, wnhi=wnhi, temp=temp, nlayer=nlayer

Checks the input dimensions against the DISORT_VARS block, sets the thermal
inputs, compresses the layer stack (adjacent layers with identical SSALB and
PMOM are merged, which is what happens to the Rayleigh-only layers) for the
non-Planck calls, and calls DISORT.  The DISORT Fortran itself is the
preserved production DISORT2 source reached through
oraclut.radiative_transfer.legacy_disort.disort (the Python equivalent of the
IDL DISORT DLM call).
"""

import numpy as np

from ..radiative_transfer.legacy_disort import disort


def call_disort(disort_vars, dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot,
                plank=False, wnlo=None, wnhi=None, temp=None, nlayer=None):
   """Return (rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed) from DISORT."""

   d = disort_vars
   dtauc = np.asarray(dtauc, dtype=np.float32)
   ssalb = np.asarray(ssalb, dtype=np.float32)
   pmom = np.asarray(pmom, dtype=np.float32)
   utau = np.asarray(utau, dtype=np.float32)
   umu = np.asarray(umu, dtype=np.float32)
   phi = np.asarray(phi, dtype=np.float32)

   # If the nlayer keyword has been set, we are only running on a sub-set of the
   # atmospheric layers...
   if nlayer is not None:
      nlyr1 = int(nlayer)
   else:
      nlyr1 = d.nlyr

   # Do some array size checks, just to make sure all is as expected
   if dtauc.size != nlyr1:
      raise ValueError("DTAUC dimension missmatch")
   if ssalb.size != nlyr1:
      raise ValueError("SSALB dimension missmatch")
   if pmom.shape[0] != d.nmom + 1:
      raise ValueError("PMOM dimension 1 missmatch")
   if pmom.shape[1] != nlyr1:
      raise ValueError("PMOM dimension 2 missmatch")
   if utau.size != d.ntau:
      raise ValueError("UTAU dimension missmatch")
   if umu.size != d.numu:
      raise ValueError("UMU dimension missmatch")
   if phi.size != d.nphi:
      raise ValueError("PHI dimension missmatch")

   # Set thermal related inputs based on whether the plank keyword was set or not
   if plank:
      if wnlo is None or wnhi is None:
         raise ValueError("Wavenumber must be specified for emission calculations!")
      if temp is None:
         raise ValueError("Temperature must be specified for emission calculations!")
      dtauc1 = dtauc
      ssalb1 = ssalb
      pmom1 = pmom
      wvnmlo = float(wnlo)
      wvnmhi = float(wnhi)
      temper = np.asarray(temp, dtype=np.float32)
      if temper.size != nlyr1 + 1:
         raise ValueError("TEMP dimension missmatch")
   else:
      # Compress the layer stack by combining layers with identical values for
      # SSALB and PMOM
      dtauc1 = np.zeros(nlyr1, dtype=np.float32)
      ssalb1 = np.zeros(nlyr1, dtype=np.float32)
      pmom1 = np.zeros((d.nmom + 1, nlyr1), dtype=np.float32, order="F")
      ii = 0
      dtauc1[0] = dtauc[0]
      ssalb1[0] = ssalb[0]
      pmom1[:, 0] = pmom[:, 0]
      for i in range(1, nlyr1):
         if ssalb1[ii] == ssalb[i] and np.array_equal(pmom1[:, ii], pmom[:, i]):
            dtauc1[ii] += dtauc[i]
         else:
            ii += 1
            dtauc1[ii] = dtauc[i]
            ssalb1[ii] = ssalb[i]
            pmom1[:, ii] = pmom[:, i]
      nlyr1 = ii + 1
      dtauc1 = dtauc1[:nlyr1]
      ssalb1 = ssalb1[:nlyr1]
      pmom1 = np.asfortranarray(pmom1[:, :nlyr1])
      wvnmlo = 0.0
      wvnmhi = 0.0
      temper = np.zeros(nlyr1, dtype=np.float32)

   # Call DISORT (NSTR, and the fixed flags of setup_disort, come from disort_vars / the adapter)
   result = disort(dtauc1, ssalb1, pmom1, utau, umu, phi, float(fbeam), float(umu0), float(fisot),
                   nstreams=d.nstr, plank=bool(plank), wavenumber_low=wvnmlo, wavenumber_high=wvnmhi,
                   temperature=temper)
   return (result["rfldir"], result["rfldn"], result["flup"], result["dfdt"], result["uavg"],
           result["uu"], result["albmed"], result["trnmed"])
