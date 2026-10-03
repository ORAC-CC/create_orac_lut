"""Python counterparts of generate_scattering_properties.pro and create_range.pro.

IDL:  pro generate_scattering_properties, srfstrarr, scatoffset, nwvl_max, inststr, mmstr, lutstr,
         nmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom,
         tmatrix_path=tmatrix_path, no_screen=no_screen
      pro create_range, Rat, MR, Spd, Radii, Nrat_o, NMR_o

Sequence, as in the IDL:
   1. interpolate each component's refractive index onto 0.55 microns and the
      channel wavelengths (AerM550, AerM);
   2. work out the mixing ratios / mode radii giving each LUT effective radius
      (create_range for log-normal components, the effective radius itself for
      a modified-gamma component);
   3. Gauss-Legendre points QV for the phase function (quadrature);
   4. for every component and effective radius call the scattering code
      (create_bwgp) at 0.55 microns and at the channel wavelengths;
   5. combine the components into the class properties weighted by mixing
      ratio and extinction, expand the phase functions in Legendre moments
      (legpexp) and form BextRat = Bext / Bext550.

Array index order follows the IDL: bext[m, l, r] is (SRF point, channel,
effective radius) and amom[p, m, l, r] adds the moment index first.

Deliberate difference from the IDL (2026-10, microphysical integration): for
liquid-water modified-gamma components the size integration stops at the
first node of the legacy radius lattice at or beyond 3.5 x effective radius
(beyond the IDL's fixed 100 um when necessary) instead of at 100 um; see
radius_upper_factor and create_bwgp.lattice_upper_radius.  Other components
(ice spheres, log-normal aerosols, Baum, T-matrix) are integrated as in the
IDL.

Deliberate difference from the IDL (V24, radius grid): for liquid-water and
ice-sphere modified-gamma components every interval of the legacy radius
trapezoid is divided into 2^k equal parts, with k chosen per wavelength and
effective radius so that the size-parameter step is at most
LIQUID_REFINED_XRES or ICE_SPHERE_REFINED_XRES; the limits and the legacy
nodes are unchanged.  See refined_xres and create_bwgp.radius_refinement_level.
Log-normal (aerosol) components keep the legacy grid.

Deliberate difference from the IDL (2026-10, Legendre expansion): for a class
whose components are all Mie size distributions, steps 3 and 5 no longer use
the IDL's fixed NMom (= 1000) for the quadrature order, the coefficient count
and the DISORT moments.  Each wavelength gets its own Gauss-Legendre order
Nq, above the polynomial-degree bound of its averaged phase function (from
the size integration actually performed), and every class phase function
(one per wavelength and effective radius) keeps the number of coefficients L
that King's criterion selects from its own expansion; see
legendre_expansion.py.  Baum and T-matrix (tabulated) classes keep the fixed
nmom treatment unchanged: their Legendre convergence has not been solved.
"""

import math

import numpy as np

from .create_bwgp import create_bwgp, legpexp, mie_integration_limits, quadrature, size_integration_limits
from .baum import BaumTable
from .legendre_expansion import (
   CHECK_THETA, KING_THRESHOLD, MAX_QUADRATURE_ATTEMPTS, MIE_MAX_ANGLES, QUADRATURE_NOISE_LIMIT,
   RECONSTRUCTION_TOLERANCE, expansion_length,
   gauss_legendre, initial_quadrature_order, legendre_coefficients, mie_phase_degree,
)

# Upper limit of the liquid-water modified-gamma size integration, as a
# multiple of effective radius, and the largest effective variance for which
# it has been validated (validation/size_distribution_limits/ and
# validation/REPORT_lut_numerics_development.md).  The fraction of the
# distribution beyond k re depends only on the effective variance and on the
# weighting: at 0.1111111 and k = 3.5 it is 6.6e-7 area-weighted (extinction
# by large particles) and 6.5e-5 r^6-weighted (scattering by particles small
# compared with the wavelength, the worst case), smaller for narrower
# distributions.  k = 3, the exploratory candidate, left up to 5.8e-4 in
# scattering and 2.5e-4 in g at r_e = 1-3 um in the thermal infrared.
LIQUID_UPPER_RADIUS_FACTOR = 3.5
LIQUID_MAX_EFFECTIVE_VARIANCE = 0.1111111

# Largest size-parameter step of the refined radius integration (V24), per
# class: the coarsest nested refinement of the legacy 0.4 step for which the
# radius-quadrature error of compact MODIS / dual-view SLSTR LUTs, against a
# 64x finer grid, is within 1.5e-3 in R_0v at every node (glory included;
# 0.3 x the 0.005 SLSTR noise-equivalent reflectance) and 2e-4 in extinction
# (relative) and g (validation/radius_grid/).  Ice spheres (m ~ 1.31) have
# weaker narrow resonances than liquid water (m ~ 1.33) and meet it at twice
# the step.
LIQUID_REFINED_XRES = 0.025
ICE_SPHERE_REFINED_XRES = 0.05


def radius_upper_factor(mmstr, c):
   """create_bwgp radius_upper_factor for component c (None = the IDL limits 0.001-100 um)."""

   if mmstr.substance.lower() != "liquid-water" or mmstr.distname[c] != "modified_gamma":
      return None
   if mmstr.s[c] > LIQUID_MAX_EFFECTIVE_VARIANCE * (1.0 + 1e-6):
      raise ValueError(f"{mmstr.compname[c]}: the {LIQUID_UPPER_RADIUS_FACTOR:g} x effective-radius upper limit is "
                       f"validated for effective variance <= {LIQUID_MAX_EFFECTIVE_VARIANCE}, not {mmstr.s[c]}")
   return LIQUID_UPPER_RADIUS_FACTOR


def refined_xres(mmstr, c):
   """create_bwgp refined_xres for component c (None = the IDL radius grid).

   Liquid-water and water-ice (sphere) modified-gamma Mie components only;
   other components keep the legacy grid (not validated here).
   """

   if mmstr.distname[c] != "modified_gamma" or str(mmstr.comp[c].code).lower() != "mie":
      return None
   return {"liquid-water": LIQUID_REFINED_XRES, "water-ice": ICE_SPHERE_REFINED_XRES}.get(mmstr.substance.lower())


def _interpol_complex(cm, wl, wvl):
   """IDL: interpol(mmstr.(cc).Cm, mmstr.(cc).wl, wvl) for the complex refractive index.

   The IDL interpolates a single-precision complex array.  The validated Python
   interpolates the real and imaginary parts in double precision and rounds each
   to single precision; that is retained (np.interp: no extrapolation is needed
   within the tabulated wavelength range).

   IDL INTERPOL accepts a monotonically decreasing abscissa, but np.interp
   requires an increasing one and silently returns an endpoint value otherwise.
   read_ri deliberately preserves the tabulated ordering of read_ri.pro, so a
   '#FORMAT = WAVN' table with ascending wavenumber arrives here as descending
   wavl = 1e4/wavn.  The abscissa and the refractive index are therefore sorted
   together into ascending wavelength at this interpolation boundary, which is
   where the Python requirement differs from the IDL; read_ri itself is left
   alone.  For an already ascending table the sort is the identity, so the
   previously validated results are unchanged bit for bit.

   No .ri table in the repository has duplicate wavelengths (checked by
   validation/refractive_index_ordering.py), so no duplicate abscissae are
   discarded or collapsed here.
   """

   wvl = np.asarray(wvl, dtype=np.float64)
   wl = np.asarray(wl, dtype=np.float64)
   cm = np.asarray(cm)
   order = np.argsort(wl, kind="stable")          # IDL: INTERPOL needs no such reordering
   wl = wl[order]
   cm = cm[order]
   real = np.interp(wvl, wl, np.real(cm)).astype(np.float32)
   imag = np.interp(wvl, wl, np.imag(cm)).astype(np.float32)
   return (real + 1j * imag).astype(np.complex64)


def create_range(rat, mr, spd, radii):
   """Mixing ratios and mode radii of the components giving each target effective radius.

   This is the IDL ``create_range.pro`` algorithm for one to four unique size
   modes. Components with the same effective radius are kept in a single mode;
   their within-mode proportions are restored when the mode result is expanded
   back to the component arrays.
   Returns (nrat_o, nmr_o) of shape (ncomponent, nradius).
   """

   rat = np.asarray(rat, dtype=np.float64)
   mr = np.asarray(mr, dtype=np.float64)
   spd = np.asarray(spd, dtype=np.float64)
   radii = np.asarray(radii, dtype=np.float64)

   if rat.ndim != 1 or mr.ndim != 1 or spd.ndim != 1:
      raise ValueError("create_range: Rat, MR and Spd must be one-dimensional")
   if not (rat.size == mr.size == spd.size) or rat.size == 0:
      raise ValueError("create_range: Rat, MR and Spd must have the same non-zero length")
   if radii.ndim != 1 or radii.size == 0:
      raise ValueError("create_range: Radii must be a non-empty one-dimensional array")
   if np.any(spd <= 0.0) or np.any(mr <= 0.0):
      raise ValueError("create_range: mode radii and spreads must be positive")
   if np.any(rat < 0.0) or not np.any(rat > 0.0):
      raise ValueError("create_range: Rat must be non-negative with a positive total")

   # IDL: Er = MR^3 * exp(4.5*alog(Spd)^2) /
   #           (MR^2 * exp(2.0*alog(Spd)^2))
   # IDL: Er = MR^3 * exp(4.5*alog(Spd)^2) / (Mr^2 * exp(2.0*alog(Spd)^2))
   er = mr**3 * np.exp(4.5 * np.log(spd)**2) / (mr**2 * np.exp(2.0 * np.log(spd)**2))

   # IDL: srtEr = sort(Er); unqEr = uniq(Er[srtEr]).  IDL UNIQ returns the
   # final index of each equal-valued run, which is also how the group slices
   # below are written in the original routine.
   srt_er = np.argsort(er, kind="stable")
   sorted_er = er[srt_er]
   unq_er = np.concatenate((np.flatnonzero(sorted_er[1:] != sorted_er[:-1]), [er.size - 1]))
   nunq = int(unq_er.size)
   if nunq > 4:
      raise ValueError("create_range: more than four unique size modes are not supported")

   erq = er[srt_er[unq_er]]
   spq = spd[srt_er[unq_er]]
   mrq = mr[srt_er[unq_er]]

   # IDL: Ratq is the total mixing ratio of each unique size mode and
   # Mode_Ratq retains each component's proportion within that mode.
   ratq = np.zeros(nunq, dtype=np.float64)
   max_mode_count = int(np.max(np.diff(np.concatenate(([-1], unq_er)))))
   mode_ratq = np.zeros((max_mode_count, nunq), dtype=np.float64)
   group_starts = np.concatenate(([0], unq_er[:-1] + 1))
   for mode, (start, end) in enumerate(zip(group_starts, unq_er)):
      members = srt_er[start:end + 1]
      ratq[mode] = np.sum(rat[members])
      if ratq[mode] <= 0.0:
         raise ValueError("create_range: every size mode must have a positive total mixing ratio")
      mode_ratq[:members.size, mode] = rat[members] / ratq[mode]

   # IDL: Nrat/NMR hold the unique-mode result before it is expanded to all
   # original components.  g and f are the log-normal moment factors.
   nrat = np.zeros((nunq, radii.size), dtype=np.float64)
   nmr = np.zeros((nunq, radii.size), dtype=np.float64)
   g = np.exp(2.0 * np.log(spq)**2)
   f = np.exp(4.5 * np.log(spq)**2)
   ea = np.sum(rat * mr**3 * np.exp(4.5 * np.log(spd)**2)) / np.sum(
      rat * mr**2 * np.exp(2.0 * np.log(spd)**2)
   )

   for i, radius in enumerate(radii):
      # IDL: Nmr[*,I] = Mrq; this is overwritten in the boundary branches.
      nmr[:, i] = mrq
      if nunq == 1:
         nrat[:, i] = 1.0
         nmr[0, i] = mrq[0] * radius / erq[0]
      elif nunq == 2:
         if radius <= erq[0]:
            nrat[:, i] = [1.0, 0.0]
            nmr[:, i] = [mrq[0] * radius / erq[0], mrq[1]]
         elif radius < erq[1]:
            den1 = nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])
            den2 = -nmr[1, i]**2 * (radius * g[1] - nmr[1, i] * f[1])
            nrat[0, i] = den2 / (den1 + den2)
            nrat[1, i] = 1.0 - nrat[0, i]
         else:
            nrat[:, i] = [0.0, 1.0]
            nmr[:, i] = [mrq[0], mrq[1] * radius / erq[1]]
      elif nunq == 3:
         if radius <= erq[0]:
            nrat[:, i] = [1.0, 0.0, 0.0]
            nmr[:, i] = [mrq[0] * radius / erq[0], mrq[1], mrq[2]]
         elif radius < ea:
            den1 = (ratq[2] / ratq[1]) * nmr[2, i]**2 * (radius * g[2] - nmr[2, i] * f[2])
            den2 = nmr[1, i]**2 * (radius * g[1] - nmr[1, i] * f[1])
            den3 = -(nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])) * (ratq[2] / ratq[1] + 1.0)
            nrat[1, i] = (-nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])) / (den1 + den2 + den3)
            nrat[0, i] = 1.0 - nrat[1, i] * (ratq[2] / ratq[1] + 1.0)
            nrat[2, i] = nrat[1, i] * ratq[2] / ratq[1]
         elif radius < erq[2]:
            den1 = (ratq[0] / ratq[1]) * nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])
            den2 = nmr[1, i]**2 * (radius * g[1] - nmr[1, i] * f[1])
            den3 = -(nmr[2, i]**2 * (radius * g[2] - nmr[2, i] * f[2])) * (ratq[0] / ratq[1] + 1.0)
            nrat[1, i] = (-nmr[2, i]**2 * (radius * g[2] - nmr[2, i] * f[2])) / (den1 + den2 + den3)
            nrat[0, i] = nrat[1, i] * ratq[0] / ratq[1]
            nrat[2, i] = 1.0 - nrat[1, i] * (ratq[0] / ratq[1] + 1.0)
         else:
            nrat[:, i] = [0.0, 0.0, 1.0]
            nmr[:, i] = [mrq[0], mrq[1], mrq[2] * radius / erq[2]]
      else:  # nunq == 4; direct transcription of the two interior IDL cases.
         if radius <= erq[0]:
            nrat[:, i] = [1.0, 0.0, 0.0, 0.0]
            nmr[:, i] = [mrq[0] * radius / erq[0], mrq[1], mrq[2], mrq[3]]
         elif radius < ea:
            den1 = nmr[3, i]**2 * (radius * g[3] - nmr[3, i] * f[3])
            den2 = (ratq[1] / ratq[3]) * nmr[1, i]**2 * (radius * g[1] - nmr[1, i] * f[1])
            den3 = (ratq[2] / ratq[3]) * nmr[2, i]**2 * (radius * g[2] - nmr[2, i] * f[2])
            den4 = -(1.0 + ratq[1] / ratq[3] + ratq[2] / ratq[3]) * nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])
            nrat[3, i] = -(nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])) / (den1 + den2 + den3 + den4)
            nrat[1, i] = ratq[1] / ratq[3] * nrat[3, i]
            nrat[2, i] = ratq[2] / ratq[3] * nrat[3, i]
            nrat[0, i] = 1.0 - nrat[3, i] - nrat[2, i] - nrat[1, i]
         elif radius < erq[3]:
            den1 = nmr[0, i]**2 * (radius * g[0] - nmr[0, i] * f[0])
            den2 = (ratq[1] / ratq[0]) * nmr[1, i]**2 * (radius * g[1] - nmr[1, i] * f[1])
            den3 = (ratq[2] / ratq[0]) * nmr[2, i]**2 * (radius * g[2] - nmr[2, i] * f[2])
            den4 = -(1.0 + ratq[1] / ratq[0] + ratq[2] / ratq[0]) * nmr[3, i]**2 * (radius * g[3] - nmr[3, i] * f[3])
            nrat[0, i] = -(nmr[3, i]**2 * (radius * g[3] - nmr[3, i] * f[3])) / (den1 + den2 + den3 + den4)
            nrat[1, i] = ratq[1] / ratq[0] * nrat[0, i]
            nrat[2, i] = ratq[2] / ratq[0] * nrat[0, i]
            nrat[3, i] = 1.0 - nrat[0, i] - nrat[2, i] - nrat[1, i]
         else:
            nrat[:, i] = [0.0, 0.0, 0.0, 1.0]
            nmr[:, i] = [mrq[0], mrq[1], mrq[2], mrq[3] * radius / erq[3]]

   # IDL: expand unique-mode outputs to all original components, restoring
   # Mode_Ratq and the original order recorded by srtEr.
   nrat_o = np.zeros((rat.size, radii.size), dtype=np.float64)
   nmr_o = np.zeros((rat.size, radii.size), dtype=np.float64)
   for mode, (start, end) in enumerate(zip(group_starts, unq_er)):
      members = srt_er[start:end + 1]
      nrat_o[members, :] = nrat[mode, :][None, :] * mode_ratq[:members.size, mode, None]
      nmr_o[members, :] = nmr[mode, :]
   return nrat_o, nmr_o


def uses_adaptive_legendre(mmstr):
   """True when every component is a Mie size distribution (adaptive Legendre expansion).

   Baum crystals and Dubovik T-matrix components are interpolated tables
   rather than band-limited polynomials, so King's criterion cannot
   terminate their series; those classes keep the fixed nmom treatment.
   create_bwgp uses Mie unless the scattering code is 'tmatrix'.
   """

   for c in range(mmstr.ncomp):
      if mmstr.comptype[c].lower() not in ("opac", "user"):
         return False
      if str(mmstr.comp[c].code).lower() == "tmatrix":
         return False
   return True


def generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, nmom, tmatrix_path=None):
   """Return (lmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom).

   ``lmom[m, l, r]`` is the number of Legendre coefficients amom[0:lmom, m, l, r]
   to pass to DISORT; amom is zero beyond it.

   Mie classes (uses_adaptive_legendre): ``nmom`` must be None.  The
   quadrature order and the expansion length are determined from each
   size-distribution-averaged phase function (legendre_expansion.py), and phs
   / phs550 hold the phase function on each wavelength's own Gauss-Legendre
   nodes, zero beyond that wavelength's quadrature order.

   Baum and T-matrix classes: ``nmom`` is the fixed number of quadrature
   points and Legendre moments (the IDL n_theta keyword, 1000 in the IDL),
   and lmom is nmom everywhere.
   """

   nchan = inststr.number_of_nadir_channels
   ncomp = mmstr.ncomp
   nefr = lutstr.efr_n

   # Wavelengths of the SRF quadrature points: IDL srfstrarr[*].wvl[*] is a
   # (nwvl_max, nchan) array indexed [m, l].
   wvl = np.stack([srfstrarr[l].wvl for l in range(nchan)], axis=1).astype(np.float64)

   # Interpolate the components' refractive index values onto the 0.55 micron
   # reference wavelength and the instrument channels.
   aerm550 = np.zeros(ncomp, dtype=np.complex64)
   aerm = np.zeros((nwvl_max, nchan, ncomp), dtype=np.complex64)
   if mmstr.comptype[0].lower() in ("opac", "user"):
      for c in range(ncomp):
         aerm550[c] = _interpol_complex(mmstr.comp[c].cm, mmstr.comp[c].wl, [0.55])[0]
         aerm[:, :, c] = _interpol_complex(mmstr.comp[c].cm, mmstr.comp[c].wl, wvl)
   # (force_n / force_k refractive-index overrides: not ported, never used by the validated runs)

   # Generate the range of component mixing ratios / mode radii required to
   # provide the required effective radii
   print(mmstr.comptype[0])
   component_types = [component_type.lower() for component_type in mmstr.comptype]
   if all(component_type == "baum" for component_type in component_types):
      # Baum properties are already tabulated against effective radius; the
      # legacy path therefore uses the LUT radius directly and does not call
      # create_range or create_bwgp.
      lut_mrat = np.repeat(mmstr.mrat[:, None], nefr, axis=1)
      lut_rm = np.repeat(lutstr.efr[None, :], ncomp, axis=0)
   elif mmstr.comptype[0] in ("opac", "user"):
      print(mmstr.distname[0])
      if mmstr.distname[0] == "log_normal":
         lut_mrat, lut_rm = create_range(mmstr.mrat, mmstr.rm, mmstr.s, lutstr.efr)
      elif mmstr.distname[0] == "modified_gamma":
         if ncomp != 1:
            raise NotImplementedError("modified_gamma is ported for a single component only")
         lut_mrat = np.zeros((ncomp, nefr), dtype=np.float64)
         lut_rm = np.zeros((ncomp, nefr), dtype=np.float64)
         lut_mrat[0, :] = mmstr.mrat[0]
         lut_rm[0, :] = lutstr.efr                     # the mode "radius" is the effective radius itself
      else:
         raise ValueError("Unknown size distribution: " + str(mmstr.distname[0]))
   else:
      raise NotImplementedError("Baran ice-crystal components are not ported to Python")

   # **** Mie classes: the quadrature order and Legendre expansion length are
   #      determined from each size-distribution-averaged phase function
   #      (no IDL counterpart; replaces the fixed NMom)
   if uses_adaptive_legendre(mmstr):
      if nmom is not None:
         raise ValueError("nmom is obsolete for Mie size distributions: the Legendre expansion length is determined "
                          "from each averaged phase function (legendre_expansion.py); remove nmom from the run file")
      return _adaptive_mie_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, wvl, aerm550, aerm, lut_mrat, lut_rm)
   if nmom is None:
      raise ValueError("nmom is required for Baum and T-matrix (tabulated) phase functions, whose Legendre "
                       "convergence is not covered by the adaptive Mie expansion")

   # **** Generate the quadrature points for the scattering phase function
   #      (IDL: quadrature, 'g', NMom, Abscissas, Weights; QV = cos(scattering angle))
   abscissas, weights = quadrature("g", nmom)
   qv0 = 1.0                                           # theta = 0
   qv1 = -1.0                                          # theta = 180
   qv = ((qv1 - qv0) * abscissas + (qv0 + qv1)) / 2.0

   # **** Call the scattering code for the required range of components and mode radii
   vavg_c = np.zeros((ncomp, nefr))                    # average volume per particle
   # at the reference wavelength
   bext550_c = np.zeros((ncomp, nefr))                 # extinction coefficient
   w550_c = np.zeros((ncomp, nefr))                    # single scattering albedo
   g550_c = np.zeros((ncomp, nefr))                    # asymmetry parameter
   phs550_c = np.zeros((nmom, ncomp, nefr))            # phase function
   # at the individual channels
   bext_c = np.zeros((nwvl_max, nchan, ncomp, nefr))
   w_c = np.zeros((nwvl_max, nchan, ncomp, nefr))
   g_c = np.zeros((nwvl_max, nchan, ncomp, nefr))
   phs_c = np.zeros((nmom, nwvl_max, nchan, ncomp, nefr))

   for c in range(ncomp):
      if component_types[c] == "baum":
         table = BaumTable(mmstr.compname[c])
         theta = np.degrees(np.arccos(qv))
         bext1, w1, g1, phs1 = table.interpolate([0.55], lutstr.efr, theta)
         bext550_c[c, :] = bext1[0, :]
         w550_c[c, :] = w1[0, :]
         g550_c[c, :] = g1[0, :]
         phs550_c[:, c, :] = phs1[:, 0, :]

         # IDL normalizes each interpolated Baum phase function with the
         # Gauss-Legendre weights before expanding it into moments.
         norm550 = np.sum(phs550_c[:, c, :] * weights[:, None], axis=0) / 2.0
         if np.any(norm550 == 0.0) or not np.all(np.isfinite(norm550)):
            raise ValueError(f"{mmstr.compname[c]}: invalid Baum phase-function normalization")
         phs550_c[:, c, :] /= norm550[None, :]

         wavelength_values = wvl.ravel(order="F")
         bext1, w1, g1, phs1 = table.interpolate(wavelength_values, lutstr.efr, theta)
         bext_c[:, :, c, :] = bext1.reshape((nwvl_max, nchan, nefr), order="F")
         w_c[:, :, c, :] = w1.reshape((nwvl_max, nchan, nefr), order="F")
         g_c[:, :, c, :] = g1.reshape((nwvl_max, nchan, nefr), order="F")
         phs_c[:, :, :, c, :] = phs1.reshape((nmom, nwvl_max, nchan, nefr), order="F")

         norm = np.sum(phs_c[:, :, :, c, :] * weights[:, None, None, None], axis=0) / 2.0
         if np.any(norm == 0.0) or not np.all(np.isfinite(norm)):
            raise ValueError(f"{mmstr.compname[c]}: invalid Baum phase-function normalization")
         phs_c[:, :, :, c, :] /= norm[None, :, :, :]
         continue

      for r in range(nefr):
         print(" Performing scattering calculations for component " + mmstr.compname[c] + ", EfR: ", lutstr.efr[r])
         # The mode radius of a component does not always change from one
         # effective radius to the next; only call the scattering code if it has.
         calculated = False
         if r > 0:
            if lut_rm[c, r] == lut_rm[c, r - 1] and bext_c[0, 0, c, r - 1] != 0:
               calculated = True
         if calculated:
            vavg_c[c, r] = vavg_c[c, r - 1]
            bext550_c[c, r] = bext550_c[c, r - 1]
            w550_c[c, r] = w550_c[c, r - 1]
            g550_c[c, r] = g550_c[c, r - 1]
            phs550_c[:, c, r] = phs550_c[:, c, r - 1]
            bext_c[:, :, c, r] = bext_c[:, :, c, r - 1]
            w_c[:, :, c, r] = w_c[:, :, c, r - 1]
            g_c[:, :, c, r] = g_c[:, :, c, r - 1]
            phs_c[:, :, :, c, r] = phs_c[:, :, :, c, r - 1]
         elif lut_mrat[c, r] > 0:
            scode = mmstr.comp[c].code
            eps = getattr(mmstr.comp[c], "eps", None)
            neps = getattr(mmstr.comp[c], "neps", None)
            factor = radius_upper_factor(mmstr, c)              # IDL: always 0.001-100 um
            refine = refined_xres(mmstr, c)                    # IDL: always the xres = 0.4 grid
            # Calculate Bext at 550 nm (the reference wavelength) and Vavg
            bext1, w1, g1, phs1, vavg1 = create_bwgp(mmstr.distname[c], lut_rm[c, r], mmstr.s[c], aerm550[c], 0.55, qv,
                                                     scode=scode, tmatrix_path=tmatrix_path, eps=eps, neps=neps,
                                                     radius_upper_factor=factor, refined_xres=refine)
            vavg_c[c, r] = vavg1
            bext550_c[c, r] = bext1[0]
            w550_c[c, r] = w1[0]
            g550_c[c, r] = g1[0]
            phs550_c[:, c, r] = phs1[:, 0]
            # Calculate bext, w, g and phs for each instrument channel.  The IDL
            # passes the (nwvl_max, nchan) arrays flattened column-major and
            # reforms the results; the same ordering is used here.
            bext1, w1, g1, phs1, _ = create_bwgp(mmstr.distname[c], lut_rm[c, r], mmstr.s[c],
                                                 aerm[:, :, c].ravel(order="F"), wvl.ravel(order="F"), qv,
                                                 scode=scode, tmatrix_path=tmatrix_path, eps=eps, neps=neps,
                                                 radius_upper_factor=factor, refined_xres=refine)
            bext_c[:, :, c, r] = bext1.reshape((nwvl_max, nchan), order="F")
            w_c[:, :, c, r] = w1.reshape((nwvl_max, nchan), order="F")
            g_c[:, :, c, r] = g1.reshape((nwvl_max, nchan), order="F")
            phs_c[:, :, :, c, r] = phs1.reshape((nmom, nwvl_max, nchan), order="F")

   print("")
   print("All scattering calculations completed for each component")

   # **** Calculate the scattering parameters of the class for each required
   #      effective radius (weighted by mixing ratio and extinction)
   vavg = np.zeros(nefr)
   # At the reference wavelength
   bext550 = np.zeros(nefr)                            # Extinction coefficient
   w550 = np.zeros(nefr)                               # Single scattering albedo
   g550 = np.zeros(nefr)                               # Asymmetry parameter
   phs550 = np.zeros((nmom, nefr))                     # Phase function
   amom550 = np.zeros((nmom, nefr))                    # Legendre moments
   # For each channel
   bextrat = np.zeros((nwvl_max, nchan, nefr))         # Ratio of Bext with that at the reference wavelength
   bext = np.zeros((nwvl_max, nchan, nefr))            # Extinction coefficient
   w = np.zeros((nwvl_max, nchan, nefr))               # Single scattering albedo
   g = np.zeros((nwvl_max, nchan, nefr))               # Asymmetry parameter
   phs = np.zeros((nmom, nwvl_max, nchan, nefr))       # Phase function
   amom = np.zeros((nmom, nwvl_max, nchan, nefr))      # Legendre moments

   vavg[:] = vavg_c[0, :]
   two_n_plus_one = 2.0 * np.arange(nmom, dtype=np.float64) + 1.0

   for l in range(nchan):
      for m in range(nwvl_max):
         # QM=2 stores each channel in a common rectangular array and pads
         # shorter SRFs with zero entries.  Those entries are not scientific
         # quadrature points and must remain the initialized zero outputs;
         # processing them would recreate the invalid 0/0 divisions reported
         # by the EarthCARE QM=2 preflight.
         if m >= srfstrarr[l].nwvl:
            continue
         for r in range(nefr):
            # Calculate the 550 nm quantities (once, with the first channel)
            if l == 0:
               mratbext = lut_mrat[:, r] * bext550_c[:, r]
               tmratbext = np.sum(mratbext)
               mratbextw = mratbext * w550_c[:, r]
               tmratbextw = np.sum(mratbextw)
               bext550[r] = tmratbext / np.sum(lut_mrat[:, r])
               w550[r] = tmratbextw / tmratbext
               g550[r] = np.sum(mratbextw * g550_c[:, r]) / tmratbextw
               phs550[:, r] = np.sum(mratbextw[None, :] * phs550_c[:, :, r], axis=1) / tmratbextw
               # Calculate the Legendre moments for the phase function
               inlc, alc = legpexp(nmom, qv, weights, phs550[:, r])
               amom550[:, r] = alc / two_n_plus_one

            mratbext = lut_mrat[:, r] * bext_c[m, l, :, r]
            tmratbext = np.sum(mratbext)
            mratbextw = mratbext * w_c[m, l, :, r]
            tmratbextw = np.sum(mratbextw)
            bext[m, l, r] = tmratbext / np.sum(lut_mrat[:, r])
            w[m, l, r] = tmratbextw / tmratbext
            g[m, l, r] = np.sum(mratbextw * g_c[m, l, :, r]) / tmratbextw
            phs[:, m, l, r] = np.sum(mratbextw[None, :] * phs_c[:, m, l, :, r], axis=1) / tmratbextw

            # Calculate the Legendre moments for the phase function
            inlc, alc = legpexp(nmom, qv, weights, phs[:, m, l, r])
            amom[:, m, l, r] = alc / two_n_plus_one

            # Ratio of the extinction coefficient at the current channel and at
            # 0.55 microns, relating the spectral optical depth to the reference
            bextrat[m, l, r] = bext[m, l, r] / bext550[r]

   # The IDL holds all of these in FLTARR (single precision); round here, after
   # the double-precision combination, exactly as the validated Python did.
   f32 = np.float32
   lmom = np.full((nwvl_max, nchan, nefr), nmom, dtype=np.int64)    # fixed expansion length
   return (lmom, bext550.astype(f32), w550.astype(f32), g550.astype(f32), phs550.astype(f32),
           amom550.astype(f32), bextrat.astype(f32), bext.astype(f32), w.astype(f32), g.astype(f32),
           vavg.astype(f32), phs.astype(f32), amom.astype(f32))


def _adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, ri, wl, label):
   """Mie scattering at one wavelength, with an adaptive Legendre expansion of each class phase function.

   ri[c] is the refractive index of component c at wavelength wl (microns).
   Each effective radius r is treated separately.  Its degree bound D_r comes
   from the largest size parameter of its own size integration (all
   components with a non-zero mixing ratio at r), so its coefficients beyond
   D_r are exactly zero; its Gauss-Legendre order Nq_r is chosen above D_r
   (legendre_expansion.initial_quadrature_order), so none is aliased, and
   small radii are not sampled at the order a large radius needs.  The
   expansion length L_r comes from King's criterion applied to the averaged
   phase function up to D_r (the rounding noise is measured just above D_r),
   and is accepted only if the series reproduces the directly calculated
   phase function at the CHECK_THETA angles to six significant figures and
   the noise is below QUADRATURE_NOISE_LIMIT.  Otherwise Nq_r is increased and
   the radius recalculated, and the run fails after MAX_QUADRATURE_ATTEMPTS.

   Returns (bext_c, w_c, g_c, vavg_c, phs, omega, lmom):
      bext_c, w_c, g_c, vavg_c   (ncomp, nefr) component properties, as create_bwgp
      phs                        (max Nq_r, nefr) class phase function at each radius's
                                 Gauss-Legendre nodes, zero beyond Nq_r
      omega                      (max Nq_r, nefr) its Legendre coefficients, omega_l
                                 (2l+1 included), zero beyond Nq_r
      lmom                       (nefr,) number of coefficients L_r retained
   """

   ncomp, nefr = lut_mrat.shape
   wavenumber = 1.0 / wl
   factors = [radius_upper_factor(mmstr, c) for c in range(ncomp)]
   refine = [refined_xres(mmstr, c) for c in range(ncomp)]
   check_mu = np.cos(np.deg2rad(CHECK_THETA))

   bext_c = np.zeros((ncomp, nefr))
   w_c = np.zeros((ncomp, nefr))
   g_c = np.zeros((ncomp, nefr))
   vavg_c = np.zeros((ncomp, nefr))
   g_class = np.zeros(nefr)
   lmom = np.zeros(nefr, dtype=np.int64)
   noise = np.zeros(nefr)
   reconstruction_error = np.zeros(nefr)
   degree_r = np.zeros(nefr, dtype=np.int64)
   nq_r = np.zeros(nefr, dtype=np.int64)
   phases = []
   coefficients = []
   previous = {}                                       # component -> (mode radius, Nq, results at the previous radius)

   for r in range(nefr):
      # Largest size parameter reached by this radius's size integration, with
      # the limits create_bwgp actually uses (the liquid-water lattice limit,
      # 100 um, or the log-normal quantile limits)
      x_max = 0.0
      for c in range(ncomp):
         if lut_mrat[c, r] > 0:
            params, npts = mie_integration_limits(mmstr.distname[c], lut_rm[c, r], mmstr.s[c], wavenumber, factors[c])
            rl, ru, truncated = size_integration_limits(mmstr.distname[c], params, wavenumber)
            x_max = max(x_max, 2.0 * np.pi * ru * wavenumber)
      degree = mie_phase_degree(x_max)
      nq = initial_quadrature_order(degree)

      for attempt in range(MAX_QUADRATURE_ATTEMPTS):
         if nq + check_mu.size > MIE_MAX_ANGLES:
            raise ValueError(f"{label}: the required quadrature order {nq} exceeds the Mie kernel limit of "
                             f"{MIE_MAX_ANGLES} angles (size parameter {x_max:.0f})")
         # Gauss-Legendre points for this radius; the check angles are
         # calculated by the same Mie call (IDL: QV = cos(scattering angle))
         abscissas, weights = gauss_legendre(nq)
         qv = -abscissas                               # IDL: QV = ((QV1-QV0)*Abscissas + (QV0+QV1))/2, QV0 = 1, QV1 = -1
         dqv = np.concatenate((qv, check_mu))
         phs_c = np.zeros((dqv.size, ncomp))
         for c in range(ncomp):
            # Only call the scattering code if the mode radius (and order) has changed
            if c in previous and previous[c][0] == lut_rm[c, r] and previous[c][1] == nq:
               bext_c[c, r], w_c[c, r], g_c[c, r], vavg_c[c, r], phs_c[:, c] = previous[c][2]
            elif lut_mrat[c, r] > 0:
               bext1, w1, g1, phs1, vavg1 = create_bwgp(mmstr.distname[c], lut_rm[c, r], mmstr.s[c], np.asarray([ri[c]]),
                                                        np.asarray([wl]), dqv, scode=mmstr.comp[c].code,
                                                        radius_upper_factor=factors[c], refined_xres=refine[c])
               bext_c[c, r], w_c[c, r], g_c[c, r], vavg_c[c, r] = bext1[0], w1[0], g1[0], vavg1
               phs_c[:, c] = phs1[:, 0]
               previous[c] = (lut_rm[c, r], nq, (bext_c[c, r], w_c[c, r], g_c[c, r], vavg_c[c, r], phs_c[:, c].copy()))

         # The class phase function and asymmetry parameter, weighted by mixing
         # ratio and extinction as in the class combination below
         mratbextw = lut_mrat[:, r] * bext_c[:, r] * w_c[:, r]
         tmratbextw = np.sum(mratbextw)
         phs = np.sum(mratbextw[None, :] * phs_c, axis=1) / tmratbextw
         g_class[r] = np.sum(mratbextw * g_c[:, r]) / tmratbextw

         # Legendre coefficients, King's criterion and the reconstruction test.
         # The quadrature is adequate when the coefficients beyond the degree
         # bound are at the rounding-noise level and the double-precision
         # series reproduces the directly calculated phase function.
         omega = legendre_coefficients(qv, weights, phs[:nq])
         lmom[r], noise[r], reconstruction_error[r] = expansion_length(omega, degree, check_mu, phs[nq:])
         if reconstruction_error[r] <= RECONSTRUCTION_TOLERANCE and noise[r] <= QUADRATURE_NOISE_LIMIT:
            break
         print(f" Legendre expansion {label}, radius {lut_rm[0, r]:g}: quadrature adequacy not demonstrated with "
               f"Nq = {nq} (reconstruction error {reconstruction_error[r]:.1e}, noise beyond D {noise[r]:.1e}); "
               "increasing the quadrature order")
         nq = int(math.ceil(1.5 * nq))
      else:
         raise RuntimeError(f"{label}: quadrature adequacy not demonstrated after {MAX_QUADRATURE_ATTEMPTS} orders: "
                            f"reconstruction error {reconstruction_error[r]:.1e} (limit {RECONSTRUCTION_TOLERANCE:g}), "
                            f"noise beyond D {noise[r]:.1e} (limit {QUADRATURE_NOISE_LIMIT:g})")
      degree_r[r] = degree
      nq_r[r] = nq
      phases.append(phs[:nq])
      coefficients.append(omega)

   nq_max = int(np.max(nq_r))
   phs_out = np.zeros((nq_max, nefr))
   omega_out = np.zeros((nq_max, nefr))
   for r in range(nefr):
      phs_out[:nq_r[r], r] = phases[r]
      omega_out[:nq_r[r], r] = coefficients[r]

   # Thesis diagnostics (Grainger 1990 section 4.5): omega_0 = 1 and omega_1/3 = g.
   # With Gauss-Legendre in mu they can pass while higher moments are aliased,
   # so they are reported, not used for acceptance.  "noise-limited" counts
   # the radii whose rounding noise exceeds King's 1e-9 (L conservative, <= D + 1).
   print(f" Legendre expansion {label}: degree bounds D_r = {degree_r.min()}-{degree_r.max()}, "
         f"Nq_r = {nq_r.min()}-{nq_r.max()}, L = {lmom.min()}-{lmom.max()}, noise beyond D <= {np.max(noise):.1e} "
         f"(noise-limited: {int(np.sum(noise >= KING_THRESHOLD))} of {nefr}), "
         f"reconstruction error <= {np.max(reconstruction_error):.1e}, "
         f"|omega_0 - 1| <= {np.max(np.abs(omega_out[0, :] - 1.0)):.1e}, "
         f"|omega_1/3 - g| <= {np.max(np.abs(omega_out[1, :] / 3.0 - g_class)):.1e}")
   return bext_c, w_c, g_c, vavg_c, phs_out, omega_out, lmom


def _adaptive_mie_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, wvl, aerm550, aerm, lut_mrat, lut_rm):
   """generate_scattering_properties for a class of Mie components, with adaptive Legendre expansions.

   The bulk properties are combined exactly as in the fixed-nmom path (and are
   numerically identical to it); only the phase-function sampling and the
   Legendre moments differ.
   """

   nchan = inststr.number_of_nadir_channels
   ncomp, nefr = lut_mrat.shape

   # **** Call the scattering code at 0.55 microns and at each SRF point of
   #      each channel, each with its own quadrature order
   print("")
   (bext550_c, w550_c, g550_c, vavg_c, phs550_n, omega550, lmom550) = \
      _adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, aerm550, 0.55, "0.55 um reference")
   channel_results = {}
   for l in range(nchan):
      for m in range(srfstrarr[l].nwvl):
         label = f"channel {inststr.channelid[l]} point {m} ({wvl[m, l]:.4f} um)"
         channel_results[(m, l)] = _adaptive_mie_wavelength(mmstr, lut_mrat, lut_rm, aerm[m, l, :], wvl[m, l], label)
   print("")
   print("All scattering calculations completed for each component")

   # **** Calculate the scattering parameters of the class for each required
   #      effective radius (weighted by mixing ratio and extinction)
   lmax = max(int(np.max(result[6])) for result in channel_results.values())
   nqmax = max(result[4].shape[0] for result in channel_results.values())
   vavg = np.zeros(nefr)
   # At the reference wavelength
   bext550 = np.zeros(nefr)                            # Extinction coefficient
   w550 = np.zeros(nefr)                               # Single scattering albedo
   g550 = np.zeros(nefr)                               # Asymmetry parameter
   phs550 = np.zeros((phs550_n.shape[0], nefr))        # Phase function (0.55 um Gauss-Legendre nodes)
   amom550 = np.zeros((int(np.max(lmom550)), nefr))    # Legendre moments
   # For each channel
   bextrat = np.zeros((nwvl_max, nchan, nefr))         # Ratio of Bext with that at the reference wavelength
   bext = np.zeros((nwvl_max, nchan, nefr))            # Extinction coefficient
   w = np.zeros((nwvl_max, nchan, nefr))               # Single scattering albedo
   g = np.zeros((nwvl_max, nchan, nefr))               # Asymmetry parameter
   phs = np.zeros((nqmax, nwvl_max, nchan, nefr))      # Phase function (each wavelength's own nodes)
   amom = np.zeros((lmax, nwvl_max, nchan, nefr))      # Legendre moments
   lmom = np.zeros((nwvl_max, nchan, nefr), dtype=np.int64)   # Number of moments retained

   vavg[:] = vavg_c[0, :]

   for l in range(nchan):
      for m in range(nwvl_max):
         # QM=2 zero-padded SRF entries are not quadrature points (see the fixed path)
         if m >= srfstrarr[l].nwvl:
            continue
         bext_c, w_c, g_c, _, phs_n, omega, lmom_ml = channel_results[(m, l)]
         for r in range(nefr):
            # Calculate the 550 nm quantities (once, with the first channel)
            if l == 0:
               mratbext = lut_mrat[:, r] * bext550_c[:, r]
               tmratbext = np.sum(mratbext)
               mratbextw = mratbext * w550_c[:, r]
               tmratbextw = np.sum(mratbextw)
               bext550[r] = tmratbext / np.sum(lut_mrat[:, r])
               w550[r] = tmratbextw / tmratbext
               g550[r] = np.sum(mratbextw * g550_c[:, r]) / tmratbextw
               phs550[:, r] = phs550_n[:, r]
               # The Legendre moments of the phase function (DISORT convention chi_l = omega_l / (2l+1))
               amom550[:lmom550[r], r] = omega550[:lmom550[r], r] / (2.0 * np.arange(lmom550[r]) + 1.0)

            mratbext = lut_mrat[:, r] * bext_c[:, r]
            tmratbext = np.sum(mratbext)
            mratbextw = mratbext * w_c[:, r]
            tmratbextw = np.sum(mratbextw)
            bext[m, l, r] = tmratbext / np.sum(lut_mrat[:, r])
            w[m, l, r] = tmratbextw / tmratbext
            g[m, l, r] = np.sum(mratbextw * g_c[:, r]) / tmratbextw
            phs[:phs_n.shape[0], m, l, r] = phs_n[:, r]

            # The Legendre moments of the phase function (DISORT convention chi_l = omega_l / (2l+1))
            lmom[m, l, r] = lmom_ml[r]
            amom[:lmom_ml[r], m, l, r] = omega[:lmom_ml[r], r] / (2.0 * np.arange(lmom_ml[r]) + 1.0)

            # Ratio of the extinction coefficient at the current channel and at
            # 0.55 microns, relating the spectral optical depth to the reference
            bextrat[m, l, r] = bext[m, l, r] / bext550[r]

   # Single precision, as the fixed path (IDL FLTARR)
   f32 = np.float32
   return (lmom, bext550.astype(f32), w550.astype(f32), g550.astype(f32), phs550.astype(f32),
           amom550.astype(f32), bextrat.astype(f32), bext.astype(f32), w.astype(f32), g.astype(f32),
           vavg.astype(f32), phs.astype(f32), amom.astype(f32))
