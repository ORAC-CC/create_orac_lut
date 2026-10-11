"""Python counterpart of setup_disort.pro.

IDL:  pro setup_disort, NStreams, NLayers, NSatzen, NRelAzi, NMoments, alb=alb

The IDL fills the common block DISORT_VARS with the DISORT arguments that are
the same for every call.  Here the same names are returned in a namespace,
``disort_vars``, which is passed to call_disort().  The DISORT C adapter
(src/oraclut/radiative_transfer/legacy_disort.c) fixes USRTAU, USRANG, IBCND,
LAMBER, ALBEDO = 0, BTEMP, TTEMP, TEMIS, ONLYFL, ACCUR = 1e-8 and PHI0 at the
values listed below, so those are recorded here for reference only.
"""

from types import SimpleNamespace


def setup_disort(nstreams, nlayers, nsatzen, nrelazi, nmoments, alb=0.0):
   """Return the DISORT_VARS common block as a namespace."""

   if alb != 0.0:
      raise NotImplementedError("A non-zero surface albedo is fixed at 0 in the DISORT adapter")
   disort_vars = SimpleNamespace(
      # DISORT requires most array sizes twice: the "N***" values are the numbers
      # of values passed or requested, the "MAX***" values size DISORT's arrays.
      nlyr=nlayers,                  # Actual number of layers
      nmom=nmoments - 1,             # Actual number of phase moments
      ntau=2,                        # Actual number of user output layers
      numu=nsatzen * 2,              # Actual number of output view zeniths
      nphi=nrelazi,                  # Actual number of output azimuths
      nstr=nstreams,                 # Number of streams to calculate
      maxcly=nlayers,                # Maximum number of layers
      maxulv=2,                      # Maximum number of output layers
      maxumu=nsatzen * 2,            # Maximum number of output view zeniths
      maxphi=nrelazi,                # Maximum number of output azimuths
      maxmom=nmoments - 1,           # Maximum number of phase moments
      usrtau=1,                      # Flag for user defined output levels
      usrang=1,                      # Flag for user defined output angles
      phi0=0.0,                      # Azimuth angle of incident beam
      ibcnd=0,                       # Boundary condition flag
      lamber=1,                      # Flag for Lambertian bottom boundary
      albedo=alb,                    # Albedo of bottom boundary
      onlyfl=0,                      # Flag for disabling intensity output
      accur=1e-8,                    # Convergence criteria
      prnt=[0, 0, 0, 0, 0],          # Printing flags
      header=" " * 127,              # Header string placeholder
   )
   return disort_vars
