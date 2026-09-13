      REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )

c      Supplies surface bi-directional reflectivity.
c
c      This routine replaces the "stub" version supplied with the
c      standard distribution. It calls the MODIS AMBRALS surface BRDF
c      model via the C function ambrals-fortran.
c
c      Written by Gareth Thomas, 02 Feb 2007.
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c
c   Called by- DREF, SURFAC

c +-------------------------------------------------------------------+
c
c     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
c     .. Internal variables ..
      REAL *4   PI
      Parameter (PI = 3.14159)
      REAL *4   sz, vz, phi, wvlo, wvhi, bdr
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, AmbralsFortran

c     Make sure the parameters being passed into the C code have the right
c     size.
      sz  = 180.0 * acos(MUP) / PI
      vz  = 180.0 * acos(MU)  / PI
      phi = DPHI
      wvlo = WVNMLO
      wvhi = WVNMHI

      CALL AmbralsFortran(wvlo,wvhi,sz,vz,phi,bdr)

      BDREF = bdr

      RETURN
      END
