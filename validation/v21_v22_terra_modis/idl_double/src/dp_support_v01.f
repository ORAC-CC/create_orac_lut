c     Support routines for the diagnostic double-precision DISORT build.
c
c     DPREALBYTES reports the storage size of default REAL, so the DLM can
c     refuse to load unless the kernel really was compiled with 8-byte REAL.
c
c     BDREF replaces BDREF.f in this build.  The original calls the C routine
c     ambralsfortran_ with float arguments, which would be mis-passed once
c     REAL is promoted; ORAC only uses LAMBER = .TRUE., so any call is fatal.

      INTEGER FUNCTION DPREALBYTES()
      REAL X
      DPREALBYTES = STORAGE_SIZE(X) / 8
      END

      REAL FUNCTION BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )
      REAL WVNMLO, WVNMHI, MU, MUP, DPHI
      BDREF = 0.0
      CALL ERRMSG( 'BDREF--BRDF not available in double-precision '//
     &             'diagnostic build', .TRUE. )
      END
