      SUBROUTINE  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT,
     &                    UU, TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU,
     &                    MAXTAU, MAXMU, MAXAZ )

c        Print DISORT results and, directly beneath them, their
c        ratios to the correct answers;  print number of non-unit
c        ratios that occur but try to count just the cases where
c        there is a real disagreement and not those where flux or
c        intensity are down at their noise level (defined as 10^(-6)
c        times their maximum value).  d(flux)/d(tau) is treated the
c        same as fluxes in this noise estimation even though it
c        is a different type of quantity (although with flux units).

c     INPUT :   TSTFIR  correct direct flux
c               TSTFDN  correct diffuse down flux
c               TSTFUP  correct diffuse up flux
c               TSTDFD  correct d(flux)/d(optical depth)
c               TSTUU   correct intensity
c               (remaining input = DISORT I/O variables)

c --------------------------------------------------------------------

c     .. Parameters ..

      INTEGER   MAXRAT
      PARAMETER ( MAXRAT = 100 )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   ONLYFL
      INTEGER   MAXAZ, MAXMU, MAXTAU, MAXULV, MAXUMU, NPHI, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      DFDT( * ), FLUP( * ), PHI( * ), RFLDIR( * ), RFLDN( * ),
     &          TSTDFD( * ), TSTFDN( * ), TSTFIR( * ), TSTFUP( * ),
     &          TSTUU( MAXTAU, MAXMU, MAXAZ ), UMU( * ), UTAU( * ),
     &          UU( MAXUMU, MAXULV, * )
c     ..
c     .. Local Scalars ..

      INTEGER  IU, J, LU, NUMBAD
      REAL     FLXMAX, FNOISE, RAT, RAT1, RAT2, RAT3, RAT4, UMAX, UNOISE
c     ..
c     .. Local Arrays ..

      REAL      RATV( MAXRAT )
c     ..
c     .. External Functions ..

      REAL      RATIO
      EXTERNAL  RATIO
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Statement Functions ..

      LOGICAL   BADRAT
c     ..
c     .. Statement Function definitions ..

      BADRAT( RAT ) = (RAT.LT.0.99) .OR. (RAT.GT.1.01)
c     ..


      IF ( NTAU.GT.MAXTAU .OR. NUMU.GT.MAXMU .OR. NPHI.GT.MAXAZ )  CALL
     &   ERRMSG( 'PRTFIN--out of bounds in comparator arrays', .TRUE.)

      FLXMAX = 0.0
      DO 5  LU = 1, NTAU
         FLXMAX = MAX( FLXMAX, TSTFIR(LU), TSTFDN(LU), TSTFUP(LU) )
 5    CONTINUE
      FNOISE = 1.E-6 * FLXMAX
      IF( FLXMAX.LE.0.0 )
     &    CALL ERRMSG( 'PRTFIN--all fluxes zero or negative', .FALSE.)
      IF( FNOISE.LE.0.0 )
     &    CALL ERRMSG( 'PRTFIN--all fluxes near underflowing', .FALSE.)

      NUMBAD = 0

      WRITE(*,'(//,A,/,A,/,A)')
     &  '                  <-------------- FLUXES -------------->',
     &  '    Optical       Downward       Downward         Upward'//
     &  '    d(Net Flux)',
     &  '      Depth         Direct        Diffuse        Diffuse'//
     &  '    / d(Op Dep)'

      DO 10  LU = 1, NTAU

         WRITE( *,'(0P,F11.4,1P,4E15.4)')  UTAU(LU), RFLDIR(LU),
     &          RFLDN(LU), FLUP(LU), DFDT(LU)
         RAT1 = RATIO( RFLDIR(LU), TSTFIR(LU) )
         RAT2 = RATIO( RFLDN(LU),  TSTFDN(LU) )
         RAT3 = RATIO(  FLUP(LU),  TSTFUP(LU) )
         RAT4 = RATIO(  DFDT(LU),  TSTDFD(LU) )
         WRITE( *,'(11X,4( ''    ('',F9.4,'')''))')
     &          RAT1, RAT2, RAT3, RAT4

         IF( BADRAT(RAT1) .AND. ABS(RFLDIR(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT2) .AND. ABS(RFLDN(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT3) .AND. ABS(FLUP(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT4) .AND. ABS(DFDT(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1

10    CONTINUE


      IF ( ONLYFL )  GO TO 100

      IF ( NUMU.GT.MAXRAT .OR. NPHI.GT.MAXRAT )
     &     CALL ERRMSG( 'PRTFIN--increase parameter MAXRAT', .TRUE.)


c                                       ** Print intensities

      IF ( NPHI.GT.8 ) CALL ERRMSG
     &      ( 'PRTFIN--intensity FORMATs inadequate',.FALSE.)

      UMAX = 0.0
      DO 36  LU = 1, NTAU
         DO 35  IU = 1, NUMU
            DO 34  J = 1, NPHI
               UMAX = MAX( UMAX, TSTUU(LU,IU,J) )
 34         CONTINUE
 35      CONTINUE
 36   CONTINUE
      UNOISE = 1.E-6 * UMAX
      IF( UMAX.LE.0.0 )  CALL ERRMSG
     &     ( 'PRTFIN--all intensities zero or negative',.FALSE.)
      IF( UNOISE.LE.0.0 ) CALL ERRMSG
     &     ( 'PRTFIN--all intensities near underflowing',.FALSE.)

      WRITE( *,'(//,A,//,A,/,A,/,A,8(F10.1,4X))' )
     &            ' ********  I N T E N S I T I E S  *********',
     &            '             Polar   Azimuthal Angles (Degrees)',
     &            '   Optical   Angle',
     &            '     Depth  Cosine', ( PHI(J), J = 1, NPHI )

      DO 60  LU = 1, NTAU

         DO 50  IU = 1, NUMU

            IF( IU.EQ.1 ) WRITE( *,'(/,0P,F10.3,F8.3,1P,8E14.4)')
     &          UTAU(LU), UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

            IF( IU.GT.1 ) WRITE( *,'(10X,0P,F8.3, 1P,8E14.4)')
     &                    UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

            DO 40  J = 1, NPHI
               RATV(J) = RATIO( UU(IU,LU,J), TSTUU(LU,IU,J) )
               IF( BADRAT(RATV(J)) .AND. ABS(UU(IU,LU,J)).GT.UNOISE )
     &                  NUMBAD = NUMBAD + 1
   40       CONTINUE

            WRITE( *,'(18X, 8(:,''   ('',F9.4,'')''))')
     &           ( RATV(J), J = 1, NPHI )

   50    CONTINUE
   60 CONTINUE

  100 CONTINUE
      IF( NUMBAD.GT.0 )  WRITE( *,300)  ' ====  ', NUMBAD,
     &    '  SERIOUSLY NON-UNIT RATIOS    ===='

      RETURN

300   FORMAT( //,1X,45('='),/,A,I4,A,/,1X,45('=') )
      END

      BLOCK DATA  CHEKDO

c       Correct answers to test problems ( as produced by DISORT
c       running entirely in double precision (14 significant digits)
c       on a Digital Alpha Workstation computer ).

c     .. Parameters ..

      INTEGER   MXPROB, MXCASE, MAXTAU, MAXMU, MAXAZ
      PARAMETER ( MXPROB = 9, MXCASE = 8, MAXTAU = 5, MAXMU = 10,
     &            MAXAZ = 3 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, J, K
c     ..
c     .. Common blocks ..

      COMMON    / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU

      REAL      TSTDFD( MAXTAU, MXCASE, MXPROB ),
     &          TSTFDN( MAXTAU, MXCASE, MXPROB ),
     &          TSTFIR( MAXTAU, MXCASE, MXPROB ),
     &          TSTFUP( MAXTAU, MXCASE, MXPROB ),
     &          TSTUU( MAXTAU, MAXMU, MAXAZ, MXCASE, MXPROB )
c     ..

c ********************* Test Case 1A *********************************

      DATA (TSTFIR(I,1,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,1,1), I = 1, 2) / 0.0, 7.94108E-02 /
      DATA (TSTFUP(I,1,1), I = 1, 2) / 7.99451E-02, 0.0 /
      DATA (TSTDFD(I,1,1), I = 1, 2) / 2.54067E+01, 1.86531E+01 /
      DATA ((TSTUU(I,J,1,1,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.17771E-01, 2.64170E-02, 1.34041E-02,
     &    1.33826E-02, 2.63324E-02, 1.15898E-01, 3*0.0 /

c ********************* Test Case 1B *********************************

      DATA (TSTFIR(I,2,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,2,1), I = 1, 2) / 0.0, 4.20233E-01 /
      DATA (TSTFUP(I,2,1), I = 1, 2) / 4.22922E-01, 0.0 /
      DATA (TSTDFD(I,2,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 6.22884E-01, 1.39763E-01, 7.09192E-02,
     &    7.08109E-02, 1.39337E-01, 6.13458E-01, 3*0.0 /

c ********************* Test Case 1C *********************************

      DATA (TSTFIR(I,3,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,3,1), I = 1, 2) / 3.14159E+00, 3.04897E+00 /
      DATA (TSTFUP(I,3,1), I = 1, 2) / 9.06556E-02, 0.0 /
      DATA (TSTDFD(I,3,1), I = 1, 2) / 6.66870E-02, 5.88936E-02 /
      DATA ((TSTUU(I,J,1,3,1), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 1.33177E-01, 2.99879E-02, 1.52233E-02,
     &    9.84447E-01, 9.69363E-01, 8.63946E-01, 3*0.0 /

c ********************* Test Case 1D *********************************

      DATA (TSTFIR(I,4,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,4,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,4,1), I = 1, 2) / 2.59686E-01, 0.0 /
      DATA (TSTDFD(I,4,1), I = 1, 2) / 2.57766E+01, 0.0 /
      DATA ((TSTUU(I,J,1,4,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 2.62972E-01, 9.06967E-02, 5.02853E-02,
     &    1.22980E-15, 1.30698E-17, 6.88840E-18, 3*0.0 /

c ********************* Test Case 1E *********************************

      DATA (TSTFIR(I,5,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,5,1), I = 1, 2) / 0.0, 6.76954E-02 /
      DATA (TSTFUP(I,5,1), I = 1, 2) / 3.07390E+00, 0.0 /
      DATA (TSTDFD(I,5,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,5,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.93321E+00, 1.02732E+00, 7.97199E-01,
     &    2.71316E-02, 1.87805E-02, 1.16385E-02, 3*0.0 /

c ********************* Test Case 1F *********************************

      DATA (TSTFIR(I,6,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,6,1), I = 1, 2) / 3.14159E+00, 4.60048E-03 /
      DATA (TSTFUP(I,6,1), I = 1, 2) / 2.49618E+00, 0.0 /
      DATA (TSTDFD(I,6,1), I = 1, 2) / 1.14239E-01, 7.93633E-05 /
      DATA ((TSTUU(I,J,1,6,1), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 8.77510E-01, 8.15136E-01, 7.52715E-01,
     &    1.86840E-03, 1.26492E-03, 7.79280E-04, 3*0.0 /


c ********************* Test Case 2A *********************************

      DATA (TSTFIR(I,1,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,1,2), I = 1, 2) / 0.0, 4.41791E-02 /
      DATA (TSTFUP(I,1,2), I = 1, 2) / 5.35063E-02, 0.0 /
      DATA (TSTDFD(I,1,2), I = 1, 2) / 1.66570E+00, 1.89848E-01 /
      DATA ((TSTUU(I,J,1,1,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.61796E-01, 2.11501E-02, 7.86713E-03,
     &    7.71897E-03, 2.00778E-02, 2.57685E-02, 3*0.0 /

c ********************* Test Case 2B *********************************

      DATA (TSTFIR(I,2,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,2,2), I = 1, 2) / 0.0, 1.06123E-01 /
      DATA (TSTFUP(I,2,2), I = 1, 2) / 1.25561E-01, 0.0 /
      DATA (TSTDFD(I,2,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 3.47678E-01, 4.87120E-02, 1.89387E-02,
     &    1.86027E-02, 4.64061E-02, 6.77603E-02, 3*0.0 /

c ********************* Test Case 2C *********************************

      DATA (TSTFIR(I,3,2), I = 1, 2) / 2.52716E-01, 2.56077E-28 /
      DATA (TSTFDN(I,3,2), I = 1, 2) / 0.0, 2.51683E-04 /
      DATA (TSTFUP(I,3,2), I = 1, 2) / 6.24730E-02, 0.0 /
      DATA (TSTDFD(I,3,2), I = 1, 2) / 1.67462E+00, 1.75464E-04 /
      DATA ((TSTUU(I,J,1,3,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.62566E-01, 2.45786E-02, 1.01498E-02,
     &    1.70004E-04, 3.97168E-05, 1.32472E-05, 3*0.0 /

c ********************* Test Case 2D *********************************

      DATA (TSTFIR(I,4,2), I = 1, 2) / 2.52716E-01, 0.0 /
      DATA (TSTFDN(I,4,2), I = 1, 2) / 0.0, 2.68008E-02 /
      DATA (TSTFUP(I,4,2), I = 1, 2) / 2.25915E-01, 0.0 /
      DATA (TSTDFD(I,4,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,4,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 3.64010E-01, 8.26993E-02, 4.92370E-02,
     &    1.05950E-02, 7.69337E-03, 3.79276E-03, 3*0.0 /


c ********************* Test Case 3A *********************************

      DATA (TSTFIR(I,1,3), I = 1, 2) / 3.14159E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,3), I = 1, 2) / 0.0, 1.73849E+00 /
      DATA (TSTFUP(I,1,3), I = 1, 2) / 2.47374E-01, 0.0 /
      DATA (TSTDFD(I,1,3), I = 1, 2) / 0.0, 0.0 /
      DATA ((TSTUU(I,J,1,1,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 1.51159E-01, 1.01103E-01, 3.95460E-02,
     &    3.05855E+00, 2.66648E-01, 2.13750E-01, 3*0.0 /

c ********************* Test Case 3B *********************************

      DATA (TSTFIR(I,2,3), I = 1, 2) / 3.14159E+00, 1.05389E-03 /
      DATA (TSTFDN(I,2,3), I = 1, 2) / 0.0, 1.54958E+00 /
      DATA (TSTFUP(I,2,3), I = 1, 2) / 1.59096E+00, 0.0 /
      DATA (TSTDFD(I,2,3), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 3.79740E-01, 5.19598E-01, 4.93302E-01,
     &    6.69581E-01, 4.22350E-01, 2.36362E-01, 3*0.0 /


c ********************* Test Case 4A *********************************

      DATA (TSTFIR(I,1,4), I = 1, 3)
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,4), I = 1, 3)
     &  / 0.0, 1.17401E+00, 1.81264E+00 /
      DATA (TSTFUP(I,1,4), I = 1, 3)
     &  / 1.73223E-01, 1.11113E-01, 0.0 /
      DATA (TSTDFD(I,1,4), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,4), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 9.26837E-02,
     &    6.59569E-02, 3.64755E-02, 2.51608E+00, 1.19287E-01,
     &    1.34962E-01, 1.23887E-01, 4.02058E-02, 1.77746E-02,
     &    3.37302E+00, 2.19835E-01, 1.56893E-01, 3*0.0 /

c ********************* Test Case 4B *********************************

      DATA (TSTFIR(I,2,4), I = 1, 3)
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,2,4), I = 1, 3)
     &  / 0.0, 1.01517E+00, 1.51554E+00 /
      DATA (TSTFUP(I,2,4), I = 1, 3)
     &  / 1.23665E-01, 7.88690E-02, 0.0 /
      DATA (TSTDFD(I,2,4), I = 1, 3)
     &  / 3.43724E-01, 3.52390E-01, 3.19450E-01 /
      DATA ((TSTUU(I,J,1,2,4), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 6.53056E-02,
     &    4.55144E-02, 2.82693E-02, 2.24258E+00, 9.66049E-02,
     &    9.61335E-02, 8.43278E-02, 2.79473E-02, 1.38835E-02,
     &    2.97057E+00, 1.67698E-01, 1.08115E-01, 3*0.0 /

c ********************* Test Case 4C *********************************

      DATA (TSTFIR(I,3,4), I = 1, 3)
     &  / 1.57080E+00, 5.77864E-01, 2.12584E-01 /
      DATA (TSTFDN(I,3,4), I = 1, 3)
     &  / 0.0, 7.02764E-01, 8.03294E-01 /
      DATA (TSTFUP(I,3,4), I = 1, 3)
     &  / 2.25487E-01, 1.23848E-01, 0.0 /
      DATA (TSTDFD(I,3,4), I = 1, 3)
     &  / 3.85003E-01, 3.37317E-01, 2.16403E-01 /
      DATA (((TSTUU(I,J,K,3,4), J = 1, 6), I = 1, 3), K = 1, 3)
     &  / 3*0.0, 8.70812E-01,
     &    2.24960E-01, 2.27572E-02, 4.77016E-02, 3.02631E+00,
     &    1.41195E+00, 6.97692E-01, 1.09130E-01, 9.32861E-03,
     &    8.38488E-02, 2.70538E+00, 8.76523E-01, 6*0.0,
     &    8.88117E-02, 5.77411E-02, 2.27572E-02,
     &    4.77016E-02, 5.80971E-02, 1.04502E-01, 9.16071E-02,
     &    2.95842E-02, 9.32861E-03, 8.38488E-02, 9.42187E-02,
     &    8.95457E-02, 6*0.0, 6.98247E-02,
     &    5.02877E-02, 2.27572E-02, 4.77016E-02, 2.58544E-02,
     &    6.25954E-02, 5.91273E-02, 2.47702E-02, 9.32861E-03,
     &    8.38488E-02, 3.99383E-02, 4.67155E-02, 3*0.0 /


c ********************* Test Case 5A *********************************

      DATA (TSTFIR(I,1,5), I = 1, 3)
     &  / 3.14159E+00, 3.97856E-14, 5.03852E-28 /
      DATA (TSTFDN(I,1,5), I = 1, 3)
     &  / 0.0, 2.24768E+00, 4.79851E-01 /
      DATA (TSTFUP(I,1,5), I = 1, 3)
     &  / 2.66174E+00, 1.76783E+00, 0.0 /
      DATA (TSTDFD(I,1,5), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,5), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 4.58927E-01,
     &    7.72983E-01, 1.07196E+00, 7.53662E-01, 6.96362E-01,
     &    6.50541E-01, 6.27631E-01, 5.81809E-01, 5.24532E-01,
     &    1.95230E-01, 1.31990E-01, 7.20655E-02, 3*0.0 /

c ********************* Test Case 5B *********************************

      DATA (TSTFIR(I,2,5), I = 1, 3)
     &  / 1.28058E-01, 8.67322E-06, 4.47729E-21 /
      DATA (TSTFDN(I,2,5), I = 1, 3)
     &  / 1.74767E+00, 2.33975E-01, 6.38345E-05 /
      DATA (TSTFUP(I,2,5), I = 1, 3)
     &  / 2.70485E-01, 3.74252E-02, 1.02904E-05 /
      DATA (TSTDFD(I,2,5), I = 1, 3)
     &  / 3.10129E-01, 4.52671E-02, 1.25021E-05 /
      DATA ((TSTUU(I,J,1,2,5), J = 1, 6), I = 1, 3)
     &  / 6.79623E+01, 2.21027E-01, 1.36619E-01, 1.14084E-01,
     &    8.73870E-02, 8.81626E-02, 2.05706E-01, 4.92736E-02,
     &    2.65449E-02, 2.02154E-02, 1.29661E-02, 9.51334E-03,
     &    3.41286E-05, 1.39916E-05, 7.47039E-06, 5.65602E-06,
     &    3.58245E-06, 2.57858E-06 /

c ********************* Test Case 6A *********************************

      DATA (TSTFIR(I,1,6), I = 1, 2) / 2*100.0 /
      DATA (TSTFDN(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTDFD(I,1,6), I = 1, 2) / 2*200.0 /
      DATA ((TSTUU(I,J,1,1,6), J = 1, 4), I = 1, 2) / 8*0.0 /

c ********************* Test Case 6B *********************************

      DATA (TSTFIR(I,2,6), I = 1, 3)
     &  / 1.000000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTDFD(I,2,6), I = 1, 3)
     &  / 2.00000E+02, 7.35759E+01, 2.70671E+01 /
      DATA ((TSTUU(I,J,1,2,6), J = 1, 4), I = 1, 3) / 12*0.0 /

c ********************* Test Case 6C *********************************

      DATA (TSTFIR(I,3,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,3,6), I = 1, 3)
     &  / 1.48450E+00, 2.99914E+00, 6.76676E+00 /
      DATA (TSTDFD(I,3,6), I = 1, 3)
     &  / 2.02010E+02, 7.79962E+01, 4.06006E+01 /
      DATA ((TSTUU(I,J,1,3,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 9.77882E-05, 7.92386E-01,
     &    2*0.0, 1.45131E-02, 1.30642E+00,
     &    2*0.0, 2.15393E+00, 2.15393E+00 /

c ********************* Test Case 6D *********************************

      DATA (TSTFIR(I,4,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,4,6), I = 1, 3)
     &  / 6.70783E-01, 1.39084E+00, 3.31655E+00 /
      DATA (TSTDFD(I,4,6), I = 1, 3)
     &  / 2.00936E+02, 7.57187E+01, 3.45317E+01 /
      DATA ((TSTUU(I,J,1,4,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 6.80068E-05, 3.15441E-01,
     &    2*0.0, 1.00931E-02, 5.20074E-01,
     &    2*0.0, 1.49795E+00, 8.57458E-01 /

c ********************* Test Case 6E *********************************

      DATA (TSTFIR(I,5,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,5,6), I = 1, 3)
     &  / 7.95458E+01, 1.59902E+02, 3.56410E+02 /
      DATA (TSTDFD(I,5,6), I = 1, 3)
     &  / 3.07079E+02, 3.07108E+02, 7.17467E+02 /
      DATA ((TSTUU(I,J,1,5,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 4.53789E-03, 4.33773E+01,
     &    2*0.0, 6.73483E-01, 7.15170E+01,
     &    2*0.0, 9.99537E+01, 1.17912E+02 /

c ********************* Test Case 6F *********************************

      DATA (TSTFIR(I,6,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,6,6), I = 1, 3)
     &  / 3.21497E+02, 1.42493E+02, 7.05305E+01 /
      DATA (TSTFUP(I,6,6), I = 1, 3)
     &  / 8.27917E+01, 1.66532E+02, 3.71743E+02 /
      DATA (TSTDFD(I,6,6), I = 1, 3)
     &  / 9.54523E+02, 5.27085E+02, 8.45341E+02 /
      DATA ((TSTUU(I,J,1,6,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 4.80531E-03, 4.50168E+01,
     &    6.20697E+01, 6.89532E-01, 7.13172E-01, 7.42191E+01,
     &    3.76472E+01, 4.64603E-03, 1.05844E+02, 1.22368E+02 /

c ********************* Test Case 6G *********************************

      DATA (TSTFIR(I,7,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,7,6), I = 1, 3)
     &  / 3.21497E+02, 3.04775E+02, 3.63632E+02 /
      DATA (TSTFUP(I,7,6), I = 1, 3)
     &  / 3.35292E+02, 4.12540E+02, 4.41125E+02 /
      DATA (TSTDFD(I,7,6), I = 1, 3)
     &  / 5.80394E+02, 1.27117E+02, -1.68003E+02 /
      DATA ((TSTUU(I,J,1,7,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 7.80733E+01, 1.16430E+02,
     &    9.78748E+01, 1.01048E+02, 1.15819E+02, 1.34966E+02,
     &    1.10061E+02, 1.38631E+02, 1.38695E+02, 1.40974E+02 /

c ********************* Test Case 6H *********************************

      DATA (TSTFIR(I,8,6), I = 1, 3)
     &  / 1.00000E+02, 1.35335E+01, 2.06115E-07 /
      DATA ( TSTFDN(I,8,6), I = 1, 3)
     &  / 3.21497E+02, 2.55455E+02, 4.43444E+02 /
      DATA (TSTFUP(I,8,6), I = 1, 3)
     &  / 2.37350E+02, 2.61130E+02, 4.55861E+02 /
      DATA (TSTDFD(I,8,6), I = 1, 3)
     &  / 4.23780E+02, 6.19828E+01, -3.11658E+01 /
      DATA ((TSTUU(I,J,1,8,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 7.12616E+01, 7.80736E+01,
     &    8.49992E+01, 7.73186E+01, 7.88310E+01, 8.56423E+01,
     &    1.38631E+02, 1.45441E+02, 1.44792E+02, 1.45163E+02 /


c ********************* Test Case 7A *********************************

      DATA (TSTFIR(I,1,7), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,1,7), I = 1, 2) / 0.0, 1.21204E+02 /
      DATA (TSTFUP(I,1,7), I = 1, 2) / 8.62936E+01,  0.0 /
      DATA (TSTDFD(I,1,7), I = 1, 2) /-5.13731E+01,-5.41036E+02 /

c ********************* Test Case 7B *********************************

      DATA (TSTFIR(I,2,7), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,2,7), I = 1, 2) / 0.0, 2.07786E-05 /
      DATA (TSTFUP(I,2,7), I = 1, 2) / 1.10949E-06,  0.0 /
      DATA (TSTDFD(I,2,7), I = 1, 2) / 8.23219E-08, -5.06461E-06 /
      DATA ((TSTUU(I,J,1,2,7), J = 1, 2), I = 1, 2) /
     &    0.00000E+00, 4.65744E-07, 7.52311E-06, 0.00000E+00 /

c ********************* Test Case 7C *********************************

      DATA (TSTFIR(I,3,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,7), I = 1, 3)
     &  / 3.19830E+02, 3.54099E+02, 3.01334E+02 /
      DATA (TSTFUP(I,3,7), I = 1, 3)
     &  / 4.29572E+02, 4.47018E+02, 5.94576E+02 /
      DATA (TSTDFD(I,3,7), I = 1, 3)
     &  /-8.04270E+01, 2.51589E+02, 7.15964E+02 /
      DATA (((TSTUU(I,J,K,3,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.46775E+02, 1.49033E+02,
     &    1.06583E+02, 1.28565E+02, 1.04464E+02, 1.59054E+02,
     &    9.66519E+01, 8.65854E+01, 1.89259E+02, 1.89259E+02,
     &    1.01805E+02, 1.01805E+02, 1.29641E+02, 1.49033E+02,
     &    1.06583E+02, 1.06408E+02, 9.48418E+01, 1.59054E+02,
     &    9.66519E+01, 7.49310E+01, 1.89259E+02, 1.89259E+02 /

c ********************* Test Case 7D *********************************

      DATA (TSTFIR(I,4,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,7), I = 1, 3)
     &  / 3.19830E+02, 3.50555E+02, 2.92063E+02 /
      DATA (TSTFUP(I,4,7), I = 1, 3)
     &  / 3.12563E+02, 2.68126E+02, 3.05596E+02 /
      DATA (TSTDFD(I,4,7), I = 1, 3)
     &  /-1.68356E+02, 1.01251E+02, 4.09326E+02 /
      DATA (((TSTUU(I,J,K,4,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.40977E+02, 9.62764E+01,
     &    1.06203E+02, 1.23126E+02, 9.19545E+01, 8.89528E+01,
     &    9.56010E+01, 7.25576E+01, 9.72743E+01, 9.72743E+01,
     &    1.01805E+02, 1.01805E+02, 1.23843E+02, 9.62764E+01,
     &    1.06203E+02, 1.00969E+02, 8.23318E+01, 8.89528E+01,
     &    9.56010E+01, 6.09031E+01, 9.72743E+01, 9.72743E+01 /

c ********************* Test Case 7E *********************************

      DATA (TSTFIR(I,5,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,7), I = 1, 3)
     &  / 3.19830E+02, 3.53275E+02, 2.99002E+02 /
      DATA (TSTFUP(I,5,7), I = 1, 3)
     &  / 4.04300E+02, 4.07843E+02, 5.29248E+02 /
      DATA (TSTDFD(I,5,7), I = 1, 3)
     &  /-9.98568E+01, 2.17387E+02, 6.38461E+02 /
      DATA (((TSTUU(I,J,K,5,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.45448E+02, 1.38554E+02,
     &    1.06496E+02, 1.27296E+02, 1.01395E+02, 1.45229E+02,
     &    9.63993E+01, 8.29009E+01, 1.60734E+02, 1.71307E+02,
     &    1.01805E+02, 1.01805E+02, 1.28281E+02, 1.38554E+02,
     &    1.06496E+02, 1.05111E+02, 9.16726E+01, 1.45229E+02,
     &    9.63993E+01, 7.11248E+01, 1.59286E+02, 1.71307E+02 /

c ********************* Test Case 8A *********************************

      DATA (TSTFIR(I,1,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,1,8), I = 1, 3) /
     &  1.00000E+00, 7.22235E-01, 5.13132E-01 /
      DATA (TSTFUP(I,1,8), I = 1, 3) / 9.29633E-02, 2.78952E-02, 0.0 /
      DATA (TSTDFD(I,1,8), I = 1, 3) /
     &  1.12474E+00, 6.51821E-01, 5.63361E-01 /
      DATA ((TSTUU(I,J,1,1,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 5.62566E-02, 1.94423E-02,
     &  2.62711E-01, 1.36952E-01, 1.84909E-02, 5.52188E-03,
     &  2.10014E-01, 5.60376E-02, 2*0.0 /

c ********************* Test Case 8B *********************************

      DATA (TSTFIR(I,2,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,2,8), I = 1, 3) /
     &  1.00000E+00, 7.95332E-01, 6.50417E-01 /
      DATA (TSTFUP(I,2,8), I = 1, 3) / 2.25136E-01, 1.26349E-01, 0.0 /
      DATA (TSTDFD(I,2,8), I = 1, 3) /
     &  5.12692E-01, 3.56655E-01, 5.68095E-02 /
      DATA ((TSTUU(I,J,1,2,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 1.23687E-01, 4.95581E-02,
     &  2.77499E-01, 1.83950E-01, 8.35695E-02, 2.50575E-02,
     &  2.40731E-01, 1.29291E-01, 2*0.0 /

c ********************* Test Case 8C *********************************

      DATA (TSTFIR(I,3,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,3,8), I = 1, 3) /
     &  1.00000E+00, 4.86157E-01, 1.59984E-01 /
      DATA (TSTFUP(I,3,8), I = 1, 3) / 3.78578E-01, 2.43397E-01, 0.0 /
      DATA (TSTDFD(I,3,8), I = 1, 3) /
     &  5.65095E-01, 2.76697E-01, 1.35679E-02 /
      DATA ((TSTUU(I,J,1,3,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 1.49335E-01, 1.04766E-01,
     &  1.89020E-01, 9.88158E-02, 9.65192E-02, 6.54445E-02,
     &  6.84762E-02, 2.96698E-02, 2*0.0 /


c ********************* Test Case 9A *********************************

      DATA (TSTFIR(I,1,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,1,9), I = 1, 5) /
     &  1.00000E+00, 3.55151E-01, 1.44265E-01, 6.71445E-03, 6.16968E-07/
      DATA (TSTFUP(I,1,9), I = 1, 5) /
     &  2.27973E-01, 8.75098E-02, 3.61819E-02, 2.19291E-03, 0.0 /
      DATA (TSTDFD(I,1,9), I = 1, 5) /
     &  8.82116E-01, 2.32366E-01, 9.33443E-02, 3.92782E-03, 1.02500E-07/
      DATA ((TSTUU(I,J,1,1,9), J = 1, 4), I = 1, 5) /
     &  2*3.18310E-01, 9.98915E-02, 5.91345E-02,
     &  1.53507E-01, 5.09531E-02, 3.67006E-02, 2.31903E-02,
     &  7.06614E-02, 2.09119E-02, 1.48545E-02, 9.72307E-03,
     &  3.72784E-03, 1.08815E-03, 8.83316E-04, 5.94743E-04,
     &  2.87656E-07, 1.05921E-07, 2*0.0 /

c ********************* Test Case 9B *********************************

      DATA (TSTFIR(I,2,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,2,9), I = 1, 5) /
     &  1.00000E+00, 4.52357E-01, 2.36473E-01, 2.76475E-02, 7.41853E-05/
      DATA (TSTFUP(I,2,9), I = 1, 5) /
     &  1.00079E-01, 4.52014E-02, 2.41941E-02, 4.16016E-03, 0.0 /
      DATA (TSTDFD(I,2,9), I = 1, 5) /
     &  8.04577E-01, 2.55330E-01, 1.30976E-01, 1.36227E-02, 1.22022E-05/
      DATA ((TSTUU(I,J,1,2,9), J = 1, 4), I = 1, 5) /
     &  2*3.18310E-01, 7.39198E-02, 1.32768E-02,
     &  1.96609E-01, 5.92369E-02, 3.00230E-02, 7.05566E-03,
     &  1.15478E-01, 3.01809E-02, 1.52672E-02, 4.06932E-03,
     &  1.46177E-02, 3.85590E-03, 2.38301E-03, 7.77890E-04,
     &  3.37742E-05, 1.20858E-05, 2*0.0 /

c ********************* Test Case 9C *********************************

      DATA (TSTFIR(I,3,9), I = 1, 5) /
     &  1.57080E+00, 1.92354E-01, 2.35550E-02, 9.65131E-06, 9.03133E-19/
      DATA (TSTFDN(I,3,9), I = 1, 5 ) /
     &  6.09217E+00, 4.97279E+00, 4.46616E+00, 4.22731E+00, 4.73767E+00/
      DATA (TSTFUP(I,3,9), I = 1, 5) /
     &  4.68414E+00, 4.24381E+00, 4.16941E+00, 4.30667E+00, 5.11524E+00/
      DATA (TSTDFD(I,3,9), I = 1, 5) /
     &  3.49563E+00, 8.81206E-01, 3.50053E-01, 1.93471E-02, 7.15349E-02/
      DATA (((TSTUU(I,J,K,3,9), J = 1, 4), I = 1, 5), K = 1, 3)
     & / 1.93920E+00, 1.93920E+00, 1.61855E+00, 1.43872E+00,
     &   1.66764E+00, 1.44453E+00, 1.38339E+00, 1.33890E+00,
     &   1.48511E+00, 1.35009E+00, 1.33079E+00, 1.32794E+00,
     &   1.34514E+00, 1.35131E+00, 1.35980E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00,
     &   1.93920E+00, 1.93920E+00, 1.57895E+00, 1.43872E+00,
     &   1.66764E+00, 1.42925E+00, 1.37317E+00, 1.33890E+00,
     &   1.48511E+00, 1.34587E+00, 1.32921E+00, 1.32794E+00,
     &   1.34514E+00, 1.35129E+00, 1.35979E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00,
     &   1.93920E+00, 1.93920E+00, 1.56559E+00, 1.43872E+00,
     &   1.66764E+00, 1.42444E+00, 1.37034E+00, 1.33890E+00,
     &   1.48511E+00, 1.34469E+00, 1.32873E+00, 1.32794E+00,
     &   1.34514E+00, 1.35128E+00, 1.35979E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 2*1.62823E+00 /

      END
