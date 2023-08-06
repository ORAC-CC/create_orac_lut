      SUBROUTINE  PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU, U0U )          PAI 0001
                                                                        PAI 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            PAI 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                PAI 0004
C        PRINT AZIMUTHALLY AVERAGED INTENSITIES AT USER ANGLES          PAI 0005
                                                                        PAI 0006
      REAL*8     UMU(*), UTAU(*), U0U( MAXUMU,* )                       PAI 0007
                                                                        PAI 0008
                                                                        PAI 0009
      WRITE ( *, '(//,A)' )                                             PAI 0010
     $         ' *********  AZIMUTHALLY AVERAGED INTENSITIES '          PAI 0011
     $       // '(USER POLAR ANGLES)  *********'                        PAI 0012
      LENFMT = 8                                                        PAI 0013
      NPASS = 1 + NUMU / LENFMT                                         PAI 0014
      IF ( MOD(NUMU,LENFMT) .EQ. 0 )  NPASS = NPASS - 1                 PAI 0015
      DO 10  NP = 1, NPASS                                              PAI 0016
         IUMIN = 1 + LENFMT * (NP-1)                                    PAI 0017
         IUMAX = MIN0( LENFMT*NP, NUMU )                                PAI 0018
         WRITE ( *,101 )  ( UMU(IU), IU = IUMIN, IUMAX )                PAI 0019
         DO 10  LU = 1, NTAU                                            PAI 0020
            WRITE( *,102 ) UTAU(LU), ( U0U(IU,LU), IU=IUMIN,IUMAX)      PAI 0021
 10   CONTINUE                                                          PAI 0022
                                                                        PAI 0023
      RETURN                                                            PAI 0024
                                                                        PAI 0025
101   FORMAT( /, 3X,'OPTICAL   POLAR ANGLE COSINES',                    PAI 0026
     $        /, 3X,'  DEPTH', 8F14.5 )                                 PAI 0027
102   FORMAT( 0P,F10.4, 1P,8E14.4 )                                     PAI 0028
      END                                                               PAI 0029
