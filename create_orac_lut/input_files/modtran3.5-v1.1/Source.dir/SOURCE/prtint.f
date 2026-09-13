      SUBROUTINE  PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI,         INT 0001
     $                    MAXULV, MAXUMU )                              INT 0002
                                                                        INT 0003
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            INT 0004
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                INT 0005
C         PRINTS THE INTENSITY AT USER POLAR AND AZIMUTHAL ANGLES       INT 0006
                                                                        INT 0007
C     ALL ARGUMENTS ARE DISORT INPUT OR OUTPUT VARIABLES                INT 0008
                                                                        INT 0009
C+---------------------------------------------------------------------+INT 0010
      REAL*8   PHI(*), UMU(*), UTAU(*), UU( MAXUMU, MAXULV, * )         INT 0011
                                                                        INT 0012
                                                                        INT 0013
      WRITE ( *, '(//,A)' )                                             INT 0014
     $         ' *********  I N T E N S I T I E S  *********'           INT 0015
      LENFMT = 10                                                       INT 0016
      NPASS = 1 + NPHI / LENFMT                                         INT 0017
      IF ( MOD(NPHI,LENFMT) .EQ. 0 )  NPASS = NPASS - 1                 INT 0018
      DO 10  LU = 1, NTAU                                               INT 0019
         DO 10  NP = 1, NPASS                                           INT 0020
            JMIN = 1 + LENFMT * (NP-1)                                  INT 0021
            JMAX = MIN0( LENFMT*NP, NPHI )                              INT 0022
            WRITE( *,101 )  ( PHI(J), J = JMIN, JMAX )                  INT 0023
            DO 10  IU = 1, NUMU                                         INT 0024
               IF( IU.EQ.1 )  WRITE( *,102 )  UTAU(LU), UMU(IU),        INT 0025
     $           ( UU( IU,LU,J ), J = JMIN, JMAX )                      INT 0026
               IF( IU.GT.1 )  WRITE( *,103 )  UMU(IU),                  INT 0027
     $           ( UU( IU,LU,J ), J = JMIN, JMAX )                      INT 0028
10    CONTINUE                                                          INT 0029
                                                                        INT 0030
      RETURN                                                            INT 0031
                                                                        INT 0032
101   FORMAT( /, 3X,'          POLAR   AZIMUTH ANGLES (DEGREES)',       INT 0033
     $        /, 3X,'OPTICAL   ANGLE',                                  INT 0034
     $        /, 3X,' DEPTH   COSINE', 10F11.2 )                        INT 0035
102   FORMAT( F10.4, F8.4, 1P,10E11.3 )                                 INT 0036
103   FORMAT( 10X,   F8.4, 1P,10E11.3 )                                 INT 0037
      END                                                               INT 0038
