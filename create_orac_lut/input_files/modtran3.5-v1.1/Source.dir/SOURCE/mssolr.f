        SUBROUTINE  MSSOLR( CMU, CWT, FBEAM, GC, GU, KK, LAYRU, LL,     MSS 0001
     $                      LYRCUT, MAXUMU, MXCMU, MXUMU, NCUT, NN,     MSS 0002
     $                      NSTR, NTAU, NUMU, PI, TAUCPR, UMU0, UTAU,   MSS 0003
     $                      UTAUPR, BEAMMS, ZZ, FDNSRT, S0CMS )         MSS 0004
                                                                        MSS 0005
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            MSS 0006
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                MSS 0007
C    I N P U T     V A R I A B L E S:                                   MSS 0008
                                                                        MSS 0009
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE    MSS 0010
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE      MSS 0011
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)      MSS 0012
C       GU       :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES      MSS 0013
C                   (I.E., -G- IN EQ. SC(1))                            MSS 0014
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)           MSS 0015
C       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-                   MSS 0016
C       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED     MSS 0017
C                   BY SOLVING SCALED VERSION OF EQ. SC(5);             MSS 0018
C                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED         MSS 0019
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER        MSS 0020
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)           MSS 0021
C       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION      MSS 0022
C                     OPTICAL DEPTH EXCEEDS -ABSCUT-                    MSS 0023
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)           MSS 0024
C       UTAU     :  OPTICAL DEPTHS OF USER OUTPUT LEVELS                MSS 0025
C       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M     MSS 0026
C                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M      MSS 0027
C       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)                   MSS 0028
                                                                        MSS 0029
C   I N T E R N A L       V A R I A B L E S:                            MSS 0030
                                                                        MSS 0031
C       FLDIR    :  DIRECT-BEAM FLUX (DELTA-M SCALED)                   MSS 0032
C       FACT     :  EXP( - UTAUPR / UMU0 )                              MSS 0033
C       ZINT     :  INTENSITY OF M = 0 CASE, IN EQ. SC(1)               MSS 0034
C      RFLDIR    :  DIRECT-BEAM FLUX (NOT DELTA-M SCALED)               MSS 0035
                                                                        MSS 0036
C   O U T P U T    V A R I A B L E S:                                   MSS 0037
                                                                        MSS 0038
C       S0CMS    :  MULTIPLE SCATTERING SOLAR SOURCE FUNCTION           MSS 0039
C      FDNSRT    :  DOWNWARD DIFFUSE SOLAR FLUX AT SURFACE              MSS 0040
                                                                        MSS 0041
      LOGICAL LYRCUT                                                    MSS 0042
      REAL*8    S0CMS( MAXUMU,* )                                       MSS 0043
      INTEGER LAYRU( * )                                                MSS 0044
      REAL*8    CMU( * ), CWT( * ), GC( MXCMU,MXCMU,* ),                MSS 0045
     $        GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),        MSS 0046
     $        TAUCPR( 0:* ), UTAU( * ), UTAUPR( * ),                    MSS 0047
     $        BEAMMS( MXUMU,* ), ZZ( MXCMU,* )                          MSS 0048
                                                                        MSS 0049
C                                               ** LOOP OVER USER LEVELSMSS 0050
      DO 30 LU = 1, NTAU                                                MSS 0051
                                                                        MSS 0052
         LYU = LAYRU(LU)                                                MSS 0053
         IF ( LYRCUT .AND. LYU.GT.NCUT ) GOTO 30                        MSS 0054
                                                                        MSS 0055
C                                     ** NO RADIATION REACHES THIS LEVELMSS 0056
                                                                        MSS 0057
         FACT  = DEXP( - UTAUPR(LU) / UMU0 )                            MSS 0058
C                                               ** LOOP OVER USER ANGLESMSS 0059
         DO 20 IU = 1, NUMU                                             MSS 0060
            ZINT = 0.0                                                  MSS 0061
            DO 10 JQ = 1, NN                                            MSS 0062
               ZINT = ZINT + GU(IU,JQ,LYU) * LL(JQ,LYU) *               MSS 0063
     $                DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) ) MSS 0064
10          CONTINUE                                                    MSS 0065
            DO 11 JQ = NN+1, NSTR                                       MSS 0066
               ZINT = ZINT + GU(IU,JQ,LYU) * LL(JQ,LYU) *               MSS 0067
     $         DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )      MSS 0068
11          CONTINUE                                                    MSS 0069
                                                                        MSS 0070
C                   **  MS SOURCE FUNCTIONS CALCULATED STW(30) M.S. TERMMSS 0071
                                                                        MSS 0072
            S0CMS(IU,LU) = ZINT + BEAMMS(IU,LYU) * FACT                 MSS 0073
20       CONTINUE                                                       MSS 0074
30    CONTINUE                                                          MSS 0075
                                                                        MSS 0076
C                                  **  LAYER AVERAGE OF S0CMS AS MODTRANMSS 0077
      DO 40 IU = 1, NUMU                                                MSS 0078
         DO 40 LU = 1, NTAU-1                                           MSS 0079
         S0CMS(IU,LU) = ( S0CMS(IU,LU) + S0CMS(IU,LU+1) ) / 2.          MSS 0080
40    CONTINUE                                                          MSS 0081
                                                                        MSS 0082
      LU = NTAU                                                         MSS 0083
      LYU = LAYRU(LU)                                                   MSS 0084
      FDNSRT = 0.                                                       MSS 0085
                                                                        MSS 0086
      IF ( .NOT.( LYRCUT .AND. LYU.GT.NCUT ) )  THEN                    MSS 0087
C                                           ** RADIATION REACHES SURFACEMSS 0088
         FACT  = DEXP( - UTAUPR(LU) / UMU0 )                            MSS 0089
         FLDIR = UMU0 * ( FBEAM * FACT )                                MSS 0090
         RFLDIR= UMU0 * FBEAM * EXP( - UTAU( LU ) / UMU0 )              MSS 0091
                                                                        MSS 0092
         DO 60  IQ = 1, NN                                              MSS 0093
            ZINT = 0.0                                                  MSS 0094
            DO 50  JQ = 1, NN                                           MSS 0095
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               MSS 0096
     $                DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) ) MSS 0097
50          CONTINUE                                                    MSS 0098
            DO 51  JQ = NN+1, NSTR                                      MSS 0099
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               MSS 0100
     $          DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )     MSS 0101
51          CONTINUE                                                    MSS 0102
            FDNSRT = FDNSRT + CWT(NN+1-IQ) * CMU(NN+1-IQ) *             MSS 0103
     $                        ( ZINT + ZZ(IQ,LYU)*FACT )                MSS 0104
60       CONTINUE                                                       MSS 0105
                                                                        MSS 0106
         FDNSRT = 2.0 * PI * FDNSRT + FLDIR - RFLDIR                    MSS 0107
         IF ( FDNSRT.LT.0. )  FDNSRT = 0.                               MSS 0108
      END IF                                                            MSS 0109
                                                                        MSS 0110
      RETURN                                                            MSS 0111
      END                                                               MSS 0112
