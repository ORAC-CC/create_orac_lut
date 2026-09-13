        SUBROUTINE  MSTHML( CMU, CWT, GC, GU, KK, LAYRU, LL, LYRCUT,    MST 0001
     $                      MAXUMU, MXCMU, MXUMU, NCUT, NN, NSTR, NTAU, MST 0002
     $                      NUMU, OPRIM, PI, TAUCPR, UMU0, UTAUPR,      MST 0003
     $                      XR0, XR1, Z0UMS, Z1UMS, ZPLK0, ZPLK1,       MST 0004
     $                      FDNTRT, T0CMS )                             MST 0005
                                                                        MST 0006
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            MST 0007
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                MST 0008
C    I N P U T     V A R I A B L E S:                                   MST 0009
                                                                        MST 0010
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE    MST 0011
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE      MST 0012
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)      MST 0013
C       GU       :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES      MST 0014
C                   (I.E., -G- IN EQ. SC(1))                            MST 0015
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)           MST 0016
C       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-                   MST 0017
C       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED     MST 0018
C                   BY SOLVING SCALED VERSION OF EQ. SC(5);             MST 0019
C                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED         MST 0020
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER        MST 0021
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)           MST 0022
C       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION      MST 0023
C                     OPTICAL DEPTH EXCEEDS -ABSCUT-                    MST 0024
C       OPRIM    :  SINGLE SCATTERING ALBEDO                            MST 0025
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)           MST 0026
C       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M     MST 0027
C                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M      MST 0028
C       XR0      :  EXPANSION OF THERMAL SOURCE FUNCTION                MST 0029
C       XR1      :  EXPANSION OF THERMAL SOURCE FUNCTION EQS.SS(14-16)  MST 0030
C       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)  MST 0031
C       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)  MST 0032
                                                                        MST 0033
C   I N T E R N A L       V A R I A B L E S:                            MST 0034
                                                                        MST 0035
C       ZINT     :  INTENSITY OF M = 0 CASE, IN EQ. SC(1)               MST 0036
                                                                        MST 0037
C   O U T P U T    V A R I A B L E S:                                   MST 0038
                                                                        MST 0039
C       T0CMS    :  MULTIPLE SCATTERING THERMAL SOURCE FUNCTION         MST 0040
C      FDNTRT    :  DOWNWARD DIFFUSE THERMAL FLUX AT SURFACE            MST 0041
                                                                        MST 0042
      LOGICAL LYRCUT                                                    MST 0043
      REAL*8    T0CMS( MAXUMU,* )                                       MST 0044
      INTEGER LAYRU( * )                                                MST 0045
      REAL*8    T0C,CMU( * ), CWT( * ), GC( MXCMU,MXCMU,* ),            MST 0046
     $        GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),        MST 0047
     $        OPRIM( * ), TAUCPR( 0:* ), UTAUPR( * ), XR0( * ),         MST 0048
     $        XR1( * ), Z0UMS( MXUMU,* ), Z1UMS( MXUMU,* ),             MST 0049
     $        FDNTRT,ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* )                 MST 0050
C                                               ** LOOP OVER USER LEVELSMST 0051
      DO 30 LU = 1, NTAU                                                MST 0052
         LYU = LAYRU(LU)                                                MST 0053
         IF ( LYRCUT .AND. LYU.GT.NCUT ) GOTO 30                        MST 0054
                                                                        MST 0055
C                                     ** NO RADIATION REACHES THIS LEVELMST 0056
                                                                        MST 0057
         T0C = (1. - OPRIM(LYU)) * ( XR0(LYU) + XR1(LYU)*UTAUPR(LU) )   MST 0058
C                                ** LOOP OVER USER ANGLES               MST 0059
         DO 20 IU = 1, NUMU                                             MST 0060
            ZINT = 0.0                                                  MST 0061
            DO 10 JQ = 1, NN                                            MST 0062
               ZINT = ZINT + GU(IU,JQ,LYU) * LL(JQ,LYU) *               MST 0063
     $                DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) ) MST 0064
10          CONTINUE                                                    MST 0065
            DO 11 JQ = NN+1, NSTR                                       MST 0066
               ZINT = ZINT + GU(IU,JQ,LYU) * LL(JQ,LYU) *               MST 0067
     $         DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )      MST 0068
11          CONTINUE                                                    MST 0069
                                                                        MST 0070
C                   **  MS SOURCE FUNCTIONS CALCULATED STW(30) M.S. TERMMST 0071
                                                                        MST 0072
            T0CMS(IU,LU) = T0C + ZINT + Z0UMS(IU,LYU) +                 MST 0073
     $                                  Z1UMS(IU,LYU) * UTAUPR(LU)      MST 0074
20       CONTINUE                                                       MST 0075
30    CONTINUE                                                          MST 0076
C                                  **  LAYER AVERAGE OF T0CMS AS MODTRANMST 0077
      DO 40 IU = 1, NUMU                                                MST 0078
         DO 40 LU = 1, NTAU-1                                           MST 0079
                                                                        MST 0080
                                                                        MST 0081
         T0CMS(IU,LU) = ( T0CMS(IU,LU) + T0CMS(IU,LU+1) ) / 2.          MST 0082
40    CONTINUE                                                          MST 0083
      LU = NTAU                                                         MST 0084
                                                                        MST 0085
      LYU = LAYRU(LU)                                                   MST 0086
                                                                        MST 0087
C      FDNTRT = 0.0                                                     MST 0088
      FDNTRT = 0.D0                                                     MST 0089
                                                                        MST 0090
      DO 60  IQ = 1, NN                                                 MST 0091
                                                                        MST 0092
                                                                        MST 0093
         ZINT = 0.0                                                     MST 0094
         DO 50  JQ = 1, NN                                              MST 0095
            ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *                  MST 0096
     $             DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )    MST 0097
50       CONTINUE                                                       MST 0098
         DO 51  JQ = NN+1, NSTR                                         MST 0099
            ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *                  MST 0100
     $             DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )  MST 0101
51       CONTINUE                                                       MST 0102
         FDNTRT = FDNTRT + CWT(NN+1-IQ) * CMU(NN+1-IQ) * ( ZINT +       MST 0103
     $                      ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU) )    MST 0104
60       CONTINUE                                                       MST 0105
                                                                        MST 0106
         FDNTRT = 2.0 * PI * FDNTRT                                     MST 0107
         IF ( FDNTRT.LT.0. )  FDNTRT = 0.                               MST 0108
      RETURN                                                            MST 0109
      END                                                               MST 0110
