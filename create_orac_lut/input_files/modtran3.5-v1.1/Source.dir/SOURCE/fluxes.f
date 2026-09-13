      SUBROUTINE  FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,   FXS 0001
     $                    MXCMU, MXULV, NCUT, NN, NSTR, NTAU, PI,       FXS 0002
     $                    PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,      FXS 0003
     $                    XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,       FXS 0004
     $                    FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C,        FXS 0005
     $                    MAXULV )                                      FXS 0006
                                                                        FXS 0007
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            FXS 0008
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                FXS 0009
C       CALCULATES THE RADIATIVE FLUXES, MEAN INTENSITY, AND FLUX       FXS 0010
C       DERIVATIVE WITH RESPECT TO OPTICAL DEPTH FROM THE M=0 INTENSITY FXS 0011
C       COMPONENTS (THE AZIMUTHALLY-AVERAGED INTENSITY)                 FXS 0012
                                                                        FXS 0013
C    I N P U T     V A R I A B L E S:                                   FXS 0014
                                                                        FXS 0015
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE    FXS 0016
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE      FXS 0017
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)      FXS 0018
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)           FXS 0019
C       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-                   FXS 0020
C       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED     FXS 0021
C                   BY SOLVING SCALED VERSION OF EQ. SC(5);             FXS 0022
C                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED         FXS 0023
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER        FXS 0024
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)           FXS 0025
C       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION      FXS 0026
C                     OPTICAL DEPTH EXCEEDS -ABSCUT-                    FXS 0027
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)           FXS 0028
C       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M     FXS 0029
C                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M      FXS 0030
C       XR0      :  EXPANSION OF THERMAL SOURCE FUNCTION IN EQ. SS(14)  FXS 0031
C       XR1      :  EXPANSION OF THERMAL SOURCE FUNCTION EQS. SS(16)    FXS 0032
C       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)                   FXS 0033
C       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)  FXS 0034
C       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)  FXS 0035
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        FXS 0036
                                                                        FXS 0037
C   O U T P U T     V A R I A B L E S:                                  FXS 0038
                                                                        FXS 0039
C       U0C      :  AZIMUTHALLY AVERAGED INTENSITIES                    FXS 0040
C                   ( AT POLAR QUADRATURE ANGLES )                      FXS 0041
C       (RFLDIR, RFLDN, FLUP, DFDT, UAVG ARE 'DISORT' OUTPUT VARIABLES) FXS 0042
                                                                        FXS 0043
C   I N T E R N A L       V A R I A B L E S:                            FXS 0044
                                                                        FXS 0045
C       DIRINT   :  DIRECT INTENSITY ATTENUATED                         FXS 0046
C       FDNTOT   :  TOTAL DOWNWARD FLUX (DIRECT + DIFFUSE)              FXS 0047
C       FLDIR    :  DIRECT-BEAM FLUX (DELTA-M SCALED)                   FXS 0048
C       FLDN     :  DIFFUSE DOWN-FLUX (DELTA-M SCALED)                  FXS 0049
C       FNET     :  NET FLUX (TOTAL-DOWN - DIFFUSE-UP)                  FXS 0050
C       FACT     :  EXP( - UTAUPR / UMU0 )                              FXS 0051
C       PLSORC   :  PLANCK SOURCE FUNCTION (THERMAL)                    FXS 0052
C       ZINT     :  INTENSITY OF M = 0 CASE, IN EQ. SC(1)               FXS 0053
C+---------------------------------------------------------------------+FXS 0054
                                                                        FXS 0055
      LOGICAL LYRCUT, PRNT(*)                                           FXS 0056
      REAL*8 DFDT(*),FLUP(*),FLDIR(*),FLDN(*),RFLDIR(*),RFLDN(* ),      FXS 0057
     $        U0C( MXCMU,MXULV ), UAVG(*)                               FXS 0058
      INTEGER LAYRU(*)                                                  FXS 0059
      REAL*8    CMU(*), CWT(*), GC( MXCMU,MXCMU,* ), KK( MXCMU,* ),     FXS 0060
     $        LL( MXCMU,* ), SSALB(*), TAUCPR( 0:* ),                   FXS 0061
     $        UTAU(*), UTAUPR(*), XR0(*), XR1(*), ZZ( MXCMU,* ),        FXS 0062
     $        ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* )                        FXS 0063
                                                                        FXS 0064
                                                                        FXS 0065
      IF ( PRNT(2) )  WRITE( *,1010 )                                   FXS 0066
C                                          ** ZERO DISORT OUTPUT ARRAYS FXS 0067
      CALL  ZEROIT( U0C, MXULV*MXCMU )                                  FXS 0068
      CALL  ZEROIT( RFLDIR, MAXULV )                                    FXS 0069
      CALL  ZEROIT( FLDIR,  MXULV )                                     FXS 0070
      CALL  ZEROIT( RFLDN,  MAXULV )                                    FXS 0071
      CALL  ZEROIT( FLDN,   MXULV )                                     FXS 0072
      CALL  ZEROIT( FLUP,   MAXULV )                                    FXS 0073
      CALL  ZEROIT( UAVG,   MAXULV )                                    FXS 0074
      CALL  ZEROIT( DFDT,   MAXULV )                                    FXS 0075
C                                        ** LOOP OVER USER LEVELS       FXS 0076
      DO 100  LU = 1, NTAU                                              FXS 0077
                                                                        FXS 0078
         LYU = LAYRU(LU)                                                FXS 0079
                                                                        FXS 0080
         IF ( LYRCUT .AND. LYU.GT.NCUT ) THEN                           FXS 0081
C                                                ** NO RADIATION REACHESFXS 0082
C                                                ** THIS LEVEL          FXS 0083
            FDNTOT = 0.0                                                FXS 0084
            FNET   = 0.0                                                FXS 0085
            PLSORC = 0.0                                                FXS 0086
            GO TO 90                                                    FXS 0087
         END IF                                                         FXS 0088
                                                                        FXS 0089
         IF ( FBEAM.GT.0.0 )  THEN                                      FXS 0090
            FACT  = DEXP( - UTAUPR(LU) / UMU0 )                         FXS 0091
            DIRINT = FBEAM * FACT                                       FXS 0092
            FLDIR(  LU ) = UMU0 * ( FBEAM * FACT )                      FXS 0093
            RFLDIR( LU ) = UMU0 * FBEAM * DEXP( - UTAU( LU ) / UMU0 )   FXS 0094
         ELSE                                                           FXS 0095
            DIRINT = 0.0                                                FXS 0096
            FLDIR(  LU ) = 0.0                                          FXS 0097
            RFLDIR( LU ) = 0.0                                          FXS 0098
         END IF                                                         FXS 0099
                                                                        FXS 0100
         DO 20  IQ = 1, NN                                              FXS 0101
                                                                        FXS 0102
            ZINT = 0.0                                                  FXS 0103
            DO 10  JQ = 1, NN                                           FXS 0104
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               FXS 0105
     $                DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) ) FXS 0106
10          CONTINUE                                                    FXS 0107
            DO 11  JQ = NN+1, NSTR                                      FXS 0108
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               FXS 0109
     $           DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )    FXS 0110
11          CONTINUE                                                    FXS 0111
                                                                        FXS 0112
            U0C( IQ,LU ) = ZINT                                         FXS 0113
            IF ( FBEAM.GT.0.0 )  U0C( IQ,LU ) = ZINT + ZZ(IQ,LYU) * FACTFXS 0114
            U0C( IQ,LU ) = U0C( IQ,LU ) + ZPLK0(IQ,LYU)                 FXS 0115
     $                     + ZPLK1(IQ,LYU) * UTAUPR(LU)                 FXS 0116
            UAVG(LU) = UAVG(LU) + CWT(NN+1-IQ) * U0C( IQ,LU )           FXS 0117
            FLDN(LU) = FLDN(LU) + CWT(NN+1-IQ)*CMU(NN+1-IQ) * U0C(IQ,LU)FXS 0118
20       CONTINUE                                                       FXS 0119
                                                                        FXS 0120
         DO 40  IQ = NN+1, NSTR                                         FXS 0121
                                                                        FXS 0122
            ZINT = 0.0                                                  FXS 0123
            DO 30  JQ = 1, NN                                           FXS 0124
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               FXS 0125
     $                DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) ) FXS 0126
30          CONTINUE                                                    FXS 0127
            DO 31  JQ = NN+1, NSTR                                      FXS 0128
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               FXS 0129
     $           DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )    FXS 0130
31          CONTINUE                                                    FXS 0131
                                                                        FXS 0132
            U0C( IQ,LU ) = ZINT                                         FXS 0133
            IF ( FBEAM.GT.0.0 )  U0C( IQ,LU ) = ZINT + ZZ(IQ,LYU) * FACTFXS 0134
            U0C( IQ,LU ) = U0C( IQ,LU ) + ZPLK0(IQ,LYU)                 FXS 0135
     $                     + ZPLK1(IQ,LYU) * UTAUPR(LU)                 FXS 0136
            UAVG(LU) = UAVG(LU) + CWT(IQ-NN) * U0C( IQ,LU )             FXS 0137
            FLUP(LU) = FLUP(LU) + CWT(IQ-NN) * CMU(IQ-NN) * U0C( IQ,LU )FXS 0138
40       CONTINUE                                                       FXS 0139
                                                                        FXS 0140
         FLUP( LU )  = 2.0 * PI * FLUP( LU )                            FXS 0141
         FLDN( LU )  = 2.0 * PI * FLDN( LU )                            FXS 0142
         FDNTOT = FLDN( LU ) + FLDIR( LU )                              FXS 0143
         FNET   = FDNTOT - FLUP( LU )                                   FXS 0144
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )                            FXS 0145
         UAVG( LU ) = ( 2.0 * PI * UAVG(LU) + DIRINT ) / ( 4.*PI )      FXS 0146
         PLSORC =  XR0(LYU) + XR1(LYU) * UTAUPR(LU)                     FXS 0147
         DFDT( LU ) = ( 1.0-SSALB(LYU) ) * 4.*PI* ( UAVG(LU) - PLSORC ) FXS 0148
 90      IF( PRNT(2) )  WRITE( *,1020 ) UTAU(LU), LYU, RFLDIR(LU),      FXS 0149
     $                                 RFLDN(LU), FDNTOT, FLUP(LU),     FXS 0150
     $                                 FNET, UAVG(LU), PLSORC, DFDT(LU) FXS 0151
100   CONTINUE                                                          FXS 0152
                                                                        FXS 0153
      IF ( PRNT(3) )  THEN                                              FXS 0154
         WRITE ( *,1100 )                                               FXS 0155
         DO 200  LU = 1, NTAU                                           FXS 0156
            WRITE( *,1110 )  UTAU( LU )                                 FXS 0157
            DO  200  IQ = 1, NN                                         FXS 0158
               ANG1 = 180./PI * ACOS( CMU(2*NN-IQ+1) )                  FXS 0159
               ANG2 = 180./PI * ACOS( CMU(IQ) )                         FXS 0160
               WRITE( *,1120 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),        FXS 0161
     $                         ANG2, CMU(IQ),        U0C(IQ+NN,LU)      FXS 0162
200      CONTINUE                                                       FXS 0163
      END IF                                                            FXS 0164
                                                                        FXS 0165
1010  FORMAT( //, 21X,                                                  FXS 0166
     $ '<----------------------- FLUXES ----------------------->', /,   FXS 0167
     $ '   OPTICAL  COMPU    DOWNWARD    DOWNWARD    DOWNWARD     ',    FXS 0168
     $ ' UPWARD                    MEAN      PLANCK   D(NET FLUX)', /,  FXS 0169
     $ '     DEPTH  LAYER      DIRECT     DIFFUSE       TOTAL     ',    FXS 0170
     $ 'DIFFUSE         NET   INTENSITY      SOURCE   / D(OP DEP)', / ) FXS 0171
1020  FORMAT( F10.4, I7, 1P,7E12.3, E14.3 )                             FXS 0172
1100  FORMAT( //, ' ******** AZIMUTHALLY AVERAGED INTENSITIES',         FXS 0173
     $      ' ( AT POLAR QUADRATURE ANGLES ) *******' )                 FXS 0174
1110  FORMAT( /, ' OPTICAL DEPTH =', F10.4, //,                         FXS 0175
     $  '     ANGLE (DEG)   COS(ANGLE)     INTENSITY',                  FXS 0176
     $  '     ANGLE (DEG)   COS(ANGLE)     INTENSITY' )                 FXS 0177
1120  FORMAT( 2( 0P,F16.4, F13.5, 1P,E14.3 ) )                          FXS 0178
                                                                        FXS 0179
      RETURN                                                            FXS 0180
      END                                                               FXS 0181
