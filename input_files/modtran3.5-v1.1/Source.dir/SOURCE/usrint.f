      SUBROUTINE  USRINT( BPLANK, CMU, CWT, DELM0, EMU, EXPBEA,         UIN 0001
     $                    FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL,  UIN 0002
     $                    LYRCUT, MAZIM, MXCMU, MXULV, MXUMU, NCUT,     UIN 0003
     $                    NLYR, NN, NSTR, PLANK, NUMU, NTAU, PI, RMU,   UIN 0004
     $                    TAUCPR, TPLANK, UMU, UMU0, UTAUPR, WK,        UIN 0005
     $                    ZBEAM, Z0U, Z1U, ZZ, ZPLK0, ZPLK1, UUM )      UIN 0006
                                                                        UIN 0007
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            UIN 0008
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                UIN 0009
C       COMPUTES INTENSITY COMPONENTS AT USER OUTPUT ANGLES             UIN 0010
C       FOR AZIMUTHAL EXPANSION TERMS IN EQ. SD(2)                      UIN 0011
                                                                        UIN 0012
C   I N P U T    V A R I A B L E S:                                     UIN 0013
                                                                        UIN 0014
C       BPLANK :  INTEGRATED PLANCK FUNCTION FOR EMISSION FROM          UIN 0015
C                 BOTTOM BOUNDARY                                       UIN 0016
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE      UIN 0017
C       CWT    :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE        UIN 0018
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0                         UIN 0019
C       EMU    :  SURFACE DIRECTIONAL EMISSIVITY (USER ANGLES)          UIN 0020
C       EXPBEA :  TRANSMISSION OF INCIDENT BEAM, EXP(-TAUCPR/UMU0)      UIN 0021
C       GC     :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)        UIN 0022
C       GU     :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES        UIN 0023
C                 (I.E., G IN EQ. SC(1) )                               UIN 0024
C       KK     :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)             UIN 0025
C       LAYRU  :  LAYER NUMBER OF USER LEVEL -UTAU-                     UIN 0026
C       LL     :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED       UIN 0027
C                 BY SOLVING SCALED VERSION OF EQ. SC(5);               UIN 0028
C                 EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED           UIN 0029
C       LYRCUT :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER          UIN 0030
C       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT                          UIN 0031
C       NCUT   :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED       UIN 0032
C       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)             UIN 0033
C       RMU    :  SURFACE BIDIRECTIONAL REFLECTIVITY (USER ANGLES)      UIN 0034
C       TAUCPR :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)             UIN 0035
C       TPLANK :  INTEGRATED PLANCK FUNCTION FOR EMISSION FROM          UIN 0036
C                 TOP BOUNDARY                                          UIN 0037
C       UTAUPR :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M       UIN 0038
C                    COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M       UIN 0039
C       Z0U    :  Z-SUB-ZERO IN EQ. SS(16) INTERPOLATED TO USER         UIN 0040
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)           UIN 0041
C       Z1U    :  Z-SUB-ONE IN EQ. SS(16) INTERPOLATED TO USER          UIN 0042
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)           UIN 0043
C       ZZ     :  BEAM SOURCE VECTORS IN EQ. SS(19)                     UIN 0044
C       ZPLK0  :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)    UIN 0045
C       ZPLK1  :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)    UIN 0046
C       ZBEAM  :  INCIDENT-BEAM SOURCE VECTORS                          UIN 0047
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        UIN 0048
                                                                        UIN 0049
C   O U T P U T    V A R I A B L E S:                                   UIN 0050
                                                                        UIN 0051
C       UUM  :  AZIMUTHAL COMPONENTS OF THE INTENSITY IN EQ. STWJ(5)    UIN 0052
                                                                        UIN 0053
C   I N T E R N A L    V A R I A B L E S:                               UIN 0054
                                                                        UIN 0055
C       BNDDIR :  DIRECT INTENSITY DOWN AT THE BOTTOM BOUNDARY          UIN 0056
C       BNDDFU :  DIFFUSE INTENSITY DOWN AT THE BOTTOM BOUNDARY         UIN 0057
C       BNDINT :  INTENSITY ATTENUATED AT BOTH BOUNDARIES, STWJ(25-6)   UIN 0058
C       DTAU   :  OPTICAL DEPTH OF A COMPUTATIONAL LAYER                UIN 0059
C       LYREND :  END LAYER OF INTEGRATION                              UIN 0060
C       LYRSTR :  START LAYER OF INTEGRATION                            UIN 0061
C       PALINT :  INTENSITY COMPONENT FROM PARALLEL BEAM                UIN 0062
C       PLKINT :  INTENSITY COMPONENT FROM PLANCK SOURCE                UIN 0063
C       WK     :  SCRATCH VECTOR FOR SAVING 'EXP' EVALUATIONS           UIN 0064
C       ALL THE EXPONENTIAL FACTORS ( EXP1, EXPN,... ETC.)              UIN 0065
C       COME FROM THE SUBSTITUTION OF CONSTANTS OF INTEGRATION IN       UIN 0066
C       EQ. SC(12) INTO EQS. S1(8-9).  THEY ALL HAVE NEGATIVE           UIN 0067
C       ARGUMENTS SO THERE SHOULD NEVER BE OVERFLOW PROBLEMS.           UIN 0068
C+---------------------------------------------------------------------+UIN 0069
                                                                        UIN 0070
      LOGICAL  LAMBER, LYRCUT, PLANK, NEGUMU                            UIN 0071
      INTEGER  LAYRU(*)                                                 UIN 0072
      REAL*8 CMU(*),CWT(*),EMU(*),EXPBEA(0:*),GC( MXCMU,MXCMU,* ),      UIN 0073
     $         GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),       UIN 0074
     $         RMU( MXUMU,0:* ), TAUCPR( 0:* ), UUM( MXUMU,MXULV,0:* ), UIN 0075
     $         UMU(*), UTAUPR(*), WK(*), Z0U( MXUMU,* ), Z1U( MXUMU,* ),UIN 0076
     $         ZBEAM( MXUMU,* ), ZZ( MXCMU,* ), ZPLK0( MXCMU,* ),       UIN 0077
     $         ZPLK1( MXCMU,* )                                         UIN 0078
                                                                        UIN 0079
                                                                        UIN 0080
      CALL  ZEROIT( UUM, MXUMU*MXULV*(MXCMU+1) )                        UIN 0081
                                                                        UIN 0082
C                          ** INCORPORATE CONSTANTS OF INTEGRATION INTO UIN 0083
C                          ** INTERPOLATED EIGENVECTORS                 UIN 0084
      DO 10  LC = 1, NCUT                                               UIN 0085
         DO  10  IQ = 1, NSTR                                           UIN 0086
            DO 10  IU = 1, NUMU                                         UIN 0087
               GU(IU,IQ,LC) = GU(IU,IQ,LC) * LL(IQ,LC)                  UIN 0088
10    CONTINUE                                                          UIN 0089
C                           ** LOOP OVER LEVELS AT WHICH INTENSITIES    UIN 0090
C                           ** ARE DESIRED ('USER OUTPUT LEVELS')       UIN 0091
      DO 200  LU = 1, NTAU                                              UIN 0092
                                                                        UIN 0093
         IF ( FBEAM .GT. 0.0 )  EXP0 = DEXP( - UTAUPR(LU) / UMU0 )      UIN 0094
         LYU = LAYRU(LU)                                                UIN 0095
C                              ** LOOP OVER POLAR ANGLES AT WHICH       UIN 0096
C                              ** INTENSITIES ARE DESIRED               UIN 0097
         DO 100  IU = 1, NUMU                                           UIN 0098
            IF ( LYRCUT .AND. LYU.GT.NCUT )  GO TO 100                  UIN 0099
            NEGUMU = UMU(IU).LT.0.0                                     UIN 0100
            IF( NEGUMU )  THEN                                          UIN 0101
               LYRSTR = 1                                               UIN 0102
               LYREND = LYU - 1                                         UIN 0103
               SGN = - 1.0                                              UIN 0104
            ELSE                                                        UIN 0105
               LYRSTR = LYU + 1                                         UIN 0106
               LYREND = NCUT                                            UIN 0107
               SGN = 1.0                                                UIN 0108
            END IF                                                      UIN 0109
C                          ** FOR DOWNWARD INTENSITY, INTEGRATE FROM TOPUIN 0110
C                          ** TO 'LYU-1' IN EQ. S1(8); FOR UPWARD,      UIN 0111
C                          ** INTEGRATE FROM BOTTOM TO 'LYU+1' IN S1(9) UIN 0112
            PALINT = 0.0                                                UIN 0113
            PLKINT = 0.0                                                UIN 0114
            DO 30  LC = LYRSTR, LYREND                                  UIN 0115
                                                                        UIN 0116
               DTAU = TAUCPR(LC) - TAUCPR(LC-1)                         UIN 0117
               EXP1 =  DEXP( (UTAUPR(LU) - TAUCPR(LC-1)) / UMU(IU) )    UIN 0118
               EXP2 =  DEXP( (UTAUPR(LU) - TAUCPR( LC )) / UMU(IU) )    UIN 0119
                                                                        UIN 0120
               IF ( PLANK .AND. MAZIM.EQ.0 )                            UIN 0121
     $           PLKINT = PLKINT + SGN * ( Z0U(IU,LC) * (EXP1 - EXP2) + UIN 0122
     $                    Z1U(IU,LC) * ( (TAUCPR(LC-1) + UMU(IU))*EXP1 -UIN 0123
     $                                   (TAUCPR(LC) + UMU(IU))*EXP2 ) )UIN 0124
                                                                        UIN 0125
               IF ( FBEAM.GT.0.0 )  THEN                                UIN 0126
                  DENOM = 1.0 + UMU(IU) / UMU0                          UIN 0127
                  IF ( DABS(DENOM).LT.0.0001 ) THEN                     UIN 0128
C                                                   ** L'HOSPITAL LIMIT UIN 0129
                     EXPN = ( DTAU / UMU0 ) * EXP0                      UIN 0130
                  ELSE                                                  UIN 0131
                     EXPN = ( EXP1 * EXPBEA(LC-1) - EXP2 * EXPBEA(LC) ) UIN 0132
     $                      * SGN / DENOM                               UIN 0133
                  END IF                                                UIN 0134
                  PALINT = PALINT + ZBEAM(IU,LC) * EXPN                 UIN 0135
               ENDIF                                                    UIN 0136
C                                                   ** -KK- IS NEGATIVE UIN 0137
               DO 20  IQ = 1, NN                                        UIN 0138
                  WK(IQ) = DEXP( KK(IQ,LC) * DTAU )                     UIN 0139
                  DENOM = 1.0 + UMU(IU) * KK(IQ,LC)                     UIN 0140
                  IF ( DABS(DENOM).LT.0.0001 ) THEN                     UIN 0141
C                                                   ** L'HOSPITAL LIMIT UIN 0142
                     EXPN = DTAU / UMU(IU) * EXP2                       UIN 0143
                  ELSE                                                  UIN 0144
                     EXPN = SGN * ( EXP1 * WK(IQ) - EXP2 ) / DENOM      UIN 0145
                  END IF                                                UIN 0146
                  PALINT = PALINT + GU(IU,IQ,LC) * EXPN                 UIN 0147
20             CONTINUE                                                 UIN 0148
C                                                   ** -KK- IS POSITIVE UIN 0149
               DO 21  IQ = NN+1, NSTR                                   UIN 0150
                  DENOM = 1.0 + UMU(IU) * KK(IQ,LC)                     UIN 0151
                  IF ( DABS(DENOM).LT.0.0001 ) THEN                     UIN 0152
C                                                   ** L'HOSPITAL LIMIT UIN 0153
                     EXPN = - DTAU / UMU(IU) * EXP1                     UIN 0154
                  ELSE                                                  UIN 0155
                     EXPN = SGN *( EXP1 - EXP2 * WK(NSTR+1-IQ) ) / DENOMUIN 0156
                  END IF                                                UIN 0157
                  PALINT = PALINT + GU(IU,IQ,LC) * EXPN                 UIN 0158
21             CONTINUE                                                 UIN 0159
                                                                        UIN 0160
30          CONTINUE                                                    UIN 0161
C                           ** CALCULATE CONTRIBUTION FROM USER         UIN 0162
C                           ** OUTPUT LEVEL TO NEXT COMPUTATIONAL LEVEL UIN 0163
                                                                        UIN 0164
            DTAU1 = UTAUPR(LU) - TAUCPR(LYU-1)                          UIN 0165
            DTAU2 = UTAUPR(LU) - TAUCPR(LYU)                            UIN 0166
            IF( DABS(DTAU1).LT.1.E-6 .AND. NEGUMU )  GO TO 50           UIN 0167
            IF( DABS(DTAU2).LT.1.E-6 .AND. (.NOT.NEGUMU) )  GO TO 50    UIN 0168
            IF( NEGUMU ) EXP1 = DEXP( DTAU1 / UMU(IU) )                 UIN 0169
            IF( .NOT.NEGUMU ) EXP2 = DEXP( DTAU2 / UMU(IU) )            UIN 0170
                                                                        UIN 0171
            IF ( FBEAM.GT.0.0 )  THEN                                   UIN 0172
               DENOM = 1.0 + UMU(IU) / UMU0                             UIN 0173
               IF ( DABS(DENOM).LT.0.0001 ) THEN                        UIN 0174
                  EXPN =  ( DTAU1 / UMU0 ) * EXP0                       UIN 0175
               ELSE IF ( NEGUMU ) THEN                                  UIN 0176
                  EXPN = ( EXP0 - EXPBEA(LYU-1) * EXP1 ) / DENOM        UIN 0177
               ELSE                                                     UIN 0178
                  EXPN = ( EXP0 - EXPBEA(LYU) * EXP2 ) / DENOM          UIN 0179
               END IF                                                   UIN 0180
               PALINT = PALINT + ZBEAM(IU,LYU) * EXPN                   UIN 0181
            ENDIF                                                       UIN 0182
C                                                   ** -KK- IS NEGATIVE UIN 0183
            DTAU = TAUCPR(LYU) - TAUCPR(LYU-1)                          UIN 0184
            DO 40  IQ = 1, NN                                           UIN 0185
               DENOM = 1.0 + UMU(IU) * KK(IQ,LYU)                       UIN 0186
               IF ( DABS(DENOM).LT.0.0001 ) THEN                        UIN 0187
                  EXPN = - DTAU2 / UMU(IU) * EXP2                       UIN 0188
               ELSE IF ( NEGUMU ) THEN                                  UIN 0189
                  EXPN = ( DEXP( - KK(IQ,LYU) * DTAU2 ) -               UIN 0190
     $                     DEXP( KK(IQ,LYU) * DTAU ) * EXP1 ) / DENOM   UIN 0191
               ELSE                                                     UIN 0192
                  EXPN = ( DEXP( - KK(IQ,LYU) * DTAU2 ) - EXP2 ) / DENOMUIN 0193
               END IF                                                   UIN 0194
               PALINT = PALINT + GU(IU,IQ,LYU) * EXPN                   UIN 0195
40          CONTINUE                                                    UIN 0196
C                                                   ** -KK- IS POSITIVE UIN 0197
            DO 41  IQ = NN+1, NSTR                                      UIN 0198
               DENOM = 1.0 + UMU(IU) * KK(IQ,LYU)                       UIN 0199
               IF ( DABS(DENOM).LT.0.0001 ) THEN                        UIN 0200
                  EXPN = - DTAU1 / UMU(IU) * EXP1                       UIN 0201
               ELSE IF ( NEGUMU ) THEN                                  UIN 0202
                  EXPN = ( DEXP(- KK(IQ,LYU) * DTAU1 ) - EXP1 ) / DENOM UIN 0203
               ELSE                                                     UIN 0204
                  EXPN = ( DEXP( - KK(IQ,LYU) * DTAU1 ) -               UIN 0205
     $                     DEXP( - KK(IQ,LYU) * DTAU ) * EXP2 ) / DENOM UIN 0206
               END IF                                                   UIN 0207
               PALINT = PALINT + GU(IU,IQ,LYU) * EXPN                   UIN 0208
41          CONTINUE                                                    UIN 0209
                                                                        UIN 0210
            IF ( PLANK .AND. MAZIM.EQ.0 )  THEN                         UIN 0211
              IF ( NEGUMU ) THEN                                        UIN 0212
                 EXPN = EXP1                                            UIN 0213
                 FACT = TAUCPR(LYU-1) + UMU(IU)                         UIN 0214
              ELSE                                                      UIN 0215
                 EXPN = EXP2                                            UIN 0216
                 FACT = TAUCPR( LYU ) + UMU(IU)                         UIN 0217
              END IF                                                    UIN 0218
              PLKINT = PLKINT + Z0U(IU,LYU) * ( 1.- EXPN ) +            UIN 0219
     $                 Z1U(IU,LYU) *( UTAUPR(LU) + UMU(IU) - FACT*EXPN )UIN 0220
            END IF                                                      UIN 0221
C                            ** CALCULATE INTENSITY COMPONENTS          UIN 0222
C                            ** ATTENUATED AT BOTH BOUNDARIES.          UIN 0223
C                            ** NOTE:: NO AZIMUTHAL INTENSITY           UIN 0224
C                            ** COMPONENT FOR ISOTROPIC SURFACE         UIN 0225
50          BNDINT = 0.0                                                UIN 0226
            IF ( NEGUMU .AND. MAZIM.EQ.0 ) THEN                         UIN 0227
              BNDINT = ( FISOT + TPLANK ) * DEXP( UTAUPR(LU) / UMU(IU) )UIN 0228
            ELSE IF ( .NOT.NEGUMU ) THEN                                UIN 0229
              IF ( LYRCUT .OR. (LAMBER .AND. MAZIM.GT.0) )  GO TO 90    UIN 0230
              DO 60  JQ = NN+1, NSTR                                    UIN 0231
           WK(JQ) = DEXP(-KK(JQ,NLYR)*(TAUCPR(NLYR)-TAUCPR(NLYR-1)))    UIN 0232
60            CONTINUE                                                  UIN 0233
              BNDDFU = 0.0                                              UIN 0234
              DO 80  IQ = NN, 1, -1                                     UIN 0235
                 DFUINT = 0.0                                           UIN 0236
                 DO 70  JQ = 1, NN                                      UIN 0237
                    DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)      UIN 0238
70               CONTINUE                                               UIN 0239
                 DO 71  JQ = NN+1, NSTR                                 UIN 0240
                    DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)      UIN 0241
     $                                * WK(JQ)                          UIN 0242
71               CONTINUE                                               UIN 0243
                 IF ( FBEAM.GT.0.0 )                                    UIN 0244
     $                DFUINT = DFUINT + ZZ(IQ,NLYR) * EXPBEA(NLYR)      UIN 0245
                 DFUINT = DFUINT + DELM0 * ( ZPLK0(IQ,NLYR)             UIN 0246
     $                                + ZPLK1(IQ,NLYR) * TAUCPR(NLYR) ) UIN 0247
                 BNDDFU = BNDDFU + ( 1. + DELM0 ) * RMU(IU,NN+1-IQ)     UIN 0248
     $                           * CMU(NN+1-IQ) * CWT(NN+1-IQ) * DFUINT UIN 0249
80            CONTINUE                                                  UIN 0250
                                                                        UIN 0251
              BNDDIR = 0.0                                              UIN 0252
              IF (FBEAM.GT.0.0) BNDDIR = UMU0*FBEAM/PI * RMU(IU,0)      UIN 0253
     $                                   * EXPBEA(NLYR)                 UIN 0254
              BNDINT = ( BNDDFU + BNDDIR + DELM0 * EMU(IU) * BPLANK )   UIN 0255
     $                 * DEXP( (UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU) )    UIN 0256
            END IF                                                      UIN 0257
                                                                        UIN 0258
90          UUM( IU, LU, MAZIM ) = PALINT + PLKINT + BNDINT             UIN 0259
                                                                        UIN 0260
100      CONTINUE                                                       UIN 0261
200   CONTINUE                                                          UIN 0262
                                                                        UIN 0263
      RETURN                                                            UIN 0264
      END                                                               UIN 0265
