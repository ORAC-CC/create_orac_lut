      SUBROUTINE DEFALT(Z,P,T)                                          DEF 0001
C                                                                       DEF 0002
C     THIS SUBROUTINE INTERPOLATES PROFILES FROM THE 6 BUILT-IN MODEL   DEF 0003
C     ATMOSPHERIC TO ALTITUDE Z.  THE JUNIT INDICES FROM /CARD1B/       DEF 0004
C     INDICATE WHICH PROFILE SHOULD BE USED FOR EACH INTERPOLATION:     DEF 0005
C                                                                       DEF 0006
C                   JUNIT     MODEL ATMOSPHERE                          DEF 0007
C                     1          TROPICAL                               DEF 0008
C                     2          MID-LATITUDE SUMMER                    DEF 0009
C                     3          MID-LATITUDE WINTER                    DEF 0010
C                     4          HIGH-LAT SUMMER                        DEF 0011
C                     5          HIGH-LAT WINTER                        DEF 0012
C                     6          US STANDARD                            DEF 0013
C                                                                       DEF 0014
C     DECLARE ARGUMENTS:                                                DEF 0015
C       Z        INPUT ALTITUDE [KM].                                   DEF 0016
C       P        OUTPUT PRESSURE FROM MODEL ATMOSPHERE.                 DEF 0017
C       T        OUTPUT TEMPERATURE FROM MODEL ATMOSPHERE.              DEF 0018
      REAL Z,P,T                                                        DEF 0019
C                                                                       DEF 0020
C     INCLUDE PARAMETERS:                                               DEF 0021
      INCLUDE 'PARAM.LST'                                               DEF 0022
C                                                                       DEF 0023
C     LIST COMMONS:                                                     DEF 0024
      INTEGER JUNITP,JUNITT,JUNIT                                       DEF 0025
      REAL WMOL,WAIR,JLOW                                               DEF 0026
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),WAIR,JLOW          DEF 0027
      REAL ALT,PMATM,TMATM,AMOL                                         DEF 0028
      COMMON/MLATM/ALT(50),PMATM(50,6),TMATM(50,6),AMOL(50,8,6)         DEF 0029
      INTEGER LAYXMX                                                    DEF 0030
      REAL ALTX,AMOLX                                                   DEF 0031
      COMMON/MLATMX/LAYXMX,ALTX(50),AMOLX(50,35)                        DEF 0032
      REAL TRAC                                                         DEF 0033
      COMMON/TRAC/TRAC(50,21)                                           DEF 0034
      REAL CO2RAT                                                       DEF 0035
      COMMON/CO2MIX/CO2RAT                                              DEF 0036
      INTEGER JUNITX                                                    DEF 0037
      REAL WMOLX                                                        DEF 0038
      COMMON/CRD1BX/JUNITX,WMOLX(MMOLX)                                 DEF 0039
C                                                                       DEF 0040
C     DECLARE FUNCTIONS:                                                DEF 0041
      REAL EXPINT                                                       DEF 0042
C                                                                       DEF 0043
C     DECLARE LOCAL VARIABLES:                                          DEF 0044
      INTEGER I0,I1,I2,I3,K,KM7                                         DEF 0045
      LOGICAL LINIT                                                     DEF 0046
      REAL Z0,Z1,Z2,Z3,FACTOR,X0,X1,X2,X3                               DEF 0047
C                                                                       DEF 0048
C     DECLARE FUNCTIONS                                                 DEF 0049
      REAL FLG4PT                                                       DEF 0050
C                                                                       DEF 0051
C     LIST DATA:                                                        DEF 0052
C       LBOUND   NUMBER OF LAYER BOUNDARIES IN MODEL ATMOSPHERES.       DEF 0053
      INTEGER LBOUND                                                    DEF 0054
      DATA LBOUND/50/                                                   DEF 0055
C                                                                       DEF 0056
C     BEGIN CALCULATIONS:                                               DEF 0057
      IF(Z.LE.ALT(2))THEN                                               DEF 0058
C                                                                       DEF 0059
C         INPUT ALTITUDE Z IS BELOW SECOND ALTITUDE.                    DEF 0060
          I1=1                                                          DEF 0061
          I2=2                                                          DEF 0062
      ELSEIF(Z.GE.ALT(LBOUND-1))THEN                                    DEF 0063
C                                                                       DEF 0064
C         INPUT ALTITUDE Z IS ABOVE SECOND TO LAST ALTITUDE.            DEF 0065
          I1=LBOUND-1                                                   DEF 0066
          I2=LBOUND                                                     DEF 0067
      ELSE                                                              DEF 0068
C                                                                       DEF 0069
C         INTERMEDIATE ALTITUDE.                                        DEF 0070
          I2=3                                                          DEF 0071
          DO 10 I3=4,LBOUND                                             DEF 0072
              IF(Z.LE.ALT(I2))GOTO20                                    DEF 0073
   10     I2=I3                                                         DEF 0074
   20     CONTINUE                                                      DEF 0075
          I1=I2-1                                                       DEF 0076
          I0=I2-2                                                       DEF 0077
C                                                                       DEF 0078
C         LINEAR AND LAGRANGE 4-POINT INTERPOLATION COEFFICIENTS:       DEF 0079
          Z0=ALT(I0)                                                    DEF 0080
          Z1=ALT(I1)                                                    DEF 0081
          Z2=ALT(I2)                                                    DEF 0082
          Z3=ALT(I3)                                                    DEF 0083
          LINIT=.TRUE.                                                  DEF 0084
C                                                                       DEF 0085
C         TEST PRESSURE INPUT FLAG.                                     DEF 0086
          IF(JUNITP.LE.6)THEN                                           DEF 0087
              X0=LOG(PMATM(I0,JUNITP))                                  DEF 0088
              X1=LOG(PMATM(I1,JUNITP))                                  DEF 0089
              X2=LOG(PMATM(I2,JUNITP))                                  DEF 0090
              X3=LOG(PMATM(I3,JUNITP))                                  DEF 0091
              P=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)                 DEF 0092
              LINIT=.FALSE.                                             DEF 0093
              P=EXP(P)                                                  DEF 0094
          ENDIF                                                         DEF 0095
C                                                                       DEF 0096
C         TEST TEMPERATURE INPUT FLAG.                                  DEF 0097
          IF(JUNITT.LE.6)THEN                                           DEF 0098
              X0=TMATM(I0,JUNITT)                                       DEF 0099
              X1=TMATM(I1,JUNITT)                                       DEF 0100
              X2=TMATM(I2,JUNITT)                                       DEF 0101
              X3=TMATM(I3,JUNITT)                                       DEF 0102
              T=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)                 DEF 0103
              LINIT=.FALSE.                                             DEF 0104
          ENDIF                                                         DEF 0105
C                                                                       DEF 0106
C         MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:            DEF 0107
          DO 30 K=1,7                                                   DEF 0108
              IF(JUNIT(K).LE.6)THEN                                     DEF 0109
                  X0=AMOL(I0,K,JUNIT(K))                                DEF 0110
                  X1=AMOL(I1,K,JUNIT(K))                                DEF 0111
                  X2=AMOL(I2,K,JUNIT(K))                                DEF 0112
                  X3=AMOL(I3,K,JUNIT(K))                                DEF 0113
                  WMOL(K)=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)       DEF 0114
                  LINIT=.FALSE.                                         DEF 0115
C                                                                       DEF 0116
C                 ADJUST CO2 MIXING RATIO BASED ON CARD1A INPUT.        DEF 0117
                  IF(K.EQ.2)WMOL(2)=CO2RAT*WMOL(2)                      DEF 0118
                  JUNIT(K)=10                                           DEF 0119
              ENDIF                                                     DEF 0120
   30     CONTINUE                                                      DEF 0121
C                                                                       DEF 0122
C         MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:           DEF 0123
          DO 40 K=8,NSPC                                                DEF 0124
              IF(JUNIT(K).LE.6)THEN                                     DEF 0125
                  KM7=K-7                                               DEF 0126
                  X0=TRAC(I0,KM7)                                       DEF 0127
                  X1=TRAC(I1,KM7)                                       DEF 0128
                  X2=TRAC(I2,KM7)                                       DEF 0129
                  X3=TRAC(I3,KM7)                                       DEF 0130
                  WMOL(K)=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)       DEF 0131
                  LINIT=.FALSE.                                         DEF 0132
                  JUNIT(K)=10                                           DEF 0133
              ENDIF                                                     DEF 0134
   40     CONTINUE                                                      DEF 0135
          WMOL(12)=1000.*WMOL(12)                                       DEF 0136
C                                                                       DEF 0137
C         MOLECULAR DENSITIES FOR CFC SPECIES:                          DEF 0138
          IF(JUNITX.LE.6)THEN                                           DEF 0139
              DO 50 K=1,NSPECX                                          DEF 0140
                  X0=AMOLX(I0,K)                                        DEF 0141
                  X1=AMOLX(I1,K)                                        DEF 0142
                  X2=AMOLX(I2,K)                                        DEF 0143
                  X3=AMOLX(I3,K)                                        DEF 0144
                  WMOLX(K)=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)      DEF 0145
                  LINIT=.FALSE.                                         DEF 0146
   50         CONTINUE                                                  DEF 0147
          ENDIF                                                         DEF 0148
          RETURN                                                        DEF 0149
      ENDIF                                                             DEF 0150
C                                                                       DEF 0151
C     USE 2-POINT INTERPOLATION/EXTRAPOLATION NEAR ALTITUDE END POINTS: DEF 0152
      FACTOR=(Z-ALT(I1))/(ALT(I2)-ALT(I1))                              DEF 0153
      IF(JUNITP.LE.6)P=EXPINT(PMATM(I1,JUNITP),PMATM(I2,JUNITP),FACTOR) DEF 0154
      IF(JUNITT.LE.6)                                                   DEF 0155
     1  T=TMATM(I1,JUNITT)+FACTOR*(TMATM(I2,JUNITT)-TMATM(I1,JUNITT))   DEF 0156
C                                                                       DEF 0157
C     MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:                DEF 0158
      DO 60 K=1,7                                                       DEF 0159
          IF(JUNIT(K).LE.6)THEN                                         DEF 0160
              WMOL(K)=                                                  DEF 0161
     1          EXPINT(AMOL(I1,K,JUNIT(K)),AMOL(I2,K,JUNIT(K)),FACTOR)  DEF 0162
              IF(K.EQ.2)WMOL(2)=CO2RAT*WMOL(2)                          DEF 0163
              JUNIT(K)=10                                               DEF 0164
          ENDIF                                                         DEF 0165
   60 CONTINUE                                                          DEF 0166
C                                                                       DEF 0167
C     MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:               DEF 0168
      DO 70 K=8,NSPC                                                    DEF 0169
          IF(JUNIT(K).LE.6)THEN                                         DEF 0170
              WMOL(K)=EXPINT(TRAC(I1,K-7),TRAC(I2,K-7),FACTOR)          DEF 0171
              JUNIT(K)=10                                               DEF 0172
          ENDIF                                                         DEF 0173
   70 CONTINUE                                                          DEF 0174
      WMOL(12)=1000.*WMOL(12)                                           DEF 0175
C                                                                       DEF 0176
C     MOLECULAR DENSITIES FOR CFC SPECIES:                              DEF 0177
      IF(JUNITX.LE.6)THEN                                               DEF 0178
          DO 80 K=1,NSPECX                                              DEF 0179
              WMOLX(K)=EXPINT(AMOLX(I1,K),AMOLX(I2,K),FACTOR)           DEF 0180
   80     CONTINUE                                                      DEF 0181
      ENDIF                                                             DEF 0182
C                                                                       DEF 0183
C     INTERPOLATIONS COMPLETE.                                          DEF 0184
      RETURN                                                            DEF 0185
      END                                                               DEF 0186
      REAL FUNCTION FLG4PT(Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4,LINIT)             DEF 0187
C                                                                       DEF 0188
C     FUNCTION LAG4PT DOES A 4 POINT LAGRANGE INTERPOLATION FOR         DEF 0189
C     Z BETWEEN Z2 AND Z3, INCLUSIVE.  THE ABSCISSAE MUST BE            DEF 0190
C     MONOTIC (Z1<Z2<Z3<Z4 OR Z1>Z2>Z3>Z4) AND IF F(Z) HAS AN           DEF 0191
C     EXTREMUM BETWEEN Z2 AND Z3, THEN A LINEAR INTERPOLATION           DEF 0192
C     IS USED INSTEAD.  THE 4-POINT INTERPOLATION FORMULA IS:           DEF 0193
C                                                                       DEF 0194
C                                                                       DEF 0195
C     F(Z)  =  COEF1 P1     COEF2 P2     COEF3 P3     COEF4 P4          DEF 0196
C                                                                       DEF 0197
C                                                                       DEF 0198
C     WHERE                                                             DEF 0199
C                                                                       DEF 0200
C                 F1              F2              F3              F4    DEF 0201
C          P1 = ------ ;   P2 = ------ ;   P4 = ------ ;   P4 = ------  DEF 0202
C               DEMON1          DENOM2          DENOM3          DENOM4  DEF 0203
C                                                                       DEF 0204
C                                                                       DEF 0205
C          COEF1  =  (Z - Z2) (Z - Z3) (Z - Z4)                         DEF 0206
C                                                                       DEF 0207
C          COEF2  =  (Z - Z1) (Z - Z3) (Z - Z4)                         DEF 0208
C                                                                       DEF 0209
C          COEF3  =  (Z - Z1) (Z - Z2) (Z - Z4)                         DEF 0210
C                                                                       DEF 0211
C          COEF4  =  (Z - Z1) (Z - Z2) (Z - Z3)                         DEF 0212
C                                                                       DEF 0213
C                                                                       DEF 0214
C          DENOM1  =   (Z1 - Z2) (Z1 - Z3) (Z1 - Z4)                    DEF 0215
C                                                                       DEF 0216
C          DENOM2  =   (Z2 - Z1) (Z2 - Z3) (Z2 - Z4)                    DEF 0217
C                                                                       DEF 0218
C          DENOM3  =   (Z3 - Z1) (Z3 - Z2) (Z3 - Z4)                    DEF 0219
C                                                                       DEF 0220
C          DENOM4  =   (Z4 - Z1) (Z4 - Z2) (Z4 - Z3)                    DEF 0221
C                                                                       DEF 0222
C                                                                       DEF 0223
C     THE EXTREMA IN F(Z) ARE FOUND BY SOLVING                          DEF 0224
C                                                                       DEF 0225
C                                                                       DEF 0226
C           d F(Z)        2                                             DEF 0227
C           ------  =  A Z   -  2 B Z  +  C  =  0                       DEF 0228
C             dZ                                                        DEF 0229
C                                                                       DEF 0230
C                                                                       DEF 0231
C     WHERE                                                             DEF 0232
C                                                                       DEF 0233
C                                                                       DEF 0234
C        A  =  3 (P1 + P2 + P3 + P4)                                    DEF 0235
C                                                                       DEF 0236
C                                                                       DEF 0237
C        B  =  (Z2 + Z3 + Z4) P1  +  (Z1 + Z3 + Z4) P2                  DEF 0238
C                                                                       DEF 0239
C           +  (Z1 + Z2 + Z4) P3  +  (Z1 + Z2 + Z3) P4                  DEF 0240
C                                                                       DEF 0241
C                                                                       DEF 0242
C        C  =  (Z2 Z3 + Z2 Z4 + Z3 Z4) P1  +  (Z1 Z3 + Z1 Z4 + Z3 Z4) P2DEF 0243
C                                                                       DEF 0244
C           +  (Z1 Z2 + Z1 Z4 + Z2 Z4) P3  +  (Z1 Z2 + Z1 Z3 + Z2 Z3) P4DEF 0245
C                                                                       DEF 0246
C                                                                       DEF 0247
C     THE DISCRIMINANT OF THE QUADRATIC EQUATION FOR Z                  DEF 0248
C     IS EXPANDED IN QUADRATIC TERMS OF P:                              DEF 0249
C                                                                       DEF 0250
C                           _           _                               DEF 0251
C          2               \           \                                DEF 0252
C         B  - 4 A C   =    |  Pi       |       Pj  COEFij              DEF 0253
C                          /_          /_                               DEF 0254
C                           i      j=i and j>i                          DEF 0255
C                                                                       DEF 0256
C                                                                       DEF 0257
C     WHERE                                                             DEF 0258
C                                                                       DEF 0259
C                     1           2           2             2           DEF 0260
C          COEFii  =  -  { (Zj-Zk)   + (Zj-Zl)   +   (Zk-Zl)  }         DEF 0261
C                     2                                                 DEF 0262
C                                                                       DEF 0263
C                              2                                        DEF 0264
C          COEFij  =  2 (Zk-Zl)  +  (Zi-Zk) (Zj-Zl)  +  (Zi-Zl) (Zj-Zk) DEF 0265
C                                                                       DEF 0266
C                                                                       DEF 0267
C     DECLARE INPUTS                                                    DEF 0268
C       Z        ABSCISSA OF DESIRED ORDINATE.                          DEF 0269
C       Z1       FIRST ABSCISSA                                         DEF 0270
C       Z2       SECOND ABSCISSA                                        DEF 0271
C       Z3       THIRD ABSCISSA                                         DEF 0272
C       Z4       FOURTH ABSCISSA                                        DEF 0273
C       F1       FIRST ORDINATE                                         DEF 0274
C       F2       SECOND ORDINATE                                        DEF 0275
C       F3       THIRD ORDINATE                                         DEF 0276
C       F4       FOURTH ORDINATE                                        DEF 0277
C       LINIT    FLAG, FALSE IF ABSCISSAE ARE UNCHANGED.                DEF 0278
      REAL Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4                                    DEF 0279
      LOGICAL LINIT                                                     DEF 0280
C                                                                       DEF 0281
C     LIST COMMONS:                                                     DEF 0282
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               DEF 0283
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           DEF 0284
C                                                                       DEF 0285
C     SAVED VARIABLES.                                                  DEF 0286
      REAL BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  DEF 0287
     1  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, DEF 0288
     2  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4  DEF 0289
      SAVE BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  DEF 0290
     1  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, DEF 0291
     2  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4  DEF 0292
C                                                                       DEF 0293
C     LOCAL VARIABLES.                                                  DEF 0294
      REAL A,B,ZEXTRM,DSCRIM,ROOT,ZDIFF1,ZDIFF2,ZDIFF3,ZDIFF4,          DEF 0295
     1  ZDIF12,ZDIF13,ZDIF14,ZDIF23,ZDIF24,ZDIF34,P1,P2,P3,P4,          DEF 0296
     2  ZD12SQ,ZD13SQ,ZD14SQ,ZD23SQ,ZD24SQ,ZD34SQ,STORE                 DEF 0297
C     LIST DATA.                                                        DEF 0298
C       SMALL    SMALL NUMBER TOLERANCE                                 DEF 0299
      REAL SMALL                                                        DEF 0300
      DATA SMALL/1.E-30/                                                DEF 0301
C                                                                       DEF 0302
C     NEW ABSCISSAE CHECK                                               DEF 0303
      IF(LINIT)THEN                                                     DEF 0304
C                                                                       DEF 0305
C         CHECK INPUTS                                                  DEF 0306
          IF(Z2.EQ.Z3)GOTO10                                            DEF 0307
          IF(Z2.LT.Z3)THEN                                              DEF 0308
              IF(Z1.GE.Z2 .OR. Z2.GT.Z .OR. Z.GT.Z3 .OR. Z3.GE.Z4)GOTO10DEF 0309
          ELSE                                                          DEF 0310
              IF(Z1.LE.Z2 .OR. Z2.LT.Z .OR. Z.LT.Z3 .OR. Z3.LE.Z4)GOTO10DEF 0311
          ENDIF                                                         DEF 0312
C                                                                       DEF 0313
C         DEFINE CONSTANTS DEPENDENT ON Z, Z1, Z2, Z3 AND Z4 ONLY.      DEF 0314
          ZDIFF1=Z-Z1                                                   DEF 0315
          ZDIFF2=Z-Z2                                                   DEF 0316
          ZDIFF3=Z-Z3                                                   DEF 0317
          ZDIFF4=Z-Z4                                                   DEF 0318
          COEF1=ZDIFF2*ZDIFF3*ZDIFF4                                    DEF 0319
          COEF2=ZDIFF1*ZDIFF3*ZDIFF4                                    DEF 0320
          COEF3=ZDIFF1*ZDIFF2*ZDIFF4                                    DEF 0321
          COEF4=ZDIFF1*ZDIFF2*ZDIFF3                                    DEF 0322
          ZDIF12=Z1-Z2                                                  DEF 0323
          ZDIF13=Z1-Z3                                                  DEF 0324
          ZDIF14=Z1-Z4                                                  DEF 0325
          ZDIF23=Z2-Z3                                                  DEF 0326
          ZDIF24=Z2-Z4                                                  DEF 0327
          ZDIF34=Z3-Z4                                                  DEF 0328
          DENOM1= ZDIF12*ZDIF13*ZDIF14                                  DEF 0329
          DENOM2=-ZDIF12*ZDIF23*ZDIF24                                  DEF 0330
          DENOM3= ZDIF13*ZDIF23*ZDIF34                                  DEF 0331
          DENOM4=-ZDIF14*ZDIF24*ZDIF34                                  DEF 0332
C                                                                       DEF 0333
C         CONSTANTS USED IN EXTREMUM CHECK                              DEF 0334
          BCOEF1=Z2+Z3+Z4                                               DEF 0335
          BCOEF2=Z3+Z4+Z1                                               DEF 0336
          BCOEF3=Z4+Z1+Z2                                               DEF 0337
          BCOEF4=Z1+Z2+Z3                                               DEF 0338
          ZD12SQ=ZDIF12**2                                              DEF 0339
          ZD13SQ=ZDIF13**2                                              DEF 0340
          ZD14SQ=ZDIF14**2                                              DEF 0341
          ZD23SQ=ZDIF23**2                                              DEF 0342
          ZD24SQ=ZDIF24**2                                              DEF 0343
          ZD34SQ=ZDIF34**2                                              DEF 0344
          COEF11=.5*(ZD23SQ+ZD24SQ+ZD34SQ)                              DEF 0345
          COEF22=.5*(ZD13SQ+ZD14SQ+ZD34SQ)                              DEF 0346
          COEF33=.5*(ZD12SQ+ZD14SQ+ZD24SQ)                              DEF 0347
          COEF44=.5*(ZD12SQ+ZD13SQ+ZD23SQ)                              DEF 0348
          STORE=ZDIF13*ZDIF24+ZDIF14*ZDIF23                             DEF 0349
          COEF12=2*ZD34SQ+STORE                                         DEF 0350
          COEF34=2*ZD12SQ+STORE                                         DEF 0351
          STORE=ZDIF12*ZDIF34-ZDIF14*ZDIF23                             DEF 0352
          COEF13=2*ZD24SQ+STORE                                         DEF 0353
          COEF24=2*ZD13SQ+STORE                                         DEF 0354
          STORE=ZDIF12*ZDIF34+ZDIF13*ZDIF24                             DEF 0355
          COEF14=2*ZD23SQ-STORE                                         DEF 0356
          COEF23=2*ZD14SQ-STORE                                         DEF 0357
          CCOEF1=Z2*Z3+Z2*Z4+Z3*Z4                                      DEF 0358
          CCOEF2=Z1*Z3+Z1*Z4+Z3*Z4                                      DEF 0359
          CCOEF3=Z1*Z2+Z1*Z4+Z2*Z4                                      DEF 0360
          CCOEF4=Z1*Z2+Z1*Z3+Z2*Z3                                      DEF 0361
C                                                                       DEF 0362
C         LINEAR INTERPOLATION FACTOR                                   DEF 0363
          FACTOR=(Z-Z2)/(Z3-Z2)                                         DEF 0364
      ENDIF                                                             DEF 0365
C                                                                       DEF 0366
C     DETERMINE LOCATION OF EXTREMA                                     DEF 0367
      P1=F1/DENOM1                                                      DEF 0368
      P2=F2/DENOM2                                                      DEF 0369
      P3=F3/DENOM3                                                      DEF 0370
      P4=F4/DENOM4                                                      DEF 0371
      A=3*(P1+P2+P3+P4)                                                 DEF 0372
      B=P1*BCOEF1+P2*BCOEF2+P3*BCOEF3+P4*BCOEF4                         DEF 0373
      IF(ABS(A).LT.SMALL)THEN                                           DEF 0374
          IF(ABS(B).GT.SMALL)THEN                                       DEF 0375
C                                                                       DEF 0376
C             ONE EXTREMUM                                              DEF 0377
              ZEXTRM=.5*(CCOEF1*P1+CCOEF2*P2+CCOEF3*P3+CCOEF4*P4)/B     DEF 0378
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 DEF 0379
     1           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN                 DEF 0380
C                                                                       DEF 0381
C                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.DEF 0382
                  FLG4PT=F2+FACTOR*(F3-F2)                              DEF 0383
                  RETURN                                                DEF 0384
              ENDIF                                                     DEF 0385
          ENDIF                                                         DEF 0386
      ELSE                                                              DEF 0387
          DSCRIM=P1*(P1*COEF11+P2*COEF12+P3*COEF13+P4*COEF14)           DEF 0388
     1          +P2*(P2*COEF22+P3*COEF23+P4*COEF24)                     DEF 0389
     2          +P3*(P3*COEF33+P4*COEF34)+P4**2*COEF44                  DEF 0390
          IF(DSCRIM.GE.0.)THEN                                          DEF 0391
C                                                                       DEF 0392
C             TWO EXTREMUM                                              DEF 0393
              ROOT=SQRT(DSCRIM)                                         DEF 0394
              ZEXTRM=(B+ROOT)/A                                         DEF 0395
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 DEF 0396
     1           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN                 DEF 0397
C                                                                       DEF 0398
C                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.DEF 0399
                  FLG4PT=F2+FACTOR*(F3-F2)                              DEF 0400
                  RETURN                                                DEF 0401
              ENDIF                                                     DEF 0402
              ZEXTRM=(B-ROOT)/A                                         DEF 0403
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 DEF 0404
     1           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN                 DEF 0405
C                                                                       DEF 0406
C                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.DEF 0407
                  FLG4PT=F2+FACTOR*(F3-F2)                              DEF 0408
                  RETURN                                                DEF 0409
              ENDIF                                                     DEF 0410
          ENDIF                                                         DEF 0411
      ENDIF                                                             DEF 0412
C                                                                       DEF 0413
C     NO EXTREMA BETWEEN Z2 AND Z3.  USE 4-POINT LAGRANGE INTERPOLATION.DEF 0414
      FLG4PT=COEF1*P1+COEF2*P2+COEF3*P3+COEF4*P4                        DEF 0415
      RETURN                                                            DEF 0416
C                                                                       DEF 0417
C     INCORRECT INPUT ERROR.                                            DEF 0418
   10 CONTINUE                                                          DEF 0419
      WRITE(IPR,'(/A,/(18X,A,F12.4))')                                  DEF 0420
     1  ' ERROR in FLG4PT:  ABSCISSAE OUT OF ORDER.',                   DEF 0421
     2  ' Z1 =',Z1,' Z2 =',Z2,' Z  =',Z,' Z3 =',Z3,' Z4 =',Z4           DEF 0422
      STOP ' ERROR in FLG4PT:  ABSCISSAE OUT OF ORDER.'                 DEF 0423
      END                                                               DEF 0424
