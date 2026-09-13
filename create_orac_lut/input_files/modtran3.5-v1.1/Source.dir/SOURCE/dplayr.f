      SUBROUTINE DPLAYR(J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,S,BEND,RNG)   LAY 0001
      INCLUDE 'PARAM.LST'                                               LAY 0002
C                                                                       LAY 0003
C     DOUBLE PRECISION VERSION OF THE PREVIOUS ROUTINE "LAYER"          LAY 0004
C*****************************************************************      LAY 0005
C     THIS SUBROUTINE CALCULATES THE REFRACTED PATH FROM Z1 TO Z2       LAY 0006
C     WITH THE SIN OF THE INITIAL INCIDENCE ANGLE SINAI                 LAY 0007
C*****************************************************************      LAY 0008
      DOUBLE PRECISION ZP,PP,TP,RFNDXP,SP,PPSUM,TPSUM,RHOPSM,DENP,      LAY 0009
     $     AMTP,DPRE                                                    LAY 0010
      DOUBLE PRECISION HDEN(KMAX),DENA(KMAX),DENB(KMAX),DPANDX          LAY 0011
      REAL EPSILN,SPZP,SPPP,SPTP,SPDENP,SPSP,SPPPSU,SPTPSU,SPAMTP       LAY 0012
      DOUBLE PRECISION Z1,Z2,H1,R1,DHMIN,SINAI1,COSAI1,Y1,Y3,X1,SINAI,  LAY 0013
     $     COSAI,CPATH,SH,GAMMA,BEND,RATIO1,S,DSDX1,DBNDX1,PA,PB,TA,    LAY 0014
     $     TB,RHOA,RHOB,DZ,HP,HRHO,DH,H3,R3,H2,R2,SINAI2,SINAI3,RATIO2, LAY 0015
     $     RATIO3,COSAI3,X3,DX,W1,W2,W3,COSAI2,X2,D31,D32,D21,          LAY 0016
     $     DSDX2,DSDX3,DS,DBEND,DSDZ,DHHP,DHRH,H3Z1,                    LAY 0017
     $     DPRARF,DBNDX2,DBNDX3,RNG                                     LAY 0018
      INTEGER N,J,IAMT,K                                                LAY 0019
      LOGICAL LSMALL                                                    LAY 0020
C                                                                       LAY 0021
C     SSI COMMENTS ON DOUBLE PRECISION VARIABLES:                       LAY 0022
C     RFRPTH IS THE OLD COMMON BLOCK IN SINGLE PRECISION.               LAY 0023
C     DPRFRP IS THE SAME COMMON BLOCK IN DOUBLE PRECISION; IT IS NEW.   LAY 0024
C     IN THIS ROUTINE SP IS USED AS A PREFIX TO DENOTE THE              LAY 0025
C     SINGLE PRECISION VARIABLES OF RFRPTH.                             LAY 0026
C     THE FOLLOWING ARE THE EXCEPTIONS:                                 LAY 0027
C     SPRFN  STANDS FOR THE OLD SINGLE PRECISION RFNDXP                 LAY 0028
C     SPPPSU STANDS FOR THE OLD SINGLE PRECISION PPSUM                  LAY 0029
C     SPTPSU STANDS FOR THE OLD SINGLE PRECISION TPSUM                  LAY 0030
C     SPRHOP STANDS FOR THE OLD SINGLE PRECISION RHOPSM                 LAY 0031
C     THE VARIABLES OF THE DOUBLE PRECISION BLOCK DPRFRP HAVE THE       LAY 0032
C     SAME NAMES AS THOSE OF THE ORIGINAL SINGLE PRECISION BLOCK;       LAY 0033
C     THAT IS, WITHOUT ANY PREFIXES.                                    LAY 0034
C                                                                       LAY 0035
      REAL RE,ZMAX                                                      LAY 0036
      INTEGER IMAX,IMOD,IPATH                                           LAY 0037
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             LAY 0038
C                                                                       LAY 0039
C       PI       THE CONSTANT PI                                        LAY 0040
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       LAY 0041
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       LAY 0042
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         LAY 0043
      REAL PI,DEG,BIGNUM,BIGEXP                                         LAY 0044
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                LAY 0045
      COMMON /RFRPTH/ SPZP(LAYDIM+1),SPPP(LAYDIM+1),SPTP(LAYDIM+1),     LAY 0046
     $     SPRFN(LAYDIM+1),SPSP(LAYDIM+1),                              LAY 0047
     $     SPPPSU(LAYDIM+1),SPTPSU(LAYDIM+1),SPRHOP(LAYDIM+1),          LAY 0048
     $     SPDENP(KMAX,LAYDIM+1),SPAMTP(KMAX,LAYDIM+1)                  LAY 0049
      COMMON /DPRFRP/ ZP(LAYDIM+1),PP(LAYDIM+1),TP(LAYDIM+1),           LAY 0050
     $     RFNDXP(LAYDIM+1),SP(LAYDIM+1),PPSUM(LAYDIM+1),               LAY 0051
     $     TPSUM(LAYDIM+1),RHOPSM(LAYDIM+1),DENP(KMAX,LAYDIM+1),        LAY 0052
     $     AMTP(KMAX,LAYDIM+1)                                          LAY 0053
      COMMON /SMALL2/LSMALL                                             LAY 0054
C                                                                       LAY 0055
C                                                                       LAY 0056
C                                                                       LAY 0057
C     CONVENTION                                                        LAY 0058
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         LAY 0059
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            LAY 0060
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             LAY 0061
C                                                                       LAY 0062
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              LAY 0063
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       LAY 0064
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     LAY 0065
C                                                                       LAY 0066
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           LAY 0067
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      LAY 0068
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        LAY 0069
                                                                        LAY 0070
      REAL DENAX(MMOLX),DENBX(MMOLX),HDENX(MMOLX)                       LAY 0071
      COMMON /RFRPTX/ DENPX(MMOLX,LAYDIM+1),AMTPX(MMOLX,LAYDIM+1)       LAY 0072
C                                                                       LAY 0073
C       GCAIR   GAS CONSTANT FOR AIR [MB/(GM CM-3 K)].                  LAY 0074
      DOUBLE PRECISION GCAIR,DELTAS                                     LAY 0075
      DATA EPSILN/1.E-5/,GCAIR/2870.53/,DELTAS/5./                      LAY 0076
C***  INITIALIZE LOOP                                                   LAY 0077
      DPRE=DBLE(RE)                                                     LAY 0078
      N=0                                                               LAY 0079
      Z1=ZP(J)                                                          LAY 0080
      Z2=ZP(J+1)                                                        LAY 0081
      H1=Z1                                                             LAY 0082
      R1=DPRE+H1                                                        LAY 0083
      DHMIN=DBLE(.5)*DELTAS**2/R1                                       LAY 0084
      SINAI1=CPATH/(DPANDX(H1,SH,GAMMA)*R1)                             LAY 0085
      COSAI1=-SQRT((DBLE(1.)-SINAI1)*(DBLE(1.)+SINAI1))                 LAY 0086
      Y1=DBLE(.25)*COSAI1**2                                            LAY 0087
      Y1=2*Y1*(DBLE(1.)+Y1*(DBLE(1.)+Y1))                               LAY 0088
      Y3=DBLE(0.)                                                       LAY 0089
      X1=-R1*COSAI1                                                     LAY 0090
      RATIO1=R1/DPRARF(H1,SH,GAMMA)                                     LAY 0091
      DSDX1=DBLE(1.)/(DBLE(1.)-RATIO1*SINAI1**2)                        LAY 0092
      DBNDX1=DSDX1*SINAI1*RATIO1/R1                                     LAY 0093
      S=DBLE(0.)                                                        LAY 0094
      BEND=DBLE(0.)                                                     LAY 0095
      IF(IAMT.NE.2) THEN                                                LAY 0096
C***  INITIALIZE THE VARIABLES FOR THE CALCULATION OF THE               LAY 0097
C***  ABSORBER AMOUNTS                                                  LAY 0098
      PA=PP(J)                                                          LAY 0099
      PB=PP(J+1)                                                        LAY 0100
      TA=TP(J)                                                          LAY 0101
      TB=TP(J+1)                                                        LAY 0102
      RHOA=PA/(GCAIR*TA)                                                LAY 0103
      RHOB=PB/(GCAIR*TB)                                                LAY 0104
      DZ=ZP(J+1)-ZP(J)                                                  LAY 0105
      IF (LSMALL .AND. ABS(DZ) .LT. 1E-15) THEN                         LAY 0106
      DZ=1.D-20                                                         LAY 0107
      HP=1.D-20                                                         LAY 0108
      ELSEIF( .NOT. LSMALL .AND. (ABS(PB/PA) .GE. 0.99999)) THEN        LAY 0109
         HP=1.D-20                                                      LAY 0110
      ELSE                                                              LAY 0111
         HP=-DZ/LOG(PB/PA)                                              LAY 0112
      ENDIF                                                             LAY 0113
      IF(ABS(RHOB/RHOA-DBLE(1.)).GE.EPSILN)THEN                         LAY 0114
      HRHO=-DZ/LOG(RHOB/RHOA)                                           LAY 0115
      ELSE                                                              LAY 0116
         HRHO=1.D30                                                     LAY 0117
      ENDIF                                                             LAY 0118
      DO 105 K=1,KMAX+NSPECX                                            LAY 0119
         IF ( K .GT. KMAX) THEN                                         LAY 0120
            KX=K-KMAX                                                   LAY 0121
            DENAX(KX)=DENPX(KX,J)                                       LAY 0122
            DENBX(KX)=DENPX(KX,J+1)                                     LAY 0123
            IF(DENAX(KX).LE.0. .OR. DENBX(KX).LE.0.)THEN                LAY 0124
C***           USE LINEAR INTERPOLATION                                 LAY 0125
               HDENX(KX)=0.                                             LAY 0126
               GO TO 105                                                LAY 0127
            ENDIF                                                       LAY 0128
            IF(ABS(1.-DENAX(KX)/DENBX(KX)).GT.EPSILN)THEN               LAY 0129
C***           USE EXPONENTIAL INTERPOLATION                            LAY 0130
               HDENX(KX)=-REAL(DZ)/LOG(DENBX(KX)/DENAX(KX))             LAY 0131
            ELSE                                                        LAY 0132
C***           USE LINEAR INTERPOLATION                                 LAY 0133
               HDENX(KX)=0.                                             LAY 0134
            ENDIF                                                       LAY 0135
            GO TO 105                                                   LAY 0136
         ENDIF                                                          LAY 0137
C                                                                       LAY 0138
         DENA(K)=DENP(K,J)                                              LAY 0139
         DENB(K)=DENP(K,J+1)                                            LAY 0140
         HDEN(K)=DBLE(0.)                                               LAY 0141
         IF(DENA(K).GT.0. .AND. DENB(K).GT.0. .AND.                     LAY 0142
     1     K.NE.16 .AND. K.NE.66 .AND. K.NE.67 .AND.                    LAY 0143
     2     ABS(DENB(K)-DENA(K)).GT.DBLE(EPSILN)*DENB(K))                LAY 0144
     3     HDEN(K)=DZ/LOG(DENA(K)/DENB(K))                              LAY 0145
 105  CONTINUE                                                          LAY 0146
      ENDIF                                                             LAY 0147
C***                                                                    LAY 0148
C***  LOOP THROUGH PATH                                                 LAY 0149
C***  INTEGRATE PATH QUANTITIES USING QUADRATIC INTEGRATION WITH        LAY 0150
C***  UNEQUALLY SPACED POINTS                                           LAY 0151
C***                                                                    LAY 0152
 115  CONTINUE                                                          LAY 0153
      N=N+1                                                             LAY 0154
      DH=-DELTAS*COSAI1                                                 LAY 0155
      IF(DH.LT.DHMIN) DH=DHMIN                                          LAY 0156
      H3=H1+DH                                                          LAY 0157
      IF(H3.GT.Z2) H3=Z2                                                LAY 0158
      DH=H3-H1                                                          LAY 0159
      R3=DPRE+H3                                                        LAY 0160
      H2=H1+DH/2                                                        LAY 0161
      R2=DPRE+H2                                                        LAY 0162
      SINAI2=CPATH/(DPANDX(H2,SH,GAMMA)*R2)                             LAY 0163
      SINAI3=CPATH/(DPANDX(H3,SH,GAMMA)*R3)                             LAY 0164
      RATIO2=R2/DPRARF(H2,SH,GAMMA)                                     LAY 0165
      RATIO3=R3/DPRARF(H3,SH,GAMMA)                                     LAY 0166
      IF((DBLE(1.)-SINAI2).LE.EPSILN)THEN                               LAY 0167
C***  NEAR A TANGENT HEIGHT, COSAI=-SQRT(1-SINAI**2) LOSES              LAY 0168
C***  PRECISION.  USE THE FOLLOWING ALGORITHM TO GET COSAI.             LAY 0169
      Y3=Y1+(SINAI1*(DBLE(1.)-RATIO1)/R1                                LAY 0170
     1  +DBLE(4.)*SINAI2*(DBLE(1.)-RATIO2)/R2                           LAY 0171
     2  +SINAI3*(DBLE(1.)-RATIO3)/R3)*DH/DBLE(6.)                       LAY 0172
      IF(Y3.GE.0.)THEN                                                  LAY 0173
          COSAI3=-SQRT(Y3*(DBLE(2.)-Y3))                                LAY 0174
      ELSE                                                              LAY 0175
          Y3=DBLE(1.)-SINAI3                                            LAY 0176
          COSAI3=-SQRT(Y3*(DBLE(1.)+SINAI3))                            LAY 0177
      ENDIF                                                             LAY 0178
      X3=-R3*COSAI3                                                     LAY 0179
      DX=X3-X1                                                          LAY 0180
      W1=DBLE(.5)*DX                                                    LAY 0181
      W2=DBLE(.0)                                                       LAY 0182
      W3=DBLE(.5)*DX                                                    LAY 0183
      GO TO 118                                                         LAY 0184
C***                                                                    LAY 0185
      ENDIF                                                             LAY 0186
      COSAI2=-SQRT((DBLE(1.)-SINAI2)*(DBLE(1.)+SINAI2))                 LAY 0187
      COSAI3=-SQRT((DBLE(1.)-SINAI3)*(DBLE(1.)+SINAI3))                 LAY 0188
      X2=-R2*COSAI2                                                     LAY 0189
      X3=-R3*COSAI3                                                     LAY 0190
C***  CALCULATE WEIGHTS                                                 LAY 0191
      D31=X3-X1                                                         LAY 0192
      D32=X3-X2                                                         LAY 0193
      D21=X2-X1                                                         LAY 0194
      IF(D32.EQ.0. .OR. D21.EQ.0.)THEN                                  LAY 0195
      W1=DBLE(.5)*D31                                                   LAY 0196
      W2=DBLE(.0)                                                       LAY 0197
      W3=DBLE(.5)*D31                                                   LAY 0198
      ELSE                                                              LAY 0199
         W1=(2-D32/D21)*D31/6                                           LAY 0200
         W2=D31**3/(D32*D21*6)                                          LAY 0201
         W3=(2-D21/D32)*D31/6                                           LAY 0202
      ENDIF                                                             LAY 0203
C***                                                                    LAY 0204
 118  CONTINUE                                                          LAY 0205
      DSDX2=DBLE(1.)/(DBLE(1.)-RATIO2*SINAI2**2)                        LAY 0206
      DSDX3=DBLE(1.)/(DBLE(1.)-RATIO3*SINAI3**2)                        LAY 0207
      DBNDX2=DSDX2*SINAI2*RATIO2/R2                                     LAY 0208
      DBNDX3=DSDX3*SINAI3*RATIO3/R3                                     LAY 0209
C***  INTEGRATE                                                         LAY 0210
      DS=W1*DSDX1+W2*DSDX2+W3*DSDX3                                     LAY 0211
      IF(LSMALL)DS=RNG                                                  LAY 0212
      S=S+DS                                                            LAY 0213
      DBEND=W1*DBNDX1+W2*DBNDX2+W3*DBNDX3                               LAY 0214
      IF(LSMALL)DBEND=DBLE(0.)                                          LAY 0215
      BEND=BEND+DBEND                                                   LAY 0216
      IF(IAMT.NE.2)  THEN                                               LAY 0217
C***  CALCULATE AMOUNTS                                                 LAY 0218
      IF (LSMALL .AND. ABS(DZ) .LT. 1E-15 ) THEN                        LAY 0219
      DSDZ=DBLE(BIGEXP+1.)                                              LAY 0220
      DHHP=0                                                            LAY 0221
      ELSE                                                              LAY 0222
         DSDZ=DS/DH                                                     LAY 0223
         DHHP=DH/HP                                                     LAY 0224
      ENDIF                                                             LAY 0225
      IF(DHHP .LE. BIGEXP) THEN                                         LAY 0226
      PB=PA*EXP(-DHHP )                                                 LAY 0227
      ELSE                                                              LAY 0228
         PB=DBLE(0.)                                                    LAY 0229
      ENDIF                                                             LAY 0230
      DHRH=DH/HRHO                                                      LAY 0231
      IF(DHRH .LE. BIGEXP) THEN                                         LAY 0232
      RHOB=RHOA*EXP(-DHRH   )                                           LAY 0233
      ELSE                                                              LAY 0234
         RHOB=DBLE(0.)                                                  LAY 0235
      ENDIF                                                             LAY 0236
      IF((DH/HRHO).GE.EPSILN)  THEN                                     LAY 0237
      PPSUM(J)=PPSUM(J)+                                                LAY 0238
     $     DSDZ*(HP/(DBLE(1.)+HP/HRHO))*(PA*RHOA-PB*RHOB)               LAY 0239
      TPSUM(J)=TPSUM(J)+DSDZ*HP*(PA-PB)/GCAIR                           LAY 0240
      RHOPSM(J)=RHOPSM(J)+DSDZ*HRHO*(RHOA-RHOB)                         LAY 0241
      SPPPSU(J)=REAL(PPSUM(J))                                          LAY 0242
      SPTPSU(J)=REAL(TPSUM(J))                                          LAY 0243
      SPRHOP(J)=REAL(RHOPSM(J))                                         LAY 0244
      ELSE                                                              LAY 0245
         PPSUM(J)=PPSUM(J)+DS*(PA*RHOA+PB*RHOB)/2                       LAY 0246
         TPSUM(J)=TPSUM(J)+DS*(PA+PB)/GCAIR/2                           LAY 0247
         RHOPSM(J)=RHOPSM(J)+DS*(RHOA+RHOB)/2                           LAY 0248
         SPPPSU(J)=REAL(PPSUM(J))                                       LAY 0249
         SPTPSU(J)=REAL(TPSUM(J))                                       LAY 0250
         SPRHOP(J)=REAL(RHOPSM(J))                                      LAY 0251
      ENDIF                                                             LAY 0252
      DO 140 K=1,KMAX+NSPECX                                            LAY 0253
C                                                                       LAY 0254
         IF ( K .GT. KMAX) THEN                                         LAY 0255
            KX=K -KMAX                                                  LAY 0256
            IF(ABS(HDENX(KX)).EQ.0.)THEN                                LAY 0257
C***           LINEAR INTERPOLATION                                     LAY 0258
               DENBX(KX)=DENPX(KX,J)+                                   LAY 0259
     1           (DENPX(KX,J+1)-DENPX(KX,J))*REAL((H3-Z1)/DZ)           LAY 0260
               AMTPX(KX,J)=AMTPX(KX,J)+.5*(DENAX(KX)+DENBX(KX))*REAL(DS)LAY 0261
            ELSEIF(REAL(DH)/HDENX(KX).GE.EPSILN)THEN                    LAY 0262
C***           EXPONENTIAL INTERPOLATION                                LAY 0263
               H3Z1= (H3-Z1)/DBLE(HDENX(KX))                            LAY 0264
               IF(H3Z1.LE.BIGEXP)THEN                                   LAY 0265
                  DENBX(KX)=DENPX(KX,J)*EXP(-REAL(H3Z1))                LAY 0266
               ELSE                                                     LAY 0267
                  DENBX(KX)=0.                                          LAY 0268
               ENDIF                                                    LAY 0269
               AMTPX(KX,J)=AMTPX(KX,J)+                                 LAY 0270
     1           REAL(DSDZ)*HDENX(KX)*(DENAX(KX)-DENBX(KX))             LAY 0271
            ELSE                                                        LAY 0272
C***           LINEAR INTERPOLATION                                     LAY 0273
               DENBX(KX)=DENPX(KX,J)+                                   LAY 0274
     1           (DENPX(KX,J+1)-DENPX(KX,J))*REAL((H3-Z1)/DZ)           LAY 0275
               AMTPX(KX,J)=AMTPX(KX,J)+.5*(DENAX(KX)+DENBX(KX))*REAL(DS)LAY 0276
            ENDIF                                                       LAY 0277
            GO TO 140                                                   LAY 0278
         ENDIF                                                          LAY 0279
C                                                                       LAY 0280
         IF(ABS(HDEN(K)).EQ.0.)THEN                                     LAY 0281
C***        LINEAR INTERPOLATION                                        LAY 0282
            DENB(K)=DENP(K,J)+(DENP(K,J+1)-DENP(K,J))*(H3-Z1)/DZ        LAY 0283
            AMTP(K,J)=AMTP(K,J)+DBLE(.5)*(DENA(K)+DENB(K))*DS           LAY 0284
         ELSEIF(DH/HDEN(K).GE.EPSILN)THEN                               LAY 0285
C***        EXPONENTIAL INTERPOLATION                                   LAY 0286
            H3Z1= (H3-Z1)/HDEN(K)                                       LAY 0287
            IF(H3Z1 .LE.BIGEXP) THEN                                    LAY 0288
               DENB(K)=DENP(K,J)*EXP(-H3Z1)                             LAY 0289
            ELSE                                                        LAY 0290
               DENB(K)=DBLE(0.)                                         LAY 0291
            ENDIF                                                       LAY 0292
            AMTP(K,J)=AMTP(K,J)+DSDZ*HDEN(K)*(DENA(K)-DENB(K))          LAY 0293
         ELSE                                                           LAY 0294
C***        LINEAR INTERPOLATION                                        LAY 0295
            DENB(K)=DENP(K,J)+(DENP(K,J+1)-DENP(K,J))*(H3-Z1)/DZ        LAY 0296
            AMTP(K,J)=AMTP(K,J)+DBLE(.5)*(DENA(K)+DENB(K))*DS           LAY 0297
         ENDIF                                                          LAY 0298
         SPAMTP(K,J)=REAL(AMTP(K,J))                                    LAY 0299
 140     CONTINUE                                                       LAY 0300
         PA=PB                                                          LAY 0301
         RHOA=RHOB                                                      LAY 0302
         DO 145 K=1,KMAX+NSPECX                                         LAY 0303
            IF (K .LE. KMAX) THEN                                       LAY 0304
               DENA(K)=DENB(K)                                          LAY 0305
            ELSE                                                        LAY 0306
               DENAX(K-KMAX)=DENBX(K-KMAX)                              LAY 0307
            ENDIF                                                       LAY 0308
 145     CONTINUE                                                       LAY 0309
      ENDIF                                                             LAY 0310
      IF(H3.LT.Z2) THEN                                                 LAY 0311
         H1=H3                                                          LAY 0312
         R1=R3                                                          LAY 0313
         SINAI1=SINAI3                                                  LAY 0314
         RATIO1=RATIO3                                                  LAY 0315
         Y1=Y3                                                          LAY 0316
         COSAI1=COSAI3                                                  LAY 0317
         X1=X3                                                          LAY 0318
         DSDX1=DSDX3                                                    LAY 0319
         DBNDX1=DBNDX3                                                  LAY 0320
         IF (.NOT. LSMALL) GO TO 115                                    LAY 0321
      ENDIF                                                             LAY 0322
      SINAI=SINAI3                                                      LAY 0323
      COSAI=COSAI3                                                      LAY 0324
      SP(J)=S                                                           LAY 0325
      SPSP(J)=REAL(S)                                                   LAY 0326
      RETURN                                                            LAY 0327
      END                                                               LAY 0328
