      SUBROUTINE BMOD(IK,IKMX,IPATH,IV,MSOFF,MXFREQ)                    BMO 0001
C                                                                       BMO 0002
C     THIS ROUTINE RETURNS THE TRANSMITTANCE AT A SPECTRAL              BMO 0003
C     RESOLUTION OF 1 CM-1 FOR THE CFC'S AND THE "NSPC" SPECIES.        BMO 0004
C           K = ( 1, 13)      H2O LINE (CENTERS, TAILS)                 BMO 0005
C             = ( 2, 14)      CO2 LINE (CENTERS, TAILS)                 BMO 0006
C             = ( 3, 15)       O3 LINE (CENTERS, TAILS)                 BMO 0007
C             = ( 4, 16)      N2O LINE (CENTERS, TAILS)                 BMO 0008
C             = ( 5, 17)       CO LINE (CENTERS, TAILS)                 BMO 0009
C             = ( 6, 18)      CH4 LINE (CENTERS, TAILS)                 BMO 0010
C             = ( 7, 19)       O2 LINE (CENTERS, TAILS)                 BMO 0011
C             = ( 8, 20)       NO LINE (CENTERS, TAILS)                 BMO 0012
C             = ( 9, 21)      SO2 LINE (CENTERS, TAILS)                 BMO 0013
C             = (10, 22)      NO2 LINE (CENTERS, TAILS)                 BMO 0014
C             = (11, 23)      NH3 LINE (CENTERS, TAILS)                 BMO 0015
C             = (12, 24)     HNO3 LINE (CENTERS, TAILS)                 BMO 0016
C                                                                       BMO 0017
C     DECLARE INPUTS                                                    BMO 0018
      INTEGER IK,IKMX,IPATH,IV,MSOFF,MXFREQ                             BMO 0019
C                                                                       BMO 0020
C     INCLUDE PARAMETERS                                                BMO 0021
      INCLUDE 'PARAM.LST'                                               BMO 0022
C                                                                       BMO 0023
C     LIST COMMONS                                                      BMO 0024
      INTEGER IBINX,IMOLX,IALFX                                         BMO 0025
      REAL SDZX,ODZX                                                    BMO 0026
      COMMON/BMDCMX/SDZX(MXTEMP),ODZX(MXTEMP),IBINX,IMOLX,IALFX         BMO 0027
      REAL TXX,WX,WPATHX                                                BMO 0028
      COMMON/NONAME/TXX(MMOLX),WX(MMOLX),WPATHX(LAYTHR,MMOLX)           BMO 0029
      REAL WPTHSX,TBBYSX                                                BMO 0030
      COMMON/SOLSX/WPTHSX(LAYTHR,MMOLX),TBBYSX(LAYTHR,MMOLX)            BMO 0031
      INTEGER KPOINT                                                    BMO 0032
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     BMO 0033
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   BMO 0034
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   BMO 0035
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     BMO 0036
      INTEGER JTURN,LJ                                                  BMO 0037
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       BMO 0038
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  BMO 0039
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      BMO 0040
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    BMO 0041
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     BMO 0042
      REAL TBOUND,SALB                                                  BMO 0043
      LOGICAL MODTRN                                                    BMO 0044
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,                      BMO 0045
     1  NOPRNT,TBOUND,SALB,MODTRN                                       BMO 0046
      INTEGER IBNDWD,IP,ITB,NTEMP,IBIN,IMOL,IALF,JJ,JJS                 BMO 0047
      REAL TBAND,SDZ,ODZ,SD,OD,ALF0,T5,PTM75,FF,T5S,                    BMO 0048
     1  PTM75S,FFS,DOPFAC,DOP0,DOPSUM,COLSUM,ODSUM,SDSUM                BMO 0049
      COMMON/BMDCOM/IBNDWD,IP,ITB,NTEMP,TBAND(MXTEMP),SDZ(MXTEMP),      BMO 0050
     1  ODZ(MXTEMP),IBIN,IMOL,IALF,SD(MXTEMP,MMOLT2),OD(MXTEMP,MMOLT),  BMO 0051
     2  ALF0(MMOLT),T5(LAYTHR),PTM75(LAYTHR),JJ(LAYTHR),FF(LAYTHR),     BMO 0052
     3  T5S(LAYTHR,NSPC),PTM75S(LAYTHR,NSPC),JJS(LAYTHR,MMOLT),         BMO 0053
     4  FFS(LAYTHR,MMOLT),DOPFAC(MMOLT),DOP0(MMOLT),SDSUM(MMOLT2),      BMO 0054
     5  ODSUM(MMOLT),DOPSUM(MMOLT),COLSUM(MMOLT)                        BMO 0055
C                                                                       BMO 0056
C     DECLARE LOCAL FUNCTIONS                                           BMO 0057
      REAL BMTRAN,BMTRN                                                 BMO 0058
C                                                                       BMO 0059
C     DECLARE LOCAL VARIABLES                                           BMO 0060
      INTEGER ITBX,IPATH0,K,KX,IML,IT,I,LMAX,LS,                        BMO 0061
     1  NSPED,KNEW,KSD,KP,KTAIL,L,MSOFFL,J,JM1                          BMO 0062
      LOGICAL MOLTRN                                                    BMO 0063
      REAL V,COLO3,DV,F,STORE,ABSM,DINV,TAIL,                           BMO 0064
     1  TDEPTH,DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2,TRANSM                    BMO 0065
C                                                                       BMO 0066
C     SAVE LOCAL VARIABLES                                              BMO 0067
      SAVE MOLTRN,COLO3                                                 BMO 0068
C                                                                       BMO 0069
C     DOPFAC (UNITLESS) EQUALS SQRT(2 LN2 R T / M)/C WHERE T IS THE     BMO 0070
C     STANDARD TEMPERATURE (273.15K) AND M IS MOLECULAR WEIGHT.         BMO 0071
      DATA ITBX/31/,DV/1./,IPATH0/0/                                    BMO 0072
C                                                                       BMO 0073
C     IF IV IS BEYOND BAND MODEL TAPE RANGE, SET TRANSMITTANCES TO 1    BMO 0074
      IF(IV.GT.MXFREQ)THEN                                              BMO 0075
          DO 10 K=1,NSPC                                                BMO 0076
              TX(KPOINT(K))=1.                                          BMO 0077
   10     CONTINUE                                                      BMO 0078
C                                                                       BMO 0079
          DO 20 KX=1,NSPECX                                             BMO 0080
              TXX(KX)=1.                                                BMO 0081
   20     CONTINUE                                                      BMO 0082
          RETURN                                                        BMO 0083
      ENDIF                                                             BMO 0084
C                                                                       BMO 0085
C     IK=0 IS THE INITIAL CALL FOR EACH WAVENUMBER AND IS MADE PRIOR    BMO 0086
C     TO THE LOOP OVER LAYERS.  IK=-1 IS THE INITIAL CALL TO BMOD       BMO 0087
C     AFTER THE MULTIPLE SCATTERING LOOP IS COMPLETE                    BMO 0088
      IF(IK.EQ.0)THEN                                                   BMO 0089
C                                                                       BMO 0090
C         IS THERE MOLECULAR DATA FOR FREQUENCY IV?                     BMO 0091
          IF(IBIN.NE.IV)THEN                                            BMO 0092
              MOLTRN=.FALSE.                                            BMO 0093
              RETURN                                                    BMO 0094
          ENDIF                                                         BMO 0095
          MOLTRN=.TRUE.                                                 BMO 0096
C                                                                       BMO 0097
C         FOR EACH SPECIES, SET THE 'NO DATA' INDICATOR                 BMO 0098
          DO 30 IML=1,2*NSPECT                                          BMO 0099
   30     SD(1,IML)=0.                                                  BMO 0100
C                                                                       BMO 0101
C         LOAD MOLECULAR DATA FOR FREQUENCY IV                          BMO 0102
   40     CONTINUE                                                      BMO 0103
          IF(IBIN.EQ.IV)THEN                                            BMO 0104
              CALL BMLOAD                                               BMO 0105
              IF(IV.GE.MXFREQ)GOTO50                                    BMO 0106
              IP=IP+1                                                   BMO 0107
              READ(ITB,REC=IP)                                          BMO 0108
     1          IBIN,IMOL,(SDZ(I),I=1,MXTEMP),IALF,(ODZ(I),I=1,MXTEMP)  BMO 0109
              GOTO40                                                    BMO 0110
          ENDIF                                                         BMO 0111
C                                                                       BMO 0112
C         LOAD CFC DATA FOR FREQUENCY IV                                BMO 0113
   50     CONTINUE                                                      BMO 0114
          IF(IBINX.EQ.IV)THEN                                           BMO 0115
              CALL BMXLD                                                BMO 0116
              READ(ITBX,'(I6,I5,1P11E11.3)',END=60)                     BMO 0117
     1          I,IMOLX,(SDZX(IT),IT=1,NTEMP)                           BMO 0118
              IBINX=I                                                   BMO 0119
              READ(ITBX,'(6X,I5,1P11E11.3)')IALFX,(ODZX(IT),IT=1,NTEMP) BMO 0120
              GOTO50                                                    BMO 0121
          ENDIF                                                         BMO 0122
   60     CONTINUE                                                      BMO 0123
          DV=IBNDWD                                                     BMO 0124
C                                                                       BMO 0125
C         ZERO LAYER LOOP QUANTITIES; DEFINE STANDARD DOPPLER WIDTH     BMO 0126
          V=FLOAT(IV)                                                   BMO 0127
          IF(IV.EQ.0)V=.25*IBNDWD                                       BMO 0128
          DO 70 K=1,NSPECT                                              BMO 0129
              DOP0(K)=V*DOPFAC(K)                                       BMO 0130
              COLSUM(K)=0.                                              BMO 0131
              DOPSUM(K)=0.                                              BMO 0132
              ODSUM(K)=0.                                               BMO 0133
              SDSUM(K)=0.                                               BMO 0134
              SDSUM(K+NSPECT)=0.                                        BMO 0135
   70         CONTINUE                                                  BMO 0136
              COLO3=0.                                                  BMO 0137
           RETURN                                                       BMO 0138
      ELSEIF(IK.EQ.-1)THEN                                              BMO 0139
          DO 80 K=1,NSPECT                                              BMO 0140
              COLSUM(K)=0.                                              BMO 0141
              DOPSUM(K)=0.                                              BMO 0142
              ODSUM(K)=0.                                               BMO 0143
              SDSUM(K)=0.                                               BMO 0144
              SDSUM(K+NSPECT)=0.                                        BMO 0145
   80         CONTINUE                                                  BMO 0146
              COLO3=0.                                                  BMO 0147
          RETURN                                                        BMO 0148
      ENDIF                                                             BMO 0149
      IF(.NOT.MOLTRN)THEN                                               BMO 0150
          DO 90 K=1,NSPC                                                BMO 0151
              TX(KPOINT(K))=1.                                          BMO 0152
   90     CONTINUE                                                      BMO 0153
          RETURN                                                        BMO 0154
      ENDIF                                                             BMO 0155
C                                                                       BMO 0156
C     START CALCULATION OF MOLECULAR TRANSMITTANCE                      BMO 0157
      IF(IEMSCT.EQ.0 .OR. IEMSCT.EQ.3)THEN                              BMO 0158
C                                                                       BMO 0159
C         LOOP OVER ALL LAYERS FOR TRANSMITTANCE CALCULATIONS           BMO 0160
          LMAX=IKMX                                                     BMO 0161
      ELSEIF(IEMSCT.EQ.1)THEN                                           BMO 0162
C                                                                       BMO 0163
C         LOOP OVER SINGLE LAYER FOR RADIANCE CALCULATION WITHOUT SOLAR BMO 0164
          LMAX=IK                                                       BMO 0165
      ELSE                                                              BMO 0166
          GOTO(100,110,120),IPATH                                       BMO 0167
C                                                                       BMO 0168
C         SKIP LAYER LOOP COMPLETELY FOR FIRST SOLAR PATH WITH IPATH=1  BMO 0169
  100     LMAX=0                                                        BMO 0170
          LS=1                                                          BMO 0171
          GOTO130                                                       BMO 0172
C                                                                       BMO 0173
C         LOOP OVER SINGLE LAYER ONLY FOR "L PATH" WITH IPATH=2         BMO 0174
  110     LMAX=IK                                                       BMO 0175
          IF(MSOFF.GT.0)LMAX=0                                          BMO 0176
          LS=IK+1                                                       BMO 0177
          GOTO130                                                       BMO 0178
C                                                                       BMO 0179
C         LOOP OVER SINGLE LAYER IF NOT PERFORMED WHEN IPATH EQUALED 2  BMO 0180
  120     LMAX=0                                                        BMO 0181
          IF(MSOFF.GT.0 .OR. IPATH0.NE.2)LMAX=IK                        BMO 0182
  130     IPATH0=IPATH                                                  BMO 0183
      ENDIF                                                             BMO 0184
C                                                                       BMO 0185
C     START SPECIES LOOP                                                BMO 0186
      NSPED=NSPC+NSPC                                                   BMO 0187
      DO 160 KNEW=1,NSPECT                                              BMO 0188
          K=KNEW                                                        BMO 0189
          IF(KNEW.LE.NSPC)THEN                                          BMO 0190
             KSD=K                                                      BMO 0191
             KP=KPOINT(K)                                               BMO 0192
             KTAIL=K+NSPC                                               BMO 0193
          ELSE                                                          BMO 0194
             KP=KNEW-NSPC                                               BMO 0195
             KSD=KP+NSPED                                               BMO 0196
             KTAIL=KSD+NSPECX                                           BMO 0197
          ENDIF                                                         BMO 0198
C                                                                       BMO 0199
C         CHECK IF LINE CENTER CONTRIBUTES TO ABSORPTION                BMO 0200
          IF(SD(1,KSD).GT.0.)THEN                                       BMO 0201
C                                                                       BMO 0202
C             LOOP OVER LAYERS                                          BMO 0203
              DO 140 L=IK,LMAX                                          BMO 0204
C                                                                       BMO 0205
C                 INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE    BMO 0206
                  MSOFFL=MSOFF+L                                        BMO 0207
                  J=JJ(MSOFFL)                                          BMO 0208
                  F=FF(MSOFFL)                                          BMO 0209
                  JM1=J-1                                               BMO 0210
                  STORE=SD(J,KSD)                                       BMO 0211
                  ABSM=STORE+F*(SD(JM1,KSD)-STORE)                      BMO 0212
                  STORE=OD(J,K)                                         BMO 0213
                  DINV=STORE+F*(OD(JM1,K)-STORE)                        BMO 0214
                  STORE=SD(J,KTAIL)                                     BMO 0215
                  TAIL=STORE+F*(SD(JM1,KTAIL)-STORE)                    BMO 0216
C                                                                       BMO 0217
C                 PERFORM CURTIS-GODSON SUMS                            BMO 0218
                  IF(KNEW.LE.NSPC)STORE=WPATH(MSOFFL,KP)                BMO 0219
                  IF(KNEW.GT.NSPC)STORE=WPATHX(MSOFFL,KP)               BMO 0220
                  SDSUM(KTAIL)=SDSUM(KTAIL)+STORE*TAIL*PATM(MSOFFL)     BMO 0221
                  STORE=ABSM*STORE                                      BMO 0222
                  SDSUM(KSD)=SDSUM(KSD)+STORE                           BMO 0223
                  STORE=DINV*STORE                                      BMO 0224
                  ODSUM(K)=ODSUM(K)+STORE                               BMO 0225
                  DOPSUM(K)=DOPSUM(K)+STORE*T5(MSOFFL)                  BMO 0226
                  STORE=STORE*PTM75(MSOFFL)                             BMO 0227
                  COLSUM(K)=COLSUM(K)+STORE                             BMO 0228
                  IF(K.EQ.3)COLO3=COLO3+STORE*DINV*PTM75(MSOFFL)        BMO 0229
  140         CONTINUE                                                  BMO 0230
              TDEPTH=SDSUM(KTAIL)                                       BMO 0231
              DEPTH=SDSUM(KSD)                                          BMO 0232
              ODBAR=ODSUM(K)                                            BMO 0233
              ADBAR=DOPSUM(K)                                           BMO 0234
              ACBAR=COLSUM(K)                                           BMO 0235
              ACBAR2=0.                                                 BMO 0236
              IF(K.EQ.3)ACBAR2=COLO3                                    BMO 0237
C                                                                       BMO 0238
C             IF SOLAR PATH, CALCULATE ADDITIONAL LAYER                 BMO 0239
              IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN                      BMO 0240
                  MSOFFL=MSOFF+LS                                       BMO 0241
                  J=JJS(MSOFFL,K)                                       BMO 0242
                  F=FFS(MSOFFL,K)                                       BMO 0243
                  JM1=J-1                                               BMO 0244
                  STORE=SD(J,KSD)                                       BMO 0245
                  ABSM=STORE+F*(SD(JM1,KSD)-STORE)                      BMO 0246
                  STORE=OD(J,K)                                         BMO 0247
                  DINV=STORE+F*(OD(JM1,K)-STORE)                        BMO 0248
                  STORE=SD(J,KTAIL)                                     BMO 0249
                  TAIL=STORE+F*(SD(JM1,KTAIL)-STORE)                    BMO 0250
                  IF(KNEW.LE.NSPC) STORE=WPATHS(MSOFFL,KP)              BMO 0251
                  IF(KNEW.GT.NSPC) STORE=WPTHSX(MSOFFL,KP)              BMO 0252
                  IF(MSOFF.EQ.0)THEN                                    BMO 0253
C                                                                       BMO 0254
C                     L-SHAPED PATH                                     BMO 0255
                      TDEPTH=TDEPTH+STORE*TAIL*PATMS(MSOFFL,K)          BMO 0256
                      STORE=ABSM*STORE                                  BMO 0257
                      DEPTH=DEPTH+STORE                                 BMO 0258
                      STORE=DINV*STORE                                  BMO 0259
                      ODBAR=ODBAR+STORE                                 BMO 0260
                      ADBAR=ADBAR+STORE*T5S(MSOFFL,K)                   BMO 0261
                      STORE=STORE*PTM75S(MSOFFL,K)                      BMO 0262
                      ACBAR=ACBAR+STORE                                 BMO 0263
                      IF(K.EQ.3)                                        BMO 0264
     1                  ACBAR2=ACBAR2+STORE*DINV*PTM75S(MSOFFL,K)       BMO 0265
                  ELSE                                                  BMO 0266
C                                                                       BMO 0267
C                     SOLAR PATH ONLY                                   BMO 0268
                      TDEPTH=STORE*TAIL*PATMS(MSOFFL,K)                 BMO 0269
                      DEPTH=ABSM*STORE                                  BMO 0270
                      ODBAR=DINV*DEPTH                                  BMO 0271
                      ADBAR=ODBAR*T5S(MSOFFL,K)                         BMO 0272
                      ACBAR=ODBAR*PTM75S(MSOFFL,K)                      BMO 0273
                      IF(K.EQ.3)ACBAR2=ACBAR*DINV*PTM75S(MSOFFL,K)      BMO 0274
                  ENDIF                                                 BMO 0275
              ENDIF                                                     BMO 0276
C                                                                       BMO 0277
C             CHECK FOR WEAK LINE                                       BMO 0278
              IF(DEPTH.LT.0.001)THEN                                    BMO 0279
                  TRANSM=1.-DEPTH                                       BMO 0280
              ELSE                                                      BMO 0281
C                                                                       BMO 0282
C                 CALCULATE EQUIVALENT WIDTH TRANSMITTANCE              BMO 0283
                  ODBAR=ODBAR/DEPTH                                     BMO 0284
                  ADBAR=DOP0(K)*ADBAR/DEPTH                             BMO 0285
                  ACBAR=ALF0(K)*ACBAR/DEPTH                             BMO 0286
                  IF(ACBAR2.NE.0.)THEN                                  BMO 0287
                      ACBAR2=ACBAR2/DEPTH*ALF0(K)**2                    BMO 0288
                      TRANSM=BMTRAN(DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2,DV)  BMO 0289
                  ELSE                                                  BMO 0290
                      TRANSM=BMTRN(DEPTH,ODBAR,ADBAR,ACBAR,DV)          BMO 0291
                  ENDIF                                                 BMO 0292
              ENDIF                                                     BMO 0293
C                                                                       BMO 0294
C             ADD LINE TAIL CONTRIBUTIONS.                              BMO 0295
              IF(TRANSM.GT.0.)TRANSM=TRANSM*EXP(-TDEPTH)                BMO 0296
              IF(KNEW.LE.NSPC)THEN                                      BMO 0297
                  TX(KP)=TRANSM                                         BMO 0298
              ELSE                                                      BMO 0299
                  TXX(KP)=TRANSM                                        BMO 0300
              ENDIF                                                     BMO 0301
C                                                                       BMO 0302
C         CHECK IF LINE TAILS CONTRIBUTE TO ABSORPTION                  BMO 0303
          ELSEIF(SD(1,KTAIL).GT.0.)THEN                                 BMO 0304
C                                                                       BMO 0305
C             LOOP OVER LAYERS                                          BMO 0306
              DO 150 L=IK,LMAX                                          BMO 0307
C                                                                       BMO 0308
C                 INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE    BMO 0309
                  MSOFFL=MSOFF+L                                        BMO 0310
                  J=JJ(MSOFFL)                                          BMO 0311
                  STORE=SD(J,KTAIL)                                     BMO 0312
                  TAIL=STORE+FF(MSOFFL)*(SD(J-1,KTAIL)-STORE)           BMO 0313
C                                                                       BMO 0314
C                 PERFORM CURTIS-GODSON SUMS                            BMO 0315
                  IF(KNEW.LE.NSPC)                                      BMO 0316
     1                STORE=WPATH(MSOFFL,KP)*TAIL*PATM(MSOFFL)          BMO 0317
                  IF(KNEW.GT.NSPC)                                      BMO 0318
     1                STORE=WPATHX(MSOFFL,KP)*TAIL                      BMO 0319
                  SDSUM(KTAIL)=SDSUM(KTAIL)+STORE                       BMO 0320
  150         CONTINUE                                                  BMO 0321
              TDEPTH=SDSUM(KTAIL)                                       BMO 0322
C                                                                       BMO 0323
C             IF SOLAR PATH, CALCULATE ADDITIONAL LAYER.                BMO 0324
              IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN                      BMO 0325
                  MSOFFL=MSOFF+LS                                       BMO 0326
                  J=JJS(MSOFFL,K)                                       BMO 0327
                  STORE=SD(J,KTAIL)                                     BMO 0328
                  TAIL=STORE+FFS(MSOFFL,K)*(SD(J-1,KTAIL)-STORE)        BMO 0329
                  IF(MSOFF.EQ.0)THEN                                    BMO 0330
C                                                                       BMO 0331
C                     L-SHAPED PATH                                     BMO 0332
                      IF(KNEW.LE.NSPC)                                  BMO 0333
     1                    TDEPTH=TDEPTH+WPATHS(LS,KP)*TAIL*PATMS(LS,K)  BMO 0334
                      IF(KNEW.GT.NSPC)                                  BMO 0335
     1                    TDEPTH=TDEPTH+WPTHSX(LS,KP)*TAIL              BMO 0336
                  ELSE                                                  BMO 0337
C                                                                       BMO 0338
C                     SOLAR PATH ONLY                                   BMO 0339
                      IF(KNEW.LE.NSPC)                                  BMO 0340
     1                    TDEPTH=WPATHS(MSOFFL,KP)*TAIL*PATMS(MSOFFL,K) BMO 0341
                      IF(KNEW.GT.NSPC)                                  BMO 0342
     1                    TDEPTH=WPTHSX(MSOFFL,KP)*TAIL                 BMO 0343
                  ENDIF                                                 BMO 0344
              ENDIF                                                     BMO 0345
C                                                                       BMO 0346
C             CALCULATE LINE TAIL TRANSMITTANCE                         BMO 0347
              IF(KNEW.LE.NSPC)THEN                                      BMO 0348
                 TX(KP)=EXP(-TDEPTH)                                    BMO 0349
                 IF(TX(KP).GT.1.)TX(KP)=1.                              BMO 0350
              ELSE                                                      BMO 0351
                 TXX(KP)=EXP(-TDEPTH)                                   BMO 0352
                 IF(TXX(KP).GT.1.)TXX(KP)=1.                            BMO 0353
              ENDIF                                                     BMO 0354
          ELSE                                                          BMO 0355
              IF(KNEW.LE.NSPC)THEN                                      BMO 0356
                 TX(KP)=1.                                              BMO 0357
              ELSE                                                      BMO 0358
                 TXX(KP)=1.                                             BMO 0359
              ENDIF                                                     BMO 0360
          ENDIF                                                         BMO 0361
  160 CONTINUE                                                          BMO 0362
      RETURN                                                            BMO 0363
      END                                                               BMO 0364
