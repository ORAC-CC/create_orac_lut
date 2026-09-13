      SUBROUTINE CRMERG                                                 CRM 0001
C                                                                       CRM 0002
C     THIS ROUTINE MERGES TOGETHER CLOUD/RAIN AND OLD ATMOSPHERIC       CRM 0003
C     PROFILES.  WATER DROPLET DENSITIES [GM/M3], ICE PARTICLE          CRM 0004
C     DENSITIES [GM/M3] AND RAIN RATES [MM/HR] ARE STORED IN            CRM 0005
C     DENSTY(66,.), DENSTY(67,.) AND DENSTY(3,.), RESPECTIVELY.         CRM 0006
C     LAYER BOUNDARIES ARE MERGED TOGETHER IF THEY DIFFER BY LESS       CRM 0007
C     THAN "ZTOL" KM (HALF A METER).                                    CRM 0008
C                                                                       CRM 0009
C     LIST PARAMETERS                                                   CRM 0010
      INCLUDE 'PARAM.LST'                                               CRM 0011
C                                                                       CRM 0012
C     LIST COMMONS                                                      CRM 0013
      INTEGER KPOINT                                                    CRM 0014
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     CRM 0015
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   CRM 0016
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   CRM 0017
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     CRM 0018
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2                             CRM 0019
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         CRM 0020
     1  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)    CRM 0021
      REAL WNO,WSO2,WNO2,WNH3,WAIR,WHNO3                                CRM 0022
      COMMON/MDATA1/WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),              CRM 0023
     1  WNH3(LAYDIM),WAIR(LAYDIM),WHNO3(LAYDIM)                         CRM 0024
      REAL WMOLXT                                                       CRM 0025
      COMMON/MDATAX/WMOLXT(MMOLX,LAYDIM)                                CRM 0026
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               CRM 0027
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           CRM 0028
      REAL ZM,PM,TM,RFNDX,DENSTY                                        CRM 0029
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    CRM 0030
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               CRM 0031
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     CRM 0032
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                CRM 0033
      INTEGER NBND                                                      CRM 0034
      REAL ZCLDRN(NZCLD),DRPWAT(NZCLD),PRTICE(NZCLD),RNPROF(NZCLD)      CRM 0035
      COMMON/CLDRN/NBND,ZCLDRN,DRPWAT,PRTICE,RNPROF                     CRM 0036
      INTEGER NCRALT,NCRSPC                                             CRM 0037
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    CRM 0038
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      CRM 0039
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       CRM 0040
C                                                                       CRM 0041
C     DECLARE LOCAL VARIABLES                                           CRM 0042
      INTEGER NEWBND,MBND,IBND,I,MLOLD                                  CRM 0043
      LOGICAL LSET,WARNWD,WARNIP                                        CRM 0044
      REAL FAC,TRATIO,CLDHUM                                            CRM 0045
C                                                                       CRM 0046
C     DECLARE FUNCTION NAMES                                            CRM 0047
      REAL EXPINT                                                       CRM 0048
C                                                                       CRM 0049
C     LIST DATA                                                         CRM 0050
C       ZTOL     TOLERANCE FOR MERGING TOGETHER LAYERS [KM]             CRM 0051
C       TFREEZ   TEMPERATURE BELOW WHICH A WARNING IS GENERATED IF      CRM 0052
C                LIQUID WATER DROPLET DENSITY IS POSITIVE [K]           CRM 0053
C       TMELT    TEMPERATURE ABOVE WHICH A WARNING IS GENERATED IF      CRM 0054
C                ICE PARTICLE DENSITY IS POSITIVE [K]                   CRM 0055
      REAL ZTOL,TFREEZ,TMELT                                            CRM 0056
      DATA ZTOL/.0005/,TFREEZ/260./,TMELT/278./                         CRM 0057
C                                                                       CRM 0058
C     CHECK THAT CLOUD/RAIN ALTITUDES ARE BOUNDED BY ZM ALTITUDES       CRM 0059
      IF(ZCLDRN(1).LT.ZM(1) .OR. ZCLDRN(NBND).GT.ZM(ML))THEN            CRM 0060
          WRITE(IPR,'(/A,2(F10.5,A),/14X,A,2(F10.5,A))')                CRM 0061
     1      ' FATAL ERROR:  CLOUD/RAIN MODEL BOUNDING ALTITUDES (',     CRM 0062
     2      ZCLDRN(1),' AND',ZCLDRN(NBND),' KM) ARE',                   CRM 0063
     3      ' NOT BRACKETED BY ATMOSPHERE BOUNDING ALTITUDES (',        CRM 0064
     4      ZM(1),' AND', ZM(ML),' KM).'                                CRM 0065
          STOP                                                          CRM 0066
      ENDIF                                                             CRM 0067
C                                                                       CRM 0068
C     CHECK THAT RAIN RATE AND CLOUD WATER DROPLET AND ICE PARTICLE     CRM 0069
C     DENSITIES ARE INDEED ZERO AT THE TOP OF THE PROFILE.              CRM 0070
      IF(DRPWAT(NBND).NE.0. .OR. PRTICE(NBND).NE.0.                     CRM 0071
     1                      .OR. RNPROF(NBND).NE.0.)THEN                CRM 0072
          WRITE(IPR,'(/2A,/14X,A)')' FATAL ERROR:  RAIN RATE AND',      CRM 0073
     1      ' CLOUD WATER DROPLET AND ICE PARTICLE DENSITIES AND',      CRM 0074
     2      ' ARE NOT ZERO AT THE TOP OF THE CLOUD/RAIN PROFILE.'       CRM 0075
          STOP                                                          CRM 0076
      ENDIF                                                             CRM 0077
C                                                                       CRM 0078
C     INITIALIZE WARNING MESSAGE VARIABLES                              CRM 0079
C     WARNWD   A LOGICAL VARIABLE THAT IS SET TO FALSE AFTER THE USER   CRM 0080
C              HAS BEEN WARNED THAT WATER DROPLETS EXIST BELOW TFREEZ.  CRM 0081
C     WARNIP   A LOGICAL VARIABLE THAT IS SET TO FALSE AFTER THE USER   CRM 0082
C              HAS BEEN WARNED THAT ICE PARTICLES EXIST ABOVE TMELT.    CRM 0083
      WARNWD=.TRUE.                                                     CRM 0084
      WARNIP=.TRUE.                                                     CRM 0085
C                                                                       CRM 0086
C     DETERMINE THE RELATIVE HUMIDITY [%] WITHIN THE CLOUD.             CRM 0087
C     THE RELATIVE HUMIDITY MUST BE POSITIVE (A DRY ATMOSPHERE          CRM 0088
C     IS NOT ALLOWED), AND SUPER-SATURATION IS LIMITED TO 5%.           CRM 0089
      CLDHUM=100.                                                       CRM 0090
      IF(CHUMID.GT.0. .AND. CHUMID.LE.105.)CLDHUM=CHUMID                CRM 0091
C                                                                       CRM 0092
C     TEMPORARILY TRANSLATE BOUNDARY LAYER DATA TO MAXIMUM              CRM 0093
C     LAYER INDICES TO MAKE SPACE FOR ADDITIONAL LAYERS.                CRM 0094
      NEWBND=LAYDIM+1                                                   CRM 0095
      DO 20 MBND=ML,2,-1                                                CRM 0096
          NEWBND=NEWBND-1                                               CRM 0097
          ZM(NEWBND)=ZM(MBND)                                           CRM 0098
          P(NEWBND)=P(MBND)                                             CRM 0099
          T(NEWBND)=T(MBND)                                             CRM 0100
          RELHUM(NEWBND)=RELHUM(MBND)                                   CRM 0101
          WH(NEWBND)=WH(MBND)                                           CRM 0102
          WCO2(NEWBND)=WCO2(MBND)                                       CRM 0103
          WO(NEWBND)=WO(MBND)                                           CRM 0104
          WN2O(NEWBND)=WN2O(MBND)                                       CRM 0105
          WCO(NEWBND)=WCO(MBND)                                         CRM 0106
          WCH4(NEWBND)=WCH4(MBND)                                       CRM 0107
          WO2(NEWBND)=WO2(MBND)                                         CRM 0108
          WHNO3(NEWBND)=WHNO3(MBND)                                     CRM 0109
          WNO(NEWBND)=WNO(MBND)                                         CRM 0110
          WSO2(NEWBND)=WSO2(MBND)                                       CRM 0111
          WNO2(NEWBND)=WNO2(MBND)                                       CRM 0112
          WNH3(NEWBND)=WNH3(MBND)                                       CRM 0113
          WAIR(NEWBND)=WAIR(MBND)                                       CRM 0114
          DO 10 I=1,NSPECX                                              CRM 0115
              WMOLXT(I,NEWBND)=WMOLXT(I,MBND)                           CRM 0116
   10     CONTINUE                                                      CRM 0117
          DENSTY(7,NEWBND)=DENSTY(7,MBND)                               CRM 0118
          DENSTY(12,NEWBND)=DENSTY(12,MBND)                             CRM 0119
          DENSTY(13,NEWBND)=DENSTY(13,MBND)                             CRM 0120
          DENSTY(14,NEWBND)=DENSTY(14,MBND)                             CRM 0121
          DENSTY(15,NEWBND)=DENSTY(15,MBND)                             CRM 0122
          DENSTY(16,NEWBND)=DENSTY(16,MBND)                             CRM 0123
   20 CONTINUE                                                          CRM 0124
C                                                                       CRM 0125
C     INITIALIZE RAIN RATE AND WATER DROPLET AND                        CRM 0126
C     ICE PARTICLE DENSITIES AT THE GROUND.                             CRM 0127
      ML=1                                                              CRM 0128
      DENSTY(66,1)=0.                                                   CRM 0129
      DENSTY(67,1)=0.                                                   CRM 0130
      DENSTY(3,1)=0.                                                    CRM 0131
C                                                                       CRM 0132
C     LOOP OVER OLD ATMOSPHERIC LAYERS                                  CRM 0133
      IBND=1                                                            CRM 0134
      DO 60 MBND=NEWBND,LAYDIM                                          CRM 0135
          LSET=.FALSE.                                                  CRM 0136
   30     CONTINUE                                                      CRM 0137
          IF(ZCLDRN(IBND).LE.ZM(MBND))THEN                              CRM 0138
C                                                                       CRM 0139
C             CLOUD LAYER BOUNDARY WITHIN OLD ATMOSPHERIC LAYER         CRM 0140
              IF(ZCLDRN(IBND).LT.ZM(ML)+ZTOL)THEN                       CRM 0141
C                                                                       CRM 0142
C                 MERGE INTO LOWER LAYER BOUNDARY                       CRM 0143
                  DENSTY(66,ML)=DRPWAT(IBND)                            CRM 0144
                  DENSTY(67,ML)=PRTICE(IBND)                            CRM 0145
                  DENSTY(3,ML)=RNPROF(IBND)                             CRM 0146
                  IF(DRPWAT(IBND).GT.0. .OR. PRTICE(IBND).GT.0.         CRM 0147
     1                                  .OR. RNPROF(IBND).GT.0.)THEN    CRM 0148
                      RELHUM(ML)=CLDHUM                                 CRM 0149
                      TRATIO=273.15/T(ML)                               CRM 0150
                      WH(ML)=.01*CLDHUM*TRATIO*                         CRM 0151
     1                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)    CRM 0152
                      IF(WARNWD .AND. DRPWAT(IBND).GT.0.                CRM 0153
     1                          .AND. T(ML).LT.TFREEZ)THEN              CRM 0154
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   CRM 0155
     1                      ' WARNING:  AT ALTITUDE',ZM(ML),            CRM 0156
     2                      ' KM, THE LIQUID WATER DROPLET DENSITY',    CRM 0157
     3                      ' IS POSITIVE [',DRPWAT(IBND),' GM/M3]',    CRM 0158
     4                      ' EVEN THOUGH',' THE TEMPERATURE IS',       CRM 0159
     5                      T(ML)-273.15,' DEGREES',' CELSIUS.'         CRM 0160
                          WARNWD=.FALSE.                                CRM 0161
                      ELSEIF(WARNIP .AND. PRTICE(IBND).GT.0.            CRM 0162
     1                              .AND. T(ML).GT.TMELT)THEN           CRM 0163
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   CRM 0164
     1                      ' WARNING:  AT ALTITUDE',ZM(ML),            CRM 0165
     2                      ' KM, THE ICE PARITCLE DENSITY',            CRM 0166
     3                      ' IS POSITIVE [',PRTICE(IBND),' GM/M3]',    CRM 0167
     4                      ' EVEN THOUGH',' THE TEMPERATURE IS',       CRM 0168
     5                      T(ML)-273.15,' DEGREES',' CELSIUS.'         CRM 0169
                          WARNIP=.FALSE.                                CRM 0170
                      ENDIF                                             CRM 0171
                  ENDIF                                                 CRM 0172
              ELSEIF(ZCLDRN(IBND).GT.ZM(MBND)-ZTOL)THEN                 CRM 0173
C                                                                       CRM 0174
C                 MERGE INTO UPPER LAYER BOUNDARY                       CRM 0175
                  LSET=.TRUE.                                           CRM 0176
                  DENSTY(66,MBND)=DRPWAT(IBND)                          CRM 0177
                  DENSTY(67,MBND)=PRTICE(IBND)                          CRM 0178
                  DENSTY(3,MBND)=RNPROF(IBND)                           CRM 0179
                  IF(DRPWAT(IBND).GT.0. .OR. PRTICE(IBND).GT.0.         CRM 0180
     1                                  .OR. RNPROF(IBND).GT.0.)THEN    CRM 0181
                      RELHUM(MBND)=CLDHUM                               CRM 0182
                      TRATIO=273.15/T(MBND)                             CRM 0183
                      WH(MBND)=.01*CLDHUM*TRATIO*                       CRM 0184
     1                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)    CRM 0185
                  ENDIF                                                 CRM 0186
              ELSEIF(ML+1.LT.MBND)THEN                                  CRM 0187
C                                                                       CRM 0188
C                 ADD NEW ZM LAYER BOUNDARY                             CRM 0189
                  MLOLD=ML                                              CRM 0190
                  ML=ML+1                                               CRM 0191
                  ZM(ML)=ZCLDRN(IBND)                                   CRM 0192
                  DENSTY(66,ML)=DRPWAT(IBND)                            CRM 0193
                  DENSTY(67,ML)=PRTICE(IBND)                            CRM 0194
                  DENSTY(3,ML)=RNPROF(IBND)                             CRM 0195
                  FAC=(ZM(ML)-ZM(MLOLD))/(ZM(MBND)-ZM(MLOLD))           CRM 0196
                  P(ML)=EXPINT(P(MLOLD),P(MBND),FAC)                    CRM 0197
                  T(ML)=EXPINT(T(MLOLD),T(MBND),FAC)                    CRM 0198
                  TRATIO=273.15/T(ML)                                   CRM 0199
                  IF(DRPWAT(IBND).GT.0. .OR. PRTICE(IBND).GT.0.         CRM 0200
     1                                  .OR. RNPROF(IBND).GT.0.)THEN    CRM 0201
                      RELHUM(ML)=CLDHUM                                 CRM 0202
                      WH(ML)=.01*CLDHUM*TRATIO*                         CRM 0203
     1                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)    CRM 0204
                      IF(WARNWD .AND. DRPWAT(IBND).GT.0.                CRM 0205
     1                  .AND. T(ML).LT.TFREEZ)THEN                      CRM 0206
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   CRM 0207
     1                      ' WARNING:  AT ALTITUDE',ZM(ML),            CRM 0208
     2                      ' KM, THE LIQUID WATER DROPLET DENSITY',    CRM 0209
     3                      ' IS POSITIVE [',DRPWAT(IBND),' GM/M3]',    CRM 0210
     4                      ' EVEN THOUGH',' THE TEMPERATURE IS',       CRM 0211
     5                      T(ML)-273.15,' DEGREES',' CELSIUS.'         CRM 0212
                          WARNWD=.FALSE.                                CRM 0213
                      ELSEIF(WARNIP .AND. PRTICE(IBND).GT.0.            CRM 0214
     1                              .AND. T(ML).GT.TMELT)THEN           CRM 0215
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   CRM 0216
     1                      ' WARNING:  AT ALTITUDE',ZM(ML),            CRM 0217
     2                      ' KM, THE ICE PARITCLE DENSITY',            CRM 0218
     3                      ' IS POSITIVE [',PRTICE(IBND),' GM/M3]',    CRM 0219
     4                      ' EVEN THOUGH',' THE TEMPERATURE IS',       CRM 0220
     5                      T(ML)-273.15,' DEGREES',' CELSIUS.'         CRM 0221
                          WARNIP=.FALSE.                                CRM 0222
                      ENDIF                                             CRM 0223
                  ELSE                                                  CRM 0224
                      RELHUM(ML)=EXPINT(RELHUM(MLOLD),RELHUM(MBND),FAC) CRM 0225
                      WH(ML)=.01*RELHUM(ML)*TRATIO*                     CRM 0226
     1                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)    CRM 0227
                  ENDIF                                                 CRM 0228
                  WCO2(ML)=EXPINT(WCO2(MLOLD),WCO2(MBND),FAC)           CRM 0229
                  WO(ML)=EXPINT(WO(MLOLD),WO(MBND),FAC)                 CRM 0230
                  WN2O(ML)=EXPINT(WN2O(MLOLD),WN2O(MBND),FAC)           CRM 0231
                  WCO(ML)=EXPINT(WCO(MLOLD),WCO(MBND),FAC)              CRM 0232
                  WCH4(ML)=EXPINT(WCH4(MLOLD),WCH4(MBND),FAC)           CRM 0233
                  WO2(ML)=EXPINT(WO2(MLOLD),WO2(MBND),FAC)              CRM 0234
                  WHNO3(ML)=EXPINT(WHNO3(MLOLD),WHNO3(MBND),FAC)        CRM 0235
                  WNO(ML)=EXPINT(WNO(MLOLD),WNO(MBND),FAC)              CRM 0236
                  WSO2(ML)=EXPINT(WSO2(MLOLD),WSO2(MBND),FAC)           CRM 0237
                  WNO2(ML)=EXPINT(WNO2(MLOLD),WNO2(MBND),FAC)           CRM 0238
                  WNH3(ML)=EXPINT(WNH3(MLOLD),WNH3(MBND),FAC)           CRM 0239
                  WAIR(ML)=EXPINT(WAIR(MLOLD),WAIR(MBND),FAC)           CRM 0240
                  DO 40 I=1,NSPECX                                      CRM 0241
                      WMOLXT(I,ML)=                                     CRM 0242
     1                  EXPINT(WMOLXT(I,MLOLD),WMOLXT(I,MBND),FAC)      CRM 0243
   40             CONTINUE                                              CRM 0244
                  DENSTY(7,ML)=                                         CRM 0245
     1              EXPINT(DENSTY(7,MLOLD),DENSTY(7,MBND),FAC)          CRM 0246
                  DENSTY(12,ML)=                                        CRM 0247
     1              EXPINT(DENSTY(12,MLOLD),DENSTY(12,MBND),FAC)        CRM 0248
                  DENSTY(13,ML)=                                        CRM 0249
     1              EXPINT(DENSTY(13,MLOLD),DENSTY(13,MBND),FAC)        CRM 0250
                  DENSTY(14,ML)=                                        CRM 0251
     1              EXPINT(DENSTY(14,MLOLD),DENSTY(14,MBND),FAC)        CRM 0252
                  DENSTY(15,ML)=                                        CRM 0253
     1              EXPINT(DENSTY(15,MLOLD),DENSTY(15,MBND),FAC)        CRM 0254
                  DENSTY(16,ML)=                                        CRM 0255
     1              EXPINT(DENSTY(16,MLOLD),DENSTY(16,MBND),FAC)        CRM 0256
              ELSE                                                      CRM 0257
C                                                                       CRM 0258
C                 NO MORE SPACE IN ARRAYS FOR AN ADDITIONAL LAYER       CRM 0259
                  WRITE(IPR,'(/3A,/14X,A,2(A,I4))')' FATAL ERROR: ',    CRM 0260
     1              ' FILE "PARAM.LST" PARAMETER "LAYDIM" MUST',        CRM 0261
     2              ' BE INCREASED.',' IT SUFFICES TO INCREASE',        CRM 0262
     3              ' LAYDIM FROM',LAYDIM,' TO',LAYDIM+NBND-IBND+1      CRM 0263
                  STOP                                                  CRM 0264
              ENDIF                                                     CRM 0265
C                                                                       CRM 0266
C             EXIT LOOP IF ALL CLOUD BOUNDARIES HAVE BEEN INTEGRATED.   CRM 0267
              IF(IBND.GE.NBND)GOTO70                                    CRM 0268
C                                                                       CRM 0269
C             INCREMENT CLOUD BOUNDARY INDEX AND START AGAIN            CRM 0270
              IBND=IBND+1                                               CRM 0271
              GOTO30                                                    CRM 0272
          ENDIF                                                         CRM 0273
C                                                                       CRM 0274
C         LINEARLY INTERPOLATE WATER PARTICLE DENSITIES AND RAIN RATE.  CRM 0275
          IF(.NOT.LSET)THEN                                             CRM 0276
              FAC=(ZM(MBND)-ZM(ML))/(ZCLDRN(IBND)-ZM(ML))               CRM 0277
              DENSTY(66,MBND)=DENSTY(66,ML)                             CRM 0278
     1          +FAC*(DRPWAT(IBND)-DENSTY(66,ML))                       CRM 0279
              DENSTY(67,MBND)=DENSTY(67,ML)                             CRM 0280
     1          +FAC*(PRTICE(IBND)-DENSTY(67,ML))                       CRM 0281
              DENSTY(3,MBND)=DENSTY(3,ML)                               CRM 0282
     1          +FAC*(RNPROF(IBND)-DENSTY(3,ML))                        CRM 0283
              IF(DENSTY(66,MBND).GT.0. .OR. DENSTY(67,MBND).GT.0.       CRM 0284
     1                                 .OR. DENSTY( 3,MBND).GT.0.)THEN  CRM 0285
                  RELHUM(MBND)=CLDHUM                                   CRM 0286
                  TRATIO=273.15/T(MBND)                                 CRM 0287
                  WH(MBND)=.01*CLDHUM*TRATIO*                           CRM 0288
     1              EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)        CRM 0289
              ENDIF                                                     CRM 0290
          ENDIF                                                         CRM 0291
C                                                                       CRM 0292
C         TRANSLATE LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX.          CRM 0293
          ML=ML+1                                                       CRM 0294
          ZM(ML)=ZM(MBND)                                               CRM 0295
          DENSTY(66,ML)=DENSTY(66,MBND)                                 CRM 0296
          DENSTY(67,ML)=DENSTY(67,MBND)                                 CRM 0297
          DENSTY(3,ML)=DENSTY(3,MBND)                                   CRM 0298
          P(ML)=P(MBND)                                                 CRM 0299
          T(ML)=T(MBND)                                                 CRM 0300
          IF(WARNWD .AND. DENSTY(66,ML).GT.0. .AND. T(ML).LT.TFREEZ)THENCRM 0301
              WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,A))')                CRM 0302
     1          ' WARNING:  AT ALTITUDE',ZM(ML),                        CRM 0303
     2          ' KM, THE LIQUID WATER DROPLET DENSITY',                CRM 0304
     3          ' IS POSITIVE [',DENSTY(66,ML),' GM/M3] EVEN THOUGH',   CRM 0305
     4          ' THE TEMPERATURE IS',T(ML)-273.15,' DEGREES CELSIUS.'  CRM 0306
              WARNWD=.FALSE.                                            CRM 0307
          ENDIF                                                         CRM 0308
          IF(WARNIP .AND. DENSTY(67,ML).GT.0. .AND. T(ML).GT.TMELT)THEN CRM 0309
              WRITE(IPR,'(//2A,F9.4,A,3(/10X,A,F9.4,A))')' WARNING: ',  CRM 0310
     1          ' AT ALTITUDE',ZM(ML),' KM, THE ICE PARITCLE DENSITY',  CRM 0311
     2          ' IS POSITIVE [',DENSTY(67,ML),' GM/M3] EVEN THOUGH',   CRM 0312
     3          ' THE TEMPERATURE IS',T(ML)-273.15,' DEGREES CELSIUS.'  CRM 0313
              WARNIP=.FALSE.                                            CRM 0314
          ENDIF                                                         CRM 0315
          RELHUM(ML)=RELHUM(MBND)                                       CRM 0316
          WH(ML)=WH(MBND)                                               CRM 0317
          WCO2(ML)=WCO2(MBND)                                           CRM 0318
          WO(ML)=WO(MBND)                                               CRM 0319
          WN2O(ML)=WN2O(MBND)                                           CRM 0320
          WCO(ML)=WCO(MBND)                                             CRM 0321
          WCH4(ML)=WCH4(MBND)                                           CRM 0322
          WO2(ML)=WO2(MBND)                                             CRM 0323
          WHNO3(ML)=WHNO3(MBND)                                         CRM 0324
          WNO(ML)=WNO(MBND)                                             CRM 0325
          WSO2(ML)=WSO2(MBND)                                           CRM 0326
          WNO2(ML)=WNO2(MBND)                                           CRM 0327
          WNH3(ML)=WNH3(MBND)                                           CRM 0328
          WAIR(ML)=WAIR(MBND)                                           CRM 0329
          DO 50 I=1,NSPECX                                              CRM 0330
              WMOLXT(I,ML)=WMOLXT(I,MBND)                               CRM 0331
   50     CONTINUE                                                      CRM 0332
          DENSTY(7,ML)=DENSTY(7,MBND)                                   CRM 0333
          DENSTY(12,ML)=DENSTY(12,MBND)                                 CRM 0334
          DENSTY(13,ML)=DENSTY(13,MBND)                                 CRM 0335
          DENSTY(14,ML)=DENSTY(14,MBND)                                 CRM 0336
          DENSTY(15,ML)=DENSTY(15,MBND)                                 CRM 0337
          DENSTY(16,ML)=DENSTY(16,MBND)                                 CRM 0338
   60 CONTINUE                                                          CRM 0339
   70 CONTINUE                                                          CRM 0340
C                                                                       CRM 0341
C     TRANSLATE REMAINING LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX     CRM 0342
C     SETTING CLOUD PARTICLE DENSITIES AND RAIN RATES TO ZERO.          CRM 0343
      NEWBND=MBND                                                       CRM 0344
      DO 90 MBND=NEWBND,LAYDIM                                          CRM 0345
          ML=ML+1                                                       CRM 0346
          ZM(ML)=ZM(MBND)                                               CRM 0347
          DENSTY(66,ML)=0.                                              CRM 0348
          DENSTY(67,ML)=0.                                              CRM 0349
          DENSTY(3,ML)=0.                                               CRM 0350
          P(ML)=P(MBND)                                                 CRM 0351
          T(ML)=T(MBND)                                                 CRM 0352
          RELHUM(ML)=RELHUM(MBND)                                       CRM 0353
          WH(ML)=WH(MBND)                                               CRM 0354
          WCO2(ML)=WCO2(MBND)                                           CRM 0355
          WO(ML)=WO(MBND)                                               CRM 0356
          WN2O(ML)=WN2O(MBND)                                           CRM 0357
          WCO(ML)=WCO(MBND)                                             CRM 0358
          WCH4(ML)=WCH4(MBND)                                           CRM 0359
          WO2(ML)=WO2(MBND)                                             CRM 0360
          WHNO3(ML)=WHNO3(MBND)                                         CRM 0361
          WNO(ML)=WNO(MBND)                                             CRM 0362
          WSO2(ML)=WSO2(MBND)                                           CRM 0363
          WNO2(ML)=WNO2(MBND)                                           CRM 0364
          WNH3(ML)=WNH3(MBND)                                           CRM 0365
          WAIR(ML)=WAIR(MBND)                                           CRM 0366
          DO 80 I=1,NSPECX                                              CRM 0367
              WMOLXT(I,ML)=WMOLXT(I,MBND)                               CRM 0368
   80     CONTINUE                                                      CRM 0369
          DENSTY(7,ML)=DENSTY(7,MBND)                                   CRM 0370
          DENSTY(12,ML)=DENSTY(12,MBND)                                 CRM 0371
          DENSTY(13,ML)=DENSTY(13,MBND)                                 CRM 0372
          DENSTY(14,ML)=DENSTY(14,MBND)                                 CRM 0373
          DENSTY(15,ML)=DENSTY(15,MBND)                                 CRM 0374
          DENSTY(16,ML)=DENSTY(16,MBND)                                 CRM 0375
   90 CONTINUE                                                          CRM 0376
      RETURN                                                            CRM 0377
      END                                                               CRM 0378
