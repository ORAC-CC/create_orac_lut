      SUBROUTINE GEO(IERROR,BENDNG,MSOFF,ICH1)                          GEO 0001
      INCLUDE 'PARAM.LST'                                               GEO 0002
C                                                                       GEO 0003
C     ROUTINE 'GEO' SERVES AS AN INTERFACE BETWEEN ROUTINE              GEO 0004
C     'DRIVER' AND THE GEOMETRY SUBROUTINES INCLUDING 'GEOINP',         GEO 0005
C     'REDUCE', 'FDBETA', 'EXPINT', 'DPEXNT', 'DPFNMN',                 GEO 0006
C     'DPFISH', 'DPSCHT', 'DPANDX', 'DPRARF', 'DPRFPA', 'DPFILL'        GEO 0007
C     AND 'DPLAYR'.  THESE ROUTINES CALCULATE THE ABSORBER              GEO 0008
C     AMOUNTS FOR A REFRACTED PATH THROUGH THE ATMOSPHERE.              GEO 0009
      INTEGER IERROR,MSOFF,ICH1                                         GEO 0010
      REAL BENDNG                                                       GEO 0011
C                                                                       GEO 0012
C     MMOLX  = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")        GEO 0013
C     MMOL   = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")           GEO 0014
C     NSPC   = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL     GEO 0015
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     GEO 0016
C                                                                       GEO 0017
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           GEO 0018
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      GEO 0019
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        GEO 0020
      CHARACTER*8 CNAMEX                                                GEO 0021
      COMMON/NAMEX/CNAMEX(MMOLX)                                        GEO 0022
      REAL DNSTYX                                                       GEO 0023
      COMMON/MODELX/DNSTYX(MMOLX,LAYDIM)                                GEO 0024
      REAL DENPX,AMTPX                                                  GEO 0025
      COMMON/RFRPTX/DENPX(MMOLX,LAYDIM+1),AMTPX(MMOLX,LAYDIM+1)         GEO 0026
      REAL TXX,WX,WPATHX                                                GEO 0027
      COMMON/NONAME/TXX(MMOLX),WX(MMOLX),WPATHX(LAYTHR,MMOLX)           GEO 0028
      DOUBLE PRECISION DPH1,DPH2,DPANGL,DPPHI,DPHMIN,                   GEO 0029
     1  DPBETA,DPBEND,DPRANG,SMMIN,DPZP,DPPP,DPTP,                      GEO 0030
     2  DPRFN,DPSP,DPTPSU,DPRHOP,DPDENP,DPAMTP,DPPPSU                   GEO 0031
      LOGICAL LSAVE,LNOGEO                                              GEO 0032
C                                                                       GEO 0033
C     SSI COMMENTS ON DOUBLE PRECISION VARIABLES:  /RFRPTH/ IS THE OLD  GEO 0034
C     REFRACTED PATH COMMON BLOCK IN SINGLE PRECISION.  /DPRFRP/ IS     GEO 0035
C     THE SAME COMMON BLOCK IN DOUBLE PRECISION AND IS NEW.  IN THIS    GEO 0036
C     ROUTINE "DP" IS USED AS A PREFIX TO DENOTE THE DOUBLE PRECISION   GEO 0037
C     VARIABLES OF DPRFRP.  THE FOLLOWING ARE THE EXCEPTIONS:           GEO 0038
C       DPRFN    STANDS FOR THE OLD DOUBLE PRECISION RFNDXP             GEO 0039
C       DPPPSU   STANDS FOR THE OLD DOUBLE PRECISION PPSUM              GEO 0040
C       DPTPSU   STANDS FOR THE OLD DOUBLE PRECISION TPSUM              GEO 0041
C       DPRHOP   STANDS FOR THE OLD DOUBLE PRECISION RHOPSM             GEO 0042
C                                                                       GEO 0043
C     SOME OTHER VARIABLES WERE DECLARED DOUBLE PRECISION; THEY ARE ALSOGEO 0044
C     IDENTIFIED BY THE "DP" PREFIX, BUT ARE NOT COMMON BLOCK VARIABLES.GEO 0045
      INTEGER KPOINT                                                    GEO 0046
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     GEO 0047
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   GEO 0048
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   GEO 0049
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     GEO 0050
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               GEO 0051
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           GEO 0052
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     GEO 0053
      REAL TBOUND,SALB                                                  GEO 0054
      LOGICAL MODTRN                                                    GEO 0055
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,        GEO 0056
     1  SALB,MODTRN                                                     GEO 0057
      INTEGER LENN                                                      GEO 0058
      REAL H1,H2,ANGLE,RANGE,BETA,REE                                   GEO 0059
      COMMON/CARD3/H1,H2,ANGLE,RANGE,BETA,REE,LENN                      GEO 0060
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     GEO 0061
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                GEO 0062
      REAL ZM,PM,TM,RFNDX,DENSTY                                        GEO 0063
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    GEO 0064
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               GEO 0065
      REAL RE,ZMAX                                                      GEO 0066
      INTEGER IMAX,IMOD,IPATH                                           GEO 0067
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             GEO 0068
      REAL ZP,PP,TP,RFNDXP,SP,PPSUM,TPSUM,RHOPSM,DENP,AMTP              GEO 0069
      COMMON/RFRPTH/ZP(LAYDIM+1),PP(LAYDIM+1),TP(LAYDIM+1),             GEO 0070
     1  RFNDXP(LAYDIM+1),SP(LAYDIM+1),PPSUM(LAYDIM+1),TPSUM(LAYDIM+1),  GEO 0071
     2  RHOPSM(LAYDIM+1),DENP(KMAX,LAYDIM+1),AMTP(KMAX,LAYDIM+1)        GEO 0072
      COMMON/DPRFRP/DPZP(LAYDIM+1),DPPP(LAYDIM+1),                      GEO 0073
     1  DPTP(LAYDIM+1),DPRFN(LAYDIM+1),DPSP(LAYDIM+1),                  GEO 0074
     2  DPPPSU(LAYDIM+1),DPTPSU(LAYDIM+1),DPRHOP(LAYDIM+1),             GEO 0075
     3  DPDENP(KMAX,LAYDIM+1),DPAMTP(KMAX,LAYDIM+1)                     GEO 0076
      INTEGER JTURN,LJ                                                  GEO 0077
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       GEO 0078
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  GEO 0079
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      GEO 0080
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    GEO 0081
C                                                                       GEO 0082
C     /PATH/                                                            GEO 0083
C       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.               GEO 0084
C       AHT     ALTITUDES AT PATH BOUNDARIES.                           GEO 0085
C       TPH     TEMPERATURE AT PATH BOUNDARIES.                         GEO 0086
C       IMAP    MAPPING FROM PATH SEGMENTS TO LAYERS.                   GEO 0087
      INTEGER IMAP                                                      GEO 0088
      REAL QTHETA,AHT,TPH                                               GEO 0089
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),TPH(LAYTWO),IMAP(LAYTWO)   GEO 0090
      DOUBLE PRECISION DHALFR,DPRNG2                                    GEO 0091
      COMMON/SMALL1/DHALFR,DPRNG2                                       GEO 0092
      LOGICAL LSMALL                                                    GEO 0093
      COMMON/SMALL2/LSMALL                                              GEO 0094
      REAL SMALL                                                        GEO 0095
      COMMON/SMALL3/SMALL                                               GEO 0096
      LOGICAL LPRINT                                                    GEO 0097
      COMMON/CPRINT/LPRINT                                              GEO 0098
      REAL GNDALT                                                       GEO 0099
      COMMON/GRAUND/GNDALT                                              GEO 0100
      REAL EXTV,ABSV,ASYV                                               GEO 0101
      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       GEO 0102
C                                                                       GEO 0103
C     DECLARE FUNCTIONS                                                 GEO 0104
      INTEGER PSLCT                                                     GEO 0105
      REAL EXPINT                                                       GEO 0106
C                                                                       GEO 0107
C     DECLARE LOCAL VARIABLES                                           GEO 0108
      INTEGER JMAX,JMAXP1,I,J,K,L,ILO,IHI,IH1,IH2,JOFF,JOFFM1,ISLCT,    GEO 0109
     1  LENNSV,IAMT,LDEL                                                GEO 0110
      LOGICAL LOGLOS                                                    GEO 0111
      REAL TLRNCE,PZERO,PHI,FAC,WTEM,RANGSV,HMIN,PHISV,PRCNT            GEO 0112
C                                                                       GEO 0113
C     KMOL(K) IS A POINTER USED TO REORDER THE AMOUNTS WHEN PRINTING    GEO 0114
      INTEGER KMOL(17)                                                  GEO 0115
      DATA KMOL/1,2,3,11,8,5,9,10,4,6,7,12,13,14,16,15,17/              GEO 0116
      DATA TLRNCE/0.001/,PZERO/1013.25/                                 GEO 0117
      LSMALL=.FALSE.                                                    GEO 0118
      LPRINT=.TRUE.                                                     GEO 0119
      LNOGEO=.FALSE.                                                    GEO 0120
C                                                                       GEO 0121
C     INITIALIZE CONSTANTS AND CLEAR CUMULATIVE VARIABLES.              GEO 0122
      IF(GNDALT.GT.ZM(1))WRITE(6,'(/A,2F10.6,I5)')                      GEO 0123
     1  ' GNDALT IS > FIRST PROFILE ALTITUDE:',GNDALT,ZM(1),MODEL       GEO 0124
      IERROR=0                                                          GEO 0125
      RE=REE                                                            GEO 0126
      IMOD=ML                                                           GEO 0127
      IMAX=ML                                                           GEO 0128
C                                                                       GEO 0129
C     INITIALIZE CUMULATIVE VARIABLES                                   GEO 0130
      DO 30 I=1,LAYDIM+1                                                GEO 0131
          LJ(I)=0                                                       GEO 0132
          SP(I)=0.                                                      GEO 0133
          PPSUM(I)=0.                                                   GEO 0134
          TPSUM(I)=0.                                                   GEO 0135
          RHOPSM(I)=0.                                                  GEO 0136
          DPSP(I)=DBLE(0.)                                              GEO 0137
          DPPPSU(I)=DBLE(0.)                                            GEO 0138
          DPTPSU(I)=DBLE(0.)                                            GEO 0139
          DPRHOP(I)=DBLE(0.)                                            GEO 0140
          DO 10 K=1,KMAX                                                GEO 0141
              AMTP(K,I)=0.                                              GEO 0142
              DPAMTP(K,I)=DBLE(0.)                                      GEO 0143
   10     CONTINUE                                                      GEO 0144
          DO 20 K=1,NSPECX                                              GEO 0145
              AMTPX(K,I)=0.                                             GEO 0146
   20     CONTINUE                                                      GEO 0147
   30 CONTINUE                                                          GEO 0148
      ZMAX=ZM(IMAX)                                                     GEO 0149
      IF(ITYPE.LT.2)THEN                                                GEO 0150
C                                                                       GEO 0151
C         HORIZONTAL PATH, MODEL EQ 1 TO 7:  INTERPOLATE PROFILE TO H1  GEO 0152
          ZP(1)=H1                                                      GEO 0153
          DPZP(1)=DBLE(ZP(1))                                           GEO 0154
          IF(ML.EQ.1)THEN                                               GEO 0155
              PP(1)=PM(1)                                               GEO 0156
              DPPP(1)=DBLE(PM(1))                                       GEO 0157
              TP(1)=TM(1)                                               GEO 0158
              TPH(1)=TP(1)                                              GEO 0159
              DPTP(1)=DBLE(TM(1))                                       GEO 0160
              LJ(1)=1                                                   GEO 0161
              SP(1)=RANGE                                               GEO 0162
              DPSP(1)=DBLE(RANGE)                                       GEO 0163
              DO 40 K=1,KMAX                                            GEO 0164
                  DENP(K,1)=DENSTY(K,1)                                 GEO 0165
   40         CONTINUE                                                  GEO 0166
              DO 50 K=1,NSPECX                                          GEO 0167
                  DENPX(K,1)=DNSTYX(K,1)                                GEO 0168
   50         CONTINUE                                                  GEO 0169
          ELSEIF(H1.LT.ZM(1))THEN                                       GEO 0170
              WRITE(IPR,'(A,2F10.6)')                                   GEO 0171
     1          ' ERROR HORIZ PATH H1 < ZM(1)',H1,ZM(1)                 GEO 0172
              IERROR=1                                                  GEO 0173
              RETURN                                                    GEO 0174
          ELSE                                                          GEO 0175
              ILO=1                                                     GEO 0176
              DO 60 IHI=2,ML                                            GEO 0177
                  IF(H1.LT.ZM(IHI))GOTO70                               GEO 0178
   60         ILO=IHI                                                   GEO 0179
   70         FAC=(H1-ZM(ILO))/(ZM(IHI)-ZM(ILO))                        GEO 0180
              PP(1)=EXPINT(PM(ILO),PM(IHI),FAC)                         GEO 0181
              TP(1)=TM(ILO)+(TM(IHI)-TM(ILO))*FAC                       GEO 0182
              TPH(1)=TP(1)                                              GEO 0183
              DPTP(1)=DBLE(TM(ILO)+(TM(IHI)-TM(ILO))*FAC)               GEO 0184
              LJ(1)=ILO                                                 GEO 0185
              IF(FAC.GT..5)LJ(1)=IHI                                    GEO 0186
              SP(LJ(1))=RANGE                                           GEO 0187
              DPSP(LJ(1))=DBLE(RANGE)                                   GEO 0188
              DO 80 K=1,KMAX                                            GEO 0189
                  DENP(K,1)=EXPINT(DENSTY(K,ILO),DENSTY(K,IHI),FAC)     GEO 0190
   80         CONTINUE                                                  GEO 0191
              DO 90 K=1,NSPECX                                          GEO 0192
                  DENPX(K,1)=EXPINT(DNSTYX(K,ILO),DNSTYX(K,IHI),FAC)    GEO 0193
   90         CONTINUE                                                  GEO 0194
C                                                                       GEO 0195
C             USE LINEAR INTERPOLATION FOR CLOUDS                       GEO 0196
              DENP(16,1)=DENSTY(16,ILO)                                 GEO 0197
     1          +FAC*(DENSTY(16,IHI)-DENSTY(16,ILO))                    GEO 0198
              DENP(66,1)=DENSTY(66,ILO)                                 GEO 0199
     1          +FAC*(DENSTY(66,IHI)-DENSTY(66,ILO))                    GEO 0200
              DENP(67,1)=DENSTY(67,ILO)                                 GEO 0201
     1          +FAC*(DENSTY(67,IHI)-DENSTY(67,ILO))                    GEO 0202
          ENDIF                                                         GEO 0203
C                                                                       GEO 0204
C         CALCULATE ABSORBER AMOUNTS FOR A HORIZONTAL PATH              GEO 0205
          WRITE(IPR,'(2(A,F11.3),A,I4)')                                GEO 0206
     1      '0HORIZONTAL PATH AT ALTITUDE =',H1,                        GEO 0207
     2      ' KM WITH RANGE =',RANGE,' KM, MODEL =',MODEL               GEO 0208
          IKMAX=1                                                       GEO 0209
          JMAX=1                                                        GEO 0210
          IF(MODEL.EQ.0)THEN                                            GEO 0211
              TP(1)=TM(1)                                               GEO 0212
              PP(1)=PM(1)                                               GEO 0213
              DPTP(1)=DBLE(TM(1))                                       GEO 0214
              DPPP(1)=DBLE(PM(1))                                       GEO 0215
          ENDIF                                                         GEO 0216
          JOFF=MSOFF+1                                                  GEO 0217
          PATM(JOFF)=PP(1)/PZERO                                        GEO 0218
          TBBY(JOFF)=TP(1)                                              GEO 0219
          DO 100 K=1,KMAX                                               GEO 0220
              W(K)=DENP(K,1)*RANGE                                      GEO 0221
              WPATH(JOFF,K)=W(K)                                        GEO 0222
  100     CONTINUE                                                      GEO 0223
          DO 110 K=1,NSPECX                                             GEO 0224
              WX(K)=DENPX(K,1)*RANGE                                    GEO 0225
              WPATHX(JOFF,K)=WX(K)                                      GEO 0226
  110     CONTINUE                                                      GEO 0227
          W(9)=W(5)*(296.-TP(1))/(296.-260.)                            GEO 0228
          WTEM=TP(1)-273.15                                             GEO 0229
          W(59)=W(8)*.269*WTEM                                          GEO 0230
          W(60)=W(59)*WTEM                                              GEO 0231
          WPATH(JOFF,9)=W(9)                                            GEO 0232
          WPATH(JOFF,59)=W(59)                                          GEO 0233
          WPATH(JOFF,60)=W(60)                                          GEO 0234
      ELSE                                                              GEO 0235
C                                                                       GEO 0236
C         SLANT PATH SELECTED.  INTERPRET SLANT PATH PARAMETERS.        GEO 0237
C                                                                       GEO 0238
C         LOGLOS IS A LOGICAL VARIABLE THAT IS TRUE                     GEO 0239
C         ONLY THE OPTICAL PATH WITH ITYPE = 2.                         GEO 0240
          LOGLOS=.FALSE.                                                GEO 0241
          IF(ITYPE.EQ.2 .AND. MSOFF.EQ.0)THEN                           GEO 0242
              LOGLOS=.TRUE.                                             GEO 0243
              ISLCT=PSLCT(ANGLE,RANGE,BETA)                             GEO 0244
C                                                                       GEO 0245
C             SPECIAL TREATMENT EXCEPT FOR CASE 2A (ISLCT=21)           GEO 0246
              IF(ISLCT.GT.21)THEN                                       GEO 0247
C                                                                       GEO 0248
C                 IF RANGE IS SMALL, CONVERT TO CASE 2C (ISLCT=23)      GEO 0249
                  CALL SMPREP(H1,H2,ANGLE,RANGE,BETA,ISLCT)             GEO 0250
                  IF(RANGE.GT.0 .AND. RANGE.LE.SMALL)THEN               GEO 0251
                      LSMALL=.TRUE.                                     GEO 0252
                      RANGSV=RANGE                                      GEO 0253
                      ISLCT= 23                                         GEO 0254
                  ELSEIF(ISLCT.EQ.22)THEN                               GEO 0255
C                                                                       GEO 0256
C                     CASE 2B:  (H1,ANGLE,RANGE)                        GEO 0257
C                     IF PATH TYPE IS CASE 2B, CHECK THAT THE RANGE     GEO 0258
C                     USED IN THE CALCULATION EQUALS THE INPUT RANGE.   GEO 0259
C                     DETERMINE H2 AND LENN USING ROUTINE NEWH2.        GEO 0260
                      LENNSV=LENN                                       GEO 0261
                      RANGSV=RANGE                                      GEO 0262
                      DPH1=DBLE(H1)                                     GEO 0263
                      DPANGL=DBLE(ANGLE)                                GEO 0264
                      DPRANG=DBLE(RANGE)                                GEO 0265
                      CALL NEWH2(DPH1,DPH2,DPANGL,                      GEO 0266
     1                  DPRANG,DPBETA,LENN,DPHMIN,DPPHI)                GEO 0267
                      IF(LENN.EQ.0)DPHMIN=MIN(DPH2,DPH1)                GEO 0268
                      H2=REAL(DPH2)                                     GEO 0269
                      HMIN=REAL(DPHMIN)                                 GEO 0270
                      PHISV=REAL(DPPHI)                                 GEO 0271
                      LPRINT=.FALSE.                                    GEO 0272
                      IAMT=2                                            GEO 0273
                      CALL DPRFPA(DPH1,DPH2,DPANGL,DPPHI,               GEO 0274
     1                  LENN,DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)          GEO 0275
                      LPRINT=.TRUE.                                     GEO 0276
                      PRCNT=100*ABS(REAL(DPRANG)-RANGSV)/RANGSV         GEO 0277
                      WRITE(IPR,'((A))')'1',                            GEO 0278
     1                  ' SOME INTERNAL DETAILS:'                       GEO 0279
                      WRITE(IPR,'(/(A))')                               GEO 0280
     1                  ' LOS IS INTERNAL CASE 2B (H1, ANGLE, RANGE).', GEO 0281
     2                  ' USING H2 OBTAINED FROM SUBROUTINE NEWH2:'     GEO 0282
                      WRITE(IPR,'(/(A,F12.5))')                         GEO 0283
     1                  ' H1                        =',H1,              GEO 0284
     2                  ' H2                        =',H2,              GEO 0285
     3                  ' ANGLE                     =',ANGLE,           GEO 0286
     4                  ' PHI                       =',PHISV,           GEO 0287
     5                  ' BETA                      =',BETA,            GEO 0288
     6                  ' HMIN (MINIMUM ALTITUDE)   =',HMIN,            GEO 0289
     7                  ' RANGE (OUTPUT)            =',DPRANG,          GEO 0290
     8                  ' RANGE (INPUT)             =',RANGSV,          GEO 0291
     9                  ' PERCENT DIFFERENCE        =',PRCNT            GEO 0292
                      WRITE(IPR,'(A,I12,/)')                            GEO 0293
     1                  ' LENN                      =',LENN             GEO 0294
                      IF(ANGLE.GT.90. .AND.                             GEO 0295
     1                  ABS(H2-GNDALT).LT.TLRNCE)THEN                   GEO 0296
                          LNOGEO=.TRUE.                                 GEO 0297
                          PHI=PHISV                                     GEO 0298
                          ANGLE=REAL(DPANGL)                            GEO 0299
                          HMIN=REAL(DPHMIN)                             GEO 0300
                          BETA=REAL(DPBETA)                             GEO 0301
                          RANGE=REAL(DPRANG)                            GEO 0302
                          WRITE(IPR,'(/2A,/)')'***** WARNING ***** ',   GEO 0303
     1                      ' PATH HITS THE EARTH.'                     GEO 0304
                      ELSEIF(ANGLE.LE.90. .AND.                         GEO 0305
     1                  ABS(H2-ZMAX).LT.TLRNCE)THEN                     GEO 0306
                          LNOGEO=.TRUE.                                 GEO 0307
                          PHI=PHISV                                     GEO 0308
                          ANGLE=REAL(DPANGL)                            GEO 0309
                          HMIN=REAL(DPHMIN)                             GEO 0310
                          BETA=REAL(DPBETA)                             GEO 0311
                          RANGE=REAL(DPRANG)                            GEO 0312
                          WRITE(IPR,'(/2A,/)')'***** WARNING ***** ',   GEO 0313
     1                      ' PATH HITS THE UPPERMOST LAYER BOUNDARY.'  GEO 0314
                      ELSEIF(PRCNT.LE.1. .OR. H2.GE.ZMAX)THEN           GEO 0315
                          WRITE(IPR,'(/(2A))')                          GEO 0316
     1                      ' PERCENT DIFFERENCE BEING LESS THAN 1.',   GEO 0317
     2                      ' OR THE PATH TERMINATES AT THE TOP OF',    GEO 0318
     3                      ' THE ATMOSPHERE.  THESE PATH PARAMETERS',  GEO 0319
     4                      ' WILL BE USED WITHOUT CALLING GEOINP.'     GEO 0320
                          LNOGEO=.TRUE.                                 GEO 0321
                          PHI=PHISV                                     GEO 0322
                          ANGLE=REAL(DPANGL)                            GEO 0323
                          HMIN=REAL(DPHMIN)                             GEO 0324
                          BETA=REAL(DPBETA)                             GEO 0325
                          RANGE=REAL(DPRANG)                            GEO 0326
                      ELSE                                              GEO 0327
                          WRITE(IPR,'(/2A,/2A,/)')' SINCE THE PERCENT', GEO 0328
     1                      ' DIFFERENCE EXCEEDS 1,',' "EQUIVALENT"',   GEO 0329
     2                      ' INTERNAL CASE 2C WILL BE USED.'           GEO 0330
                          LNOGEO=.FALSE.                                GEO 0331
                          RANGE=RANGSV                                  GEO 0332
                          BETA=0.                                       GEO 0333
                          HMIN=0.                                       GEO 0334
                          PHI=0.                                        GEO 0335
                          LENN=LENNSV                                   GEO 0336
                          ANGLE=0.                                      GEO 0337
                          ISLCT=PSLCT(ANGLE,RANGE,BETA)                 GEO 0338
                      ENDIF                                             GEO 0339
                  ENDIF                                                 GEO 0340
              ENDIF                                                     GEO 0341
          ENDIF                                                         GEO 0342
          LSAVE=LSMALL                                                  GEO 0343
          LSMALL=.FALSE.                                                GEO 0344
          IF(.NOT.LNOGEO)CALL GEOINP(H1,H2,ANGLE,RANGE,                 GEO 0345
     1      BETA,ITYPE,LENN,HMIN,PHI,IERROR,ISLCT,LOGLOS)               GEO 0346
          LSMALL=LSAVE                                                  GEO 0347
C                                                                       GEO 0348
C         CHECK FOR IERROR                                              GEO 0349
          IF(IERROR.NE.0)THEN                                           GEO 0350
              IF(ISSGEO.NE.1)WRITE(IPR,'(2A)')'0GEO:  IERROR NE 0:',    GEO 0351
     1          ' END THIS CALCULATION AND SKIP TO THE NEXT CASE'       GEO 0352
              RETURN                                                    GEO 0353
          ELSEIF(LSMALL)THEN                                            GEO 0354
              RANGE=RANGSV                                              GEO 0355
              CALL SMGEO(ANGLE,BETA,PHI,DHALFR,DPRNG2,BENDNG,LENN,SMMIN)GEO 0356
              HMIN=REAL(SMMIN)                                          GEO 0357
          ENDIF                                                         GEO 0358
C                                                                       GEO 0359
C         CALCULATE THE PATH THROUGH THE ATMOSPHERE                     GEO 0360
          IAMT=1                                                        GEO 0361
          DPH1=DBLE(H1)                                                 GEO 0362
          DPH2=DBLE(H2)                                                 GEO 0363
          DPANGL=DBLE(ANGLE)                                            GEO 0364
          DPPHI=DBLE(PHI)                                               GEO 0365
          DPHMIN=DBLE(HMIN)                                             GEO 0366
          DPBETA=DBLE(BETA)                                             GEO 0367
          DPRANG=DBLE(RANGE)                                            GEO 0368
          DPBEND=DBLE(BENDNG)                                           GEO 0369
C                                                                       GEO 0370
C         IF ANGLE IS NEAR 90., HMIN IS NEAR H1, AND RANGSV EXCEEDS     GEO 0371
C         SMALL, THEN YOU HAVE A HALF-TANGENT PATH I.E. H1, H2 AND      GEO 0372
C         THE EARTH CENTER FORM A RIGHT TRIANGLE.  IN THIS CASE, DO     GEO 0373
C         NOT LET LENN EQUAL 1 BECAUSE THAT MAY MAKE HMIN IDENTICALLY   GEO 0374
C         EQUAL TO H1 PRODUCING PROBLEMS IN SUBROUTINE FILL.  THE       GEO 0375
C         CHECK ON THE VARIABLE SMALL SIMPLY POINTS OUT THE FACT        GEO 0376
C         SMALL PATHS ARE ALREADY DEALT WITH BY OTHER METHODS.          GEO 0377
          IF(ABS(ANGLE-90.).LT..001 .AND. ABS(HMIN-H1).LT.TLRNCE        GEO 0378
     1      .AND. RANGSV.GT.SMALL)LENN=0                                GEO 0379
          CALL DPRFPA(DPH1,DPH2,DPANGL,DPPHI,                           GEO 0380
     1      LENN,DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)                      GEO 0381
          H1=REAL(DPH1)                                                 GEO 0382
          H2=REAL(DPH2)                                                 GEO 0383
          ANGLE=REAL(DPANGL)                                            GEO 0384
          PHI=REAL(DPPHI)                                               GEO 0385
          HMIN=REAL(DPHMIN)                                             GEO 0386
          IF(HMIN.LT.GNDALT-TLRNCE .AND. (LOGLOS .OR. ITYPE.EQ.3))THEN  GEO 0387
C             LOGLOS ONLY DEALS WITH ITYPE=2.                           GEO 0388
C             THEREFORE, THE CHECK FOR ITYPE=3.                         GEO 0389
              WRITE(IPR,'(2(/A,F12.6),/(A))')'GNDALT =',GNDALT,         GEO 0390
     1          'HMIN =',HMIN,'HMIN IS LESS THAN GNDALT.',              GEO 0391
     2          'THIS RUN ABORTED, NEXT RUN ATTEMPTED',' '              GEO 0392
              IERROR=1                                                  GEO 0393
              RETURN                                                    GEO 0394
          ENDIF                                                         GEO 0395
          BETA=REAL(DPBETA)                                             GEO 0396
          RANGE=REAL(DPRANG)                                            GEO 0397
          BENDNG=REAL(DPBEND)                                           GEO 0398
C                                                                       GEO 0399
C         LOAD LAYER AMOUNTS IN AMTP INTO WPATH FROM H1 TO H2           GEO 0400
          DO 120 I=1,IPATH                                              GEO 0401
              IF(H1.EQ.ZP(I))IH1=I                                      GEO 0402
              IF(H2.EQ.ZP(I))IH2=I                                      GEO 0403
  120     CONTINUE                                                      GEO 0404
          JMAX=(IPATH-1)+LENN*(MIN0(IH1,IH2)-1)                         GEO 0405
          IF(JMAX.GT.LAYTWO)STOP 'JMAX IS GREATER THAN PARAMETER LAYTWO'GEO 0406
          IKMAX=JMAX                                                    GEO 0407
C                                                                       GEO 0408
C         DETERMINE LJ(J), WHICH IS THE LAYER NUMBER L=LJ(J)            GEO 0409
C         IN AMTP(K,L), STARTING FROM HMIN, WHICH CORRESPONDS           GEO 0410
C         TO THE LAYER J IN WPATH(J+MSOFF,K), STARTING FROM             GEO 0411
C         H1 INITIAL DIRECTION OF PATH IS DOWN.                         GEO 0412
          L=IH1                                                         GEO 0413
          LDEL=-1                                                       GEO 0414
          IF(LENN.EQ.0 .AND. H1.LE.H2)THEN                              GEO 0415
C                                                                       GEO 0416
C             INITIAL DIRECTION OF PATH IS UP                           GEO 0417
              L=0                                                       GEO 0418
              LDEL=1                                                    GEO 0419
          ENDIF                                                         GEO 0420
          JTURN=0                                                       GEO 0421
          JMAXP1=JMAX+1                                                 GEO 0422
          DO 130 J=1,JMAXP1                                             GEO 0423
C                                                                       GEO 0424
C             TEST FOR REVERSING DIRECTION OF PATH FROM DOWN TO UP      GEO 0425
              IF(L.EQ.1 .AND. LDEL.EQ.-1)THEN                           GEO 0426
                  JTURN=J                                               GEO 0427
                  L=0                                                   GEO 0428
                  LDEL=1                                                GEO 0429
              ENDIF                                                     GEO 0430
              L=L+LDEL                                                  GEO 0431
              LJ(J)=L                                                   GEO 0432
  130     CONTINUE                                                      GEO 0433
C                                                                       GEO 0434
C         LOAD TBBY = THE DENSITY WEIGHTED MEAN TEMPERATURE             GEO 0435
C         AND WPATH = THE INCREMENTAL LAYER AMOUNTS.                    GEO 0436
          IMAX=0                                                        GEO 0437
          DO 140 K=1,KMAX                                               GEO 0438
              W(K)=0.                                                   GEO 0439
  140     CONTINUE                                                      GEO 0440
          DO 150 K=1,NSPECX                                             GEO 0441
              WX(K)=0.                                                  GEO 0442
  150     CONTINUE                                                      GEO 0443
          DO 180 J=1,JMAX                                               GEO 0444
              L=LJ(J)                                                   GEO 0445
              IF(L.GE.ML)L=ML                                           GEO 0446
              IF(L.GE.IMAX)IMAX=L                                       GEO 0447
              JOFF=MSOFF+J                                              GEO 0448
              TBBY(JOFF)=TPSUM(L)/RHOPSM(L)                             GEO 0449
              PATM(JOFF)=PPSUM(L)/(PZERO*RHOPSM(L))                     GEO 0450
              AMTP(9,L)=AMTP(5,L)*(296.-TBBY(JOFF))/36.                 GEO 0451
              WTEM=TBBY(JOFF)-273.15                                    GEO 0452
              AMTP(59,L)=.269*AMTP(8,L)*WTEM                            GEO 0453
              AMTP(60,L)=AMTP(59,L)*WTEM                                GEO 0454
              DO 160 K=1,KMAX                                           GEO 0455
                  WPATH(JOFF,K)=AMTP(K,L)                               GEO 0456
                  W(K)=W(K)+AMTP(K,L)                                   GEO 0457
  160         CONTINUE                                                  GEO 0458
              DO 170 K=1,NSPECX                                         GEO 0459
                  WPATHX(JOFF,K)=AMTPX(K,L)                             GEO 0460
                  WX(K)=WX(K)+AMTPX(K,L)                                GEO 0461
  170         CONTINUE                                                  GEO 0462
  180     CONTINUE                                                      GEO 0463
          IKMAX=JMAX                                                    GEO 0464
C                                                                       GEO 0465
C         INCLUDE BOUNDARY EMISSION IF CARD 1 CONTAINS A NON-ZERO       GEO 0466
C         TBOUND OR IF THE SLANT PATH INTERSECTS THE EARTH (TBOUND      GEO 0467
C         SET TO TEMPERATURE OF LOWEST BOUNDARY IN THIS CASE).          GEO 0468
          IF(TBOUND.EQ.0. .AND. H2.LE.ZM(1))TBOUND=TM(1)                GEO 0469
C                                                                       GEO 0470
C         PRINT LAYER PATH SEGMENT ABSORBER AMOUNTS                     GEO 0471
          IF(NPR.NE.2)THEN                                              GEO 0472
              ILO=1                                                     GEO 0473
              DO 190 IHI=2,ML-1                                         GEO 0474
                  IF(ZM(ILO).GE.H1)GOTO200                              GEO 0475
  190         ILO=IHI                                                   GEO 0476
              IHI=ML                                                    GEO 0477
  200         FAC=(H1-ZM(ILO))/(ZM(IHI)-ZM(ILO))                        GEO 0478
              TPH(1)=TM(ILO)+FAC*(TM(IHI)-TM(ILO))                      GEO 0479
C                                                                       GEO 0480
C             LDEL=0 => GOING DOWN;    LDEL=1 => GOING UP.              GEO 0481
              LDEL=1                                                    GEO 0482
              IF(LENN.EQ.1 .OR. H1.GT.H2)LDEL=0                         GEO 0483
              AHT(1)=H1                                                 GEO 0484
              DO 210 J=1,JMAX                                           GEO 0485
                  L=LJ(J)                                               GEO 0486
                  IF(J.EQ.JTURN)LDEL=1                                  GEO 0487
                  AHT(J+1)=ZP(L+LDEL)                                   GEO 0488
                  TPH(J+1)=TP(L+LDEL)                                   GEO 0489
  210         CONTINUE                                                  GEO 0490
              IF(NPR.NE.1)THEN                                          GEO 0491
                  WRITE(IPR,'(////2A,//(3A))')'1  LAYER ABSORBER',      GEO 0492
     1              ' AMOUNTS FOR THE PATH SEGMENT ENDING AT Z',        GEO 0493
     2              '  J      Z        TBAR       HNO3         O3 UV ', GEO 0494
     3              '    CNTMSLF1     CNTMSLF2      CNTMFRN ',          GEO 0495
     4              '       O2         WAT DROP     ICE PART',          GEO 0496
     5              '       (KM)        (K)     (ATM CM)     (ATM CM)', GEO 0497
     6              '   (MOL CM-2)   (MOL CM-2)   (MOL CM-2)',          GEO 0498
     7              '   (MOL CM-2)   (KM GM/M3)   (KM GM/M3)'           GEO 0499
                  WRITE(IPR,'(/(I3,0PF10.4,F9.2,1P6E13.3,0P2F13.7))')   GEO 0500
     1              (J,AHT(J+1),TBBY(MSOFF+J),                          GEO 0501
     2              (WPATH(MSOFF+J,KMOL(K)),K=4,8),WPATH(MSOFF+J,58),   GEO 0502
     3              (WPATH(MSOFF+J,K),K=66,67),J=1,JMAX)                GEO 0503
                  WRITE(IPR,'(///A,T10,A,T19,A,T31,A,T47,A,T60,A,       GEO 0504
     1              T73,A,T86,A,T99,A,T110,A,T123,A,/T9,A,T109,A)')     GEO 0505
     2              '1 J','Z','N2 CONT','MOL SCAT','AER 1','AER 2',     GEO 0506
     3              'AER 3','AER 4','CIRRUS','WAT DROP','ICE PART',     GEO 0507
     4              '(KM)','(550 NM OPTICAL DEPTH)'                     GEO 0508
                  WRITE(IPR,'(/(I3,0PF10.4,9F13.7))')(J,AHT(J+1),       GEO 0509
     1              (WPATH(MSOFF+J,KMOL(K)),K=9,15),WPATH(MSOFF+J,66)   GEO 0510
     2              *EXTV(6),WPATH(MSOFF+J,67)*EXTV(7),J=1,JMAX)        GEO 0511
              ENDIF                                                     GEO 0512
          ENDIF                                                         GEO 0513
C                                                                       GEO 0514
C         PRINT PATH SUMMARY                                            GEO 0515
          IF(ISSGEO.NE.1)THEN                                           GEO 0516
              IF(NPR.LT.1)THEN                                          GEO 0517
                  WRITE(IPR,'(///3A)')                                  GEO 0518
     1              '1    J    Z       H2O       O3        CO2 ',       GEO 0519
     2              '      CO        CH4       N2O       O2  ',         GEO 0520
     3              '      NH3       NO        NO2       SO2'           GEO 0521
                  IF(MODTRN)THEN                                        GEO 0522
                      WRITE(IPR,'(8X,A,44X,A,45X,A)')                   GEO 0523
     1                  '(KM)    (          ','ATM CM',')'              GEO 0524
                  ELSE                                                  GEO 0525
                      WRITE(IPR,'(8X,A,44X,A,45X,A)')                   GEO 0526
     1                  '(KM)   (G/CM**2)  (','ATM CM',')'              GEO 0527
                  ENDIF                                                 GEO 0528
                  WRITE(IPR,'(/(I4,0PF10.4,1P11E10.2))')(J,AHT(J+1),    GEO 0529
     1              WPATH(MSOFF+J,17),WPATH(MSOFF+J,31),                GEO 0530
     2              WPATH(MSOFF+J,36),WPATH(MSOFF+J,44),                GEO 0531
     3              WPATH(MSOFF+J,46),WPATH(MSOFF+J,47),                GEO 0532
     4              WPATH(MSOFF+J,50),WPATH(MSOFF+J,52),                GEO 0533
     5              WPATH(MSOFF+J,54),WPATH(MSOFF+J,55),                GEO 0534
     6              WPATH(MSOFF+J,56),J=1,JMAX)                         GEO 0535
                  WRITE(IPR,'(///A14,13(1X,A8:),/(14X,13(1X,A8)))')     GEO 0536
     1              '1  J      Z   ',(CNAMEX(K),K=1,MIN(22,NSPECX))     GEO 0537
                  WRITE(IPR,'(8X,A,2(53X,A))')'(KM)    (','ATM CM',')'  GEO 0538
                  DO 220 J=1,JMAX                                       GEO 0539
                      JOFF=MSOFF+J                                      GEO 0540
                      WRITE(IPR,'(I4,F10.4,1P13E9.2:,/(14X,1P13E9.2:))')GEO 0541
     1                  J,AHT(J+1),(WPATHX(JOFF,K),K=1,MIN(22,NSPECX))  GEO 0542
 220              CONTINUE                                              GEO 0543
              ENDIF                                                     GEO 0544
              WRITE(IPR,'(//A,/8(/10X,A,F11.5,A),/10X,A,I11)')          GEO 0545
     1          '0SUMMARY OF THE GEOMETRY CALCULATION',                 GEO 0546
     2          'H1      =',H1,   ' KM', 'H2      =',H2,    ' KM',      GEO 0547
     3          'ANGLE   =',ANGLE,' DEG','RANGE   =',RANGE, ' KM',      GEO 0548
     4          'BETA    =',BETA, ' DEG','PHI     =',PHI,   ' DEG',     GEO 0549
     5          'HMIN    =',HMIN, ' KM', 'BENDING =',BENDNG,' DEG',     GEO 0550
     6          'LENN    =',LENN                                        GEO 0551
          ENDIF                                                         GEO 0552
      ENDIF                                                             GEO 0553
C                                                                       GEO 0554
C     CALCULATE THE AEROSOL WEIGHTED MEAN RELATIVE HUMIDITY (RH)        GEO 0555
      IF(W(7).GT.0. .AND. ICH1.LE.7)THEN                                GEO 0556
          W(15)=100.-EXP(W(15)/W(7))                                    GEO 0557
      ELSEIF(W(12).GT.0. .AND. ICH1.GT.7)THEN                           GEO 0558
          W(15)=100.-EXP(W(15)/W(12))                                   GEO 0559
      ELSE                                                              GEO 0560
          W(15)=0.                                                      GEO 0561
      ENDIF                                                             GEO 0562
C                                                                       GEO 0563
C     CONVERT MOLECULAR BAND WPATH TO CUMULATIVE                        GEO 0564
C     PATH AMOUNTS FOR LOWTRAN RUNS.                                    GEO 0565
      IF(.NOT.MODTRN)THEN                                               GEO 0566
          JOFFM1=MSOFF+1                                                GEO 0567
          DO 240 J=2,JMAX                                               GEO 0568
              JOFF=MSOFF+J                                              GEO 0569
              DO 230 K=17,57                                            GEO 0570
                  WPATH(JOFF,K)=WPATH(JOFFM1,K)+WPATH(JOFF,K)           GEO 0571
  230         CONTINUE                                                  GEO 0572
  240     JOFFM1=JOFF                                                   GEO 0573
      ENDIF                                                             GEO 0574
C                                                                       GEO 0575
C     PRINT TOTAL PATH AMOUNTS                                          GEO 0576
      IF(H1.LT.ZM(1))THEN                                               GEO 0577
          WRITE(IPR,'(/2A,F8.4,A,/14X,A,F8.4,A)')' FATAL ERROR: ',      GEO 0578
     1      ' INITIAL ALTITUDE [H1 =',H1,' KM] IS BELOW THE',           GEO 0579
     2      ' BOTTOM OF THE ATMOSPHERE [ZM(1) =',ZM(1),' KM].'          GEO 0580
          IERROR=1                                                      GEO 0581
          RETURN                                                        GEO 0582
      ELSEIF(H2.LT.ZM(1) .AND. ITYPE.NE.1)THEN                          GEO 0583
          WRITE(IPR,'(/2A,F8.4,A,/14X,A,F8.4,A)')' FATAL ERROR: ',      GEO 0584
     1      ' FINAL OR TANGENT ALTITUDE [H2 =',H2,' KM] IS BELOW',      GEO 0585
     2      ' THE BOTTOM OF THE ATMOSPHERE [ZM(1) =',ZM(1),' KM].'      GEO 0586
          IERROR=1                                                      GEO 0587
          RETURN                                                        GEO 0588
      ENDIF                                                             GEO 0589
      IF(ISSGEO.EQ.1)RETURN                                             GEO 0590
      WRITE(IPR,'(////A,//15X,2A,/15X,A,//10X,1P7E12.4,                 GEO 0591
     1  /(//15X,2A:,/73X,A,//10X,0P7F12.6,F12.2))')                     GEO 0592
     2  '1  EQUIVALENT SEA LEVEL TOTAL ABSORBER AMOUNTS',               GEO 0593
     3  '  HNO3      O3 UV      CNTMSLF1    CNTMSLF2    CNTMFRN  ',     GEO 0594
     4  '   N2 CONT     MOL SCAT',                                      GEO 0595
     5  '(ATM CM)   (ATM CM)   (MOL CM-2)  (MOL CM-2)  (MOL CM-2)',     GEO 0596
     6  (W(KMOL(K)),K=4,10),' AER 1       AER 2       AER 3      ',     GEO 0597
     7  ' AER 4      CIRRUS     WAT DROP    ICE PART   MEAN AER RH',    GEO 0598
     8  '(KM GM/M3)  (KM GM/M3)    (PRCNT)',                            GEO 0599
     9  (W(KMOL(K)),K=11,15),W(66),W(67),W(KMOL(16)),' H2O     ',       GEO 0600
     &  '    O3          CO2         CO          CH4         N2O'       GEO 0601
      IF(MODTRN)THEN                                                    GEO 0602
          WRITE(IPR,'(13X,A,2(35X,A))')'(','ATM CM',')'                 GEO 0603
      ELSE                                                              GEO 0604
          WRITE(IPR,'(13X,A,2(22X,A))')'(G/CM**2)     (','ATM CM',')'   GEO 0605
      ENDIF                                                             GEO 0606
      WRITE(IPR,'((/10X,1P6E12.4,//2(/15X,A),2(22X,A)))')               GEO 0607
     1  W(17),W(31),W(36),W(44),W(46),W(47),                            GEO 0608
     2  ' O2          NH3         NO          NO2         SO2',         GEO 0609
     3  '(','ATM CM',')',W(50),W(52),W(54),W(55),W(56)                  GEO 0610
      IHI=7                                                             GEO 0611
  250 IF(NSPECX.GT.IHI)THEN                                             GEO 0612
          WRITE(IPR,'(//8X,7(4X,A),/13X,A,2(36X,A),//10X,1P7E12.4)')    GEO 0613
     1      (CNAMEX(I),I=IHI-6,IHI),'(','ATM CM',')',(WX(I),I=IHI-6,IHI)GEO 0614
          IHI=IHI+7                                                     GEO 0615
          GOTO250                                                       GEO 0616
      ENDIF                                                             GEO 0617
      WRITE(IPR,'(//8X,7(4X,A))')(CNAMEX(I),I=IHI-6,NSPECX)             GEO 0618
      WRITE(IPR,'(13X,A,2(36X,A),//10X,1P7E12.4)')                      GEO 0619
     1  '(','ATM CM',')',(WX(I),I=IHI-6,NSPECX)                         GEO 0620
      RETURN                                                            GEO 0621
      END                                                               GEO 0622
