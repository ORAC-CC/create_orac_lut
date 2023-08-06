      SUBROUTINE AERNSM(JPRT,GNDALT,ICH)                                ANS 0001
      INTEGER ICH(4)                                                    ANS 0002
      INCLUDE 'PARAM.LST'                                               ANS 0003
C********************************************************************** ANS 0004
C     DEFINES ALTITUDE DEPENDENT VARIABLES Z,P,T,WH,WO AND HAZE         ANS 0005
C     CLD RAIN  CLDTYPE                                                 ANS 0006
C     IT ALSO DEFINES ALTITUDE DEPENDENT VARIABLES WAIR,WCO2,WCO,       ANS 0007
C     WCH4,WN2O,WO2,WNH3,WNO,WNO2, AND WSO2                             ANS 0008
C     LOADS HAZE INTO APPROPRATE LOCATION                               ANS 0009
C********************************************************************** ANS 0010
C                                                                       ANS 0011
C                                                                       ANS 0012
C                                                                       ANS 0013
C     CONVENTION                                                        ANS 0014
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         ANS 0015
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            ANS 0016
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             ANS 0017
C                                                                       ANS 0018
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              ANS 0019
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       ANS 0020
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     ANS 0021
C                                                                       ANS 0022
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           ANS 0023
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      ANS 0024
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        ANS 0025
                                                                        ANS 0026
C                                                                       ANS 0027
C                                                                       ANS 0028
C                                                                       ANS 0029
C     TRANS VARIABLES                                                   ANS 0030
C                                                                       ANS 0031
      INTEGER KPOINT                                                    ANS 0032
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     ANS 0033
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   ANS 0034
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   ANS 0035
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     ANS 0036
      COMMON /MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),        ANS 0037
     X WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)     ANS 0038
C                                                                       ANS 0039
      COMMON /MDATAX/ WMOLXT(MMOLX,LAYDIM)                              ANS 0040
C                                                                       ANS 0041
      COMMON /MDATA1/ WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),            ANS 0042
     X WNH3(LAYDIM),WAIR(LAYDIM),WHNO3(LAYDIM)                          ANS 0043
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          ANS 0044
      COMMON /CARD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB  ANS 0045
     1  ,MODTRN                                                         ANS 0046
      LOGICAL MODTRN                                                    ANS 0047
      COMMON /CARD1A/ M4,M5,M6,MDEF,IRD1,IRD2                           ANS 0048
      COMMON /CARD1B/ JUNIT(15),WMOL(12),WAIR1,JLOW                     ANS 0049
      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   ANS 0050
     1    RAINRT                                                        ANS 0051
      INTEGER NCRALT,NCRSPC                                             ANS 0052
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    ANS 0053
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      ANS 0054
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       ANS 0055
      COMMON /CARD2D/ IREG(4),ALTB(4),IREGC(4)                          ANS 0056
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     ANS 0057
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                ANS 0058
      COMMON /MART/ RHH                                                 ANS 0059
      REAL ZM,PM,TM,RFNDX,DENSTY                                        ANS 0060
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    ANS 0061
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               ANS 0062
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)            ANS 0063
      COMMON /MDLZ/HMDLZ(8)                                             ANS 0064
      CHARACTER*4  HZ         ,SEASN     ,VULCN     ,BLANK,             ANS 0065
     X            HMET        ,HMODEL     ,HTRRAD                       ANS 0066
      COMMON /TITL/ HZ(5,16),SEASN(5,2),VULCN(5,8),BLANK,               ANS 0067
     X HMET(5,2),HMODEL(5,8),HTRRAD(6,4)                                ANS 0068
      CHARACTER*4  AHOL1,AHOL2   ,AHOL3   , AHLVSA   ,AHUS              ANS 0069
      DIMENSION AHOL1(5),AHOL2(5),AHOL3(5), AHLVSA(5),AHUS(5)           ANS 0070
      DIMENSION  ITY1(LAYDIM),IH1(LAYDIM),IS1(LAYDIM),IVL1(LAYDIM),     ANS 0071
     1    ZGN(LAYDIM)                                                   ANS 0072
      CHARACTER*4  AHAHOL   ,HHOL                                       ANS 0073
      DIMENSION AHAHOL(5,13),HHOL(5)                                    ANS 0074
      DIMENSION CLDTOP(10),AHAST(LAYDIM)                                ANS 0075
C                                                                       ANS 0076
C                                                                       ANS 0077
      CHARACTER*1 JCHAR,BL,JCHARX                                       ANS 0078
      COMMON /CRD1BX/JUNITX,WMOLX(MMOLX)                                ANS 0079
      DIMENSION  JCHAR(15)                                              ANS 0080
      DATA AHLVSA/'VSA ','DEFI','NED ','    ','    '/                   ANS 0081
      DATA  AHUS /'USER',' DEF','INED','    ','    '/                   ANS 0082
      DATA AHAHOL/                                                      ANS 0083
     1 'CUMU','LUS ','    ','    ','    ',                              ANS 0084
     2 'ALTO','STRA','TUS ','    ','    ',                              ANS 0085
     3 'STRA','TUS ','    ','    ','    ',                              ANS 0086
     4 'STRA','TUS ','STRA','TO C','UM  ',                              ANS 0087
     5 'NIMB','OSTR','ATUS','    ','    ',                              ANS 0088
     6 'DRIZ','ZLE ','2.0 ','MM/H','R   ',                              ANS 0089
     7 'LT R','AIN ','5.0 ','MM/H','R   ',                              ANS 0090
     8 'MOD ','RAIN',' 12.','5 MM','/HR ',                              ANS 0091
     9 'HEAV','Y RA','IN 2','5 MM','/HR ',                              ANS 0092
     A 'EXTR','EME ','RAIN',' 75M','M/HR',                              ANS 0093
     B 'USER',' ATM','OSPH','ERE ','    ',                              ANS 0094
     C 'USER',' RAI','N NO',' CLO','UD  ',                              ANS 0095
     D 'CIRR','US C','LOUD','    ','    ' /                             ANS 0096
      DATA CLDTOP / 3.,3.,1.,2.,.66,1.,.66,.66,3.,3./                   ANS 0097
C                                                                       ANS 0098
C     F(A) IS SATURATED WATER WAPOR DENSITY AT TEMP T,A=TZERO/T         ANS 0099
C                                                                       ANS 0100
      F(A)=EXP(18.9766-14.9595*A-2.43882*A*A)*A                         ANS 0101
C                                                                       ANS 0102
C                                                                       ANS 0103
C     ZM  COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                     ANS 0104
C     ZK  ALTITUDE FOR CLOUD                                            ANS 0105
C     ZSC ALTITUDE FOR AEROSOLS  GNDALT                                 ANS 0106
C                                                                       ANS 0107
C                                                                       ANS 0108
      IREGC(1) = 0                                                      ANS 0109
      IREGC(2) = 0                                                      ANS 0110
      IREGC(3) = 0                                                      ANS 0111
      IREGC(4) = 0                                                      ANS 0112
      ICONV =1                                                          ANS 0113
      IRD0 = 1                                                          ANS 0114
      IF((MODEL .EQ. 0.) . OR. (MODEL .EQ. 7)) THEN                     ANS 0115
            IF(IM. NE. 1) RETURN                                        ANS 0116
      ENDIF                                                             ANS 0117
      IF((MODEL .GT. 0.) .AND. (MODEL .LT. 7)) IRD0 = 0                 ANS 0118
      IF((IRD0  .EQ. 1)  .AND. (IVSA.EQ.1)   ) THEN                     ANS 0119
           IRD0 = 0                                                     ANS 0120
           IRD1 = 0                                                     ANS 0121
           IRD2 = 0                                                     ANS 0122
           ICONV =0                                                     ANS 0123
           ML = ML + 10 - JLOW                                          ANS 0124
           IF(ML.GT.LAYDIM)WRITE(IPR,910)                               ANS 0125
           IF(ML.GT.LAYDIM)ML=LAYDIM                                    ANS 0126
           ZVSA(10)=ZVSA(9)+0.01                                        ANS 0127
           RHVSA(10)=0.                                                 ANS 0128
           AHVSA(10)=0.                                                 ANS 0129
           IHVSA(10)=0                                                  ANS 0130
           IF(MODEL.EQ.0)WRITE (IPR,900)                                ANS 0131
900   FORMAT('   ERROR MODEL EQ 0 AND ARMY MODEL CANNOT MIX')           ANS 0132
           IF(MODEL.EQ.0)STOP                                           ANS 0133
910   FORMAT('  ERROR ML GT 24 AND ARMY MODEL TOP LAYER TRUNCATED')     ANS 0134
      ENDIF                                                             ANS 0135
      ICL=0                                                             ANS 0136
      IDSR=0                                                            ANS 0137
      IF(ICLD.EQ.18 .OR. ICLD.EQ.19)THEN                                ANS 0138
          CALL CIRR18                                                   ANS 0139
          CLDD=.1*CTHIK                                                 ANS 0140
          CLD0=CALT-.5*CLDD                                             ANS 0141
          IF(CLD0.LE.GNDALT)CLD0=GNDALT                                 ANS 0142
          CLD1=CLD0+CLDD                                                ANS 0143
          CLD2=CLD0+CTHIK                                               ANS 0144
          CLD3=CLD1+CTHIK                                               ANS 0145
      ENDIF                                                             ANS 0146
      IF(MODEL.GE.1 .AND. MODEL.LE.6)CALL FLAYZ(ML,ICLD,GNDALT,IVSA)    ANS 0147
      JPRT=1                                                            ANS 0148
      IF(MODEL.EQ.0  .OR. MODEL.EQ.7) JPRT=0                            ANS 0149
      IF(IVSA.EQ.1) JPRT=0                                              ANS 0150
      IF(ICLD.GE.1 .AND.ICLD.LT.20) JPRT=0                              ANS 0151
      IF(RAINRT.GT.0.) JPRT=0                                           ANS 0152
      IF(GNDALT.GT.0.) JPRT=0                                           ANS 0153
C                                                                       ANS 0154
      DO 6 II = 1,4                                                     ANS 0155
6     ALTB(II) = 0.                                                     ANS 0156
      T0=273.15                                                         ANS 0157
      IC1=1                                                             ANS 0158
      N=7                                                               ANS 0159
      IF(IVULCN.LE.0) IVULCN=1                                          ANS 0160
      IF(ISEASN.LE.0) ISEASN=1                                          ANS 0161
      IF(JPRT.EQ.0) THEN                                                ANS 0162
         WRITE(IPR,925) MODEL,ICLD                                      ANS 0163
         IF(MODEL .EQ.7)WRITE(IPR,94)                                   ANS 0164
 94      FORMAT(/,10X,' MODEL 0 / 7 USER INPUT DATA ',//)               ANS 0165
C                                                                       ANS 0166
      ENDIF                                                             ANS 0167
C                                                                       ANS 0168
C                                                                       ANS 0169
      DO 100 K=1,ML                                                     ANS 0170
C                                                                       ANS 0171
C    LOOP OVER LAYERS                                                   ANS 0172
C                                                                       ANS 0173
      RH = 0.                                                           ANS 0174
      WH(K)  =0.                                                        ANS 0175
      WO(K)  =0.                                                        ANS 0176
C                                                                       ANS 0177
      IHA1=0                                                            ANS 0178
      ICLD1=0                                                           ANS 0179
      ISEA1=0                                                           ANS 0180
      IVUL1=0                                                           ANS 0181
      VIS1=0.                                                           ANS 0182
      AHAZE=0.                                                          ANS 0183
      EQLWCZ=0.                                                         ANS 0184
      RRATZ=0.                                                          ANS 0185
      ICHR = 0                                                          ANS 0186
      HAZECZ=0.                                                         ANS 0187
      DENSTY(16,K)=0.                                                   ANS 0188
C     NEW                                                               ANS 0189
C                                                                       ANS 0190
      WAIR(K) = 0                                                       ANS 0191
      WCO2(K) = 0                                                       ANS 0192
      WCO(K)  = 0                                                       ANS 0193
      WCH4(K) = 0                                                       ANS 0194
      WN2O(K) = 0                                                       ANS 0195
      WO2(K)  = 0                                                       ANS 0196
      WNH3(K) = 0                                                       ANS 0197
      WNO (K) = 0                                                       ANS 0198
      WNO2(K) = 0                                                       ANS 0199
      WSO2(K) = 0                                                       ANS 0200
      WHNO3(K)= 0                                                       ANS 0201
      DO 10 KM = 1,15                                                   ANS 0202
      JCHAR(KM) = ' '                                                   ANS 0203
      IF(KM. GT. 12) GO TO 10                                           ANS 0204
      WMOL(KM) = 0.                                                     ANS 0205
10    CONTINUE                                                          ANS 0206
      JCHARX = ' '                                                      ANS 0207
C                                                                       ANS 0208
C                                                                       ANS 0209
C        PARAMETERS - JCHAR = INPUT KEY (SEE BELOW)                     ANS 0210
C                                                                       ANS 0211
C                                                                       ANS 0212
C     ***  ROUTINE ALSO ACCEPTS VARIABLE UNITS ON PRESS AND TEMP        ANS 0213
C                                                                       ANS 0214
C          SEE INPUT KEY BELOW                                          ANS 0215
C                                                                       ANS 0216
C                                                                       ANS 0217
C                                                                       ANS 0218
C     FOR MOLECULAR SPECIES ONLY                                        ANS 0219
C                                                                       ANS 0220
C       JCHAR   JUNIT                                                   ANS 0221
C                                                                       ANS 0222
C     " ",A      10    VOLUME MIXING RATIO (PPMV)                       ANS 0223
C         B      11    NUMBER DENSITY (CM-3)                            ANS 0224
C         C      12    MASS MIXING RATIO (GM(K)/KG(AIR))                ANS 0225
C         D      13    MASS DENSITY (GM M-3)                            ANS 0226
C         E      14    PARTIAL PRESSURE (MB)                            ANS 0227
C         F      15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY           ANS 0228
C         G      16     "    "     "  (TD IN T(C)) - H2O ONLY           ANS 0229
C         H      17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY     ANS 0230
C         I      18    AVAILABLE FOR USER DEFINITION                    ANS 0231
C        1-6    1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            ANS 0232
C                                                                       ANS 0233
C     ****************************************************************  ANS 0234
C     ****************************************************************  ANS 0235
C                                                                       ANS 0236
C     ***** OTHER 'JCHAR' SPECIFICATIONS -                              ANS 0237
C                                                                       ANS 0238
C       JCHAR   JUNIT                                                   ANS 0239
C                                                                       ANS 0240
C      " ",A     10    PRESSURE IN (MB)                                 ANS 0241
C          B     11       "     "  (ATM)                                ANS 0242
C          C     12       "     "  (TORR)                               ANS 0243
C         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            ANS 0244
C                                                                       ANS 0245
C      " ",A     10    AMBIENT TEMPERATURE IN DEG(K)                    ANS 0246
C          B     11       "         "       "  " (C)                    ANS 0247
C          C     12       "         "       "  " (F)                    ANS 0248
C         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            ANS 0249
C                                                                       ANS 0250
C     ***** DEFINITION OF "DEFAULT" CHOICES FOR PROFILE SELECTION ***** ANS 0251
C                                                                       ANS 0252
C      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL          ANS 0253
C      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS   ANS 0254
C      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE                ANS 0255
C                                                                       ANS 0256
C     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)                         ANS 0257
C                                                                       ANS 0258
C      FOR MOLECULES 8-35, ONLY US STD PROFILES ARE AVIALABLE           ANS 0259
C      THEREFORE, WHEN  'JCHAR(K) = 1-5', JCHAR(K) WILL BE RESET TO 6   ANS 0260
C                                                                       ANS 0261
C                                                                       ANS 0262
      IF(IRD0 .EQ. 1) THEN                                              ANS 0263
          READ(IRD,80)ZM(K),P(K),T(K),WMOL(1),WMOL(2),WMOL(3),          ANS 0264
     X     (JCHAR(KM),KM=1,15),JCHARX                                   ANS 0265
80         FORMAT ( F10.3,5E10.3,15A1,A1)                               ANS 0266
          WRITE(IPR,81)ZM(K),P(K),T(K),WMOL(1),WMOL(2),WMOL(3),         ANS 0267
     X     JCHAR,JCHARX                                                 ANS 0268
81         FORMAT (F10.3,1P5E10.3,10X,15A1,A1)                          ANS 0269
      ENDIF                                                             ANS 0270
      IF(IRD1 .EQ. 1) THEN                                              ANS 0271
         READ(IRD,83)(WMOL(KM),KM=4,12)                                 ANS 0272
 83      FORMAT((8E10.3))                                               ANS 0273
         WRITE(IPR,84)(WMOL(KM),KM=4,12)                                ANS 0274
 84      FORMAT((1P8E10.3))                                             ANS 0275
         IF (MDEF .EQ. 2) THEN                                          ANS 0276
C           THE EXTRA SPECIES (I. E., THE WMOLX SPECIES) WILL BE READ IFANS 0277
C           MDEF .EQ. 2 AND IRD1 .EQ. 1.                                ANS 0278
C           THEY WILL BE READ 8 SPECIES AT A TIME UNTIL ALL SPECIES ARE ANS 0279
C           READ.                                                       ANS 0280
            NROWX = INT(REAL(NSPECX)/REAL(8))+1                         ANS 0281
            MROWX = MOD(NSPECX,8)                                       ANS 0282
            IF (MROWX .EQ. 0) THEN                                      ANS 0283
               MAXROW =  NROWX-1                                        ANS 0284
            ELSE                                                        ANS 0285
               MAXROW = NROWX                                           ANS 0286
            ENDIF                                                       ANS 0287
            DO 85 KROWX = 1, MAXROW                                     ANS 0288
               IBEGX = (KROWX-1)*8+1                                    ANS 0289
               IENDX = IBEGX+8-1                                        ANS 0290
               IF (KROWX .EQ. NROWX) IENDX = IBEGX + MROWX-1            ANS 0291
               READ(IRD,83)(WMOLX(KMX),KMX=IBEGX,IENDX)                 ANS 0292
               WRITE(IPR,84)(WMOLX(KMX),KMX=IBEGX,IENDX)                ANS 0293
 85         CONTINUE                                                    ANS 0294
         ENDIF                                                          ANS 0295
      ENDIF                                                             ANS 0296
C                                                                       ANS 0297
C                                                                       ANS 0298
C     AHAZE =  AEROSOL VISIBLE EXTINCTION COFF (KM-1)                   ANS 0299
C     AT A WAVELENGTH OF 0.55 MICROMETERS                               ANS 0300
C                                                                       ANS 0301
C     EQLWCZ=LIQUID WATER CONTENT (PPMV) AT ALT Z                       ANS 0302
C            FOR AEROSOL, CLOUD OR FOG MODELS                           ANS 0303
C                                                                       ANS 0304
C     RRATZ=RAIN RATE (MM/HR) AT ALT Z                                  ANS 0305
C                                                                       ANS 0306
C     IHA1 AEROSOL MODEL USED FOR SPECTRAL DEPENDENCE OF EXTINCTION     ANS 0307
C                                                                       ANS 0308
C     IVUL1 STRATOSPHERIC AEROSOL MODEL USED FOR SPECTRAL DEPENDENCE    ANS 0309
C     OF EXT AT Z                                                       ANS 0310
C                                                                       ANS 0311
C     ICLD1 CLOUD MODEL USED FOR SPECTRAL DEPENDENCE OF EXT AT Z        ANS 0312
C                                                                       ANS 0313
C     ONLY ONE OF IHA1,ICLD1  OR IVUL1 IS ALLOWED                       ANS 0314
C     IHA1 NE 0 OTHERS IGNORED                                          ANS 0315
C     IHA1 EQ 0 AND ICLD1 NE 0 USE ICLD1                                ANS 0316
C                                                                       ANS 0317
C     IF AHAZE AND EQLWCZ ARE BOUTH ZERO                                ANS 0318
C        DEFAULT PROFILE ARE LOADED FROM IHA1,ICLD1,IVUL1               ANS 0319
C     ISEA1 = AEROSOL SEASON CONTROL FOR ALTITUDE Z                     ANS 0320
C                                                                       ANS 0321
C     ICHR  CHANGE AEROSOL PROFILE REGION FOR IHA1 = 7                  ANS 0322
C                                                                       ANS 0323
      IF(IRD2 .EQ. 1) THEN                                              ANS 0324
           READ(IRD,82)    AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,   ANS 0325
     X ICHR                                                             ANS 0326
           WRITE(IPR,82)    AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,  ANS 0327
     X ICHR                                                             ANS 0328
82         FORMAT(10X,3F10.3,5I5)                                       ANS 0329
      ELSE                                                              ANS 0330
           RRATZ = RAINRT                                               ANS 0331
           IF(ZM(K).GT.6.)RRATZ = 0                                     ANS 0332
      ENDIF                                                             ANS 0333
      ICLDS = ICLD1                                                     ANS 0334
      IF( ICHR.EQ. 1) THEN                                              ANS 0335
         IF(IHA1. EQ. 0) THEN                                           ANS 0336
             IF(ICLD1. NE. 11) ICHR = 0                                 ANS 0337
         ELSE                                                           ANS 0338
             IF(IHA1 . NE.  7) ICHR = 0                                 ANS 0339
         ENDIF                                                          ANS 0340
      ENDIF                                                             ANS 0341
      IF(MODEL .EQ. 0) THEN                                             ANS 0342
           HMDLZ(1) =   ZM(K)                                           ANS 0343
           HMDLZ(2) =    P(K)                                           ANS 0344
           HMDLZ(3) =    T(K)                                           ANS 0345
           HMDLZ(4) = WMOL(1)                                           ANS 0346
           HMDLZ(5) = WMOL(2)                                           ANS 0347
           HMDLZ(6) = WMOL(3)                                           ANS 0348
           HMDLZ(7) = AHAZE                                             ANS 0349
      ENDIF                                                             ANS 0350
      DO 12 KM = 1,15                                                   ANS 0351
         JUNIT(KM) = JOU(JCHAR(KM))                                     ANS 0352
 12   CONTINUE                                                          ANS 0353
      JUNITX = JOU(JCHARX)                                              ANS 0354
      IF(IRD0 .EQ. 0) THEN                                              ANS 0355
          JUNIT(1) = M1                                                 ANS 0356
          JUNIT(2) = M1                                                 ANS 0357
          JUNIT(3) = M2                                                 ANS 0358
          JUNIT(4) =  6                                                 ANS 0359
          JUNIT(5) = M3                                                 ANS 0360
          JUNIT(6) = M5                                                 ANS 0361
          JUNIT(7) = M6                                                 ANS 0362
          JUNIT(8) = M4                                                 ANS 0363
          JUNIT(9) =  6                                                 ANS 0364
          JUNIT(10)=  6                                                 ANS 0365
          JUNIT(11)=  6                                                 ANS 0366
          JUNIT(12)=  6                                                 ANS 0367
          JUNIT(13)=  6                                                 ANS 0368
          JUNIT(14)=  6                                                 ANS 0369
          JUNIT(15)=  6                                                 ANS 0370
          JUNITX = 6                                                    ANS 0371
      ELSE                                                              ANS 0372
          BL = ' '                                                      ANS 0373
          IF((M1   .GT.0).AND.(JCHAR(1).EQ.BL))                         ANS 0374
     X    JUNIT(1) = M1                                                 ANS 0375
          IF((M1   .GT.0).AND.(JCHAR(2).EQ.BL))                         ANS 0376
     X    JUNIT(2) = M1                                                 ANS 0377
          IF((M2   .GT.0).AND.(JCHAR(3).EQ.BL))                         ANS 0378
     X    JUNIT(3) = M2                                                 ANS 0379
          IF((MDEF .GT.0).AND.(JCHAR(4).EQ.BL))                         ANS 0380
     X    JUNIT(4) = 6                                                  ANS 0381
          IF((M3   .GT.0).AND.(JCHAR(5).EQ.BL))                         ANS 0382
     X    JUNIT(5) = M3                                                 ANS 0383
          IF((M5   .GT.0).AND.(JCHAR(6).EQ.BL))                         ANS 0384
     X    JUNIT(6) = M5                                                 ANS 0385
          IF((M6   .GT.0).AND.(JCHAR(7).EQ.BL))                         ANS 0386
     X    JUNIT(7) = M6                                                 ANS 0387
          IF((M4   .GT.0).AND.(JCHAR(8).EQ.BL))                         ANS 0388
     X    JUNIT(8) = M4                                                 ANS 0389
          IF((MDEF .GT.0).AND.(JCHAR(9).EQ.BL))                         ANS 0390
     X    JUNIT(9) = 6                                                  ANS 0391
          IF((MDEF .GT.0).AND.(JCHAR(10).EQ.BL))                        ANS 0392
     X    JUNIT(10)= 6                                                  ANS 0393
          IF((MDEF .GT.0).AND.(JCHAR(11).EQ.BL))                        ANS 0394
     X    JUNIT(11)= 6                                                  ANS 0395
          IF((MDEF .GT.0).AND.(JCHAR(12).EQ.BL))                        ANS 0396
     X    JUNIT(12)= 6                                                  ANS 0397
          IF((MDEF .GT.0).AND.(JCHAR(13).EQ.BL))                        ANS 0398
     X    JUNIT(13)= 6                                                  ANS 0399
          IF((MDEF .GT.0).AND.(JCHAR(14).EQ.BL))                        ANS 0400
     X    JUNIT(14)= 6                                                  ANS 0401
          IF((MDEF .GT.0).AND.(JCHARX   .EQ.BL))                        ANS 0402
     X    JUNITX   = 6                                                  ANS 0403
      ENDIF                                                             ANS 0404
      IF(ICONV .EQ. 1) THEN                                             ANS 0405
         CALL CHECK(P(K),JUNIT(1),1)                                    ANS 0406
         CALL CHECK(T(K),JUNIT(2),2)                                    ANS 0407
         CALL DEFALT(ZM(K),P(K),T(K))                                   ANS 0408
         CALL CONVRT (P(K),T(K))                                        ANS 0409
         WH(K)    = WMOL(1)                                             ANS 0410
         WCO2(K)  = WMOL(2)                                             ANS 0411
         WO(K)    = WMOL(3)                                             ANS 0412
         WN2O(K)  = WMOL(4)                                             ANS 0413
         WCO(K)   = WMOL(5)                                             ANS 0414
         WCH4(K)  = WMOL(6)                                             ANS 0415
         WO2(K)   = WMOL(7)                                             ANS 0416
         WNO(K)   = WMOL(8)                                             ANS 0417
         WSO2(K)  = WMOL(9)                                             ANS 0418
         WNO2(K)  = WMOL(10)                                            ANS 0419
         WNH3(K)  = WMOL(11)                                            ANS 0420
         WHNO3(K) = WMOL(12)                                            ANS 0421
         WAIR(K)  = WAIR1                                               ANS 0422
         DO 500 I = 1, NSPECX                                           ANS 0423
            WMOLXT(I,K) = WMOLX(I)                                      ANS 0424
 500     CONTINUE                                                       ANS 0425
      ELSE                                                              ANS 0426
         CALL VSANSM(K,AHAZE,IHA1,ZNEW)                                 ANS 0427
         ZM(K) = ZNEW                                                   ANS 0428
      ENDIF                                                             ANS 0429
C                                                                       ANS 0430
C                                                                       ANS 0431
C     GNDALT NOT ZERO                                                   ANS 0432
C                                                                       ANS 0433
      ZSC=ZM(K)                                                         ANS 0434
      IF(GNDALT.GT.0. .AND. ZM(K).LT.6.)THEN                            ANS 0435
           ASC=6./(6.-GNDALT)                                           ANS 0436
           CON=-ASC*GNDALT                                              ANS 0437
           ZSC=ASC*ZM(K)+CON -1.E-4                                     ANS 0438
           IF(ZSC.LT.0.)ZSC=0.                                          ANS 0439
      ENDIF                                                             ANS 0440
      ZGN(K)=ZSC                                                        ANS 0441
C                                                                       ANS 0442
C                                                                       ANS 0443
      IF(ICLD1.EQ.0) ICLD1=ICLD                                         ANS 0444
      IF(ICLD1.GT.11) ICLD1=0                                           ANS 0445
      IF(IHA1.NE.0) IVUL1=0                                             ANS 0446
      IF(IHA1.NE.0) ICLD1=0                                             ANS 0447
      IF(ICLD1.NE.0) IVUL1=0                                            ANS 0448
      IF((AHAZE.NE.0.).OR.(EQLWCZ.NE.0.)) GO TO 8                       ANS 0449
C      ********        ERRATA SEPT 19                                   ANS 0450
CC    IF(RRATZ.NE.0.) GO TO 8                                           ANS 0451
C      *****           END ERRATA                                       ANS 0452
      IF((IVSA.EQ.1).AND.(ICLD1.EQ.0)) THEN                             ANS 0453
           IF(MODEL.NE.7)CALL LAYVSA(K,RH,AHAZE,IHA1,ZM)                ANS 0454
      ELSE                                                              ANS 0455
           CALL LAYCLD(K,EQLWCZ,RRATZ,ICLD1,GNDALT)                     ANS 0456
C          ***********************  ERRATA SEPT 19                      ANS 0457
           IF(RAINRT.GT.0. AND. ZM(K).LT.6.)RRATZ=RAINRT                ANS 0458
C           *********************      END ERRATA                       ANS 0459
           IF(ICLD1 . LT.  1) GO TO 8                                   ANS 0460
           IF(ICLD1 . GT. 10) GO TO 8                                   ANS 0461
C     ***************  ERRATA JUNE 21 89 NEXT CARD                      ANS 0462
           IF(ZM(K).GT.CLDTOP(ICLD1)+GNDALT)THEN                        ANS 0463
C     ***************  END ERRATA                                       ANS 0464
              RRATZ = 0.                                                ANS 0465
           ENDIF                                                        ANS 0466
      ENDIF                                                             ANS 0467
8     CONTINUE                                                          ANS 0468
      ICLDC = ICLD                                                      ANS 0469
      IF(ICLDS .NE. 0) ICLDC = ICLDS                                    ANS 0470
C                                                                       ANS 0471
      IF(ICLDS. EQ. 18 .OR. ICLDS .EQ. 19) THEN                         ANS 0472
           IF(AHAZE . GT. 0) THEN                                       ANS 0473
                DENSTY(16,K) = AHAZE                                    ANS 0474
                AHAZE    = 0.                                           ANS 0475
                GO TO 46                                                ANS 0476
           ENDIF                                                        ANS 0477
           IF(EQLWCZ .GT. 0) THEN                                       ANS 0478
                IF(ICLDS .EQ. 18) CON = 3.446E-3                        ANS 0479
                IF(ICLDS .EQ. 19) CON = 5.811E-2                        ANS 0480
                DENSTY(16,K) = EQLWCZ/CON                               ANS 0481
                EQLWCZ = 0.                                             ANS 0482
                GO TO 46                                                ANS 0483
           ENDIF                                                        ANS 0484
      ENDIF                                                             ANS 0485
      IF(ICLDC.EQ.18 .OR. ICLDC.EQ.19)THEN                              ANS 0486
          IF(ZM(K).LE.CLD0 .OR. ZM(K).GE.CLD3)THEN                      ANS 0487
              DENSTY(16,K)=0.                                           ANS 0488
          ELSEIF(ZM(K).LT.CLD1)THEN                                     ANS 0489
              DENSTY(16,K)=CEXT*(ZM(K)-CLD0)/CLDD                       ANS 0490
          ELSEIF(ZM(K).LE.CLD2)THEN                                     ANS 0491
              DENSTY(16,K)=CEXT                                         ANS 0492
          ELSE                                                          ANS 0493
              DENSTY(16,K)=CEXT*(CLD3-ZM(K))/CLDD                       ANS 0494
          ENDIF                                                         ANS 0495
      ENDIF                                                             ANS 0496
46    DENSTY(66,K)=EQLWCZ                                               ANS 0497
      IF(ICLDS.EQ.0. AND. DENSTY(66,K).EQ.0.)ICLD1=0                    ANS 0498
      DENSTY(3,K)=RRATZ                                                 ANS 0499
      IF(MODEL  .EQ. 0 .OR. MODEL .EQ. 7) THEN                          ANS 0500
C     DONT CHANGE RH                                                    ANS 0501
      ELSE                                                              ANS 0502
            IF(EQLWCZ.GT.0.0) RH=100.0                                  ANS 0503
            IF(RRATZ .GT.0.0) RH=100.0                                  ANS 0504
      ENDIF                                                             ANS 0505
      AHAST(K)=AHAZE                                                    ANS 0506
C     IHA1  IS IHAZE FOR THIS LAYER                                     ANS 0507
C     ISEA1 IS ISEASN FOR THIS LAYER                                    ANS 0508
C     IVUL1 IS IVULCN FOR THE LAYER                                     ANS 0509
      IF(ISEA1.EQ.0) ISEA1=ISEASN                                       ANS 0510
      ITYAER=IHAZE                                                      ANS 0511
      IF(IHA1.GT.0) ITYAER=IHA1                                         ANS 0512
      IF(IVUL1.GT.0) IVULCN=IVUL1                                       ANS 0513
      IF(IVUL1.LE.0) IVUL1=IVULCN                                       ANS 0514
C                                                                       ANS 0515
      IF(K.EQ.1) GO TO 23                                               ANS 0516
      IF(ICHR .EQ. 1) GO TO 20                                          ANS 0517
      IF(ICLD1.NE.IREGC(IC1))GO TO 19                                   ANS 0518
      IF(IHA1 .EQ. 0 .AND. ICLD1. EQ. 0) THEN                           ANS 0519
           IF(ZSC.GT.2.)ITYAER=6                                        ANS 0520
           IF(ZSC.GT.10.)ITYAER=IVULCN+10                               ANS 0521
           IF(ZSC.GT.30.) ITYAER=19                                     ANS 0522
           IF(ITYAER.EQ.ICH(IC1))GO TO 23                               ANS 0523
      ENDIF                                                             ANS 0524
      IF(ICLD1 .EQ. 0 .AND. IHA1.EQ. 0) GO TO 20                        ANS 0525
      N = 7                                                             ANS 0526
      IF(IC1 .GT. 1) N= IC1 + 10                                        ANS 0527
      IF(IHA1 .EQ. 0) GO TO 23                                          ANS 0528
      IF(IHA1 .NE.ICH(IC1)) GO TO 20                                    ANS 0529
      GO TO 23                                                          ANS 0530
19    IF(ICLD1 .NE. 0) THEN                                             ANS 0531
           IF(ICLD1. EQ. IREGC(1)) THEN                                 ANS 0532
               N = 7                                                    ANS 0533
               ALTB(1) = ZM(K)                                          ANS 0534
               GO TO 24                                                 ANS 0535
           ENDIF                                                        ANS 0536
           IF(IC1 .EQ. 1) GO TO 20                                      ANS 0537
           IF(ICLD1. EQ. IREGC(2)) THEN                                 ANS 0538
               N = 12                                                   ANS 0539
               ALTB(2) = ZM(K)                                          ANS 0540
               GO TO 24                                                 ANS 0541
           ENDIF                                                        ANS 0542
           IF(IC1 .EQ. 2) GO TO 20                                      ANS 0543
           IF(ICLD1. EQ. IREGC(3)) THEN                                 ANS 0544
               N = 13                                                   ANS 0545
               ALTB(3) = ZM(K)                                          ANS 0546
               GO TO 24                                                 ANS 0547
           ENDIF                                                        ANS 0548
      ELSE                                                              ANS 0549
            IF(IHA1 .EQ. 0 .AND. ICLD1. EQ. 0) THEN                     ANS 0550
                 IF(ZSC.GT.2.)ITYAER=6                                  ANS 0551
                 IF(ZSC.GT.10.)ITYAER=IVULCN+10                         ANS 0552
                 IF(ZSC.GT.30.) ITYAER=19                               ANS 0553
            ENDIF                                                       ANS 0554
           IF(ITYAER.EQ.ICH(  1))THEN                                   ANS 0555
               N = 7                                                    ANS 0556
               ALTB(1) = ZM(K)                                          ANS 0557
               GO TO 24                                                 ANS 0558
           ENDIF                                                        ANS 0559
           IF(IC1 .EQ. 1) GO TO 20                                      ANS 0560
           IF(ITYAER.EQ.ICH(  2))THEN                                   ANS 0561
               N = 12                                                   ANS 0562
               ALTB(2) = ZM(K)                                          ANS 0563
               GO TO 24                                                 ANS 0564
           ENDIF                                                        ANS 0565
           IF(IC1 .EQ. 2) GO TO 20                                      ANS 0566
           IF(ITYAER.EQ.ICH(  3))THEN                                   ANS 0567
               N = 13                                                   ANS 0568
               ALTB(3) = ZM(K)                                          ANS 0569
               GO TO 24                                                 ANS 0570
           ENDIF                                                        ANS 0571
      ENDIF                                                             ANS 0572
20    IC1=IC1+1                                                         ANS 0573
C                                                                       ANS 0574
C                                                                       ANS 0575
C                                                                       ANS 0576
      N=IC1+10                                                          ANS 0577
      IF(RH.GT.0.) RHH=RH                                               ANS 0578
      IF(IC1.LE.4) GO TO 23                                             ANS 0579
      IC1=4                                                             ANS 0580
      N=14                                                              ANS 0581
      ITYAER=ICH(IC1)                                                   ANS 0582
23    ICH(IC1)=ITYAER                                                   ANS 0583
      IREGC(IC1) = ICLD1                                                ANS 0584
      ALTB(IC1) = ZM(K)                                                 ANS 0585
C                                                                       ANS 0586
C     FOR LVSA OR CLD OR RAIN ONLY                                      ANS 0587
C                                                                       ANS 0588
24    IF (RH.GT.  0.0) THEN                                             ANS 0589
            TA=T0/T(K)                                                  ANS 0590
            WH(K)  =F(TA)*0.01*RH                                       ANS 0591
C                                                                       ANS 0592
      ENDIF                                                             ANS 0593
   40 CONTINUE                                                          ANS 0594
      DENSTY(7,K)=0.                                                    ANS 0595
      DENSTY(12,K)=0.                                                   ANS 0596
      DENSTY(13,K)=0.                                                   ANS 0597
      DENSTY(14,K)=0.                                                   ANS 0598
      DENSTY(15,K)=0.                                                   ANS 0599
      TS=273.15/T(K)                                                    ANS 0600
      WTEMP=WH(K)                                                       ANS 0601
      RELHUM(K)=0.                                                      ANS 0602
      IF(WTEMP.LE.0.) GO TO 45                                          ANS 0603
      RELHUM(K) = 100.0*WTEMP/F(TS)                                     ANS 0604
      IF(RELHUM(K).GT.100.)WRITE(IPR,930)RELHUM(K),ZM(K)                ANS 0605
      IF( RELHUM(K) .GT. 100.) RELHUM(K)=100.                           ANS 0606
      IF(RELHUM(K).LT.0.)WRITE(IPR,930)RELHUM(K),ZM(K)                  ANS 0607
930   FORMAT(' ***ERROR RELHUM ' ,E15.4,'  AT ALT  ',F12.3)             ANS 0608
      IF( RELHUM(K) .LT.   0.) RELHUM(K)=0.                             ANS 0609
45    RHH=RELHUM(K)                                                     ANS 0610
      RH=RHH                                                            ANS 0611
      IF (VIS1.LE.0.0) VIS1=VIS                                         ANS 0612
      IF (AHAZE.EQ.0.0) GO TO 47                                        ANS 0613
      DENSTY(N,K)=AHAZE                                                 ANS 0614
      IF(ITYAER.EQ.3) GO TO 47                                          ANS 0615
      IF(ITYAER.EQ.10)GO TO 47                                          ANS 0616
C     AHAZE IS IN LOWTRAN NUMBER DENSTY UNITS                           ANS 0617
      GO TO 50                                                          ANS 0618
47    CONTINUE                                                          ANS 0619
C                                                                       ANS 0620
C     AHAZE NOT INPUT OR NAVY MARITIME MODEL IS CALLED                  ANS 0621
C                                                                       ANS 0622
C     CHECK IF GNDALT NOT ZERO                                          ANS 0623
C                                                                       ANS 0624
      IF(GNDALT.GT.0. .AND. ZM(K).LT.6.)THEN                            ANS 0625
           J=IFIX(ZSC+1.0E-6)+1                                         ANS 0626
           FAC=ZSC-FLOAT(J-1)                                           ANS 0627
      ELSE                                                              ANS 0628
      J=IFIX(ZM(K)+1.0E-6)+1                                            ANS 0629
      IF(ZM(K).GE.25.)J=IFIX((ZM(K)-25.)/ 5.+26.)                       ANS 0630
      IF(ZM(K).GE.50.)J=IFIX((ZM(K)-50.)/20.+31.)                       ANS 0631
      IF(ZM(K).GE.70.)J=IFIX((ZM(K)-70.)/30.+32.)                       ANS 0632
      IF (J.GT.32) J=32                                                 ANS 0633
      FAC=ZM(K)-FLOAT(J-1)                                              ANS 0634
      IF (J.LT.26) GO TO 125                                            ANS 0635
      FAC=(ZM(K)-5.*FLOAT(J-26)-25.)/5.                                 ANS 0636
      IF(J.GE.31)FAC=(ZM(K)-50.)/20.                                    ANS 0637
      IF(J.GE.32)FAC=(ZM(K)-70.)/30.                                    ANS 0638
      ENDIF                                                             ANS 0639
125   L=J+1                                                             ANS 0640
      IF (FAC.GT.1.0) FAC=1.0                                           ANS 0641
      IF(ITYAER.EQ.3.AND.ICL.EQ.0)THEN                                  ANS 0642
           CALL MARINE(VIS1,MODEL,WSS,WHH,ICSTL,EXTC,ABSC,IC1)          ANS 0643
           IREG(IC1) = 1                                                ANS 0644
           VIS=VIS1                                                     ANS 0645
           ICL = ICL + 1                                                ANS 0646
      ENDIF                                                             ANS 0647
      IF(ITYAER.EQ.10.AND.IDSR.EQ.0)THEN                                ANS 0648
           CALL DESATT(WSS,VIS1)                                        ANS 0649
           IREG(IC1) = 1                                                ANS 0650
           VIS=VIS1                                                     ANS 0651
           IDSR = IDSR + 1                                              ANS 0652
      ENDIF                                                             ANS 0653
      IF(AHAZE.GT.0.0) GO TO 50                                         ANS 0654
      IF(IHA1.LE.0) IHA1=IHAZE                                          ANS 0655
      CALL CLDPRF(K,ICLD1,IHA1,IC1,ICH(IC1),HAZECZ)                     ANS 0656
      CALL AERPRF(J,  VIS1,HAZ1,IHA1,      ISEA1,IVUL1,NN)              ANS 0657
      CALL AERPRF(L,  VIS1,HAZ2,IHA1,      ISEA1,IVUL1,NN)              ANS 0658
      HAZE=0.                                                           ANS 0659
      IF ((HAZ1.LE.0.0).OR.(HAZ2.LE.0.0)) GO TO 48                      ANS 0660
      HAZE=HAZ1*(HAZ2/HAZ1)**FAC                                        ANS 0661
48    CONTINUE                                                          ANS 0662
      IF(DENSTY(66,K).GT.0.)HAZE=HAZECZ                                 ANS 0663
      DENSTY(N,K)=HAZE                                                  ANS 0664
50    CONTINUE                                                          ANS 0665
      ITY1(K)=ITYAER                                                    ANS 0666
      IH1(K)=IHA1                                                       ANS 0667
      IF(AHAZE.NE.0)IH1(K)=-99                                          ANS 0668
      IS1(K)=ISEA1                                                      ANS 0669
      IVL1(K)=IVUL1                                                     ANS 0670
100   CONTINUE                                                          ANS 0671
C                                                                       ANS 0672
C     END OF LOOP                                                       ANS 0673
C                                                                       ANS 0674
      IHH=ICLD                                                          ANS 0675
      IF(IHH.LE.0) IHH=12                                               ANS 0676
      IF(IHH.GT.12)IHH=12                                               ANS 0677
      IF(ICLD.GE.18)IHH=13                                              ANS 0678
      DO 105 II=1,5                                                     ANS 0679
      HHOL(II)=AHAHOL(II,IHH)                                           ANS 0680
      IF(IVSA.NE.0) HHOL(II)=AHLVSA(II)                                 ANS 0681
105   CONTINUE                                                          ANS 0682
      IF(ICLD .NE. 0) THEN                                              ANS 0683
           IF(JPRT.EQ.0)WRITE (IPR,904) HHOL                            ANS 0684
904        FORMAT(//'0 CLOUD AND OR RAIN TYPE CHOSEN IS   ',5A4)        ANS 0685
      ENDIF                                                             ANS 0686
      IF(JPRT.EQ.0)WRITE(IPR,905)                                       ANS 0687
C                                                                       ANS 0688
C 905 FORMAT(//,T7,'Z',T17,'P',T26,'T',T32,'REL H', T41,'H2O',          ANS 0689
905   FORMAT(1H1,//,T7,'Z',T17,'P',T26,'T',T32,'REL H', T41,'H2O',      ANS 0690
     1 T49,'CLD AMT',T59,'RAIN RATE', T90,'AEROSOL'/,                   ANS 0691
     2 T6,'(KM)',T16,'(MB)',T25,'(K)',T33,'(%)',T39,'(GM M-3)',T49,     ANS 0692
     3 '(GM M-3)',T59,'(MM HR-1)',T69,                                  ANS 0693
     4 'TYPE', T90,'PROFILE')                                           ANS 0694
      IF(JPRT.EQ.1) RETURN                                              ANS 0695
      DO 60 KK=1,ML                                                     ANS 0696
      DO 52 IJ=1,5                                                      ANS 0697
      AHOL1(IJ)=BLANK                                                   ANS 0698
      AHOL2(IJ)=BLANK                                                   ANS 0699
52    AHOL3(IJ)=BLANK                                                   ANS 0700
      ITYAER=ITY1(KK)                                                   ANS 0701
      IF(ITYAER.LE.0) ITYAER=1                                          ANS 0702
      IF(ITYAER . EQ.  16) ITYAER = 11                                  ANS 0703
      IF(ITYAER . EQ.  17) ITYAER = 11                                  ANS 0704
      IF(ITYAER . EQ.  18) ITYAER = 13                                  ANS 0705
C     ***************  ERRATA JUNE 21 89 NEXT CARD                      ANS 0706
      IF(ITYAER . EQ.  19) ITYAER = 11                                  ANS 0707
C     ***************  END ERRATA                                       ANS 0708
      IHA1=IH1(KK)                                                      ANS 0709
      ISEA1=IS1(KK)                                                     ANS 0710
      IVUL1=IVL1(KK)                                                    ANS 0711
      DO 54 IJ=1,5                                                      ANS 0712
      AHOL1(IJ)=  HZ(IJ,ITYAER)                                         ANS 0713
      IF(IVSA.EQ.1) AHOL1(IJ)=HHOL(IJ)                                  ANS 0714
      IF(DENSTY(66,KK).GT.0. .OR. DENSTY(3,KK).GT.0.)AHOL1(IJ)=HHOL(IJ) ANS 0715
      IF(IHAZE.EQ.0) AHOL1(IJ)=HHOL(IJ)                                 ANS 0716
      AHOL2(IJ)=AHUS(IJ)                                                ANS 0717
      IF(AHAST(KK).EQ.0) AHOL2(IJ)=AHOL1(IJ)                            ANS 0718
      IF(DENSTY(66,KK).GT.0. .OR. DENSTY(3,KK).GT.0.)AHOL2(IJ)=HHOL(IJ) ANS 0719
54    IF (ZGN(KK).GT.2.0) AHOL3(IJ)=SEASN(IJ,ISEA1)                     ANS 0720
60    WRITE(IPR,915)ZM(KK),P(KK),T(KK),RELHUM(KK),WH(KK),               ANS 0721
     1  DENSTY(66,KK),DENSTY(3,KK),AHOL1,AHOL2,AHOL3                    ANS 0722
915   FORMAT(2F10.3,2F8.2,1P3E10.3,1X,5A4,5A4,5A4)                      ANS 0723
      RETURN                                                            ANS 0724
C                                                                       ANS 0725
925   FORMAT(//,' MODEL ATMOSPHERE NO. ',I5,' ICLD =',I5,//)            ANS 0726
      END                                                               ANS 0727
