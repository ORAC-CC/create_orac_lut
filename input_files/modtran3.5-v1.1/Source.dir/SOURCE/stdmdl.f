      SUBROUTINE STDMDL(ICH1)                                           STD 0001
C                                                                       STD 0002
C     THIS SUBROUTINE LOADS ONE OF THE 6 STANDARD ATMOSPHERIC           STD 0003
C     PROFILES INTO COMMON /MODEL/ AND CALCULATES THE DENSITIES         STD 0004
C     OF THE VARIOUS ABSORBING GASES AND AEROSOLS.                      STD 0005
C                                                                       STD 0006
C     DECLARE INPUTS.                                                   STD 0007
      INTEGER ICH1                                                      STD 0008
C                                                                       STD 0009
C     INCLUDE PARAMETERS                                                STD 0010
      INCLUDE 'PARAM.LST'                                               STD 0011
C                                                                       STD 0012
C     LIST COMMONS.                                                     STD 0013
      CHARACTER*8 CNAMEX                                                STD 0014
      COMMON/NAMEX/CNAMEX(MMOLX)                                        STD 0015
      REAL WMOLXT                                                       STD 0016
      COMMON/MDATAX/WMOLXT(MMOLX,LAYDIM)                                STD 0017
      REAL DNSTYX                                                       STD 0018
      COMMON/MODELX/DNSTYX(MMOLX,LAYDIM)                                STD 0019
      INTEGER KPOINT                                                    STD 0020
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     STD 0021
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   STD 0022
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   STD 0023
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     STD 0024
      REAL ZM,PM,TM,RFNDX,DENSTY                                        STD 0025
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    STD 0026
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               STD 0027
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               STD 0028
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           STD 0029
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     STD 0030
      REAL TBOUND,SALB                                                  STD 0031
      LOGICAL MODTRN                                                    STD 0032
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,               STD 0033
     1  TBOUND,SALB,MODTRN                                              STD 0034
      INTEGER IV1,IV2,IDV,IFWHM                                         STD 0035
      COMMON/CARD4/IV1,IV2,IDV,IFWHM                                    STD 0036
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2                             STD 0037
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         STD 0038
     1  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)    STD 0039
      REAL WNO,WSO2,WNO2,WNH3,WAIR,WHNO3                                STD 0040
      COMMON/MDATA1/WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),              STD 0041
     1  WNH3(LAYDIM),WAIR(LAYDIM),WHNO3(LAYDIM)                         STD 0042
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     STD 0043
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                STD 0044
      REAL EXTV,ABSV,ASYV                                               STD 0045
      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       STD 0046
C                                                                       STD 0047
C     DECLARE LOCAL VARIABLES                                           STD 0048
      INTEGER I,IX,K,KX                                                 STD 0049
      REAL CONJOE,RHZERO,PP,TT,PSS,TSS,F1,WTEMP,RHOAIR,RHOH2O,          STD 0050
     1  RHOFRN,RELH,RHLOG,DT,WO2D,STORE,CONO2,CONH2O,CONO3,CONCO2,      STD 0051
     2  CONCO,CONCH4,CONN2O,CONNH3,CONNO,CONNO2,CONSO2,PPW,AVW          STD 0052
C                                                                       STD 0053
C     DECLARE DATA.                                                     STD 0054
      INTEGER MNAER                                                     STD 0055
      REAL V550NM,PZERO,TZERO,XLOSCH,RV,CON,A0,A1,A2,B1,B2,C0,C1        STD 0056
C                                                                       STD 0057
C     CONSTANTS FOR RETRIEVING 550 NM CLOUD EXTINCTION COEFFICIENTS.    STD 0058
C       MNAER    MINIMUM AEROSOL INDEX USED IN CALL TO AEREXT.          STD 0059
C       V550NM   FREQUENCY EQUIVALENT OF 550 NM [CM-1].                 STD 0060
      DATA MNAER/6/,V550NM/18181.818/                                   STD 0061
C     XLOSCH = LOSCHMIDT'S NUMBER,MOLECULES CM-2,KM-1                   STD 0062
      DATA PZERO/1013.25/,TZERO/273.15/,XLOSCH/2.6868E24/               STD 0063
C     RV GAS CONSTANT FOR WATER IN MB/(GM M-3 K)                        STD 0064
C     CON CONVERTS WATER VAPOR FROM GM M-3 TO MOLECULES CM-2 KM-1       STD 0065
      DATA RV/4.6152E-3/,CON/3.3429E21/                                 STD 0066
C     CONSTANTS FOR INDEX OF REFRACTION, AFTER EDLEN, 1965              STD 0067
      DATA A0/83.42/,A1/185.08/,A2/4.11/,                               STD 0068
     X     B1/1.140E5/,B2/6.24E4/,C0/43.49/,C1/1.70E4/                  STD 0069
C                                                                       STD 0070
C     CONJOE=(1/XLOSCH)*1.E5*1.E-6 WITH                                 STD 0071
C        1.E5 ARISING FROM CM TO KM CONVERSION AND                      STD 0072
C        1.E-6  "       "  PPMV                                         STD 0073
C                                                                       STD 0074
      CONJOE = 3.7194E-21                                               STD 0075
C                                                                       STD 0076
C     H20 CONTINUUM IS STORED AT 296 K RHZERO IS AIR DENSITY AT 296 K   STD 0077
C     IN UNITS OF LOSCHMIDT'S                                           STD 0078
C                                                                       STD 0079
      RHZERO=(273.15/296.0)                                             STD 0080
C                                                                       STD 0081
C     LOAD ATMOSPHERE PROFILE INTO /MODEL/                              STD 0082
        DO 25 I=1,ML                                                    STD 0083
      PM(I)=P(I)                                                        STD 0084
      TM(I)=T(I)                                                        STD 0085
      PP=PM(I)                                                          STD 0086
      TT=TM(I)                                                          STD 0087
      PSS=PP/PZERO                                                      STD 0088
      TSS=TZERO/TT                                                      STD 0089
      F1=(PP/PZERO)/(TT/TZERO)                                          STD 0090
      WTEMP=WH(I)                                                       STD 0091
C     UV OZONE                                                          STD 0092
C     THE UNIT FOR O3 HAS BEEN CHANGED FROM G/M**3 TO PPMV.             STD 0093
      DENSTY(8,I)= CONJOE   *WAIR(I)    *WO(I)                          STD 0094
C     N2 CONTINUUM (FASCODE APPROACH)                                   STD 0095
      DENSTY(4,I)=0.781*(296.*PSS/TT)**2                                STD 0096
      RHOAIR = F1                                                       STD 0097
      RHOH2O = CON *WTEMP/XLOSCH                                        STD 0098
      RHOFRN = RHOAIR - RHOH2O                                          STD 0099
      DENSTY(5,I)= XLOSCH*RHOH2O**2/RHZERO                              STD 0100
C     FOREIGN BROADENED                                                 STD 0101
      DENSTY(10,I)= XLOSCH*RHOH2O*RHOFRN/RHZERO                         STD 0102
C     MOLECULAR SCATTERING                                              STD 0103
      DENSTY(6,I) = F1                                                  STD 0104
C     RELATIVE HUMIDITY WEIGHTED BY BOUNDARY LAYER AEROSOL (0 TO 2 KM)  STD 0105
C                                                                       STD 0106
C     LOG WEIGHTING OF REL HUMIDITY                                     STD 0107
C                                                                       STD 0108
      RELH = RELHUM(I)                                                  STD 0109
      IF(RELHUM(I). GT.99.) RELH = 99.                                  STD 0110
      RHLOG = ALOG(100. - RELH)                                         STD 0111
C     DENSTY(15,I)=RELHUM(I)*DENSTY(7,I)                                STD 0112
      DENSTY(15,I)=RHLOG    *DENSTY(7,I)                                STD 0113
C     DENSITY (9,I) TEMP DEP OF WATER SET IN GEO                        STD 0114
      DENSTY(9,I)=0.                                                    STD 0115
      IF(ICH1.GT.7)DENSTY(15,I)=RHLOG*DENSTY(12,I)                      STD 0116
C     HNO3 IN ATM * CM /KM                                              STD 0117
C     NEW PROFILE IS IN UNIT OF PART PER 10**6 BY VOLUME                STD 0118
C     DENSTY(11,I)= F1* HMIX(I)*1.0E-4                                  STD 0119
                     DENSTY(11,I)=F1*WHNO3(I)*1.0E-4                    STD 0120
C                                                                       STD 0121
C      O2 TEMP DEP                                                      STD 0122
      DT = TT  - 220.                                                   STD 0123
      WO2D       = CONJOE   *WAIR(I)    *WO2(I)  * PSS                  STD 0124
C                                                                       STD 0125
C     DT CAN BE NEGATIVE                                                STD 0126
C     EFFECTIVE DT CALCULATED IN TRANS                                  STD 0127
      DENSTY(1,I)  = WO2D * TT                                          STD 0128
      DENSTY(2,I)  = WO2D * DT * DT                                     STD 0129
      DENSTY(63,I) = WO2D                                               STD 0130
C                                                                       STD 0131
C   NEW  MICROWAVE TEMP  RAIN                                           STD 0132
C                                                                       STD 0133
      DENSTY(61,I)= 0.                                                  STD 0134
      DENSTY(62,I)= 0.                                                  STD 0135
      IF(DENSTY(3,I).GT.0.)THEN                                         STD 0136
           DENSTY(3,I)=DENSTY(3,I)**.63                                 STD 0137
           DENSTY(61,I)= DENSTY(3,I) * T(I)                             STD 0138
           DENSTY(62,I)= 1.                                             STD 0139
      ENDIF                                                             STD 0140
C     MODTRAN DENSITIES. I E ACTUAL DENISITIES.  IN LOWTRAN             STD 0141
C     SOME DENSITIES ARE SCALED.  BUT WE NEED THE ACTUAL ONES FOR UV-VISSTD 0142
C     STUFF.  THEREFORE DENSITY ARRAY'S FIRST DIMENSION WAS EXTENDED BY STD 0143
C     DENSITY(64,I) WILL HOLD NO2                                       STD 0144
C     DENSITY(65,I)WILL HOLD SO2                                        STD 0145
      DENSTY(64,I)=CONJOE*WAIR(I)*WNO2(I)                               STD 0146
      DENSTY(65,I)=CONJOE*WAIR(I)*WSO2(I)                               STD 0147
C                                                                       STD 0148
C     WHEN MODERATE RESOLUTION OPTION IS USED, ACTUAL DENSITIES ARE     STD 0149
C     NEEDED.  FOR ALL SPECIES, DENSITIES ARE STORED IN AMAGATS-CM/KM   STD 0150
C     WHERE 1 AMAGAT = 1 ATM AT STP.  FOR WATER, THE G/CM2 VALUE,       STD 0151
C     .1*WH(I), IS CONVERTED TO AMAGAT-CM BY MULTIPLYING BY             STD 0152
C     RT/(MOL WT), WHICH EQUALS 1244.2 CM3-ATM/GM/K.                    STD 0153
C                                                                       STD 0154
      CONO2 = CONJOE   *WAIR(I)  *WO2 (I)                               STD 0155
      IF(MODTRN)THEN                                                    STD 0156
          DENSTY(17,I)=WH(I)*124.42                                     STD 0157
          DENSTY(18,I) =0                                               STD 0158
          DENSTY(19,I) =0                                               STD 0159
          DENSTY(20,I) =0                                               STD 0160
          DENSTY(21,I) =0                                               STD 0161
          DENSTY(22,I) =0                                               STD 0162
          DENSTY(23,I) =0                                               STD 0163
          DENSTY(24,I) =0                                               STD 0164
          DENSTY(25,I) =0                                               STD 0165
          DENSTY(26,I) =0                                               STD 0166
          DENSTY(27,I) =0                                               STD 0167
          DENSTY(28,I) =0                                               STD 0168
          DENSTY(29,I) =0                                               STD 0169
          DENSTY(30,I) =0                                               STD 0170
          STORE=CONJOE*WAIR(I)                                          STD 0171
          DENSTY(31,I)=STORE*WO(I)                                      STD 0172
          DENSTY(32,I) =0                                               STD 0173
          DENSTY(33,I) =0                                               STD 0174
          DENSTY(34,I) =0                                               STD 0175
          DENSTY(35,I) =0                                               STD 0176
          DENSTY(36,I)=STORE*WCO2(I)                                    STD 0177
C         CO2  CONTINUM                                                 STD 0178
C         CO2 CONTINUM NOW THROUGH BAND MODEL                           STD 0179
C                                                                       STD 0180
C         CO2 = DENC *WCO2(I)* 1.0E-6                                   STD 0181
CC        DENSTY(37,I)=CO2 * PSS * 1.E-15 * 296./TT                     STD 0182
          DENSTY(37,I) =0                                               STD 0183
          DENSTY(38,I) =0                                               STD 0184
          DENSTY(39,I) =0                                               STD 0185
          DENSTY(40,I) =0                                               STD 0186
          DENSTY(41,I) =0                                               STD 0187
          DENSTY(42,I) =0                                               STD 0188
          DENSTY(43,I) =0                                               STD 0189
          DENSTY(44,I)=STORE*WCO(I)                                     STD 0190
          DENSTY(45,I) =0                                               STD 0191
          DENSTY(46,I)=STORE*WCH4(I)                                    STD 0192
          DENSTY(47,I)=STORE*WN2O(I)                                    STD 0193
          DENSTY(48,I) =0                                               STD 0194
          DENSTY(49,I) =0                                               STD 0195
          DENSTY(50,I)=STORE*WO2(I)                                     STD 0196
          DENSTY(51,I)=CONO2 *PSS**0.9353*TSS**( 0.1936)                STD 0197
CC        DENSTY(51,I) =0                                               STD 0198
          DENSTY(52,I)=STORE*WNH3(I)                                    STD 0199
          DENSTY(53,I) =0                                               STD 0200
          DENSTY(54,I)=STORE*WNO(I)                                     STD 0201
          DENSTY(55,I)=STORE*WNO2(I)                                    STD 0202
          DENSTY(56,I)=STORE*WSO2(I)                                    STD 0203
          DENSTY(57,I) =0                                               STD 0204
C                                                                       STD 0205
C         STUFF THE DNSTYX ARRAY WITH THE PROFILE INFO FOR THE          STD 0206
C         EXTRA SPECIES. THAT IS THE ADDITIONAL SPECIES BEYOND          STD 0207
C         MODTRAN'S 12 REGULAR SPECIES.                                 STD 0208
          DO 400 IX = 1,NSPECX                                          STD 0209
             DNSTYX(IX,I)=STORE*WMOLXT(IX,I)                            STD 0210
 400      CONTINUE                                                      STD 0211
C                                                                       STD 0212
      ELSE                                                              STD 0213
C  --- FOR H2O -----                                                    STD 0214
      CONH2O=WH(I)  *.1                                                 STD 0215
      DENSTY(17,I)=CONH2O*PSS**0.9810*TSS**( 0.3324)                    STD 0216
      DENSTY(18,I)=CONH2O*PSS**1.1406*TSS**(-2.6343)                    STD 0217
      DENSTY(19,I)=CONH2O*PSS**0.9834*TSS**(-2.5294)                    STD 0218
      DENSTY(20,I)=CONH2O*PSS**1.0443*TSS**(-2.4359)                    STD 0219
      DENSTY(21,I)=CONH2O*PSS**0.9681*TSS**(-1.9537)                    STD 0220
      DENSTY(22,I)=CONH2O*PSS**0.9555*TSS**(-1.5378)                    STD 0221
      DENSTY(23,I)=CONH2O*PSS**0.9362*TSS**(-1.6338)                    STD 0222
      DENSTY(24,I)=CONH2O*PSS**0.9233*TSS**(-0.9398)                    STD 0223
      DENSTY(25,I)=CONH2O*PSS**0.8658*TSS**(-0.1034)                    STD 0224
      DENSTY(26,I)=CONH2O*PSS**0.8874*TSS**(-0.2576)                    STD 0225
      DENSTY(27,I)=CONH2O*PSS**0.7982*TSS**( 0.0588)                    STD 0226
      DENSTY(28,I)=CONH2O*PSS**0.8088*TSS**( 0.2816)                    STD 0227
      DENSTY(29,I)=CONH2O*PSS**0.6642*TSS**( 0.2764)                    STD 0228
      DENSTY(30,I)=CONH2O*PSS**0.6656*TSS**( 0.5061)                    STD 0229
C  --- FOR O3                                                           STD 0230
      CONO3 = CONJOE   *WAIR(I)    *WO(I)                               STD 0231
      DENSTY(31,I)=CONO3 *PSS**0.4200*TSS**( 1.3909)                    STD 0232
      DENSTY(32,I)=CONO3 *PSS**0.4221*TSS**( 0.7678)                    STD 0233
      DENSTY(33,I)=CONO3 *PSS**0.3739*TSS**( 0.1225)                    STD 0234
      DENSTY(34,I)=CONO3 *PSS**0.1770*TSS**( 0.9827)                    STD 0235
      DENSTY(35,I)=CONO3 *PSS**0.3921*TSS**( 0.1942)                    STD 0236
C  --- FOR CO2                                                          STD 0237
      CONCO2= CONJOE   *WAIR(I)  *WCO2(I)                               STD 0238
      DENSTY(36,I)=CONCO2*PSS**0.6705*TSS**(-2.2560)                    STD 0239
      DENSTY(37,I)=CONCO2*PSS**0.7038*TSS**(-5.0768)                    STD 0240
      DENSTY(38,I)=CONCO2*PSS**0.7258*TSS**(-1.6740)                    STD 0241
      DENSTY(39,I)=CONCO2*PSS**0.6982*TSS**(-1.8107)                    STD 0242
      DENSTY(40,I)=CONCO2*PSS**0.8867*TSS**(-0.5327)                    STD 0243
      DENSTY(41,I)=CONCO2*PSS**0.7883*TSS**(-1.3244)                    STD 0244
      DENSTY(42,I)=CONCO2*PSS**0.6899*TSS**(-0.8152)                    STD 0245
      DENSTY(43,I)=CONCO2*PSS**0.6035*TSS**( 0.6026)                    STD 0246
C  --- FOR CO                                                           STD 0247
      CONCO = CONJOE   *WAIR(I)  *WCO (I)                               STD 0248
      DENSTY(44,I)=CONCO *PSS**0.7589*TSS**( 0.6911)                    STD 0249
      DENSTY(45,I)=CONCO *PSS**0.9267*TSS**( 0.1716)                    STD 0250
C  --- FOR CH4                                                          STD 0251
      CONCH4= CONJOE   *WAIR(I)  *WCH4(I)                               STD 0252
      DENSTY(46,I)=CONCH4*PSS**0.7139*TSS**(-0.4185)                    STD 0253
C  --- FOR N2O                                                          STD 0254
      CONN2O= CONJOE   *WAIR(I)  *WN2O(I)                               STD 0255
      DENSTY(47,I)=CONN2O*PSS**0.3783*TSS**( 0.9399)                    STD 0256
      DENSTY(48,I)=CONN2O*PSS**0.7203*TSS**(-0.1836)                    STD 0257
      DENSTY(49,I)=CONN2O*PSS**0.7764*TSS**( 1.1931)                    STD 0258
C  --- FOR O2                                                           STD 0259
      CONO2 = CONJOE   *WAIR(I)  *WO2 (I)                               STD 0260
      DENSTY(50,I)=CONO2 *PSS**1.1879*TSS**( 2.9738)                    STD 0261
      DENSTY(51,I)=CONO2 *PSS**0.9353*TSS**( 0.1936)                    STD 0262
C  --- FOR NH3                                                          STD 0263
      CONNH3= CONJOE   *WAIR(I)  *WNH3(I)                               STD 0264
      DENSTY(52,I)=CONNH3*PSS**0.8023*TSS**(-0.9111)                    STD 0265
      DENSTY(53,I)=CONNH3*PSS**0.6968*TSS**( 0.3377)                    STD 0266
C  --- FOR NO                                                           STD 0267
      CONNO = CONJOE   *WAIR(I)  *WNO (I)                               STD 0268
      DENSTY(54,I)=CONNO *PSS**0.5265*TSS**(-0.4702)                    STD 0269
C  --- FOR NO2                                                          STD 0270
      CONNO2= CONJOE   *WAIR(I)  *WNO2(I)                               STD 0271
      DENSTY(55,I)=CONNO2*PSS**0.3956*TSS**(-0.0545)                    STD 0272
C  --- FOR SO2                                                          STD 0273
      CONSO2= CONJOE   *WAIR(I)  *WSO2(I)                               STD 0274
      DENSTY(56,I)=CONSO2*PSS**0.2943*TSS**( 1.2316)                    STD 0275
      DENSTY(57,I)=CONSO2*PSS**0.2135*TSS**( 0.0733)                    STD 0276
C                                                                       STD 0277
C     STUFF THE DNSTYX ARRAY WITH THE PROFILE INFO FOR THE              STD 0278
C     EXTRA SPECIES. THAT IS, THE ADDITIONAL SPECIES BEYOND             STD 0279
C     MODTRAN'S 12 REGULAR SPECIES.                                     STD 0280
C     DNSTYX ARRAY IS STUFFED THE SAME WAY IN LOWTRAN AS IN MODTRAN.    STD 0281
      DO 450 IX = 1, NSPECX                                             STD 0282
         DNSTYX(IX,I)=CONJOE*WAIR(I)*WMOLXT(IX,I)                       STD 0283
 450  CONTINUE                                                          STD 0284
      ENDIF                                                             STD 0285
C***********************************************************************STD 0286
C   HERZBERG CONTINUUM PRESSURE DEPENDENCE CALCULATION, SHARDANAND 1977 STD 0287
C      AND   YOSHINO ET AL 1988                                         STD 0288
C                                                                       STD 0289
C     OXYGEN                                                            STD 0290
C                                                                       STD 0291
C      ********       ERRATA JULY 25                                    STD 0292
C     DENSTY(58,I)=(1.+.73*F1)*CONO2                                    STD 0293
      DENSTY(58,I)=(1.+.83*F1)*CONO2                                    STD 0294
C       ********      END  ERRATA                                       STD 0295
      DENSTY(59,I) = 0.                                                 STD 0296
      DENSTY(60,I) = 0.                                                 STD 0297
C     THESE (DENSTY(55,I) FOR N2O AND DENSTY(56,I) FOR SO2 ARE          STD 0298
C     MODTRAN DENSITIES. I E ACTUAL DENISITIES.  IN LOWTRAN             STD 0299
C     SOME DENSITIES ARE SCALED.  BUT WE NEED THE ACTUAL ONES FOR UV-VISSTD 0300
C     STUFF.  THEREFORE DENSITY ARRAY'S FIRST DIMENSION WAS EXTENDED BY STD 0301
C     DENSITY(64,I) WILL HOLD NO2                                       STD 0302
C     DENSITY(65,I)WILL HOLD SO2                                        STD 0303
      DENSTY(64,I)=CONJOE*WAIR(I)*WNO2(I)                               STD 0304
      DENSTY(65,I)=CONJOE*WAIR(I)*WSO2(I)                               STD 0305
C                                                                       STD 0306
C     RFNDX = REFRACTIVITY 1-INDEX OF REFRACTION FROM EDLEN, 1966       STD 0307
      PPW=RV*WTEMP*TT                                                   STD 0308
      AVW=.5*(IV1+IV2)                                                  STD 0309
      RFNDX(I)=((A0+A1/(1.-(AVW/B1)**2) +A2/(1.0-(AVW/B2)**2))*         STD 0310
     X (PP/PZERO)*(TZERO+15.0)/TT-(C0-(AVW/C1)**2)*PPW/PZERO)*1.E-6     STD 0311
25    CONTINUE                                                          STD 0312
      IF(NPR.EQ.1)RETURN                                                STD 0313
      WRITE(IPR,910)                                                    STD 0314
  910 FORMAT('1',/'  ATMOSPHERIC PROFILES',                             STD 0315
     1  //3X,'I',T10,'Z',T18,'P',T26,'T',T35,'N2',T44,'CNTMSLF',        STD 0316
     2  T54,'CNTMFRN',T62,'MOL SCAT',T75,'N-1',T83,'O3 (UV)',           STD 0317
     3  T93,'O2 (UV)',T103,'WAT DROP  ICE PART  RAIN RATE',             STD 0318
     4  /T9,'(KM)',T17,'(MB)',T25,'(K)',T43,'(  MOL/CM2 KM  )',         STD 0319
     5  T65,'(-)',T75,'(-)',T82,'(  ATM CM/KM  )',                      STD 0320
     6  T104,'(GM/M3)   (GM/M3)   (MM/HR)')                             STD 0321
      WRITE(IPR,'(/(I4,0PF9.4,F9.3,F7.1,1X,1P7E10.3,0P3F10.3))')        STD 0322
     1  (I,ZM(I),PM(I),TM(I),DENSTY(4,I),DENSTY(5,I),DENSTY(10,I),      STD 0323
     2  DENSTY(6,I),RFNDX(I),DENSTY(8,I),DENSTY(58,I),                  STD 0324
     3  (DENSTY(K,I),K=66,67),DENSTY(3,I)**1.5873,I=1,ML)               STD 0325
C                                                                       STD 0326
C     CALCULATE 550 NM CLOUD EXTINCTION COEFFICIENTS (KM-1 M3/GM).      STD 0327
      CALL AEREXT(V550NM,MNAER)                                         STD 0328
      WRITE (IPR,915)                                                   STD 0329
  915 FORMAT('1',/'  ATMOSPHERIC PROFILES',                             STD 0330
     1  //3X,'I',T10,'Z',T18,'P',T26,'T',T33,'AEROSOL 1',T43,           STD 0331
     2  'AEROSOL 2',T53,'AEROSOL 3',T63,'AEROSOL 4',T74,'AER1*RH',      STD 0332
     3  T84,'CIRRUS',T96,'RH',T103,'WAT DROP  ICE PART',                STD 0333
     4  /T9,'(KM)',T17,'(MB)',T25,'(K)',T35,'(-)',T45,'(-)',            STD 0334
     5  T55,'(-)',T65,'(-)',T76,'(-)',T85,'(-)',                        STD 0335
     6  T93,'(PERCNT)',T103,'(550nm VIS [KM-1])')                       STD 0336
      WRITE(IPR,'(/(I4,0PF9.4,F9.3,F7.1,1X,1P6E10.3,0P3F10.5))')        STD 0337
     1  (I,ZM(I),PM(I),TM(I),DENSTY(7,I),(DENSTY(K,I),K=12,16),         STD 0338
     2  RELHUM(I),EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=1,ML)     STD 0339
      IF(MODTRN)THEN                                                    STD 0340
          WRITE(IPR,'(1H1,/22H  ATMOSPHERIC PROFILES,//11H   I      Z,  STD 0341
     1      44H       P       H2O      O3       CO2      CO,            STD 0342
     2      55H       CH4      N2O      O2       NH3      NO       NO2, STD 0343
     3      19H      SO2      HNO3,/23H         (KM)    (MB)  ,         STD 0344
     4      55H(                                            ATM CM/KM , STD 0345
     5      53H                                                    ))') STD 0346
      ELSE                                                              STD 0347
          WRITE(IPR,'(A,/A,//2A,//(3A))')'1','  ATMOSPHERIC PROFILES',  STD 0348
     1      '  (IF A MOLECULE HAS MORE THAN ONE BAND, THEN',            STD 0349
     2      ' THE DATA FOR THE FIRST BAND ARE SHOWN.)',                 STD 0350
     3            '   I      Z       P       H2O      O3  ',            STD 0351
     4      '     CO2      CO       CH4      N2O      O2  ',            STD 0352
     5      '     NH3      NO       NO2      SO2      HNO3',            STD 0353
     6            '         (KM)    (MB) G/CM**2/KM  (    ',            STD 0354
     7      '                                    ATM CM/KM',            STD 0355
     8      '                                            )'             STD 0356
      ENDIF                                                             STD 0357
      WRITE(IPR,'((I4,0PF9.4,F9.3,1X,1P12E9.2))')(I,ZM(I),PM(I),        STD 0358
     1  DENSTY(17,I),DENSTY(31,I),DENSTY(36,I),DENSTY(44,I),            STD 0359
     2  DENSTY(46,I),DENSTY(47,I),DENSTY(50,I),DENSTY(52,I),            STD 0360
     3  DENSTY(54,I),DENSTY(55,I),DENSTY(56,I),DENSTY(11,I),I=1,ML)     STD 0361
      WRITE(IPR,'(1H1,/22H  ATMOSPHERIC PROFILES)')                     STD 0362
      WRITE(IPR,'(/A,13(1X,A8:),/(14X,13(1X,A8:)))')                    STD 0363
     1  '   I      Z   ',(CNAMEX(KX),KX=1,NSPECX)                       STD 0364
      WRITE(IPR,'(9X,A,50X,A,50X,A)')'(KM)   (','ATM CM/KM',')'         STD 0365
      DO 37 I = 1,ML                                                    STD 0366
          WRITE(IPR,'(I4,0PF9.4,1X,1P13E9.2:,/(14X,13E9.2:))')          STD 0367
     1      I,ZM(I),(DNSTYX(IX,I),IX=1,NSPECX)                          STD 0368
   37 CONTINUE                                                          STD 0369
      RETURN                                                            STD 0370
      END                                                               STD 0371
