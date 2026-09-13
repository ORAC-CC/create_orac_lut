      SUBROUTINE DRIVER                                                 DRV 0001
C                                                                       DRV 0002
C     INCLUDE PARAMETERS                                                DRV 0003
      INCLUDE 'PARAM.LST'                                               DRV 0004
C                                                                       DRV 0005
C     LIST COMMONS                                                      DRV 0006
      INTEGER KPOINT                                                    DRV 0007
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     DRV 0008
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   DRV 0009
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   DRV 0010
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     DRV 0011
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               DRV 0012
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           DRV 0013
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     DRV 0014
      REAL TBOUND,SALB                                                  DRV 0015
      LOGICAL MODTRN                                                    DRV 0016
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   DRV 0017
     1  MODTRN                                                          DRV 0018
      INTEGER M4,M5,M6,MDEF,IRD1,IRD2                                   DRV 0019
      COMMON/CARD1A/M4,M5,M6,MDEF,IRD1,IRD2                             DRV 0020
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA                       DRV 0021
      REAL VIS,WSS,WHH,RAINRT                                           DRV 0022
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,     DRV 0023
     1  RAINRT                                                          DRV 0024
      INTEGER NCRALT,NCRSPC                                             DRV 0025
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    DRV 0026
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      DRV 0027
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       DRV 0028
      INTEGER IREG,IREGC                                                DRV 0029
      REAL ALTB                                                         DRV 0030
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)                            DRV 0031
      INTEGER LENN                                                      DRV 0032
      REAL H1,H2,ANGLE,RANGE,BETA,RE                                    DRV 0033
      COMMON/CARD3/H1,H2,ANGLE,RANGE,BETA,RE,LENN                       DRV 0034
      INTEGER IPH                                                       DRV 0035
      REAL G                                                            DRV 0036
      COMMON/CARD3A/IPH,G                                               DRV 0037
      INTEGER IV1,IV2,IDV,IFWHM                                         DRV 0038
      COMMON/CARD4/IV1,IV2,IDV,IFWHM                                    DRV 0039
C                                                                       DRV 0040
C       PI       THE CONSTANT PI                                        DRV 0041
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       DRV 0042
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       DRV 0043
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         DRV 0044
      REAL PI,DEG,BIGNUM,BIGEXP                                         DRV 0045
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                DRV 0046
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     DRV 0047
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                DRV 0048
      REAL ZM,PM,TM,RFNDX,DENSTY                                        DRV 0049
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    DRV 0050
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               DRV 0051
      INTEGER JTURN,LJ                                                  DRV 0052
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       DRV 0053
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  DRV 0054
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      DRV 0055
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    DRV 0056
      REAL RHH                                                          DRV 0057
      COMMON/MART/RHH                                                   DRV 0058
      INTEGER NANGLS                                                    DRV 0059
      REAL ANGF,F                                                       DRV 0060
      COMMON/USRDTA/NANGLS,ANGF(50),F(4,50)                             DRV 0061
      REAL HMDLZ                                                        DRV 0062
      COMMON/MDLZ/HMDLZ(8)                                              DRV 0063
      INTEGER IHVSA                                                     DRV 0064
      REAL ZVSA,RHVSA,AHVSA                                             DRV 0065
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)              DRV 0066
      CHARACTER*4 HHAZE,HSEASN,HVULCN,BLANK,HMET,HMODEL,HTRRAD          DRV 0067
      COMMON/TITL/HHAZE(5,16),HSEASN(5,2),HVULCN(5,8),BLANK,            DRV 0068
     1  HMET(5,2),HMODEL(5,8),HTRRAD(6,4)                               DRV 0069
      REAL VSB                                                          DRV 0070
      COMMON/VSBD/VSB(10)                                               DRV 0071
C                                                                       DRV 0072
C     /PATH/                                                            DRV 0073
C       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.               DRV 0074
C       AHT     ALTITUDES AT PATH BOUNDARIES.                           DRV 0075
C       TPH     TEMPERATURE AT PATH BOUNDARIES.                         DRV 0076
C       IMAP    MAPPING FROM PATH SEGMENTS TO LAYERS.                   DRV 0077
      INTEGER IMAP                                                      DRV 0078
      REAL QTHETA,AHT,TPH                                               DRV 0079
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),TPH(LAYTWO),IMAP(LAYTWO)   DRV 0080
      REAL GNDALT                                                       DRV 0081
      COMMON/GRAUND/GNDALT                                              DRV 0082
      REAL SMALL                                                        DRV 0083
      COMMON/SMALL3/SMALL                                               DRV 0084
      LOGICAL LSAME                                                     DRV 0085
      COMMON/SOLAR/LSAME                                                DRV 0086
      REAL CO2RAT                                                       DRV 0087
      COMMON/CO2MIX/CO2RAT                                              DRV 0088
C                                                                       DRV 0089
C     DECLARE LOCAL ARRAYS                                              DRV 0090
      INTEGER ICH(4)                                                    DRV 0091
      REAL QTHETS(LAYTWO)                                               DRV 0092
C                                                                       DRV 0093
C     DECLARE LOCAL VARIABLES                                           DRV 0094
      CHARACTER CODE*1,YFLAG*1,XFLAG*1,DLIMIT*8                         DRV 0095
      LOGICAL DIS,SUN1,GROUND,LOPEN                                     DRV 0096
      DOUBLE PRECISION UANG                                             DRV 0097
      INTEGER NCORK,IRPT,NSTR,NL,I,J,JPRT,KNTRVL,LSALB,ISUN,            DRV 0098
     1  IK,ISEED,MDLSAV,ICLDSV,IPARM,IDAY,ISOURC,ITYPSV,                DRV 0099
     2  LENNSV,IPRMSV,IPHSV,IDAYSV,ICR,ISRCSV,MM1,MM2,                  DRV 0100
     3  MM3,IHVUL,IHMET,ISSGSV,MSOFF,IERROR,IMSMX,IPLOT                 DRV 0101
      REAL TBNDSV,SALBSV,CO2MX,ZCVSA,ZTVSA,ZINVSA,RAINSV,CPROB,         DRV 0102
     1  PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,RO,H1SAV,H2SAV,       DRV 0103
     2  ANGSAV,RNGSAV,BETASV,PRM1SV,PRM2SV,PRM3SV,PRM4SV,TIMESV,        DRV 0104
     3  GSV,ANGMSV,SALBV,DUMN,ALAM1,ALAM2,BENDNG,BETAH2                 DRV 0105
C@                                                                      DRV 0106
C@    HDATE AND HTIME CARRY THE DATA AND TIME AND MUST                  DRV 0107
C@    BE DOUBLE PRECISION ON A 32 BIT WORD COMPUTER                     DRV 0108
C@    DOUBLE PRECISION HDATE,HTIME                                      DRV 0109
CSSI  REAL TSEC,TSEC0                                                   DRV 0110
C                                                                       DRV 0111
C     LIST DATA                                                         DRV 0112
C       NCORK2   SAVED VALUE OF 2**NCORK, USED AS INITIALIZATION FLAG.  DRV 0113
C       LFIRST   LOGICAL FLAG, TRUE WHEN FIRST SOLAR PARAMETERS ARE     DRV 0114
C                READ IN A SERIES OF RUNS INVOLVING SOLAR PARAMETERS.   DRV 0115
C       LRDSUN   LOGICAL FLAG, SET TO TRUE IF SOLAR IRRADIANCES ARE     DRV 0116
C                USED AT A RESOLUTION OF OTHER THAN 5 CM-1.             DRV 0117
C       SUNFIL   NAME OF FILE CONTAINING SOLAR IRRADIANCES.             DRV 0118
      INTEGER NCORK2                                                    DRV 0119
      LOGICAL LFIRST,LRDSUN                                             DRV 0120
      CHARACTER*60 SUNFIL                                               DRV 0121
      DATA NCORK2/0/,LFIRST/.TRUE./,IRPT/0/,LRDSUN/.FALSE./,SUNFIL/     DRV 0122
     1  'DATA/sun_kur.dat                                            '/ DRV 0123
C                                                                       DRV 0124
C     IRD, IPR, IPU AND IPR1 ARE UNIT NUMBERS FOR INPUT, STANDARD       DRV 0125
C     OUTPUT, PLOT OUTPUT AND EXTRA OUTPUT FILES, RESPECTIVELY.         DRV 0126
      IRD=1                                                             DRV 0127
      IPR=2                                                             DRV 0128
      IPU=7                                                             DRV 0129
      IPR1=8                                                            DRV 0130
      ISCRCH=10                                                         DRV 0131
      IPLOT=12                                                          DRV 0132
      ICR=13                                                            DRV 0133
      OPEN(IRD,FILE='tape5',STATUS='OLD')                               DRV 0134
      OPEN(IPR,FILE='tape6',STATUS='UNKNOWN')                           DRV 0135
      CLOSE(IPR,STATUS='DELETE')                                        DRV 0136
      OPEN(IPR,FILE='tape6',STATUS='NEW')                               DRV 0137
      OPEN(IPU,FILE='tape7',STATUS='UNKNOWN')                           DRV 0138
      CLOSE(IPU,STATUS='DELETE')                                        DRV 0139
      OPEN(IPU,FILE='tape7',STATUS='NEW')                               DRV 0140
      OPEN(IPR1,FILE='tape8',STATUS='UNKNOWN')                          DRV 0141
      CLOSE(IPR1,STATUS='DELETE')                                       DRV 0142
      OPEN(IPR1,FILE='tape8',STATUS='NEW')                              DRV 0143
      OPEN(ISCRCH,STATUS='SCRATCH',FORM='UNFORMATTED')                  DRV 0144
      OPEN(28,FILE='DATA/refbkg',STATUS='OLD')                          DRV 0145
C                                                                       DRV 0146
C     DEFINE CONSTANTS                                                  DRV 0147
C       BIGNUM  CUTOFF FOR AVOIDING OVERFLOWS                           DRV 0148
C       BIGEXP  ARGUMENT CUTOFF FOR AVOIDING OVERFLOWS IN EXP FUNCTION  DRV 0149
C       NL      NUMBER OF BOUNDARIES IN THE STANDARD MODELS 1 TO 6      DRV 0150
      PI=2.*ASIN(1.)                                                    DRV 0151
      DEG=180./PI                                                       DRV 0152
      RANGE=0.                                                          DRV 0153
      SMALL=2.                                                          DRV 0154
      BIGNUM=1.0E35                                                     DRV 0155
      BIGEXP=87.0                                                       DRV 0156
      NL=36                                                             DRV 0157
C@                                                                      DRV 0158
C@    TIME AND DATE                                                     DRV 0159
C@      THE USER MAY WISH TO INCLUDE SUBROUTINES FDATE AND FCLOCK WHICH DRV 0160
C@      RETURN THE DATE AND TIME IN MM/DD/YY AND HH.MM.SS FORMATS,      DRV 0161
C@      RESPECTIVELY. THE REQUIRED ROUTINES FOR A CDC 6600 ARE INCLUDED DRV 0162
C@      AT THE END OF DRIVER IN COMMENT CARDS STARTING WITH "C@".       DRV 0163
C@    CALL FDATE(HDATE)                                                 DRV 0164
C@    CALL FCLOCK(HTIME)                                                DRV 0165
CSSI  CALL SECOND(TSEC0)                                                DRV 0166
C                                                                       DRV 0167
C     START CALCULATION                                                 DRV 0168
   10 CONTINUE                                                          DRV 0169
      IREG(1)=0                                                         DRV 0170
      IREG(2)=0                                                         DRV 0171
      IREG(3)=0                                                         DRV 0172
      IREG(4)=0                                                         DRV 0173
      WRITE(IPR,'(A,20X,A)')                                            DRV 0174
     1  '1','*****  MODTRAN 3.5 Version 1.1   Jan 97  *****'            DRV 0175
C@    WRITE(IPR,'(22X,A,3X,A)HDATE,HTIME                                DRV 0176
      DO 30 I=1,NAER                                                    DRV 0177
          DO 20 J=1,MXWVLN                                              DRV 0178
              EXTC(I,J)=0.                                              DRV 0179
              ABSC(I,J)=0.                                              DRV 0180
              ASYM(I,J)=0.                                              DRV 0181
   20     CONTINUE                                                      DRV 0182
   30 CONTINUE                                                          DRV 0183
      JPRT=0                                                            DRV 0184
C                                                                       DRV 0185
C     CARD 1                                                            DRV 0186
      NCORK=0                                                           DRV 0187
      READ(IRD,'(A1,I4,12I5,F8.3,F7.2)')CODE,MODEL,ITYPE,               DRV 0188
     1  IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB       DRV 0189
      IF(IMULT.NE.0 .AND. IEMSCT.NE.1 .AND. IEMSCT.NE.2)THEN            DRV 0190
          WRITE(IPR,'(A)')'0   MULTIPLE SCATTERING HAS BEEN TURNED OFF' DRV 0191
          IMULT=0                                                       DRV 0192
      ENDIF                                                             DRV 0193
      WRITE(IPR,'(2A,I1,I3,12I5,F8.3,F7.2)')                            DRV 0194
     1  '0 CARD 1  *****',CODE,NCORK,MODEL,ITYPE,IEMSCT,                DRV 0195
     2  IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB              DRV 0196
      MODTRN=.TRUE.                                                     DRV 0197
CORK  IF((IEMSCT.EQ.1 .OR. IEMSCT.EQ.2) .AND.                           DRV 0198
CORK 1  (CODE.EQ.'C' .OR. CODE.EQ.'c' .OR.                              DRV 0199
CORK 2   CODE.EQ.'K' .OR. CODE.EQ.'k'))THEN                             DRV 0200
CORK      IF(NCORK.LE.0 .OR. NCORK.GT.4)THEN                            DRV 0201
CORK          NCORK=1                                                   DRV 0202
CORK      ELSE                                                          DRV 0203
CORK          NCORK=2**NCORK                                            DRV 0204
CORK      ENDIF                                                         DRV 0205
CORK      IF(NCORK.NE.NCORK2)THEN                                       DRV 0206
CORK          NCORK2=NCORK                                              DRV 0207
CORK          CALL RDCORK(NCORK2,KNTRVL)                                DRV 0208
CORK      ENDIF                                                         DRV 0209
CORK  ELSE                                                              DRV 0210
          KNTRVL=1                                                      DRV 0211
          NCORK2=0                                                      DRV 0212
          IF(CODE.EQ.'L' .OR. CODE.EQ.'l' .OR.                          DRV 0213
     1      CODE.EQ.'F' .OR. CODE.EQ.'f')MODTRN=.FALSE.                 DRV 0214
CORK  ENDIF                                                             DRV 0215
C                                                                       DRV 0216
C     OPEN COOLING RATE FILE IF MULTIPLE SCATTERING                     DRV 0217
C     IS ON AND THE FILE IS NOT ALREADY OPEN.                           DRV 0218
CORK  IF(IMULT.NE.0 .AND. NOPRNT.NE.1)THEN                              DRV 0219
CORK      INQUIRE(IPLOT,OPENED=LOPEN)                                   DRV 0220
CORK      IF(.NOT.LOPEN)THEN                                            DRV 0221
CORK          OPEN(ICR,FILE='clrates',STATUS='UNKNOWN')                 DRV 0222
CORK          CLOSE(ICR,STATUS='DELETE')                                DRV 0223
CORK          OPEN(ICR,FILE='clrates',STATUS='NEW')                     DRV 0224
CORK      ENDIF                                                         DRV 0225
CORK  ENDIF                                                             DRV 0226
      TBNDSV=TBOUND                                                     DRV 0227
      LSALB=0                                                           DRV 0228
      SALBSV=SALB                                                       DRV 0229
      IF(SALB.LT.0)LSALB=-INT(SALB)                                     DRV 0230
C                                                                       DRV 0231
C     SET THE NUMBER OF SPECIES TREATED WITH THE 1 CM-1 BAND MODEL.     DRV 0232
C     ALSO, FOR EACH SPECIES, SET THE POINTER WHICH MAPS THE HITRAN     DRV 0233
C     NUMERICAL LABEL TO THE LOWTRAN NUMERICAL LABEL.                   DRV 0234
      KPOINT( 1)=17                                                     DRV 0235
      KPOINT( 2)=36                                                     DRV 0236
      KPOINT( 3)=31                                                     DRV 0237
      KPOINT( 4)=47                                                     DRV 0238
      KPOINT( 5)=44                                                     DRV 0239
      KPOINT( 6)=46                                                     DRV 0240
      KPOINT( 7)=50                                                     DRV 0241
      KPOINT( 8)=54                                                     DRV 0242
      KPOINT( 9)=56                                                     DRV 0243
      KPOINT(10)=55                                                     DRV 0244
      KPOINT(11)=52                                                     DRV 0245
      KPOINT(12)=11                                                     DRV 0246
      IRD1=0                                                            DRV 0247
      IRD2=0                                                            DRV 0248
      IF(MODEL.EQ.0)LENN=0                                              DRV 0249
      IF(MODEL.GE.1 .AND. MODEL.LE.6)THEN                               DRV 0250
          IF(M1.EQ.0)M1=MODEL                                           DRV 0251
          IF(M2.EQ.0)M2=MODEL                                           DRV 0252
          IF(M3.EQ.0)M3=MODEL                                           DRV 0253
          IF(M4.EQ.0)M4=MODEL                                           DRV 0254
          IF(M5.EQ.0)M5=MODEL                                           DRV 0255
          IF(M6.EQ.0)M6=MODEL                                           DRV 0256
          IF(MDEF.EQ.0)MDEF=1                                           DRV 0257
      ENDIF                                                             DRV 0258
      NPR=NOPRNT                                                        DRV 0259
C                                                                       DRV 0260
C     CARD 1B                                                           DRV 0261
      READ(IRD,'(2(L1,I4),F10.3)')DIS,NSTR,SUN1,ISUN,CO2MX              DRV 0262
      IF(CO2MX.LE.0.)CO2MX=330.                                         DRV 0263
      IF(IMULT.EQ.0)DIS=.FALSE.                                         DRV 0264
      WRITE(IPR,'(A,2(L1,I4),F10.3)')                                   DRV 0265
     1  '0 CARD 1B *****',DIS,NSTR,SUN1,ISUN,CO2MX                      DRV 0266
      CO2RAT=CO2MX/330.                                                 DRV 0267
      IF(IEMSCT.GE.2)THEN                                               DRV 0268
C                                                                       DRV 0269
C         SUN TO BE USED.                                               DRV 0270
          IF(SUN1)THEN                                                  DRV 0271
C                                                                       DRV 0272
C             READ SOLAR IRRADIANCES.                                   DRV 0273
              CALL RDSUN(ISUN,SUNFIL)                                   DRV 0274
              IF(ISUN.NE.5)LRDSUN=.TRUE.                                DRV 0275
          ELSEIF(LRDSUN)THEN                                            DRV 0276
C                                                                       DRV 0277
C             SOLAR IRRADIANCE BLOCK DATA WAS OVER-WRITTEN,             DRV 0278
C             BUT IS NOW NEEDED SO IT MUST BE REGENERATED.              DRV 0279
              ISUN=5                                                    DRV 0280
              CALL RDSUN(ISUN,SUNFIL)                                   DRV 0281
              LRDSUN=.FALSE.                                            DRV 0282
          ENDIF                                                         DRV 0283
      ENDIF                                                             DRV 0284
C                                                                       DRV 0285
C     CARD 2 AEROSOL MODEL                                              DRV 0286
      READ(IRD,'(6I5,5F10.5)')                                          DRV 0287
     1  IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT   DRV 0288
      WRITE(IPR,'(A,6I5,5F10.5)')'0 CARD 2  *****',                     DRV 0289
     1  IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT   DRV 0290
C                                                                       DRV 0291
C     CHECK IF IHAZE OR ICLD NEED TO BE RESET.                          DRV 0292
C       IF IHAZE < 0, THEN NO AEROSOLS BUT CLOUDS IF 0 < ICLD < 11      DRV 0293
C       IF IHAZE = 0, THEN NO AEROSOLS AND NO CLOUDS                    DRV 0294
C       IF IHAZE > 0, THEN AEROSOLS AND, IF ICLD > 0, CLOUDS            DRV 0295
      IF(IHAZE.EQ.0 .AND. ICLD.NE.0)THEN                                DRV 0296
C                                                                       DRV 0297
C         RESET ICLD TO ZERO - NO CLOUDS                                DRV 0298
          WRITE(IPR,'(/2A,I3,A)')' WARNING:  INPUT ICLD IS BEING',      DRV 0299
     1      ' RESET FROM',ICLD,' TO 0 SINCE IHAZE EQUALS 0.'            DRV 0300
          ICLD=0                                                        DRV 0301
      ELSEIF(IHAZE.LT.0)THEN                                            DRV 0302
C                                                                       DRV 0303
C         FOR INTERNAL USE, SET IHAZE TO ZERO (CLOUDS WILL BE           DRV 0304
C         INCLUDED IF ICLD IS BETWEEN 1 AND 10, INCLUSIVE).             DRV 0305
          IHAZE=0                                                       DRV 0306
      ENDIF                                                             DRV 0307
C                                                                       DRV 0308
C     CHECK GROUND ALTITUDE.                                            DRV 0309
      IF(GNDALT.NE.0.)WRITE(IPR,'(A,F10.5)')'0   GNDALT =',GNDALT       DRV 0310
      IF(GNDALT.GE.6.)THEN                                              DRV 0311
          WRITE(IPR,'(A,F10.5)')                                        DRV 0312
     1      '0   GNDALT > 6 KM RESET TO 0 KM; GNDALT WAS',GNDALT        DRV 0313
          GNDALT=0.                                                     DRV 0314
      ENDIF                                                             DRV 0315
      IF(VIS.LE.0. .AND. IHAZE.GT.0)VIS=VSB(IHAZE)                      DRV 0316
      RHH=0.                                                            DRV 0317
      IF(MODEL.GE.1 .AND. MODEL.LE.6)THEN                               DRV 0318
          ML=NL                                                         DRV 0319
          IF((MODEL.EQ.3 .OR. MODEL.EQ.5) .AND. ISEASN.EQ.0)ISEASN=2    DRV 0320
          IF(IVSA.EQ.1 .AND. IHAZE.EQ.3)                                DRV 0321
     1      CALL MARINE(VIS,MODEL,WSS,WHH,ICSTL,EXTC,ABSC,1)            DRV 0322
          ICH(1)=IHAZE                                                  DRV 0323
          ICH(2)=6                                                      DRV 0324
          ICH(3)=9+IVULCN                                               DRV 0325
      ENDIF                                                             DRV 0326
      ICH(4)=18                                                         DRV 0327
      IF(ICH(1).LE.0)ICH(1)=1                                           DRV 0328
      IF(ICH(3).LE.9)ICH(3)=10                                          DRV 0329
      IF(ICLD.EQ.11)THEN                                                DRV 0330
          ICH(4)=ICH(3)                                                 DRV 0331
          ICH(3)=ICH(2)                                                 DRV 0332
          ICH(2)=ICLD                                                   DRV 0333
      ENDIF                                                             DRV 0334
      IF(RAINRT.NE.0.)WRITE(IPR,'(A,F9.3,A)')                           DRV 0335
     1  '0   RAIN MODEL CALLED, RAIN RATE = ',RAINRT,' MM/HR'           DRV 0336
      CTHIK=-99.                                                        DRV 0337
      CALT=-99.                                                         DRV 0338
      CEXT=-99.                                                         DRV 0339
      ISEED=-99                                                         DRV 0340
      NCRALT=-99                                                        DRV 0341
      NCRSPC=-99                                                        DRV 0342
      CWAVLN=-99.                                                       DRV 0343
      CCOLWD=-99.                                                       DRV 0344
      CCOLIP=-99.                                                       DRV 0345
      CHUMID=-99.                                                       DRV 0346
      ASYMWD=-99.                                                       DRV 0347
      ASYMIP=-99.                                                       DRV 0348
      IF(ICLD.GE.18)THEN                                                DRV 0349
C                                                                       DRV 0350
C         CARD 2A MODEL CIRRUS                                          DRV 0351
          READ(IRD,'(3F8.3,I4)')CTHIK,CALT,CEXT,ISEED                   DRV 0352
          IF(CTHIK.LT.0.)CTHIK=0.                                       DRV 0353
          IF(CALT.LT.0.)CALT=0.                                         DRV 0354
          IF(CEXT.LT.0.)CEXT=0.                                         DRV 0355
          WRITE(IPR,'(A,3F8.3,I4)')                                     DRV 0356
     1      '0 CARD 2A *****',CTHIK,CALT,CEXT,ISEED                     DRV 0357
      ELSEIF(ICLD.GE.1 .AND. ICLD.LE.10)THEN                            DRV 0358
C                                                                       DRV 0359
C         CARD 2A MODEL CLOUDS                                          DRV 0360
          READ(IRD,'(3F8.3,2I4,6F8.3)')CTHIK,CALT,CEXT,NCRALT,          DRV 0361
     1      NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP            DRV 0362
          WRITE(IPR,'(A,3F8.3,2I4,6F8.3)')'0 CARD 2A *****',CTHIK,CALT, DRV 0363
     1      CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIPDRV 0364
      ENDIF                                                             DRV 0365
C                                                                       DRV 0366
C     CARD 2B VERTICAL STRUCTURE ALGORITHM                              DRV 0367
      ZCVSA=-99.                                                        DRV 0368
      ZTVSA=-99.                                                        DRV 0369
      ZINVSA=-99.                                                       DRV 0370
      IF(IVSA.EQ.1)THEN                                                 DRV 0371
          READ(IRD,'(3F10.5)')ZCVSA,ZTVSA,ZINVSA                        DRV 0372
          WRITE(IPR,'(A,3F10.5)')'0 CARD 2B *****',ZCVSA,ZTVSA,ZINVSA   DRV 0373
          CALL VSA(IHAZE,VIS,ZCVSA,ZTVSA,ZINVSA,ZVSA,RHVSA,AHVSA,IHVSA) DRV 0374
      ENDIF                                                             DRV 0375
      MDLSAV=MODEL                                                      DRV 0376
      IF(MODEL.EQ.0)MDLSAV=8                                            DRV 0377
      DO 70 I=1,5                                                       DRV 0378
          HMODEL(I,7)=HMODEL(I,MDLSAV)                                  DRV 0379
   70 CONTINUE                                                          DRV 0380
      IF(IM.EQ.1)THEN                                                   DRV 0381
          IF(MODEL.EQ.7 .OR. MODEL.EQ.0)THEN                            DRV 0382
C                                                                       DRV 0383
C             CARD 2C:  USER SUPPLIED ATMOSPHERIC PROFILE               DRV 0384
              READ(IRD,'(3I5,18A4)')ML,IRD1,IRD2,(HMODEL(I,7),I=1,5)    DRV 0385
              WRITE(IPR,'(A,3I5,18A4)')                                 DRV 0386
     1          '0 CARD 2C *****',ML,IRD1,IRD2,(HMODEL(I,7),I=1,5)      DRV 0387
              IF(IVSA.EQ.1)CALL RDNSM                                   DRV 0388
          ENDIF                                                         DRV 0389
      ENDIF                                                             DRV 0390
      IF(ICLD.GE.1 .AND. ICLD.LE.10)THEN                                DRV 0391
C                                                                       DRV 0392
C         CLOUD/RAIN MODELS 1-10 ARE NOW SET UP IN ROUTINE CRDRIV, NOT  DRV 0393
C         ROUTINE AERNSM; TEMPORARILY SET ICLD AND RAINRT TO ZERO.      DRV 0394
          ICLDSV=ICLD                                                   DRV 0395
          RAINSV=RAINRT                                                 DRV 0396
          ICLD=0                                                        DRV 0397
          RAINRT=0.                                                     DRV 0398
          CALL AERNSM(JPRT,GNDALT,ICH)                                  DRV 0399
          ICLD=ICLDSV                                                   DRV 0400
          RAINRT=RAINSV                                                 DRV 0401
          CALL CRDRIV                                                   DRV 0402
      ELSE                                                              DRV 0403
          CALL AERNSM(JPRT,GNDALT,ICH)                                  DRV 0404
      ENDIF                                                             DRV 0405
C                                                                       DRV 0406
C     CHECK GROUND ALTITUDE.                                            DRV 0407
      IF(GNDALT.LT.ZM(1))THEN                                           DRV 0408
C                                                                       DRV 0409
C         RAISE GROUND ALTITUDE TO THE BOTTOM OF ATMOSPHERE.            DRV 0410
          WRITE(IPR,'(2A,F8.4,A,/10X,A,F8.4,A)')' WARNING: ',           DRV 0411
     1      ' THE INPUT GROUND ALTITUDE (',GNDALT,' KM) IS BEING',      DRV 0412
     2      ' RAISED TO THE BOTTOM OF THE ATMOSPHERE,',ZM(1),' KM.'     DRV 0413
          GNDALT=ZM(1)                                                  DRV 0414
      ENDIF                                                             DRV 0415
      IF(ICLD.GE.20)THEN                                                DRV 0416
C                                                                       DRV 0417
C         SET UP CIRRUS MODEL                                           DRV 0418
          CALL CIRRUS(CTHIK,CALT,ISEED,CPROB,CEXT)                      DRV 0419
          WRITE(IPR,'(15X,A)')                                          DRV 0420
     1      'CIRRUS ATTENUATION INCLUDED (N O A A CIRRUS)'              DRV 0421
          IF(ISEED.EQ.0)THEN                                            DRV 0422
              WRITE(IPR,'((15X,2A,F10.5,A))')' CIRRUS THICKNESS',       DRV 0423
     1          ' DEFAULTED TO MEAN VALUE OF',CTHIK,'KM','CIRRUS',      DRV 0424
     2          ' BASE ALTITUDE DEFAULTED TO MEAN VALUE OF',CALT,'KM'   DRV 0425
          ELSE                                                          DRV 0426
              IF(CTHIK.NE.0.)THEN                                       DRV 0427
                  WRITE(IPR,'(15X,2A,F10.5,A)')' CIRRUS THICKNESS',     DRV 0428
     1              ' USER DETERMINED TO BE',CTHIK,'KM'                 DRV 0429
              ELSE                                                      DRV 0430
                  WRITE(IPR,'(15X,2A,F10.5,A)')'CIRRUS ATTENUATION',    DRV 0431
     1              ' STATISTICALLY DETERMINED TO BE',CTHIK,'KM'        DRV 0432
              ENDIF                                                     DRV 0433
              IF(CALT.NE.0)THEN                                         DRV 0434
                  WRITE(IPR,'(15X,2A,F10.5,A)')'CIRRUS BASE ALTITUDE',  DRV 0435
     1              ' USER DETERMINED TO BE',CALT,'KM'                  DRV 0436
              ELSE                                                      DRV 0437
                  WRITE(IPR,'(15X,2A,F10.5,A)')'CIRRUS BASE ALTITUDE',  DRV 0438
     1              ' STATISTICALLY DETERMINED TO BE',CALT,'KM'         DRV 0439
              ENDIF                                                     DRV 0440
          ENDIF                                                         DRV 0441
          WRITE(IPR,'(15X,A,F7.1,A)')                                   DRV 0442
     1      'PROBABILITY OF CLOUD OCCURRING IS',CPROB,' PERCENT'        DRV 0443
      ENDIF                                                             DRV 0444
C                                                                       DRV 0445
C*****CARD 2E:  USER-SUPPLIED AEROSOL EXTINCTION,                       DRV 0446
C               ABSORPTION, AND ASYMMETRY                               DRV 0447
      IF(IHAZE.EQ.7 .OR. ICLD.EQ.11)CALL RDEXA                          DRV 0448
   80 CONTINUE                                                          DRV 0449
      IPARM =-99                                                        DRV 0450
      IPH   =-99                                                        DRV 0451
      IDAY  =-99                                                        DRV 0452
      ISOURC=-99                                                        DRV 0453
      PARM1 =-99.                                                       DRV 0454
      PARM2 =-99.                                                       DRV 0455
      PARM3 =-99.                                                       DRV 0456
      PARM4 =-99.                                                       DRV 0457
      TIME  =-99.                                                       DRV 0458
      PSIPO =-99.                                                       DRV 0459
      ANGLEM=-99.                                                       DRV 0460
      G     =-99.                                                       DRV 0461
C                                                                       DRV 0462
C*****CARD 3 GEOMETRY PARAMETERS                                        DRV 0463
      IF(IEMSCT.LT.3)THEN                                               DRV 0464
          READ(IRD,'(6F10.5,I5)')H1,H2,ANGLE,RANGE,BETA,RO,LENN         DRV 0465
          WRITE(IPR,'(A,6F10.5,I5)')                                    DRV 0466
     1      '0 CARD 3  *****',H1,H2,ANGLE,RANGE,BETA,RO,LENN            DRV 0467
      ELSE                                                              DRV 0468
C                                                                       DRV 0469
C         CARD 3 FOR DIRECTLY TRANSMITTED SOLAR IRRADIANCE (IEMSCT=3)   DRV 0470
          READ(IRD,'(2(3F10.5,I5,F15.5,I5))')                           DRV 0471
     1      H1,H2,ANGLE,IDAY,RO,ISOURC,ANGLEM                           DRV 0472
          WRITE(IPR,'(A,2(3F10.5,I5,F15.5,I5))')                        DRV 0473
     1      '0 CARD 3   *****',H1,H2,ANGLE,IDAY,RO,ISOURC,ANGLEM        DRV 0474
          ITYPE=3                                                       DRV 0475
          RANGE=0.                                                      DRV 0476
          BETA=0.                                                       DRV 0477
          LENN=0                                                        DRV 0478
      ENDIF                                                             DRV 0479
C                                                                       DRV 0480
C     RE IS THE RADIUS OF THE EARTH USED BY MODTRAN                     DRV 0481
      RE=6371.23                                                        DRV 0482
      IF(MODEL.EQ.0)THEN                                                DRV 0483
          RO=RE                                                         DRV 0484
      ELSEIF(RO.GT.0.)THEN                                              DRV 0485
          RE=RO                                                         DRV 0486
      ELSEIF(MODEL.EQ.1)THEN                                            DRV 0487
          RE=6378.39                                                    DRV 0488
      ELSEIF(MODEL.EQ.4)THEN                                            DRV 0489
          RE=6356.91                                                    DRV 0490
      ELSEIF(MODEL.EQ.5)THEN                                            DRV 0491
          RE=6356.91                                                    DRV 0492
      ENDIF                                                             DRV 0493
      IF(H1.LT.GNDALT)THEN                                              DRV 0494
          WRITE(IPR,'((2(A,F10.5),A))')'  OBSERVER ALTITUDE (H1 =',     DRV 0495
     1      H1,' KM) BELOW GROUND (GNDALT =',                           DRV 0496
     2      GNDALT,' KM.','    H1 WAS RESET TO GNDALT.'                 DRV 0497
          H1=GNDALT                                                     DRV 0498
      ENDIF                                                             DRV 0499
      H1SAV=H1                                                          DRV 0500
      H2SAV=H2                                                          DRV 0501
      ANGSAV=ANGLE                                                      DRV 0502
      RNGSAV=RANGE                                                      DRV 0503
      BETASV=BETA                                                       DRV 0504
      ITYPSV=ITYPE                                                      DRV 0505
      LENNSV=LENN                                                       DRV 0506
      IF(IEMSCT.EQ.2)THEN                                               DRV 0507
C                                                                       DRV 0508
C         CARD 3A1                                                      DRV 0509
          READ(IRD,'(4I5)')IPARM,IPH,IDAY,ISOURC                        DRV 0510
          WRITE(IPR,'(A,4I5)')'0 CARD 3A1*****',IPARM,IPH,IDAY,ISOURC   DRV 0511
C                                                                       DRV 0512
C         CARD 3A2                                                      DRV 0513
          READ(IRD,'(8F10.5)')                                          DRV 0514
     1      PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G                 DRV 0515
          WRITE(IPR,'(A,8F10.5)')'0 CARD 3A2*****',                     DRV 0516
     1      PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G                 DRV 0517
          REWIND(ISCRCH)                                                DRV 0518
          IF(LFIRST .AND. IMULT.EQ.1)THEN                               DRV 0519
C                                                                       DRV 0520
C             SAVE SOLAR PARAMETERS FOR COMPARING LATER.                DRV 0521
              LFIRST=.FALSE.                                            DRV 0522
              CALL SVSOLA(IPARM,IPH,IDAY,ISOURC,PARM1,PARM2,PARM3,      DRV 0523
     1          PARM4,TIME,G,ANGLEM,IPRMSV,IPHSV,IDAYSV,ISRCSV,         DRV 0524
     2          PRM1SV,PRM2SV,PRM3SV,PRM4SV,TIMESV,GSV,ANGMSV)          DRV 0525
              LSAME=.FALSE.                                             DRV 0526
C                                                                       DRV 0527
C         PRESENTLY, FLUXES ARE NOT SAVED IF THE CORRELATED-K METHOD    DRV 0528
C         IS USED SINCE THE SCRATCH FILE CAN GROW TOO LARGE.            DRV 0529
          ELSEIF(IMULT.EQ.1 .AND. IRPT.EQ.3 .AND.                       DRV 0530
     1      .NOT.DIS .AND. KNTRVL.EQ.1)THEN                             DRV 0531
C                                                                       DRV 0532
C             NOW COMPARE SOLAR PARAMETERS; LSAME IS TRUE IF THEY MATCH.DRV 0533
              CALL COMPAR(IPARM,IPH,IDAY,ISOURC,PARM1,PARM2,PARM3,      DRV 0534
     1          PARM4,TIME,G,ANGLEM,IPRMSV,IPHSV,IDAYSV,ISRCSV,         DRV 0535
     2          PRM1SV,PRM2SV,PRM3SV,PRM4SV,TIMESV,GSV,ANGMSV,LSAME)    DRV 0536
              CALL SVSOLA(IPARM,IPH,IDAY,ISOURC,PARM1,PARM2,PARM3,      DRV 0537
     1          PARM4,TIME,G,ANGLEM,IPRMSV,IPHSV,IDAYSV,ISRCSV,         DRV 0538
     2          PRM1SV,PRM2SV,PRM3SV,PRM4SV,TIMESV,GSV,ANGMSV)          DRV 0539
          ELSE                                                          DRV 0540
C                                                                       DRV 0541
C             PREPARE FOR ADDITIONAL MULTIPLE SOLAR SCATTERING RUNS.    DRV 0542
              LFIRST=.TRUE.                                             DRV 0543
              LSAME=.FALSE.                                             DRV 0544
          ENDIF                                                         DRV 0545
          IF(IPH.EQ.0 .AND. ABS(G).GE.1.)THEN                           DRV 0546
              G=.9999                                                   DRV 0547
              IF(G.LT.0.)G=-.9999                                       DRV 0548
          ELSEIF(IPH.EQ.1)THEN                                          DRV 0549
C                                                                       DRV 0550
C             CARD 3B1 USER DEFINED PHASE FUNCTION                      DRV 0551
              READ(IRD,'(I5)')NANGLS                                    DRV 0552
              WRITE(IPR,'(A,I5)')' CARD 3B1*****',NANGLS                DRV 0553
C                                                                       DRV 0554
C             CARD 3B2                                                  DRV 0555
              READ(IRD,'((5E10.3))')                                    DRV 0556
     1          (ANGF(I),F(1,I),F(2,I),F(3,I),F(4,I),I=1,NANGLS)        DRV 0557
              WRITE(IPR,'(A,5E10.3,/(15X,5E10.3))')'0 CARD 3B2*****',   DRV 0558
     1          (ANGF(I),F(1,I),F(2,I),F(3,I),F(4,I),I=1,NANGLS)        DRV 0559
          ENDIF                                                         DRV 0560
      ENDIF                                                             DRV 0561
   90 CONTINUE                                                          DRV 0562
      IF(IRPT.EQ.3)THEN                                                 DRV 0563
          IF(IPARM.EQ.1)CALL SUBSOL(PARM3,PARM4,TIME,IDAY)              DRV 0564
      ELSE                                                              DRV 0565
C                                                                       DRV 0566
C         CARD 4 WAVENUMBER                                             DRV 0567
          READ(IRD,'(4I10,2A1,A8)')IV1,IV2,IDV,IFWHM,YFLAG,XFLAG,DLIMIT DRV 0568
C                                                                       DRV 0569
C         CHECK plot.dat FILE FLAGS                                     DRV 0570
          IF(YFLAG.EQ.'R' .OR. YFLAG.EQ.'r')THEN                        DRV 0571
C                                                                       DRV 0572
C             WRITE SPECTRAL RADIANCES (TRANSMITTANCES IF IEMSCT=0      DRV 0573
C             OR TRANSMITTED SOLAR IRRADIANCES IF IEMSCT=3)             DRV 0574
              YFLAG='R'                                                 DRV 0575
              IF(IEMSCT.EQ.0)YFLAG='T'                                  DRV 0576
          ELSEIF(YFLAG.EQ.'T' .OR. YFLAG.EQ.'t')THEN                    DRV 0577
C                                                                       DRV 0578
C             WRITE TRANSMITTANCES                                      DRV 0579
              YFLAG='T'                                                 DRV 0580
          ELSE                                                          DRV 0581
C                                                                       DRV 0582
C             DO NOT WRITE TO THE plot.dat FILE.                        DRV 0583
              YFLAG='N'                                                 DRV 0584
          ENDIF                                                         DRV 0585
          IF(YFLAG.NE.'N')THEN                                          DRV 0586
              IF(XFLAG.EQ.'N' .OR. XFLAG.EQ.'n')THEN                    DRV 0587
                  XFLAG='N'                                             DRV 0588
              ELSEIF(XFLAG.EQ.'M' .OR. XFLAG.EQ.'m')THEN                DRV 0589
                  XFLAG='M'                                             DRV 0590
              ELSE                                                      DRV 0591
                  XFLAG='W'                                             DRV 0592
              ENDIF                                                     DRV 0593
              INQUIRE(IPLOT,OPENED=LOPEN)                               DRV 0594
              IF(.NOT.LOPEN)THEN                                        DRV 0595
C                                                                       DRV 0596
C                 OPEN plot.dat FILE                                    DRV 0597
                  OPEN(IPLOT,FILE='plot.dat',STATUS='UNKNOWN')          DRV 0598
                  CLOSE(IPLOT,STATUS='DELETE')                          DRV 0599
                  OPEN(IPLOT,FILE='plot.dat',STATUS='NEW')              DRV 0600
              ENDIF                                                     DRV 0601
          ENDIF                                                         DRV 0602
          IF(DIS .AND. IV2.LE.20)DIS=.FALSE.                            DRV 0603
          IF(DIS .AND. IV1.LE.10)IV1=10                                 DRV 0604
          WRITE(IPR,'(A,4I10)')'0 CARD 4  *****',IV1,IV2,IDV,IFWHM      DRV 0605
          SALB=SALBSV                                                   DRV 0606
          IF(SALB.LT.0)CALL RHOEPS(-1,LSALB,SALBV,DUMN)                 DRV 0607
          IF(IDV.LE.0)THEN                                              DRV 0608
              WRITE(IPR,'(/2A,I10,A)')' WARNING:  IDV IS BEING',        DRV 0609
     1          ' RESET FROM',IDV,' CM-1 TO 1 CM-1.'                    DRV 0610
              IDV=1                                                     DRV 0611
          ENDIF                                                         DRV 0612
          IF(IFWHM.LE.0)THEN                                            DRV 0613
              WRITE(IPR,'(/2A,I10,A)')' WARNING:  IFWHM IS BEING ',     DRV 0614
     1          ' RESET FROM',IFWHM,' CM-1 TO 2 CM-1'                   DRV 0615
              IFWHM=2                                                   DRV 0616
          ENDIF                                                         DRV 0617
          IF(IHAZE.EQ.3 .AND. IV1.LT.250)THEN                           DRV 0618
              IHAZE=4                                                   DRV 0619
              WRITE(IPR,'(//2A,/12X,2A)')' **WARNING** NAVY HAZE MODEL',DRV 0620
     1          ' CAN NOT BE USED BELOW 250 CM-1.',' PROGRAM WILL',     DRV 0621
     2          ' SWITCH TO LOWTRAN 5 MARITIME HAZE MODEL (IHAZE=4).'   DRV 0622
          ENDIF                                                         DRV 0623
          IF(IRPT.LE.1)THEN                                             DRV 0624
              WRITE(IPR,'(7A)')                                         DRV 0625
     1          '0 PROGRAM WILL COMPUTE ',(HTRRAD(I,IEMSCT+1),I=1,6)    DRV 0626
              IF(ISOURC.EQ.1)WRITE(IPR,'(A)')'   LUNAR SOURCE ONLY'     DRV 0627
              IF(IMULT.NE.0)WRITE(IPR,'(A)')                            DRV 0628
     1          '0 CALCULATIONS WILL BE DONE USING MULTIPLE SCATTERING' DRV 0629
              IF(MODEL.GT.0)THEN                                        DRV 0630
                  MM1=M1                                                DRV 0631
                  IF(MM1.EQ.0)MM1=MODEL                                 DRV 0632
                  MM2=M2                                                DRV 0633
                  IF(MM2.EQ.0)MM2=MODEL                                 DRV 0634
                  MM3=M3                                                DRV 0635
                  IF(MM3.EQ.0)MM3=MODEL                                 DRV 0636
                  WRITE(IPR,'(A,/(10X,A,I5,6A))')'0 ATMOSPHERIC MODEL', DRV 0637
     1              'TEMPERATURE =',MM1,'     ',(HMODEL(I,MM1),I=1,5),  DRV 0638
     2              'WATER VAPOR =',MM2,'     ',(HMODEL(I,MM2),I=1,5),  DRV 0639
     3              'OZONE       =',MM3,'     ',(HMODEL(I,MM3),I=1,5)   DRV 0640
                  WRITE(IPR,'(22X,4(A,I6))')                            DRV 0641
     1              'M4 =',M4,' M5 =',M5,' M6 =',M6,' MDEF =' ,MDEF     DRV 0642
              ENDIF                                                     DRV 0643
              IF(IHAZE.GT.0 .AND. JPRT.NE.0)THEN                        DRV 0644
                  IF(ISEASN.EQ.0)ISEASN=1                               DRV 0645
                  IF(IVULCN.LE.0)IVULCN=1                               DRV 0646
                  IHVUL=IVULCN+10                                       DRV 0647
                  IF(IVULCN.EQ.6)IHVUL=11                               DRV 0648
                  IF(IVULCN.EQ.7)IHVUL=11                               DRV 0649
                  IF(IVULCN.EQ.8)IHVUL=13                               DRV 0650
                  IHMET=1                                               DRV 0651
                  IF(IVULCN.GT.1)IHMET=2                                DRV 0652
                  WRITE(IPR,'(A,/10X,A,T35,A,T60,A,T85,A,               DRV 0653
     1              //10X,A,T35,5A,T60,F5.1,A,                          DRV 0654
     2              /(10X,A,T35,5A,T60,5A,T85,5A))')'0 AEROSOL MODEL',  DRV 0655
     3              'REGIME','AEROSOL TYPE','PROFILE','SEASON',         DRV 0656
     4              'BOUNDARY LAYER (0-2 KM)',(HHAZE(I,IHAZE),I=1,5),   DRV 0657
     5              VIS,' KM VIS AT SEA LEVEL',                         DRV 0658
     6              'TROPOSPHERE  (2-10KM)',(HHAZE(I,6),I=1,5),         DRV 0659
     7              (HHAZE(I,6),I=1,5),(HSEASN(I,ISEASN),I=1,5),        DRV 0660
     8              'STRATOSPHERE (10-30KM)',(HHAZE(I,IHVUL),I=1,5),    DRV 0661
     9              (HVULCN(I,IVULCN),I=1,5),(HSEASN(I,ISEASN),I=1,5),  DRV 0662
     &              'UPPER ATMOS (30-100KM)',(HHAZE(I,16),I=1,5),       DRV 0663
     1              (HMET(I,IHMET),I=1,5)                               DRV 0664
              ENDIF                                                     DRV 0665
              IF(ITYPE.EQ.1)THEN                                        DRV 0666
                  WRITE(IPR,'(A,/(10X,A,F11.5,A))')'0 HORIZONTAL PATH', DRV 0667
     1              'ALTITUDE =',H1,' KM','RANGE    =',RANGE,' KM'      DRV 0668
              ELSEIF(ITYPE.EQ.2)THEN                                    DRV 0669
                  WRITE(IPR,'(A,5(/10X,A,F11.5,A),/10X,A,I7)')          DRV 0670
     1              '0 SLANT PATH, H1 TO H2',                           DRV 0671
     2              'H1    =',H1   ,' KM' ,'H2    =',H2   ,' KM',       DRV 0672
     3              'ANGLE =',ANGLE,' DEG','RANGE =',RANGE,' KM',       DRV 0673
     4              'BETA  =',BETA ,' DEG','LENN  =',LENN               DRV 0674
              ELSE                                                      DRV 0675
                  WRITE(IPR,'(A,/(10X,A,F11.5,A))')                     DRV 0676
     1              '0 SLANT PATH TO SPACE','H1    =',H1,' KM',         DRV 0677
     2              'HMIN  =',H2,' KM','ANGLE =',ANGLE,' DEG'           DRV 0678
              ENDIF                                                     DRV 0679
              IF(IEMSCT.EQ.2)THEN                                       DRV 0680
C                                                                       DRV 0681
C                 INTREPRET SOLAR SCATTERING PARAMETERS                 DRV 0682
                  WRITE(IPR,'(A)')                                      DRV 0683
     1              '0 SINGLE SCATTERING CONTROL PARAMETERS SUMMARY'    DRV 0684
                  IF(IPARM.EQ.2)THEN                                    DRV 0685
                      WRITE(IPR,'(/(10X,A,T35,F10.4,A))')               DRV 0686
     1                  'RELATIVE AZIMUTH =',PARM1,' DEG EAST OF NORTH',DRV 0687
     2                  'SOLAR ZENITH ='    ,PARM2,' DEG'               DRV 0688
                  ELSE                                                  DRV 0689
                      IF(IPARM.EQ.1)CALL SUBSOL(PARM3,PARM4,TIME,IDAY)  DRV 0690
                      WRITE(IPR,'(/(10X,A,T35,F10.4,A))')               DRV 0691
     1                  'OBSERVER LATITUDE =',PARM1,                    DRV 0692
     2                  ' DEG NORTH OF EQUATOR',                        DRV 0693
     3                  'OBSERVER LONGITUDE =',PARM2,                   DRV 0694
     4                  ' DEG WEST OF GREENWICH',                       DRV 0695
     5                  'SUBSOLAR LATITUDE =',PARM3,                    DRV 0696
     6                  ' DEG NORTH OF EQUATOR',                        DRV 0697
     7                  'SUBSOLAR LONGITUDE =',PARM4,                   DRV 0698
     8                  ' DEG WEST OF GREENWICH'                        DRV 0699
                  ENDIF                                                 DRV 0700
                  WRITE(IPR,'(10X,2(A,T35,F10.4,A,/10X),A,T35,I10)')    DRV 0701
     1              'TIME (<0 UNDEF) =' ,TIME ,' GREENWICH TIME',       DRV 0702
     2              'PATH AZIMUTH ='    ,PSIPO,' DEG EAST OF NORTH',    DRV 0703
     3              'DAY OF THE YEAR =' ,IDAY                           DRV 0704
                  IF(ISOURC.EQ.0)WRITE(IPR,'(A)')                       DRV 0705
     1              '0 EXTRATERRESTIAL SOURCE IS THE SUN'               DRV 0706
                  IF(ISOURC.EQ.1)WRITE(IPR,'(2A,F10.3,A)')              DRV 0707
     1              '0 EXTRATERRESTIAL SOURCE IS THE MOON,',            DRV 0708
     2              ' MOON PHASE ANGLE =',ANGLEM,' DEG'                 DRV 0709
                  IF(IPH.EQ.0)WRITE(IPR,'(A,F10.5)')                    DRV 0710
     1              '0 H-G PHASE FUNCTION ,G=',G                        DRV 0711
                  IF(IPH.EQ.1)WRITE(IPR,'(A)')                          DRV 0712
     1              '0 USER SUPPLIED PHASE FUNCTION'                    DRV 0713
                  IF(IPH.EQ.2)WRITE (IPR,'(A)')                         DRV 0714
     1              '0 PHASE FUNCTION FROM MIE DATA BASE'               DRV 0715
              ENDIF                                                     DRV 0716
          ENDIF                                                         DRV 0717
          IF(.NOT.MODTRN)THEN                                           DRV 0718
              IV1=5*(IV1/5)                                             DRV 0719
              IV2=5*((IV2+4)/5)                                         DRV 0720
              IDV=5+5*((IDV-5)/5)                                       DRV 0721
          ELSEIF(IV1.GE.22681)THEN                                      DRV 0722
              IV1=5*(IV1/5)                                             DRV 0723
              IV2=5*((IV2+4)/5)                                         DRV 0724
              IDV=5+5*((IDV-5)/5)                                       DRV 0725
              WRITE(IPR,'(/A)')' IDV RESET BEYOND MODTRAN BAND MODEL.'  DRV 0726
          ENDIF                                                         DRV 0727
          IF(IV2.LT.IV1+IDV)THEN                                        DRV 0728
              WRITE(IPR,'(/A)')                                         DRV 0729
     1          ' IV2 WAS LESS THAN IV1 + IDV AND HAS BEEN RESET.'      DRV 0730
              IV2=IV1+IDV                                               DRV 0731
          ENDIF                                                         DRV 0732
          ALAM1= 99999.98                                               DRV 0733
          IF(IV1.NE.0)ALAM1=10000./IV1                                  DRV 0734
          ALAM2=10000./IV2                                              DRV 0735
          IF(IFWHM.LT.1)IFWHM=1                                         DRV 0736
          IF(IFWHM.GT.50)IFWHM=50                                       DRV 0737
          WRITE(IPR,'(A,2(/13X,A,I10,A,F10.2,A),/(11X,A,I10,A))')       DRV 0738
     1      ' FREQUENCY RANGE','IV1 =',IV1,' CM-1  (',ALAM1,            DRV 0739
     2      ' MICRONS)','IV2 =',IV2,' CM-1  (',ALAM2,' MICRONS)',       DRV 0740
     3      '  IDV =',IDV,' CM-1','IFWHM =',IFWHM,' CM-1'               DRV 0741
C                                                                       DRV 0742
C         LOAD ATMOSPHERIC PROFILE INTO /MODEL/                         DRV 0743
          CALL XPROFL                                                   DRV 0744
          CALL STDMDL(ICH(1))                                           DRV 0745
      ENDIF                                                             DRV 0746
C                                                                       DRV 0747
C     INITIALIZE COOLING RATE FITTING COEFFICIENTS                      DRV 0748
CORK  IF(IMULT.NE.0)CALL COOL0(ICR,IRPT)                                DRV 0749
      TBOUND=TBNDSV                                                     DRV 0750
      DO 110 I=1,LAYTHR                                                 DRV 0751
          DO 100 J=1,KMAX                                               DRV 0752
              WPATH(I,J)=0.                                             DRV 0753
              WPATHS(I,J)=0.                                            DRV 0754
  100 CONTINUE                                                          DRV 0755
  110 CONTINUE                                                          DRV 0756
C                                                                       DRV 0757
C     ADJUST GEOMETRY INPUTS FOR SOLAR/LUNAR CALCULATIONS               DRV 0758
C     WITH H1 ABOVE THE TOP OF THE ATMOSPHERE.                          DRV 0759
      IF(IEMSCT.EQ.2)CALL NEWSRC(H1SAV,H2SAV,ANGSAV,RNGSAV,             DRV 0760
     1  BETASV,LENNSV,IPARM,PARM1,PARM2,PSIPO,BETAH2)                   DRV 0761
      IF(IMULT.NE.0)THEN                                                DRV 0762
          H1=GNDALT                                                     DRV 0763
          H2=ZM(ML)                                                     DRV 0764
          ITYPE=2                                                       DRV 0765
          ANGLE=0.                                                      DRV 0766
          BETA=0.                                                       DRV 0767
          RANGE =0.                                                     DRV 0768
          ISSGSV=ISSGEO                                                 DRV 0769
          ISSGEO=0                                                      DRV 0770
          MSOFF=LAYTWO                                                  DRV 0771
          CALL GEO(IERROR,BENDNG,MSOFF,ICH(1))                          DRV 0772
          ISSGEO=ISSGSV                                                 DRV 0773
          IMSMX=IKMAX                                                   DRV 0774
          IF(IEMSCT.EQ.2)THEN                                           DRV 0775
              IF(IMULT.EQ.1 .OR. BETAH2.LE.0.)THEN                      DRV 0776
C                                                                       DRV 0777
C                 SOLAR CONTRIBUTIONS TO MULTIPLE SCATTERING AT H1      DRV 0778
                  CALL SSGEO(IERROR,IPH,IPARM,PARM1,PARM2,              DRV 0779
     1              PARM3,PARM4,PSIPO,G,MSOFF,ICH(1),KNTRVL)            DRV 0780
              ELSE                                                      DRV 0781
C                                                                       DRV 0782
C                 SOLAR CONTRIBUTIONS TO MULTIPLE SCATTERING AT H2      DRV 0783
                  CALL H2SRC(BETAH2,IPH,IPARM,PARM1,PARM2,              DRV 0784
     1              PARM3,PARM4,PSIPO,G,MSOFF,ICH(1),KNTRVL)            DRV 0785
              ENDIF                                                     DRV 0786
          ENDIF                                                         DRV 0787
      ENDIF                                                             DRV 0788
      H1=H1SAV                                                          DRV 0789
      H2=H2SAV                                                          DRV 0790
      ANGLE=ANGSAV                                                      DRV 0791
      RANGE=RNGSAV                                                      DRV 0792
      BETA=BETASV                                                       DRV 0793
      ITYPE=ITYPSV                                                      DRV 0794
      LENN=LENNSV                                                       DRV 0795
C                                                                       DRV 0796
C     TRACE PATH THROUGH THE ATMOSPHERE AND CALCULATE ABSORBER AMOUNTS  DRV 0797
      ISSGEO=0                                                          DRV 0798
      MSOFF=0                                                           DRV 0799
      CALL GEO(IERROR,BENDNG,MSOFF,ICH(1))                              DRV 0800
      CALL AERTMP                                                       DRV 0801
C                                                                       DRV 0802
C     SAVE ZENITH ANGLE DATA.                                           DRV 0803
      UANG=DBLE(ANGLE)                                                  DRV 0804
      IF(IMULT.NE.0)THEN                                                DRV 0805
          DO 120 IK=1,IKMAX                                             DRV 0806
              QTHETS(IK)=QTHETA(IK)                                     DRV 0807
  120     CONTINUE                                                      DRV 0808
      ENDIF                                                             DRV 0809
      IF(IERROR.GT.0)GOTO140                                            DRV 0810
      IF(IEMSCT.EQ.3 .AND. IERROR.EQ. -5)THEN                           DRV 0811
          WRITE(IPR,'(/2A)')' DIRECT PATH TO SUN INTERSECTS',           DRV 0812
     1      ' THE EARTH: SKIP TO NEXT CASE'                             DRV 0813
          GOTO140                                                       DRV 0814
      ENDIF                                                             DRV 0815
C                                                                       DRV 0816
C     THE SECOND CALL TO SSGEO IS TO GET THE CORRECT ANGLES             DRV 0817
C     FOR PHASE FUNCTIONS AND TO SAVE SOLAR PATH INFORMATION.           DRV 0818
      IF(IEMSCT.EQ.2)CALL SSGEO(IERROR,IPH,IPARM,PARM1,                 DRV 0819
     1  PARM2,PARM3,PARM4,PSIPO,G,MSOFF,ICH(1),KNTRVL)                  DRV 0820
      IF(IERROR.GT.0)GOTO140                                            DRV 0821
      IF(IMULT.NE.0)THEN                                                DRV 0822
          DO 130 IK=1,IKMAX                                             DRV 0823
              QTHETA(IK)=QTHETS(IK)                                     DRV 0824
  130     CONTINUE                                                      DRV 0825
      ENDIF                                                             DRV 0826
C                                                                       DRV 0827
C     LOAD AEROSOL EXTINCTION, ABSORPTION, AND ASYMMETRY COEFFICIENTS   DRV 0828
      CALL EXABIN(ICH)                                                  DRV 0829
C                                                                       DRV 0830
C     WRITE HEADER DATA TO TAPE 7                                       DRV 0831
      WRITE(IPU ,'(A1,I4,12I5,F8.3,F7.2)')CODE,MODEL,ITYPE,IEMSCT,      DRV 0832
     1  IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB              DRV 0833
      WRITE(IPR1,'(A1,I4,12I5,F8.3,F7.2)')CODE,MODEL,ITYPE,IEMSCT,      DRV 0834
     1  IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB              DRV 0835
      WRITE(IPU ,'(6I5,5F10.5)')                                        DRV 0836
     1  IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT   DRV 0837
      WRITE(IPR1,'(6I5,5F10.5)')                                        DRV 0838
     1  IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT   DRV 0839
      IF(ICLD.LT.1 .OR. ICLD.GT.10)THEN                                 DRV 0840
          WRITE(IPU ,'(3F8.3,I4)')CTHIK,CALT,CEXT,ISEED                 DRV 0841
          WRITE(IPR1,'(3F8.3,I4)')CTHIK,CALT,CEXT,ISEED                 DRV 0842
      ELSE                                                              DRV 0843
          WRITE(IPU ,'(3F8.3,2I4,6F8.3)')CTHIK,CALT,CEXT,               DRV 0844
     1      NCRALT,NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP     DRV 0845
          WRITE(IPR1,'(3F8.3,2I4,6F8.3)')CTHIK,CALT,CEXT,               DRV 0846
     1      NCRALT,NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP     DRV 0847
      ENDIF                                                             DRV 0848
      WRITE(IPU ,'(3F10.5)')ZCVSA,ZTVSA,ZINVSA                          DRV 0849
      WRITE(IPR1,'(3F10.5)')ZCVSA,ZTVSA,ZINVSA                          DRV 0850
      WRITE(IPU ,'(I5,18A4)')ML,(HMODEL(I,7),I=1,5)                     DRV 0851
      WRITE(IPR1,'(I5,18A4)')ML,(HMODEL(I,7),I=1,5)                     DRV 0852
      IF(MODEL.EQ.0)THEN                                                DRV 0853
          HMDLZ(8)=RANGE                                                DRV 0854
          WRITE(IPU ,'(3F10.5,1P5E10.3)')(HMDLZ(I),I=1,8)               DRV 0855
          WRITE(IPR1,'(3F10.5,1P5E10.3)')(HMDLZ(I),I=1,8)               DRV 0856
      ELSE                                                              DRV 0857
          WRITE(IPU ,'(6F10.5,I5)')H1,H2,ANGLE,RANGE,BETA,RO,LENN       DRV 0858
          WRITE(IPR1,'(6F10.5,I5)')H1,H2,ANGLE,RANGE,BETA,RO,LENN       DRV 0859
      ENDIF                                                             DRV 0860
      WRITE(IPU ,'(4I5)')IPARM,IPH,IDAY,ISOURC                          DRV 0861
      WRITE(IPR1,'(4I5)')IPARM,IPH,IDAY,ISOURC                          DRV 0862
      WRITE(IPU ,'(8F10.5)')PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G DRV 0863
      WRITE(IPR1,'(8F10.5)')PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G DRV 0864
      WRITE(IPU ,'(4I10)')IV1,IV2,IDV,IFWHM                             DRV 0865
      WRITE(IPR1,'(4I10)')IV1,IV2,IDV,IFWHM                             DRV 0866
      REWIND 28                                                         DRV 0867
C                                                                       DRV 0868
C     CARD 5                                                            DRV 0869
      READ(IRD,'(I5)')IRPT                                              DRV 0870
      WRITE(IPU,'(I5)')IRPT                                             DRV 0871
      WRITE(IPR1,'(I5)')IRPT                                            DRV 0872
      GROUND=.FALSE.                                                    DRV 0873
      IF(H2.LE.GNDALT)GROUND=.TRUE.                                     DRV 0874
      CALL TRANS(DIS,NSTR,UANG,IPH,ISOURC,IDAY,ANGLEM,                  DRV 0875
     1  GROUND,LSALB,IMSMX,KNTRVL,YFLAG,XFLAG,DLIMIT,IPLOT)             DRV 0876
  140 CONTINUE                                                          DRV 0877
C                                                                       DRV 0878
C     WRITE END OF FILE LABEL FOR TAPE7 AND TAPE8                       DRV 0879
      WRITE(IPU ,'(A)')' -9999.'                                        DRV 0880
      WRITE(IPR1,'(A)')' -9999.'                                        DRV 0881
      IF(IERROR.GT.0)THEN                                               DRV 0882
          READ(IRD,'(I5)',END=150)IRPT                                  DRV 0883
          WRITE(IPU ,'(I5)')IRPT                                        DRV 0884
          WRITE(IPR1,'(I5)')IRPT                                        DRV 0885
      ENDIF                                                             DRV 0886
      WRITE(IPR,'(A,I5)')'0 CARD 5 *****',IRPT                          DRV 0887
CSSI  CALL SECOND(TSEC)                                                 DRV 0888
CSSI  WRITE(IPR,'(/F10.3,A)')TSEC-TSEC0,                                DRV 0889
CSSI 1  ' SECONDS WERE REQUIRED TO COMPLETE THIS RUN.'                  DRV 0890
CSSI  TSEC0=TSEC                                                        DRV 0891
      IF(IRPT.LE.0 .OR. IRPT.EQ.2 .OR. IRPT.GE.5)RETURN                 DRV 0892
      IF(IRPT.EQ.4)GOTO90                                               DRV 0893
      SALB=SALBSV                                                       DRV 0894
      IF(IRPT.EQ.3)GOTO80                                               DRV 0895
C                                                                       DRV 0896
C     IRPT=1:  START NEXT FULL CALCULATION                              DRV 0897
      GOTO10                                                            DRV 0898
  150 CONTINUE                                                          DRV 0899
C                                                                       DRV 0900
C     RETURN TO MAIN                                                    DRV 0901
      RETURN                                                            DRV 0902
      END                                                               DRV 0903
