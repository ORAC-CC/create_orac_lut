      SUBROUTINE SSGEO(IERROR,IPH,IPARM,PARM1,PARM2,PARM3,PARM4,        SSG 0001
     1  PSIPO,G,MSOFF,ICH1,KNTRVL)                                      SSG 0002
C                                                                       SSG 0003
C     THIS ROUTINE CALLS THE LOWTRAN GEOMETRY ROUTINES REPEATEDLY       SSG 0004
C     TO OBTAIN THE ABSORBER AMOUNTS FROM THE SCATTERING POINTS ON      SSG 0005
C     THE OPTICAL PATH TO THE EXTRATERRESTRIAL SOURCE, AS IS NECESSARY  SSG 0006
C     FOR THE LAYER BY LAYER SINGLE SCATTERING RADIANCE CALCULATION.    SSG 0007
C                                                                       SSG 0008
C     DECLARE INPUTS                                                    SSG 0009
C       IERROR   ERROR FLAG (=0 FOR SUCCESSFUL CALLS)                   SSG 0010
C       IPH      PHASE FUNCTION FLAG (=0 FOR HENYEY-GREENSTEIN)         SSG 0011
C                                    (=1 FOR USER-DEFINED)              SSG 0012
C                                    (=2 FOR MIE-GENERATED)             SSG 0013
C       IPARM    SOLAR INPUT GEOMETRY FLAG                              SSG 0014
C       PARM1    OBSERVER LATITUDE (IPARM=0 OR IPARM=1) [DEG NORTH]     SSG 0015
C                SOLAR AZIMUTH ANGLE (IPARM=2) [DEG]                    SSG 0016
C       PARM2    OBSERVER LONGITUDE (IPARM=0 OR IPARM=1) [DEG WEST]     SSG 0017
C                SOLAR ZENITH ANGLE (IPARM=2) [DEG]                     SSG 0018
C       PARM3    SOURCE LATITUDE (IPARM=0 OR IPARM=1) [DEG NORTH]       SSG 0019
C       PARM4    SOURCE LONGITUDE (IPARM=0 OR IPARM=1) [DEG WEST]       SSG 0020
C       PSIPO    PATH AZIMUTH ANGLE [DEG EAST OF NORTH]                 SSG 0021
C       G        HENYEY-GREENSTEIN ASYMMETRY FACTOR (IPH=0)             SSG 0022
C       MSOFF    MULTIPLE SCATTERING LAYER OFFSET (EQUALS 0 FOR         SSG 0023
C                OPTICAL PATH AND EQUALS LAYTWO FOR VERTICAL PATH)      SSG 0024
C       ICH1     HAZE MODEL FLAG                                        SSG 0025
C       KNTRVL   NUMBER OF CORRELATED-K SUB-INTERVALS                   SSG 0026
C                (=1 IF CORRELATED-K APPROACH IS NOT USED)              SSG 0027
      INTEGER IERROR,IPH,IPARM,MSOFF,ICH1,KNTRVL                        SSG 0028
      REAL PARM1,PARM2,PARM3,PARM4,PSIPO,G                              SSG 0029
C                                                                       SSG 0030
C     LIST PARAMETERS                                                   SSG 0031
      INCLUDE 'PARAM.LST'                                               SSG 0032
C                                                                       SSG 0033
C     LIST COMMONS                                                      SSG 0034
      REAL TXX,WX,WPATHX                                                SSG 0035
      COMMON/NONAME/TXX(MMOLX),WX(MMOLX),WPATHX(LAYTHR,MMOLX)           SSG 0036
      REAL WPTHSX,TBBYSX                                                SSG 0037
      COMMON/SOLSX/WPTHSX(LAYTHR,MMOLX),TBBYSX(LAYTHR,MMOLX)            SSG 0038
      INTEGER KPOINT                                                    SSG 0039
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     SSG 0040
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   SSG 0041
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   SSG 0042
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     SSG 0043
      INTEGER JTURN,LJ                                                  SSG 0044
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       SSG 0045
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  SSG 0046
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      SSG 0047
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    SSG 0048
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               SSG 0049
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           SSG 0050
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     SSG 0051
      REAL TBOUND,SALB                                                  SSG 0052
      LOGICAL MODTRN                                                    SSG 0053
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   SSG 0054
     1  MODTRN                                                          SSG 0055
      INTEGER NCRALT,NCRSPC                                             SSG 0056
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    SSG 0057
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      SSG 0058
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       SSG 0059
      INTEGER LENN                                                      SSG 0060
      REAL H1,H2,ANGLE,RANGE,BETA,REE                                   SSG 0061
      COMMON/CARD3/H1,H2,ANGLE,RANGE,BETA,REE,LENN                      SSG 0062
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     SSG 0063
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                SSG 0064
      REAL RE,ZMAX                                                      SSG 0065
      INTEGER IMAX,IMOD,IPATH                                           SSG 0066
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             SSG 0067
C                                                                       SSG 0068
C       PI       THE CONSTANT PI                                        SSG 0069
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       SSG 0070
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       SSG 0071
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         SSG 0072
      REAL PI,DEG,BIGNUM,BIGEXP                                         SSG 0073
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                SSG 0074
      REAL ZP,PP,TP,RFNDXP,SP,PPSUM,TPSUM,RHOPSM,DENP,AMTP              SSG 0075
      COMMON/RFRPTH/ZP(LAYDIM+1),PP(LAYDIM+1),TP(LAYDIM+1),             SSG 0076
     1  RFNDXP(LAYDIM+1),SP(LAYDIM+1),PPSUM(LAYDIM+1),TPSUM(LAYDIM+1),  SSG 0077
     2  RHOPSM(LAYDIM+1),DENP(KMAX,LAYDIM+1),AMTP(KMAX,LAYDIM+1)        SSG 0078
      INTEGER NANGLS                                                    SSG 0079
      REAL ANGF,F                                                       SSG 0080
      COMMON/USRDTA/NANGLS,ANGF(50),F(4,50)                             SSG 0081
C                                                                       SSG 0082
C     DECLARE LOCAL VARIABLES                                           SSG 0083
      CHARACTER*48 MESSAG                                               SSG 0084
      INTEGER MSOFFX,JTRNSV,IKMAX1,LENNSV,ITYPSV,J,MSOFFJ,K,IARBO,      SSG 0085
     1  IARB,LM1,L,LJ0,JITER,MM,LLM1,LL,MSOFFL,KP,IK,MSOFFK             SSG 0086
      REAL X,THETAO,PHIO,THETAS,PHIS,PSIO,DELO,H1SAV,H2SAV,ANGSAV,      SSG 0087
     1  RNGSAV,BETASV,PSIPOS,BETAST,PSIST,ANGL0,RELH,THTST,ANGERR,ANGMX,SSG 0088
     2  BENDNG,HTOP,ANGSV,ANGMN,SANGLE,COSANG,PNUM,TNUM,DENOM,WPTH      SSG 0089
C                                                                       SSG 0090
C     DECLARE LOCAL ARRAYS                                              SSG 0091
      INTEGER LJSAV(LAYTWO)                                             SSG 0092
      REAL WPSAV(LAYTWO,KMAX),TBSAV(LAYTWO),WSAV(KMAX),ZPSAV(LAYDIM+1), SSG 0093
     1  PSAV(LAYTWO),RHSAV(LAYDIM+1),WPSUM(KMAX),WPSAVX(LAYTHR,MMOLX)   SSG 0094
C                                                                       SSG 0095
C     DECLARE LOCAL FUNCTIONS                                           SSG 0096
      REAL PFMOL,PSI,DEL,SCTANG,HENGNS                                  SSG 0097
C                                                                       SSG 0098
C     MOLECULAR PHASE FUNCTION [STER-1]; X=COS(SCATTERING ANGLE)        SSG 0099
      PFMOL(X)=.06050402+.0572197*X**2                                  SSG 0100
      MSOFFX=0                                                          SSG 0101
      NPR=2                                                             SSG 0102
      ISSGEO=1                                                          SSG 0103
C                                                                       SSG 0104
C     SPECIFY THE GEOMETRICAL CONFIGURATION                             SSG 0105
      IF(IPARM.EQ.2)THEN                                                SSG 0106
          PSIO=PARM1                                                    SSG 0107
          DELO=PARM2                                                    SSG 0108
      ELSE                                                              SSG 0109
          THETAO=PARM1                                                  SSG 0110
          PHIO=PARM2                                                    SSG 0111
          THETAS=PARM3                                                  SSG 0112
          PHIS=PARM4                                                    SSG 0113
          IF(THETAO.LT.-89.5)THEN                                       SSG 0114
C                                                                       SSG 0115
C             OBSERVER IS AT OR NEAR THE SOUTH POLE.                    SSG 0116
C             REMAP TO THE EQUATOR.                                     SSG 0117
              WRITE(IPR,'(/2A)')                                        SSG 0118
     1          '  THETAO < -89.5 (OBSERVER AT OR NEAR THE SOUTH',      SSG 0119
     2          ' POLE).  PROBLEM HAS BEEN REMAPPED TO THE EQUATOR.'    SSG 0120
              PSIPO=PSIPO-PHIS                                          SSG 0121
              THETAO=0.                                                 SSG 0122
              PHIO=0.                                                   SSG 0123
              THETAS=0.                                                 SSG 0124
              PHIS=90.+THETAS                                           SSG 0125
          ELSEIF(THETAO.GT.89.5)THEN                                    SSG 0126
C                                                                       SSG 0127
C             OBSERVER IS AT OR NEAR THE NORTH POLE.  REMAP TO EQUATOR. SSG 0128
              WRITE(IPR,'(/2A)')                                        SSG 0129
     1          '  THETAO > +89.5 (OBSERVER AT OR NEAR THE NORTH',      SSG 0130
     2          ' POLE).  PROBLEM HAS BEEN REMAPPED TO THE EQUATOR.'    SSG 0131
              PSIPO=PHIS-PSIPO                                          SSG 0132
              THETAO=0.                                                 SSG 0133
              PHIO=0.                                                   SSG 0134
              THETAS=0.                                                 SSG 0135
              PHIS=90.-THETAS                                           SSG 0136
          ENDIF                                                         SSG 0137
      ENDIF                                                             SSG 0138
      WRITE(IPR,'(//A)')' SINGLE SCATTERING POINT TO SOURCE PATHS '     SSG 0139
C                                                                       SSG 0140
C     SAVE OPTICAL PATH PARAMETERS AND AMOUNTS                          SSG 0141
      JTRNSV=JTURN                                                      SSG 0142
      IKMAX1=IKMAX+1                                                    SSG 0143
      H1SAV=H1                                                          SSG 0144
      H2SAV=H2                                                          SSG 0145
      ANGSAV=ANGLE                                                      SSG 0146
      RNGSAV=RANGE                                                      SSG 0147
      BETASV=BETA                                                       SSG 0148
      BETA=0.                                                           SSG 0149
      LENNSV=LENN                                                       SSG 0150
      ITYPSV=ITYPE                                                      SSG 0151
      DO 10 J=1,ML                                                      SSG 0152
          ZPSAV(J)=ZP(J)                                                SSG 0153
          RHSAV(J)=RELHUM(J)                                            SSG 0154
   10 CONTINUE                                                          SSG 0155
      DO 40 J=1,IKMAX1                                                  SSG 0156
          MSOFFJ=MSOFF+J                                                SSG 0157
          TBSAV(J)=TBBY(MSOFFJ)                                         SSG 0158
          PSAV(J)=PATM(MSOFFJ)                                          SSG 0159
          LJSAV(J)=LJ(J)                                                SSG 0160
          IF(LJSAV(J).GT.ML)LJSAV(J)=ML                                 SSG 0161
          DO 20 K=1,KMAX                                                SSG 0162
              WPSAV(J,K)=WPATH(MSOFFJ,K)                                SSG 0163
   20     CONTINUE                                                      SSG 0164
          DO 30 K=1,NSPECX                                              SSG 0165
              WPSAVX(J,K)=WPATHX(MSOFFJ,K)                              SSG 0166
   30     CONTINUE                                                      SSG 0167
   40 CONTINUE                                                          SSG 0168
      DO 50 K=1,KMAX                                                    SSG 0169
          WSAV(K)=W(K)                                                  SSG 0170
          WPSUM(K)=0.                                                   SSG 0171
   50 CONTINUE                                                          SSG 0172
      IMAX=ML-1                                                         SSG 0173
C                                                                       SSG 0174
C     ESTABLISH PSIO AND DELO                                           SSG 0175
      IARBO=0                                                           SSG 0176
      IF(ANGLE.LT..01 .OR. ANGLE.GT.179.99)IARBO=1                      SSG 0177
      IARB=IARBO                                                        SSG 0178
      BETAST=0.                                                         SSG 0179
      IF(IPARM.NE.2)THEN                                                SSG 0180
          CALL PSIECA(THETAO,PHIO,THETAS,PHIS,PSIPOS,DELO)              SSG 0181
          PSIO=PSIPOS-PSIPO                                             SSG 0182
          IF(PSIO.GT.180.)PSIO=PSIO-360.                                SSG 0183
          IF(PSIO.LT.-180.)PSIO=PSIO+360.                               SSG 0184
      ENDIF                                                             SSG 0185
      PSIST=PSIO                                                        SSG 0186
      ANGL0=DELO                                                        SSG 0187
C                                                                       SSG 0188
C     LOOP OVER THE POINT TO SUN PATHS TO OBTAIN AMOUNTS                SSG 0189
      WRITE(IPR,'((2A))')' SCTTR SCTTR SUBTENDED  SOLAR   PATH',        SSG 0190
     1  '  RELATIVE  SCTTR   MOLECULAR',' POINT  ALT    ANGLE',         SSG 0191
     2  '   ZENITH  ZENITH AZIMUTH   ANGLE   PHASE F'                   SSG 0192
      LM1=0                                                             SSG 0193
      DO 190 L=1,IKMAX1                                                 SSG 0194
          IF(LENNSV.EQ.0 .AND. JTRNSV.EQ.0)THEN                         SSG 0195
C                                                                       SSG 0196
C             SHORT PATH, UP                                            SSG 0197
              IF(L.GE.2)BETAST=BETAST+ADBETA(LM1)                       SSG 0198
              H1=ZPSAV(L)                                               SSG 0199
              RELH=RHSAV(L)                                             SSG 0200
              THTST=ATHETA(L)                                           SSG 0201
          ELSE                                                          SSG 0202
C                                                                       SSG 0203
C             LONG PATH, OR SHORT PATH DOWN                             SSG 0204
              IF(L.GE.2)BETAST=BETAST+ADBETA(LJSAV(LM1))                SSG 0205
              LJ0=LJSAV(L)                                              SSG 0206
              IF(L.LT.JTRNSV)LJ0=LJ0+1                                  SSG 0207
              THTST=ATHETA(LJ0)                                         SSG 0208
              IF(L.LE.JTRNSV)THTST=180.-THTST                           SSG 0209
              H1=ZPSAV(LJ0)                                             SSG 0210
              RELH=RHSAV(LJ0)                                           SSG 0211
          ENDIF                                                         SSG 0212
          AH1(L)=H1                                                     SSG 0213
          ARH(L)=RELH                                                   SSG 0214
          IF(L.GE.2)THEN                                                SSG 0215
              PSIST=PSI(PSIO,DELO,BETAST,IARB,IARBO)                    SSG 0216
              ANGL0=DEL(PSIO,DELO,BETAST,IARBO)                         SSG 0217
          ENDIF                                                         SSG 0218
          ANGERR=0.                                                     SSG 0219
          ANGMX=ANGL0                                                   SSG 0220
          BENDNG=0.                                                     SSG 0221
          ANGLE=ANGL0                                                   SSG 0222
C                                                                       SSG 0223
C         RANGE=UNKNOWN                                                 SSG 0224
          ITYPE=3                                                       SSG 0225
          DO 80 JITER=1,12                                              SSG 0226
              MESSAG=' '                                                SSG 0227
C                                                                       SSG 0228
C             SET H2 TO ZERO TO INSURE THAT ANGLE                       SSG 0229
C             IS USED TO DEFINE SOLAR PATH.                             SSG 0230
              H2=0.                                                     SSG 0231
              LENN=0                                                    SSG 0232
              IF(ANGLE.GT.90.)THEN                                      SSG 0233
                  IF(ANGLE.LE.90.0001)THEN                              SSG 0234
C                                                                       SSG 0235
C                     MULTIPLE SCATTERING CORRECTION FOR                SSG 0236
C                     SCATTERING POINT TO SUN PATHS.                    SSG 0237
                      ANGLE=90.                                         SSG 0238
                  ELSE                                                  SSG 0239
                      LENN=1                                            SSG 0240
                      MESSAG(1:23)='THIS SOLAR PATH PASSES '            SSG 0241
                      MESSAG(24:48)='THROUGH A TANGENT HEIGHT.'         SSG 0242
                  ENDIF                                                 SSG 0243
              ENDIF                                                     SSG 0244
              HTOP=ZMAX                                                 SSG 0245
              IF(H1.GE.HTOP .AND. LENN.EQ.0)THEN                        SSG 0246
C                                                                       SSG 0247
C                 SCATTERING POINT IS AT OR ABOVE HTOP                  SSG 0248
C                 AND LENN=0;  SET W(K)=0.0 AND CONTINUE.               SSG 0249
                  DO 60 K=1,KMAX                                        SSG 0250
                      W(K)=0.                                           SSG 0251
   60             CONTINUE                                              SSG 0252
                  DO 70 K=1,NSPECX                                      SSG 0253
                      WX(K)=0.                                          SSG 0254
   70             CONTINUE                                              SSG 0255
                  GOTO90                                                SSG 0256
              ENDIF                                                     SSG 0257
              ANGSV=ANGLE                                               SSG 0258
              CALL GEO(IERROR,BENDNG,MSOFFX,ICH1)                       SSG 0259
              ANGLE=ANGSV                                               SSG 0260
C                                                                       SSG 0261
C             IERROR=-5 IF SCATTERING POINT IS SHADED; SET W(36)=-5.    SSG 0262
              IF(IERROR.EQ.-5)THEN                                      SSG 0263
                  MESSAG='THIS SCATTERING POINT IS IN THE SHADE.'       SSG 0264
                  W(36)=-5.                                             SSG 0265
                  IERROR=0                                              SSG 0266
                  GOTO90                                                SSG 0267
              ENDIF                                                     SSG 0268
C                                                                       SSG 0269
C             SOLAR ZENITH ERROR                                        SSG 0270
              ANGERR=ANGLE+BENDNG-ANGL0                                 SSG 0271
              IF(ABS(ANGERR).LT..001)GOTO90                             SSG 0272
C                                                                       SSG 0273
C             CALCULATE ANGLE USING THE BISECTION METHOD.  ANGL0 IS     SSG 0274
C             A MAXIMUM AND ANGL0 MINUS ITS BENDING IS A MINIMUM.       SSG 0275
              IF(JITER.EQ.1)THEN                                        SSG 0276
                  ANGMN=ANGL0-BENDNG                                    SSG 0277
                  ANGLE=ANGMN                                           SSG 0278
              ELSE                                                      SSG 0279
                  IF(ANGERR.GT.0)THEN                                   SSG 0280
                      ANGMX=ANGLE                                       SSG 0281
                  ELSE                                                  SSG 0282
                      ANGMN=ANGLE                                       SSG 0283
                  ENDIF                                                 SSG 0284
                  ANGLE=.5*(ANGMX+ANGMN)                                SSG 0285
              ENDIF                                                     SSG 0286
   80     CONTINUE                                                      SSG 0287
          MESSAG='THE SOLAR ZENITH ANGLE DID NOT CONVERGE.'             SSG 0288
          WRITE(IPR,'(2A,F12.5,A)')                                     SSG 0289
     1      ' AFTER 12 ITERATIONS, THE SOLAR ZENITH EXITING',           SSG 0290
     2      ' THE ATMOSPHERE IS STILL IN ERROR BY',ANGERR,'DEG'         SSG 0291
   90     CONTINUE                                                      SSG 0292
          SANGLE=SCTANG(ANGLE,THTST,PSIST,IARB)                         SSG 0293
          COSANG=COS(SANGLE/DEG)                                        SSG 0294
C                                                                       SSG 0295
C         LOAD MOLECULAR PHASE FUNCTION ARRAY                           SSG 0296
          PHASFN(L,1)=PFMOL(COSANG)                                     SSG 0297
C                                                                       SSG 0298
C         LOAD AEROSOL PHASE FUNCTION ARRAY                             SSG 0299
C         HENYEY-GREENSTEIN                                             SSG 0300
          IF(IPH.EQ.0)THEN                                              SSG 0301
              PHASFN(L,2)=HENGNS(G,COSANG)                              SSG 0302
          ELSEIF(IPH.EQ.1)THEN                                          SSG 0303
C                                                                       SSG 0304
C             USER SUPPLIED PHASE FUNCTION                              SSG 0305
C             DETERMINE ALTITUDE AND ANGLE INDICES                      SSG 0306
              MM=4                                                      SSG 0307
              IF(H1.LE.30.)MM=3                                         SSG 0308
              IF(H1.LE.10.)MM=2                                         SSG 0309
              IF(H1.LE.2.)MM=1                                          SSG 0310
              LLM1=1                                                    SSG 0311
              DO 100 LL=2,NANGLS-1                                      SSG 0312
                  IF(ANGF(LL).GT.SANGLE)GOTO110                         SSG 0313
  100         LLM1=LL                                                   SSG 0314
              LL=NANGLS                                                 SSG 0315
  110         CONTINUE                                                  SSG 0316
              CALL INTERP(2,SANGLE,ANGF(LLM1),ANGF(LL),                 SSG 0317
     1          PHASFN(L,2),F(MM,LLM1),F(MM,LL))                        SSG 0318
          ELSE                                                          SSG 0319
C                                                                       SSG 0320
C             V DEPENDENT MIE DATA BASE, SAVE SCATTERING ANGLE INSTEAD  SSG 0321
              PHASFN(L,2)=SANGLE                                        SSG 0322
          ENDIF                                                         SSG 0323
C                                                                       SSG 0324
C         LOAD WATER DROPLET CLOUD PHASE FUNCTION INFORMATION           SSG 0325
          IF(ABS(ASYMWD).LT.1.)THEN                                     SSG 0326
C                                                                       SSG 0327
C             SPECTRALLY INDEPENDENT ASYMMETRY                          SSG 0328
C             FACTOR; STORE PHASE FUNCTION.                             SSG 0329
              PHASFN(L,3)=HENGNS(ASYMWD,COSANG)                         SSG 0330
          ELSE                                                          SSG 0331
C                                                                       SSG 0332
C             SPECTRALLY DEPENDENT ASYMMETRY                            SSG 0333
C             FACTOR; STORE SCATTERING COSINE.                          SSG 0334
              PHASFN(L,3)=COSANG                                        SSG 0335
          ENDIF                                                         SSG 0336
C                                                                       SSG 0337
C         LOAD ICE PARTICLE CLOUD PHASE FUNCTION INFORMATION            SSG 0338
          IF(ABS(ASYMIP).LT.1.)THEN                                     SSG 0339
C                                                                       SSG 0340
C             SPECTRALLY INDEPENDENT ASYMMETRY                          SSG 0341
C             FACTOR; STORE PHASE FUNCTION.                             SSG 0342
              PHASFN(L,4)=HENGNS(ASYMIP,COSANG)                         SSG 0343
          ELSE                                                          SSG 0344
C                                                                       SSG 0345
C             SPECTRALLY DEPENDENT ASYMMETRY                            SSG 0346
C             FACTOR; STORE SCATTERING ANGLE.                           SSG 0347
              PHASFN(L,4)=COSANG                                        SSG 0348
          ENDIF                                                         SSG 0349
C                                                                       SSG 0350
C         STORE SOLAR ZENITH DATA FOR MULTIPLE SCATTERING.              SSG 0351
          IF(MSOFF.GT.0)CALL SOLZEN(L,IKMAX1,REE,DEG,ANGLE)             SSG 0352
C                                                                       SSG 0353
C         WRITE SOLAR SCATTER PATH DATA.                                SSG 0354
          WRITE(IPR,'(I4,6F8.2,E11.3,3X,A)')                            SSG 0355
     1      L,H1,BETAST,ANGLE,THTST,PSIST,SANGLE,PHASFN(L,1),MESSAG     SSG 0356
C                                                                       SSG 0357
C         LOAD AMOUNTS FROM W(K) INTO WPATHS(L,K)                       SSG 0358
          MSOFFL=MSOFF+L                                                SSG 0359
          IF(MSOFF.EQ.0 .AND. L.GT.1 .AND.                              SSG 0360
     1      W(36).GE.0. .AND. KNTRVL.EQ.1)THEN                          SSG 0361
C                                                                       SSG 0362
C             ADD OBSERVER TO SCATTERING POINT PATH AMOUNTS TO          SSG 0363
C             W FOR ALL BUT THE MODTRAN BAND MODEL ABSORBERS.           SSG 0364
              DO 120 K=1,KMAX                                           SSG 0365
                  WPSUM(K)=WPSUM(K)+WPSAV(LM1,K)                        SSG 0366
                  WPATHS(MSOFFL,K)=W(K)+WPSUM(K)                        SSG 0367
  120         CONTINUE                                                  SSG 0368
              IF(MODTRN)THEN                                            SSG 0369
                  DO 130 K=1,NSPC                                       SSG 0370
                      KP=KPOINT(K)                                      SSG 0371
                      WPATHS(MSOFFL,KP)=W(KP)                           SSG 0372
  130             CONTINUE                                              SSG 0373
              ENDIF                                                     SSG 0374
          ELSE                                                          SSG 0375
              DO 140 K=1,KMAX                                           SSG 0376
                  WPATHS(MSOFFL,K)=W(K)                                 SSG 0377
  140         CONTINUE                                                  SSG 0378
          ENDIF                                                         SSG 0379
          DO 150 K=1,NSPECX                                             SSG 0380
              WPTHSX(MSOFFL,K)=WX(K)                                    SSG 0381
  150     CONTINUE                                                      SSG 0382
C                                                                       SSG 0383
C         WHEN THE MODERATE RESOLUTION OPTION IS USED, A CURTIS-GODSON  SSG 0384
C         AVERAGE PRESSURE, PATMS, AND TEMPERATURE, TBBYS, IS DEFINED   SSG 0385
C         FOR THE SCATTERING POINT TO EXTRATERRESTRIAL SOURCE "LAYER"   SSG 0386
          IF(MODTRN)THEN                                                SSG 0387
              DO 180 K=1,NSPECT                                         SSG 0388
                  TNUM=0.                                               SSG 0389
                  DENOM=0.                                              SSG 0390
                  IF(K.LE.NSPC)THEN                                     SSG 0391
                      PNUM=0.                                           SSG 0392
                      KP=KPOINT(K)                                      SSG 0393
                      DO 160 IK=1,IKMAX                                 SSG 0394
                          MSOFFK=MSOFFX+IK                              SSG 0395
                          WPTH=WPATH(MSOFFK,KP)                         SSG 0396
                          IF(WPTH.LE.0)GOTO160                          SSG 0397
                          PNUM=PNUM+PATM(MSOFFK)*WPTH                   SSG 0398
                          TNUM=TNUM+TBBY(MSOFFK)*WPTH                   SSG 0399
                          DENOM=DENOM+WPTH                              SSG 0400
  160                 CONTINUE                                          SSG 0401
                      IF(DENOM.NE.0.)THEN                               SSG 0402
                          PATMS(MSOFFL,K)=PNUM/DENOM                    SSG 0403
                          TBBYS(MSOFFL,K)=TNUM/DENOM                    SSG 0404
                      ENDIF                                             SSG 0405
                  ELSE                                                  SSG 0406
                      DO 170 IK=1,IKMAX                                 SSG 0407
                          MSOFFK=MSOFFX+IK                              SSG 0408
                          WPTH=WPATHX(MSOFFK,K-NSPC)                    SSG 0409
                          IF(WPTH.LE.0)GOTO170                          SSG 0410
                          TNUM=TNUM+TBBY(MSOFFK)*WPTH                   SSG 0411
                          DENOM=DENOM+WPTH                              SSG 0412
  170                 CONTINUE                                          SSG 0413
                      IF(DENOM.NE.0.)TBBYSX(MSOFFL,K-NSPC)=TNUM/DENOM   SSG 0414
                  ENDIF                                                 SSG 0415
  180         CONTINUE                                                  SSG 0416
          ENDIF                                                         SSG 0417
          LM1=L                                                         SSG 0418
  190 CONTINUE                                                          SSG 0419
      ANGSUN=ANGLE                                                      SSG 0420
C                                                                       SSG 0421
C     RESTORE OPTICAL PATH AMOUNTS                                      SSG 0422
      IKMAX=IKMAX1-1                                                    SSG 0423
      H1=H1SAV                                                          SSG 0424
      H2=H2SAV                                                          SSG 0425
      ANGLE=ANGSAV                                                      SSG 0426
      RANGE=RNGSAV                                                      SSG 0427
      BETA=BETASV                                                       SSG 0428
      LENN=LENNSV                                                       SSG 0429
      ITYPE=ITYPSV                                                      SSG 0430
      NPR=NOPRNT                                                        SSG 0431
      DO 220 J=1,IKMAX1                                                 SSG 0432
          MSOFFJ=MSOFF+J                                                SSG 0433
          TBBY(MSOFFJ)=TBSAV(J)                                         SSG 0434
          PATM(MSOFFJ)=PSAV(J)                                          SSG 0435
          LJ(J)=LJSAV(J)                                                SSG 0436
          DO 200 K=1,KMAX                                               SSG 0437
              WPATH(MSOFFJ,K)=WPSAV(J,K)                                SSG 0438
  200     CONTINUE                                                      SSG 0439
          DO 210 K=1,NSPECX                                             SSG 0440
              WPATHX(MSOFFJ,K)=WPSAVX(J,K)                              SSG 0441
  210     CONTINUE                                                      SSG 0442
  220 CONTINUE                                                          SSG 0443
      DO 230 K=1,KMAX                                                   SSG 0444
          W(K)=WSAV(K)                                                  SSG 0445
  230 CONTINUE                                                          SSG 0446
C                                                                       SSG 0447
C     RETURN TO DRIVER                                                  SSG 0448
      RETURN                                                            SSG 0449
      END                                                               SSG 0450
