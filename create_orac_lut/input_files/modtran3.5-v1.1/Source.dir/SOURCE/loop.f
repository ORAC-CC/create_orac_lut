C
C This routine modified by RS to output aerosol extinction and
C absorption in each layer.
C File aerosol.out is created, containing the data
C
C R.S. 13/11/97
C
C
      SUBROUTINE LOOP(LOOP0,IV,IVX,IKMX,MXFREQ,SUMTMS,SUMMS,TRANSM,     LOP 0001
     1  IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,RADCUM,S0,FRAC,UANG,        LOP 0002
     2  DIS,NSTR,GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL)             LOP 0003
C                                                                       LOP 0004
C     THIS ROUTINE PERFORMS THE LOOP OVER LAYERS FOR EACH FREQUENCY     LOP 0005
      INCLUDE 'PARAM.LST'                                               LOP 0006
      integer aerlun,iptx
      data iptx/0/
      data aerlun/77/
      INTEGER IV,IVX,IKMX,MXFREQ,IPH,NSTR,KNTRVL                        LOP 0007
      REAL SUMTMS,SUMMS,SUMSSS,UNIF,TRACE,TRANSX,                       LOP 0008
     1  RADCUM,S0,FRAC,TSNOBS,TSNREF,FDNTRT,FDNSRT                      LOP 0009
      DOUBLE PRECISION UANG                                             LOP 0010
      LOGICAL LOOP0,TRANSM,IVTEST,DIS,GROUND                            LOP 0011
C                                                                       LOP 0012
C     CONVENTION                                                        LOP 0013
C     MMOLX=MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")           LOP 0014
C     MMOL =MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")              LOP 0015
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             LOP 0016
C                                                                       LOP 0017
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              LOP 0018
C     NSPEC=ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL        LOP 0019
C     NSPECX=ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX       LOP 0020
C                                                                       LOP 0021
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           LOP 0022
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      LOP 0023
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        LOP 0024
      REAL TXX,WX,WPATHX                                                LOP 0025
      COMMON/NONAME/TXX(MMOLX),WX(MMOLX),WPATHX(LAYTHR,MMOLX)           LOP 0026
      INTEGER KPOINT                                                    LOP 0027
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     LOP 0028
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   LOP 0029
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   LOP 0030
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     LOP 0031
C                                                                       LOP 0032
C       PI       THE CONSTANT PI.                                       LOP 0033
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       LOP 0034
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       LOP 0035
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         LOP 0036
      REAL PI,DEG,BIGNUM,BIGEXP                                         LOP 0037
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                LOP 0038
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               LOP 0039
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           LOP 0040
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     LOP 0041
      REAL TBOUND,SALB                                                  LOP 0042
      LOGICAL MODTRN                                                    LOP 0043
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   LOP 0044
     1  MODTRN                                                          LOP 0045
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,     DRV 0023
     1  RAINRT
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     LOP 0046
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                LOP 0047
C                                                                       LOP 0048
C     /PATH/                                                            LOP 0049
C       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.               LOP 0050
C       AHT     ALTITUDES AT PATH BOUNDARIES.                           LOP 0051
C       TPH     TEMPERATURE AT PATH BOUNDARIES.                         LOP 0052
C       IMAP    MAPPING FROM PATH SEGMENTS TO LAYERS.                   LOP 0053
      INTEGER IMAP                                                      LOP 0054
      REAL QTHETA,AHT,TPH                                               LOP 0055
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),TPH(LAYTWO),IMAP(LAYTWO)   LOP 0056
      INTEGER JTURN,LJ                                                  LOP 0057
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       LOP 0058
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  LOP 0059
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      LOP 0060
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    LOP 0061
C                                                                       LOP 0062
C     COMMON /MSRD/                                                     LOP 0063
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            LOP 0064
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             LOP 0065
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               LOP 0066
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    LOP 0067
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  LOP 0068
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. LOP 0069
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               LOP 0070
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         LOP 0071
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          LOP 0072
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              LOP 0073
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  LOP 0074
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               LOP 0075
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           LOP 0076
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     LOP 0077
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   LOP 0078
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       LOP 0079
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     LOP 0080
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              LOP 0081
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 LOP 0082
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          LOP 0083
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             LOP 0084
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                LOP 0085
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      LOP 0086
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         LOP 0087
      REAL TXNEW,TXOLD                                                  LOP 0088
      COMMON/LAY5/TXNEW(20,LAYTHR,3),TXOLD(20,LAYTHR,3)                 LOP 0089
      SAVE /LAY5/                                                       LOP 0090
      LOGICAL LSAME                                                     LOP 0091
      COMMON/SOLAR/LSAME                                                LOP 0092
      REAL ZM,PM,TM,RFNDX,DENSTY                                        LOP 0093
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    LOP 0094
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               LOP 0095
C                                                                       LOP 0096
C     COMMON /CORKDT/                                                   LOP 0097
C       WTKSUB   SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.           LOP 0098
C       DEPLAY   INCREMENTAL EXTINCTION OPTICAL DEPTHS                  LOP 0099
C       TRNLAY   INCREMENTAL TRANSMITTANCES                             LOP 0100
C       TRNCUM   CUMULATIVE TRANSMITTANCES                              LOP 0101
      REAL WTKSUB,DEPLAY,TRNLAY,TRNCUM                                  LOP 0102
      COMMON/CORKDT/WTKSUB(NKSUB),DEPLAY(NKSUB),                        LOP 0103
     1  TRNLAY(NKSUB),TRNCUM(NKSUB)                                     LOP 0104
      SAVE /CORKDT/                                                     LOP 0105
C                                                                       LOP 0106
C       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.       LOP 0107
C       UPFLX    LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].    LOP 0108
C       DNFLX    LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].  LOP 0109
C       UPFLXS   BOUNDARY UPWARD SCATTERED SOLAR FLUX [W CM-2 / CM-1].  LOP 0110
C       DNFLXS   BOUNDARY DOWNWARD SCATTERED SOLAR FLUX [W CM-2 / CM-1].LOP 0111
C       NTFLX    LAYER BOUNDARY NET (THERMAL PLUS SCATTERED SOLAR       LOP 0112
C                PLUS DIRECT SOLAR) UPWARD FLUX [W CM-2 / CM-1].        LOP 0113
      REAL SUBINT,UPFLX,DNFLX,UPFLXS,DNFLXS,NTFLX                       LOP 0114
      COMMON/NETFLX/SUBINT(NKSUB),UPFLX(LAYDIM),DNFLX(LAYDIM),          LOP 0115
     1  UPFLXS(LAYDIM),DNFLXS(LAYDIM),NTFLX(LAYDIM)                     LOP 0116
C                                                                       LOP 0117
C     DECLARE FUNCTION NAMES                                            LOP 0118
      REAL BETABS,BBFN                                                  LOP 0119
C                                                                       LOP 0120
C     DECLARE LOCAL VARIABLES                                           LOP 0121
      INTEGER I,J,K,IPATH,MSOFF,IK,IKP1,N,IKOFF,INDEXX,                 LOP 0122
     1  INTRVL,NFRNT,NBACK,NP1                                          LOP 0123
      REAL V,STORE,BLAYER,BBOUND,DTAU,THMLAY,SOLLAY,COEF,WTTRAN,        LOP 0124
     1  UP,DN,THMBND,SOLBND,THMDIF,SOLDIF,BNCOEF,DFCOEF,THMSCT,         LOP 0125
     2  TRNOLD,TRNNEW,TX9LAY,TX9CUM,OPTDEP,EMSLAY,B0,OMEGA,TMOL,        LOP 0126
     3  FACTOR,SMTRNL,SMTHML,SMSOLL                                     LOP 0127
      if(iptx.eq.0)
     1     open(aerlun,file='aerosol.out',status='unknown')
C                                                                       LOP 0128
C     LOOP0 IS .TRUE. FOR FIRST CALL FROM ROUTINE TRANS                 LOP 0129
      IF(LOOP0)THEN                                                     LOP 0130
C                                                                       LOP 0131
C         DEFINE NUMBER OF CORRELATED-K METHOD SPECTRAL                 LOP 0132
C         BINS (=1 IF METHOD IS NOT BEING USED).                        LOP 0133
          SUBINT(1)=1.                                                  LOP 0134
          IF(KNTRVL.GT.1)THEN                                           LOP 0135
              DO 10 INTRVL=1,KNTRVL                                     LOP 0136
                  SUBINT(INTRVL)=WTKSUB(INTRVL)                         LOP 0137
   10         CONTINUE                                                  LOP 0138
          ENDIF                                                         LOP 0139
C                                                                       LOP 0140
C         INITIALIZE TRANSMISSION ARRAYS                                LOP 0141
          DO 20 I=1,20                                                  LOP 0142
              DO 20 J=1,LAYTHR                                          LOP 0143
                  DO 20 K=1,3                                           LOP 0144
                      TXOLD(I,J,K)=0.                                   LOP 0145
                      TXNEW(I,J,K)=0.                                   LOP 0146
   20     CONTINUE                                                      LOP 0147
          DO 30 K=1,NSPECX                                              LOP 0148
              TXX(K)=1.                                                 LOP 0149
   30     CONTINUE                                                      LOP 0150
C                                                                       LOP 0151
C         INITIALIZATION FOR MULTIPLE SCATTERING CALCULATIONS           LOP 0152
          IF(IMULT.NE.0)CALL MAPMS(ML,IKMAX)                            LOP 0153
C                                                                       LOP 0154
C         INITIALIZATION FOR LOWTRAN CALCULATIONS                       LOP 0155
          IF(.NOT.MODTRN)CALL FRQ5DT(LOOP0,IV)                          LOP 0156
          RETURN                                                        LOP 0157
      ENDIF                                                             LOP 0158
C                                                                       LOP 0159
C     DEFINE THE LAYER INDEPENDENT 5 CM-1 DATA                          LOP 0160
      IF(IVTEST)CALL FRQ5DT(LOOP0,IV)                                   LOP 0161
      V=IVX                                                             LOP 0162
      IPATH=1                                                           LOP 0163
C                                                                       LOP 0164
C     MSOFF IS EQUAL LAYTWO FOR THE MULTIPLE SCATTERING LAYER LOOP.     LOP 0165
      IF(.NOT.LSAME)THEN                                                LOP 0166
          MSOFF=ABS(IMULT)*LAYTWO                                       LOP 0167
      ELSE                                                              LOP 0168
C                                                                       LOP 0169
C         IF (LSAME) FLUX DATA ARE UNCHANGED FROM THE PREVIOUS RUN.     LOP 0170
C         SET MSOFF TO 0 TO SKIP THE MULTIPLE SCATTERING CALCULATION.   LOP 0171
          MSOFF=0                                                       LOP 0172
          READ(ISCRCH)(COSBAR(IK),TSCAT(IK),                            LOP 0173
     1      (TAUT(INTRVL,IK),INTRVL=1,KNTRVL),IK=1,ML-1),               LOP 0174
     2      ((UPF(INTRVL,IK),DNF(INTRVL,IK),                            LOP 0175
     3      UPFS(INTRVL,IK),DNFS(INTRVL,IK),INTRVL=1,KNTRVL),IK=1,ML)   LOP 0176
      ENDIF                                                             LOP 0177
C                                                                       LOP 0178
C     FOR EACH WAVENUMBER, AN INITIAL CALL (IK=0) TO BMOD IS REQUIRED   LOP 0179
      IK=0                                                              LOP 0180
   40 IF(MODTRN)CALL BMOD(IK,IKMX,IPATH,IVX,MSOFF,MXFREQ)               LOP 0181
C                                                                       LOP 0182
C     INITIALIZE PARAMETERS                                             LOP 0183
      IPATH=1                                                           LOP 0184
      TRNOLD=1.                                                         LOP 0185
      TX9CUM=1.                                                         LOP 0186
      IF(.NOT.TRANSM)THEN                                               LOP 0187
          SUMMS=0.                                                      LOP 0188
          SUMTMS=0.                                                     LOP 0189
          DO 50 INTRVL=1,KNTRVL                                         LOP 0190
              TRNCUM(INTRVL)=1.                                         LOP 0191
   50     CONTINUE                                                      LOP 0192
      ENDIF                                                             LOP 0193
C                                                                       LOP 0194
C     BEGINNING OF LAYER LOOP                                           LOP 0195
      IF(MSOFF.GT.0)THEN                                                LOP 0196
          BBNDRY(1)=PI*BBFN(TM(1),V)                                    LOP 0197
          IF(GROUND .AND. TBOUND.GT.0.)THEN                             LOP 0198
              BBGRND=PI*BBFN(TBOUND,V)                                  LOP 0199
          ELSE                                                          LOP 0200
              BBGRND=BBNDRY(1)                                          LOP 0201
          ENDIF                                                         LOP 0202
      ENDIF                                                             LOP 0203
      DO 180 IK=1,IKMAX                                                 LOP 0204
C                                                                       LOP 0205
C         MSOFF IS A LAYER OFFSET PARAMETER EQUAL TO 0 FOR              LOP 0206
C         THE LINE-OF-SIGHT PATH AND EQUAL TO LAYTWO FOR THE            LOP 0207
C         MULTIPLE SCATTERING VERTICAL PATH FROM GROUND TO SPACE.       LOP 0208
          IKOFF=IK+MSOFF                                                LOP 0209
C                                                                       LOP 0210
C         FOR TRANSMISSION RUNS, W(K) WAS DEFINED IN ROUTINE GEO.       LOP 0211
          IF(TRANSM)GOTO110                                             LOP 0212
C                                                                       LOP 0213
C         LOAD APPROPRIATE ABSORBER AMOUNTS INTO W(K)                   LOP 0214
          IF(IEMSCT.EQ.1)GOTO90                                         LOP 0215
C                                                                       LOP 0216
C         FOR FIRST LAYER (IK=1) CALCULATE SUN-TO-OBSERVER TRANSMITTANCELOP 0217
          IF(IPATH.EQ.1)THEN                                            LOP 0218
              INDEXX=MSOFF+1                                            LOP 0219
              IF(WPATHS(INDEXX,36).GE.0.)THEN                           LOP 0220
C                                                                       LOP 0221
C                 LOAD W WITH SUN-TO-OBSERVER PATH ABSORBER AMOUNTS     LOP 0222
                  DO 60 K=1,KMAX                                        LOP 0223
   60             W(K)=WPATHS(INDEXX,K)                                 LOP 0224
                  GOTO110                                               LOP 0225
              ENDIF                                                     LOP 0226
C                                                                       LOP 0227
C             IN SSGEO, SCATTERING POINT INDEXX WAS DETERMINED TO       LOP 0228
C             BE IN THE SHADE AND WPATHS(INDEXX,36) WAS SET TO -5.      LOP 0229
              CALL SHADE(IPH,IK,MSOFF,IPATH,                            LOP 0230
     1          KNTRVL,V,TSNOBS,TSNREF,SUMSSS)                          LOP 0231
          ENDIF                                                         LOP 0232
   70     IPATH=2                                                       LOP 0233
          INDEXX=MSOFF+IK+1                                             LOP 0234
          IF(WPATHS(INDEXX,36).GE.0.)THEN                               LOP 0235
C                                                                       LOP 0236
C             LOAD W WITH SUN-TO-SCATTERING POINT ABSORBER AMOUNTS.     LOP 0237
C             FOR THE OPTICAL PATH (MSOFF=0) WHEN THE CORRELATED-K      LOP 0238
C             APPROACH IS NOT USED (KNTRVL=1), MOST OF THE WPATHS ALSO  LOP 0239
C             INCLUDE THE OBSERVER-TO-SCATTERING POINT ABSORBER AMOUNTS.LOP 0240
              DO 80 K=1,KMAX                                            LOP 0241
                  W(K)=WPATHS(INDEXX,K)                                 LOP 0242
   80         CONTINUE                                                  LOP 0243
              GOTO110                                                   LOP 0244
          ENDIF                                                         LOP 0245
C                                                                       LOP 0246
C         IN SSGEO, THE W'S WERE SET TO -5 IF THE SUN WAS BLOCKED       LOP 0247
          CALL SHADE(IPH,IK,MSOFF,IPATH,KNTRVL,V,TSNOBS,TSNREF,SUMSSS)  LOP 0248
C                                                                       LOP 0249
C         LOAD W(K) WITH OPTICAL PATH ABSORBER AMOUNT                   LOP 0250
   90     IPATH=3                                                       LOP 0251
          DO 100 K=1,KMAX                                               LOP 0252
              W(K)=WPATH(IKOFF,K)                                       LOP 0253
  100     CONTINUE                                                      LOP 0254
C                                                                       LOP 0255
C         DEFINE THE LAYER DEPENDENT 5 CM-1 DATA                        LOP 0256
  110     CONTINUE                                                      LOP 0257
          IF(IVTEST)CALL LAY5DT(V,IKOFF,IPATH)                          LOP 0258
          if(msoff.eq.0.and.ipath.eq.3.and.iptx.eq.0) then
             if(ik.eq.1) then 
                write(aerlun,*)';IHAZE: ',IHAZE
		write(aerlun,*)';ISEASN: ',ISEASN
		write(aerlun,*)';IVULCN: ',IVULCN
		write(aerlun,*)';ICSTL: ',ICSTL
		write(aerlun,*)';ICLD: ',ICLD
		write(aerlun,*)';IVSA: ',IVSA
		write(aerlun,*)';VIS: ',VIS
		write(aerlun,*)';WSS: ',WSS
		write(aerlun,*)';WHH: ',WHH
		write(aerlun,*)';RAINRT: ',RAINRT
                write(aerlun,*)';V/CM-1: ',V
                write(aerlun,*)
     1 ';LAYER BOUNDARIES(KM), INCREMENTAL AEROSOL SCATTERING AND ',
     2 ';EXTINCTION OPTICAL DEPTH'
             endif
             write(aerlun,*)zm(imap(ik)),zm(imap(ik)+1),
     1             txnew(2,ikoff,ipath),
     2             txnew(13,ikoff,ipath)
          endif
C                                                                       LOP 0259
C         DEFINE TX ARRAY                                               LOP 0260
C           1   ASYMMETRY PARAMETER WEIGHTED BY SCATTERING DEPTH        LOP 0261
C           2   INCREMENTAL AEROSOL SCATTERING OPTICAL DEPTH            LOP 0262
C           3   TOTAL O2 CONTINUUM TRANSMITTANCE                        LOP 0263
C           4   N2 CONTINUUM TRANSMITTANCE                              LOP 0264
C           5   TOTAL H2O CONTINUUM TRANSMITTANCE                       LOP 0265
C           6   RAYLEIGH MOLECULAR SCATTERED TRANSMITTANCE              LOP 0266
C           7   AEROSOL EXTINCTION                                      LOP 0267
C           8   TOTAL OZONE CONTINUUM TRANSMITTANCE                     LOP 0268
C           9   PRODUCT OF ALL CONTINUUM TRANSMITTANCES EXCEPT O2 & HNO3LOP 0269
C          10   AEROSOL ABSORPTION                                      LOP 0270
C          11   HNO3 TRANSMITTANCE                                      LOP 0271
C          12   MOLECULAR CONTINUUM OPTICAL DEPTH                       LOP 0272
C          13   INCREMENTAL AEROSOL+CLOUD EXTINCTION OPTICAL DEPTH      LOP 0273
C          14   TOTAL CONTINUUM OPTICAL DEPTH                           LOP 0274
C          15   LAYER RAYLEIGH MOLECULAR SCATTERING OPTICAL DEPTH       LOP 0275
C          16   CIRRUS CLOUD TRANSMITTANCE (ICLD=20 ONLY)               LOP 0276
C          64   UV/VIS NO2 TRANSMITTANCE                                LOP 0277
C          65   UV/VIS SO2 TRANSMITTANCE                                LOP 0278
C          66   INCREMENTAL WATER DROPLET SCATTERING OPTICAL DEPTH      LOP 0279
C          67   INCREMENTAL ICE PARTICLE SCATTERING OPTICAL DEPTH       LOP 0280
C          ---  MOLECULAR LINE CENTER TRANSMITTANCE  ---                LOP 0281
C          17=H2O  36=CO2  31=O3   47=N2O  44=CO   46=CH4               LOP 0282
C          50=O2   54=NO   56=SO2  55=NO2  52=NH3  11=HNO3              LOP 0283
          DO 120 K=1,16                                                 LOP 0284
              STORE=TXNEW(K,IKOFF,IPATH)                                LOP 0285
              TX(K)=STORE+FRAC*(TXOLD(K,IKOFF,IPATH)-STORE)             LOP 0286
  120     CONTINUE                                                      LOP 0287
          STORE=TXNEW(17,IKOFF,IPATH)                                   LOP 0288
          TX(64)=STORE+FRAC*(TXOLD(17,IKOFF,IPATH)-STORE)               LOP 0289
          STORE=TXNEW(18,IKOFF,IPATH)                                   LOP 0290
          TX(65)=STORE+FRAC*(TXOLD(18,IKOFF,IPATH)-STORE)               LOP 0291
          STORE=TXNEW(19,IKOFF,IPATH)                                   LOP 0292
          TX(66)=STORE+FRAC*(TXOLD(19,IKOFF,IPATH)-STORE)               LOP 0293
          STORE=TXNEW(20,IKOFF,IPATH)                                   LOP 0294
          TX(67)=STORE+FRAC*(TXOLD(20,IKOFF,IPATH)-STORE)               LOP 0295
CORK      IF(KNTRVL.GT.1)THEN                                           LOP 0296
CORK          IF(IPATH.LT.3)THEN                                        LOP 0297
C                                                                       LOP 0298
C                 SOLAR PATH                                            LOP 0299
CORK              CALL BMCRKS(KNTRVL,INDEXX)                            LOP 0300
CORK              CALL SSCORK(IPH,IK,MSOFF,                             LOP 0301
CORK 1              IPATH,KNTRVL,V,S0,TSNOBS,SUMSSS)                    LOP 0302
CORK              GOTO(70,90),IPATH                                     LOP 0303
CORK          ELSE                                                      LOP 0304
C                                                                       LOP 0305
C                 OPTICAL (LINE-OF-SIGHT) PATH                          LOP 0306
CORK              CALL BMCORK(KNTRVL,IK,MSOFF)                          LOP 0307
CORK              IF(IEMSCT.EQ.2)CALL SSCORK(IPH,IK,MSOFF,              LOP 0308
CORK 1              IPATH,KNTRVL,V,S0,TSNOBS,SUMSSS)                    LOP 0309
CORK          ENDIF                                                     LOP 0310
CORK      ELSE                                                          LOP 0311
              IF(MODTRN)CALL BMOD(IK,IKMX,IPATH,IVX,MSOFF,MXFREQ)       LOP 0312
C                                                                       LOP 0313
C             BEER'S LAW SPECIES CONTRIBUTIONS                          LOP 0314
              TRANSX=TXX(1)                                             LOP 0315
              DO 130 K=2,NSPECX                                         LOP 0316
                  TRANSX=TRANSX*TXX(K)                                  LOP 0317
  130         CONTINUE                                                  LOP 0318
C                                                                       LOP 0319
C             COMBINE TRANSMISSIONS IF NOT USING CORRELATED-K METHOD    LOP 0320
C               UNIF    UNIFORMLY MIXED GASES TRANSMITTANCE             LOP 0321
C               TRACE   TRACE GASES TRANSMITTANCE                       LOP 0322
              UNIF=TX(36)*TX(44)*TX(46)*TX(47)*TX(50)                   LOP 0323
              TRACE=TX(52)*TX(54)*TX(55)*TX(56)*TX(11)                  LOP 0324
              IF(IPATH.EQ.3)THEN                                        LOP 0325
                  IF(TRNOLD.LE.0.)THEN                                  LOP 0326
C                     TRANSMITTANCE OF THE CURRENT LAYER CANNOT         LOP 0327
C                     BE DETERMINED.  ASSUME IT IS ZERO.                LOP 0328
                      TRNLAY(1)=0.                                      LOP 0329
                      DEPLAY(1)=BIGEXP                                  LOP 0330
                  ELSE                                                  LOP 0331
                      DEPLAY(1)=0.                                      LOP 0332
                      TRNLAY(1)=1.                                      LOP 0333
                      TRNNEW=TX(17)*UNIF*TX(31)*TRACE*TRANSX            LOP 0334
                      IF(TRNNEW.LT.TRNOLD)THEN                          LOP 0335
C                                                                       LOP 0336
C                         DETERMINE DECREASE IN TRANSMITTANCE.          LOP 0337
                          TRNLAY(1)=TRNNEW/TRNOLD                       LOP 0338
                          DEPLAY(1)=BIGEXP                              LOP 0339
                          IF(TRNLAY(1).GT.0.)DEPLAY(1)=-LOG(TRNLAY(1))  LOP 0340
                          TRNOLD=TRNNEW                                 LOP 0341
                      ENDIF                                             LOP 0342
                  ENDIF                                                 LOP 0343
                  DEPLAY(1)=DEPLAY(1)+TX(14)                            LOP 0344
                  IF(MSOFF.EQ.0)THEN                                    LOP 0345
                      TX9LAY=TX(9)*TX(3)*TX(64)*TX(65)                  LOP 0346
                      TX9CUM=TX9CUM*TX9LAY                              LOP 0347
                      TX(9)=TX9CUM*TRNNEW                               LOP 0348
                      TRNLAY(1)=TRNLAY(1)*TX9LAY                        LOP 0349
                      IF(IEMSCT.EQ.2)CALL SSRAD(IPH,IK,MSOFF,           LOP 0350
     1                  IPATH,V,S0,TRNNEW,TSNOBS,TSNREF,SUMSSS)         LOP 0351
                  ENDIF                                                 LOP 0352
              ELSEIF(TRANSM)THEN                                        LOP 0353
C                                                                       LOP 0354
C                 DEFINE SPECIES DEPENDENT TRANSMITTANCES               LOP 0355
C                 FOR TRANSMISSION ONLY CALCULATIONS.                   LOP 0356
                  UNIF=UNIF*TX(3)                                       LOP 0357
                  TRACE=TRACE*TX(64)*TX(65)                             LOP 0358
                  TX(9)=TX(17)*UNIF*TX(31)*TX(9)*TRACE*TRANSX           LOP 0359
C                                                                       LOP 0360
C                 COMBINE IR AND UV/VIS CONTRIBUTIONS                   LOP 0361
C                 FOR NO2, SO2 AND O2.                                  LOP 0362
                  TX(64)=TX(64)*TX(55)                                  LOP 0363
                  TX(65)=TX(65)*TX(56)                                  LOP 0364
                  TX(50)=TX(50)*TX(3)                                   LOP 0365
C                                                                       LOP 0366
C                 OZONE CONTRIBUTIONS ARE MODELED WITH                  LOP 0367
C                 THE MODTRAN BAND MODEL BELOW 9170 CM-1                LOP 0368
C                 AND AS A CONTINUUM ABOVE 9170 CM-1.                   LOP 0369
                  IF(IV.GE.9170)TX(31)=TX(8)                            LOP 0370
C                                                                       LOP 0371
C                 IF TRANSMISSION ONLY CALCULATION, RETURN              LOP 0372
                  RETURN                                                LOP 0373
              ELSEIF(IEMSCT.EQ.2)THEN                                   LOP 0374
C                                                                       LOP 0375
C                 CALCULATE TMOL, THE MOLECULAR TRANSMITTANCES          LOP 0376
C                 TO THE SUN.  VALUES OF TMOL LESS THAN 1.E-8           LOP 0377
C                 ARE NOT ACCURATE AND ARE SET TO ZERO.                 LOP 0378
                  TMOL=TX(17)*UNIF*TX(31)*TRACE*TRANSX                  LOP 0379
                  IF(TMOL.LT.1.E-8)TMOL=0.                              LOP 0380
                  CALL SSRAD(IPH,IK,MSOFF,IPATH,                        LOP 0381
     1              V,S0,TMOL,TSNOBS,TSNREF,SUMSSS)                     LOP 0382
                  GOTO(70,90),IPATH                                     LOP 0383
              ENDIF                                                     LOP 0384
CORK      ENDIF                                                         LOP 0385
C                                                                       LOP 0386
C         INITIALIZE TOTAL TRANSMITTANCE                                LOP 0387
          TX(9)=0.                                                      LOP 0388
          IF(MSOFF.GT.0)THEN                                            LOP 0389
C                                                                       LOP 0390
C             SET UP FOR MULTIPLE SCATTERING VERTICAL PATH CALCULATION  LOP 0391
              IKP1=IK+1                                                 LOP 0392
              BBNDRY(IKP1)=PI*BBFN(TM(IKP1),V)                          LOP 0393
CORK          CALL THMFLX(KNTRVL,IK,IKP1,ML,V,                          LOP 0394
CORK 1          BBGRND,BBNDRY(IK),BBNDRY(IKP1))                         LOP 0395
              TSCAT(IK)=TX(2)+TX(66)+TX(67)+TX(15)                      LOP 0396
              COSBAR(IK)=0.                                             LOP 0397
              IF(TSCAT(IK).GT.0.)COSBAR(IK)=TX(1)/TSCAT(IK)             LOP 0398
              IF(MODTRN .OR. DIS)THEN                                   LOP 0399
                  DO 140 INTRVL=1,KNTRVL                                LOP 0400
                      TAUT(INTRVL,IK)=DEPLAY(INTRVL)                    LOP 0401
  140             CONTINUE                                              LOP 0402
              ELSE                                                      LOP 0403
                  TCONT(IK)=TX(12)+TX(13)+TX(15)                        LOP 0404
              ENDIF                                                     LOP 0405
C                                                                       LOP 0406
C             INITIALIZE NET FLUX ARRAY                                 LOP 0407
              UPFLX(IK)=0.                                              LOP 0408
              DNFLX(IK)=0.                                              LOP 0409
              UPFLXS(IK)=0.                                             LOP 0410
              DNFLXS(IK)=0.                                             LOP 0411
              NTFLX(IK)=0.                                              LOP 0412
              GOTO180                                                   LOP 0413
          ELSEIF(IMULT.EQ.0)THEN                                        LOP 0414
C                                                                       LOP 0415
C             NO MULTIPLE SCATTERING.  DEFINE BLACKBODY FUNCTIONS.      LOP 0416
              BLAYER=BBFN(TBBY(IK),V)                                   LOP 0417
              BBOUND=BBFN(TPH(IK),V)                                    LOP 0418
C                                                                       LOP 0419
C             LOOP OVER CORRELATED-K METHOD SUB-INTERVALS               LOP 0420
              SMTRNL=0.                                                 LOP 0421
              SMTHML=0.                                                 LOP 0422
              DO 150 INTRVL=1,KNTRVL                                    LOP 0423
                  OPTDEP=DEPLAY(INTRVL)                                 LOP 0424
                  IF(OPTDEP.LT..02)THEN                                 LOP 0425
                      THMLAY=OPTDEP*(BLAYER+OPTDEP*(BBOUND-4*BLAYER)/6) LOP 0426
                  ELSE                                                  LOP 0427
                      EMSLAY=1.-TRNLAY(INTRVL)                          LOP 0428
                      THMLAY=EMSLAY*BBOUND+2*(BLAYER-BBOUND)            LOP 0429
     1                  *(EMSLAY/OPTDEP-TRNLAY(INTRVL))                 LOP 0430
                  ENDIF                                                 LOP 0431
                  THMLAY=SUBINT(INTRVL)*TRNCUM(INTRVL)*THMLAY           LOP 0432
                  RADCUM=RADCUM+THMLAY                                  LOP 0433
                  TRNCUM(INTRVL)=TRNCUM(INTRVL)*TRNLAY(INTRVL)          LOP 0434
                  TX(9)=TX(9)+SUBINT(INTRVL)*TRNCUM(INTRVL)             LOP 0435
                  SMTRNL=SMTRNL+SUBINT(INTRVL)*TRNLAY(INTRVL)           LOP 0436
                  SMTHML=SMTHML+SUBINT(INTRVL)*THMLAY                   LOP 0437
  150         CONTINUE                                                  LOP 0438
C                                                                       LOP 0439
C             IF NOPRNT=-1, PRINT WEIGHTING FUNCTION DATA.              LOP 0440
              IF(NOPRNT.EQ.-1)WRITE(IPR1,                               LOP 0441
     1          '(I7,2F10.5,3X,2(1P2E11.3,0P2F11.7))')IVX,AHT(IK),      LOP 0442
     2          AHT(IK+1),BLAYER,BBOUND,TX(9),SMTRNL,SMTHML,RADCUM      LOP 0443
          ELSEIF(DIS)THEN                                               LOP 0444
C                                                                       LOP 0445
C             DISORT (DISCREET ORDINATE) MULTIPLE SCATTERING            LOP 0446
              N=ML-IMAP(IK)                                             LOP 0447
              SMTRNL=0.                                                 LOP 0448
              SMTHML=0.                                                 LOP 0449
              SMSOLL=0.                                                 LOP 0450
              DO 160 INTRVL=1,KNTRVL                                    LOP 0451
                  DTAU=1.-TRNLAY(INTRVL)                                LOP 0452
C                                                                       LOP 0453
C                 THE DISORT THERMAL (T0CMS) AND SOLAR (S0CMS)          LOP 0454
C                 SOURCE TERMS HAVE BEEN STORED IN THE UPF AND          LOP 0455
C                 UPFS ARRAYS, RESPECTIVELY.                            LOP 0456
                  THMLAY=DTAU*UPF(INTRVL,N)                             LOP 0457
                  SOLLAY=DTAU*UPFS(INTRVL,N)                            LOP 0458
                  WTTRAN=SUBINT(INTRVL)*TRNCUM(INTRVL)                  LOP 0459
                  RADCUM=RADCUM+WTTRAN*THMLAY                           LOP 0460
                  SUMMS=SUMMS+WTTRAN*SOLLAY                             LOP 0461
                  TRNCUM(INTRVL)=TRNCUM(INTRVL)*TRNLAY(INTRVL)          LOP 0462
                  TX(9)=TX(9)+SUBINT(INTRVL)*TRNCUM(INTRVL)             LOP 0463
                  SMTRNL=SMTRNL+SUBINT(INTRVL)*TRNLAY(INTRVL)           LOP 0464
                  SMTHML=SMTHML+SUBINT(INTRVL)*THMLAY                   LOP 0465
                  SMSOLL=SMSOLL+SUBINT(INTRVL)*SOLLAY                   LOP 0466
  160         CONTINUE                                                  LOP 0467
              IF(NOPRNT.EQ.-1)                                          LOP 0468
     1          WRITE(IPR1,'(I7,F8.3,3X,2F9.6,2(24X,1P2E12.5))')        LOP 0469
     2          IVX,AHT(IK+1),TX(9),SMTRNL,SMTHML,RADCUM,SMSOLL,SUMMS   LOP 0470
          ELSE                                                          LOP 0471
C                                                                       LOP 0472
C             TWO-FLUX MULTIPLE SCATTERING.  DEFINE BLACKBODY FUNCTIONS.LOP 0473
              BLAYER=BBFN(TBBY(IK),V)                                   LOP 0474
              BBOUND=BBFN(TPH(IK),V)                                    LOP 0475
              N=IMAP(IK)                                                LOP 0476
              B0=.5                                                     LOP 0477
              IF(COSBAR(N).NE.0.)B0=BETABS(QTHETA(IK),COSBAR(N))        LOP 0478
              UP=B0/PI                                                  LOP 0479
              DN=(1.-B0)/PI                                             LOP 0480
              NFRNT=N                                                   LOP 0481
              NBACK=N+1                                                 LOP 0482
              IF(QTHETA(IK).LT.0.)THEN                                  LOP 0483
                  NFRNT=NBACK                                           LOP 0484
                  NBACK=N                                               LOP 0485
              ENDIF                                                     LOP 0486
C                                                                       LOP 0487
C             LOOP OVER CORRELATED-K METHOD SUB-INTERVALS               LOP 0488
              SMTRNL=0.                                                 LOP 0489
              SMTHML=0.                                                 LOP 0490
              SMSOLL=0.                                                 LOP 0491
              DO 170 INTRVL=1,KNTRVL                                    LOP 0492
                  OMEGA=0.                                              LOP 0493
                  IF(TAUT(INTRVL,N).GT.0.)OMEGA=TSCAT(N)/TAUT(INTRVL,N) LOP 0494
C                                                                       LOP 0495
C                 THERMAL AND SOLAR SOURCE TERMS.                       LOP 0496
                  THMBND=UP*UPF(INTRVL,NFRNT)+DN*DNF(INTRVL,NFRNT)      LOP 0497
                  SOLBND=UP*UPFS(INTRVL,NFRNT)+DN*DNFS(INTRVL,NFRNT)    LOP 0498
                  THMDIF=UP*(UPF(INTRVL,NBACK)-UPF(INTRVL,NFRNT))       LOP 0499
     1                  +DN*(DNF(INTRVL,NBACK)-DNF(INTRVL,NFRNT))       LOP 0500
                  SOLDIF=UP*(UPFS(INTRVL,NBACK)-UPFS(INTRVL,NFRNT))     LOP 0501
     1                  +DN*(DNFS(INTRVL,NBACK)-DNFS(INTRVL,NFRNT))     LOP 0502
                  OPTDEP=DEPLAY(INTRVL)                                 LOP 0503
                  IF(OPTDEP.LT..02)THEN                                 LOP 0504
                      BNCOEF=OPTDEP*(1.-.5*OPTDEP)                      LOP 0505
                      DFCOEF=OPTDEP*(.5-OPTDEP/3.)                      LOP 0506
                  ELSE                                                  LOP 0507
                      BNCOEF=1.-TRNLAY(INTRVL)                          LOP 0508
                      DFCOEF=BNCOEF/OPTDEP-TRNLAY(INTRVL)               LOP 0509
                  ENDIF                                                 LOP 0510
                  THMSCT=OMEGA*(BNCOEF*THMBND+DFCOEF*THMDIF)            LOP 0511
                  THMLAY=THMSCT                                         LOP 0512
     1              +(1.-OMEGA)*(BNCOEF*BBOUND+2*DFCOEF*(BLAYER-BBOUND))LOP 0513
                  SOLLAY=OMEGA*(BNCOEF*SOLBND+DFCOEF*SOLDIF)            LOP 0514
                  WTTRAN=SUBINT(INTRVL)*TRNCUM(INTRVL)                  LOP 0515
                  RADCUM=RADCUM+WTTRAN*THMLAY                           LOP 0516
                  SUMTMS=SUMTMS+WTTRAN*THMSCT                           LOP 0517
                  SUMMS=SUMMS+WTTRAN*SOLLAY                             LOP 0518
                  TRNCUM(INTRVL)=TRNCUM(INTRVL)*TRNLAY(INTRVL)          LOP 0519
                  TX(9)=TX(9)+SUBINT(INTRVL)*TRNCUM(INTRVL)             LOP 0520
                  SMTRNL=SMTRNL+SUBINT(INTRVL)*TRNLAY(INTRVL)           LOP 0521
                  SMTHML=SMTHML+SUBINT(INTRVL)*THMLAY                   LOP 0522
                  SMSOLL=SMSOLL+SUBINT(INTRVL)*SOLLAY                   LOP 0523
  170         CONTINUE                                                  LOP 0524
              IF(NOPRNT.EQ.-1)THEN                                      LOP 0525
C                                                                       LOP 0526
C                 INTERPOLATE FLUXES IF PATH SEGMENT DOES               LOP 0527
C                 NOT TERMINATE AT A LAYER BOUNDARY.                    LOP 0528
                  FACTOR=AHT(IK+1)-ZM(NBACK)                            LOP 0529
                  IF(ABS(FACTOR).LT..0001)THEN                          LOP 0530
                      WRITE(IPR1,'(I7,F8.3,3X,2F9.6,1P8E12.5)')         LOP 0531
     1                  IVX,AHT(IK+1),TX(9),SMTRNL,                     LOP 0532
     2                  UPFLX(NBACK),DNFLX(NBACK),SMTHML,RADCUM,        LOP 0533
     3                  UPFLXS(NBACK),DNFLXS(NBACK),SMSOLL,SUMMS        LOP 0534
                  ELSE                                                  LOP 0535
                      FACTOR=FACTOR/(ZM(NFRNT)-ZM(NBACK))               LOP 0536
                      WRITE(IPR1,'(I7,F8.3,3X,2F9.6,1P8E12.5)')         LOP 0537
     1                  IVX,AHT(IK+1),TX(9),SMTRNL,                     LOP 0538
     2                  UPFLX(NBACK)+FACTOR*(UPFLX(NFRNT)-UPFLX(NBACK)),LOP 0539
     3                  DNFLX(NBACK)+FACTOR*(DNFLX(NFRNT)-DNFLX(NBACK)),LOP 0540
     4                  SMTHML,RADCUM,UPFLXS(NBACK)+FACTOR*             LOP 0541
     5                  (UPFLXS(NFRNT)-UPFLXS(NBACK)),                  LOP 0542
     6                  DNFLXS(NBACK)+FACTOR*                           LOP 0543
     7                  (DNFLXS(NFRNT)-DNFLXS(NBACK)),SMSOLL,SUMMS      LOP 0544
                  ENDIF                                                 LOP 0545
              ENDIF                                                     LOP 0546
          ENDIF                                                         LOP 0547
C                                                                       LOP 0548
C         IF TOTAL TRANSMISSION HAS DROPPED TO ZERO AND IVTEST IS       LOP 0549
C         FALSE, EXIT LAYER LOOP UNLESS THE CORRELATED-K APPROACH       LOP 0550
C         IS BEING USED AND WEIGHTING FUNCTIONS ARE BEING PRINTED.      LOP 0551
          IF(TX(9).LE.0. .AND. .NOT.IVTEST .AND.                        LOP 0552
     1      (KNTRVL.EQ.1 .OR. NOPRNT.NE.-1))GOTO190                     LOP 0553
  180 CONTINUE                                                          LOP 0554
C                                                                       LOP 0555
C     LAYER LOOP EXIT                                                   LOP 0556
  190 CONTINUE                                                          LOP 0557
      IF(MSOFF.GT.0)THEN                                                LOP 0558
C                                                                       LOP 0559
C         CALCULATE SOLAR AND THERMAL FLUXES.                           LOP 0560
          UPFLX(ML)=0.                                                  LOP 0561
          DNFLX(ML)=0.                                                  LOP 0562
          UPFLXS(ML)=0.                                                 LOP 0563
          DNFLXS(ML)=0.                                                 LOP 0564
          NTFLX(ML)=0.                                                  LOP 0565
          IF(DIS)THEN                                                   LOP 0566
              CALL MSRAD(GROUND,UANG,NSTR,V,S0,KNTRVL)                  LOP 0567
C                                                                       LOP 0568
C             CALCULATE COOLING RATES.                                  LOP 0569
CORK          CALL COOL(IVX)                                            LOP 0570
              IF(NOPRNT.EQ.-1)                                          LOP 0571
     1          WRITE(IPR1,'(/I7,0PF8.3,3X,A)')IVX,AHT(1),' 1.000000'   LOP 0572
          ELSE                                                          LOP 0573
              IF(MODTRN)THEN                                            LOP 0574
                  CALL BMFLUX(ML,KNTRVL,IEMSCT,SALB,S0)                 LOP 0575
C                                                                       LOP 0576
C                 CALCULATE COOLING RATES.                              LOP 0577
CORK              CALL COOL(IVX)                                        LOP 0578
              ELSE                                                      LOP 0579
                  CALL FLXADD(ML,IEMSCT,SALB)                           LOP 0580
              ENDIF                                                     LOP 0581
              IF(NOPRNT.EQ.-1)THEN                                      LOP 0582
                  N=IMAP(1)                                             LOP 0583
                  FACTOR=AHT(1)-ZM(N)                                   LOP 0584
                  IF(ABS(FACTOR).LE..0001)THEN                          LOP 0585
                      WRITE(IPR1,'(/(I7,0PF8.3,A,2(1P2E12.5:,24X)))')   LOP 0586
     1                  IVX,AHT(1),'    1.000000         ',             LOP 0587
     2                  UPFLX(N),DNFLX(N),UPFLXS(N),DNFLXS(N)           LOP 0588
                  ELSE                                                  LOP 0589
                      NP1=N+1                                           LOP 0590
                      FACTOR=FACTOR/(ZM(NP1)-ZM(N))                     LOP 0591
                      WRITE(IPR1,'(/(I7,0PF8.3,A,2(1P2E12.5:,24X)))')   LOP 0592
     1                  IVX,AHT(1),'    1.000000         ',             LOP 0593
     2                  UPFLX(N)+FACTOR*(UPFLX(NP1)-UPFLX(N)),          LOP 0594
     3                  DNFLX(N)+FACTOR*(DNFLX(NP1)-DNFLX(N)),          LOP 0595
     4                  UPFLXS(N)+FACTOR*(UPFLXS(NP1)-UPFLXS(N)),       LOP 0596
     5                  DNFLXS(N)+FACTOR*(DNFLXS(NP1)-DNFLXS(N))        LOP 0597
                  ENDIF                                                 LOP 0598
              ENDIF                                                     LOP 0599
          ENDIF                                                         LOP 0600
C                                                                       LOP 0601
C         VERTICAL PATH COMPLETE.  NOW PERFORM OPTICAL PATH CALCULATIONSLOP 0602
          IKMAX=IKMX                                                    LOP 0603
          MSOFF=0                                                       LOP 0604
          IK=-1                                                         LOP 0605
          GOTO40                                                        LOP 0606
      ELSEIF(GROUND)THEN                                                LOP 0607
          IF(IEMSCT.EQ.2 .AND. KNTRVL.GT.1)THEN                         LOP 0608
              IK=IKMAX+1                                                LOP 0609
              TSNREF=SUBINT(1)*S0TRN(1,IK)*TRNCUM(1)                    LOP 0610
C                                                                       LOP 0611
C             IF THE GROUND SCATTERING POINT IS IN SHADOW, THE GROUND   LOP 0612
C             REFLECTED SOLAR IRRADIANCE, TSNREF, IS ZERO.              LOP 0613
              IF(TSNREF.GT.0.)THEN                                      LOP 0614
                  DO 200 INTRVL=2,KNTRVL                                LOP 0615
                      TSNREF=TSNREF+SUBINT(INTRVL)*                     LOP 0616
     1                  S0TRN(INTRVL,IK)*TRNCUM(INTRVL)                 LOP 0617
  200             CONTINUE                                              LOP 0618
              ENDIF                                                     LOP 0619
          ENDIF                                                         LOP 0620
          IF(IMULT.NE.0)THEN                                            LOP 0621
C                                                                       LOP 0622
C             CALCULATE TRANSMITTED GROUND THERMAL AND SOLAR FLUXES     LOP 0623
              FDNTRT=0.                                                 LOP 0624
              FDNSRT=0.                                                 LOP 0625
              FACTOR=SALB/PI                                            LOP 0626
              DO 210 INTRVL=1,KNTRVL                                    LOP 0627
                  COEF=SUBINT(INTRVL)*TRNCUM(INTRVL)                    LOP 0628
                  FDNTRT=FDNTRT+COEF*DNF(INTRVL,1)                      LOP 0629
                  FDNSRT=FDNSRT+COEF*DNFS(INTRVL,1)                     LOP 0630
  210         CONTINUE                                                  LOP 0631
              FDNTRT=FACTOR*FDNTRT                                      LOP 0632
              FDNSRT=FACTOR*FDNSRT                                      LOP 0633
          ENDIF                                                         LOP 0634
      ELSEIF(IEMSCT.EQ.2)THEN                                           LOP 0635
          TSNREF=0.                                                     LOP 0636
      ENDIF                                                             LOP 0637
      if(iptx.eq.0) then
           close(aerlun)
	   iptx=1
      endif
      RETURN                                                            LOP 0638
      END                                                               LOP 0639
