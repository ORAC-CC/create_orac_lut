      SUBROUTINE TRANS(DIS,NSTR,UANG,IPH,ISOURC,IDAY,ANGLEM,            TRN 0001
     1  GROUND,LSALB,IMSMX,KNTRVL,YFLAG,XFLAG,DLIMIT,IPLOT)             TRN 0002
C                                                                       TRN 0003
C     ROUTINE TRANS CALCULATES TRANSMITTANCE AND RADIANCE VALUES        TRN 0004
C     BETWEEN IV1 AND IV2 FOR A GIVEN ATMOSPHERIC SLANT PATH.           TRN 0005
C                                                                       TRN 0006
C     DECLARE ARGUMENTS                                                 TRN 0007
C       YFLAG    Y COORDINATE FLAG FOR plot.dat FILE                    TRN 0008
C                  = "T" FOR TRANSMITTANCE                              TRN 0009
C                  = "R" FOR RADIANCE (IRRADIANCE FOR IEMSCT=3)         TRN 0010
C                  = "N" FOR NO plot.dat OUTPUT                         TRN 0011
C       XFLAG    X COORDINATE FLAG FOR plot.dat FILE                    TRN 0012
C                  = "W" FOR FREQUENCY IN WAVENUMBERS (CM-1) AND        TRN 0013
C                        RADIANCE IN W SR-1 CM-2 / CM-1                 TRN 0014
C                  = "M" FOR WAVELENGTH IN MICRONS AND                  TRN 0015
C                        RADIANCE IN W SR-1 CM-2 / MICRON               TRN 0016
C                  = "N" FOR WAVELENGTH IN MICRONS AND                  TRN 0017
C                        RADIANCE IN MICRO-WATTS SR-1 CM-2 / NANOMETER  TRN 0018
C       DLIMIT   DELIMITER CHARACTER STRING BETWEEN MODTRAN RUNS        TRN 0019
C       IPLOT    UNIT NUMBER FOR plot.dat FILE.                         TRN 0020
      LOGICAL DIS,GROUND                                                TRN 0021
      INTEGER NSTR,IPH,ISOURC,IDAY,LSALB,IMSMX,KNTRVL,IPLOT             TRN 0022
      DOUBLE PRECISION UANG                                             TRN 0023
      REAL ANGLEM                                                       TRN 0024
      CHARACTER YFLAG*1,XFLAG*1,DLIMIT*8                                TRN 0025
C                                                                       TRN 0026
C     LIST PARAMETERS                                                   TRN 0027
      INCLUDE 'PARAM.LST'                                               TRN 0028
      INTEGER NBINS,IPRINT,MAXV                                         TRN 0029
      PARAMETER(NBINS=99,IPRINT=50,MAXV=50000)                          TRN 0030
C                                                                       TRN 0031
C     INCREASE FIRST DIMENSION IN SLIT TO KMAX TO INSURE THAT ALL       TRN 0032
C     TRANSMITTANCES (TX ARRAY) ARE PASSED THROUGH THE SLIT FUNCTION.   TRN 0033
      REAL WGT(NBINS),SLIT(KMAX,NBINS)                                  TRN 0034
C                                                                       TRN 0035
C     CONVENTION                                                        TRN 0036
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         TRN 0037
C                                                                       TRN 0038
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              TRN 0039
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     TRN 0040
C                                                                       TRN 0041
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           TRN 0042
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      TRN 0043
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        TRN 0044
C                                                                       TRN 0045
      CHARACTER*8 CNAMEX                                                TRN 0046
      COMMON/NAMEX/CNAMEX(MMOLX)                                        TRN 0047
      LOGICAL IVTEST,LOOP0,TRANSM                                       TRN 0048
      INTEGER KPOINT                                                    TRN 0049
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     TRN 0050
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   TRN 0051
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   TRN 0052
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     TRN 0053
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               TRN 0054
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           TRN 0055
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     TRN 0056
      REAL TBOUND,SALB                                                  TRN 0057
      LOGICAL MODTRN                                                    TRN 0058
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   TRN 0059
     1  MODTRN                                                          TRN 0060
      INTEGER IV1,IV2,IDV,IFWHM                                         TRN 0061
      COMMON/CARD4/IV1,IV2,IDV,IFWHM                                    TRN 0062
C                                                                       TRN 0063
C       PI       THE CONSTANT PI                                        TRN 0064
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       TRN 0065
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       TRN 0066
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         TRN 0067
      REAL PI,DEG,BIGNUM,BIGEXP                                         TRN 0068
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                TRN 0069
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     TRN 0070
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                TRN 0071
      INTEGER JTURN,LJ                                                  TRN 0072
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       TRN 0073
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  TRN 0074
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      TRN 0075
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    TRN 0076
      INTEGER ICALL                                                     TRN 0077
      REAL FPHS,FALB,FORBIT                                             TRN 0078
      COMMON/ICLL/ICALL,FPHS,FALB,FORBIT                                TRN 0079
      LOGICAL LSAME                                                     TRN 0080
      COMMON/SOLAR/LSAME                                                TRN 0081
C                                                                       TRN 0082
C     DECLARE FUNCTIONS                                                 TRN 0083
      REAL BBFN,SOURCE                                                  TRN 0084
C                                                                       TRN 0085
C     DECLARE LOCAL VARIABLES                                           TRN 0086
      INTEGER I,J,IKMX,IDV5,IDVX,MXFREQ,IWIDM1,IV,IVX,IVXMAX,           TRN 0087
     1  IWRITE,IWIDTH,NWGT,NWGTM1,K,IP1,ICOUNT                          TRN 0088
      REAL YPLTMX,SALBS,RADMIN,RADMAX,EMISS,SUMA,RADSUM,SSOL,           TRN 0089
     1  BBG,SUMSSS,SUMSSR,RFLSS,RFSURF,S0,TS0,FACTOR,WNORM,             TRN 0090
     2  SUMTMS,SUMMS,UNIF,TRACE,TRANSX,SUMV,FRAC,V,CONVRT,STORE,        TRN 0091
     3  ALAM,ALTX9,SUMT,VRMAX,VRMIN,TSNOBS,TSNREF,FDNSRT,FDNTRT         TRN 0092
      CHARACTER*175 FRMT                                                TRN 0093
C                                                                       TRN 0094
C     INITIALIZE plot.dat FILE MAXIMUM.                                 TRN 0095
      YPLTMX=0.                                                         TRN 0096
C                                                                       TRN 0097
C     INITIALIZE SLIT FUNCTION ARRAY                                    TRN 0098
      SALBS=SALB                                                        TRN 0099
C                                                                       TRN 0100
C     INCREASE DIMENSION IN SLIT FUNCTION INITIALIZATION                TRN 0101
      DO 10 I=1,KMAX                                                    TRN 0102
          DO 10 J=1,NBINS                                               TRN 0103
   10 SLIT(I,J)=0.                                                      TRN 0104
C                                                                       TRN 0105
C     INITIALIZE RADIANCE MINIMUM AND MAXIMUM PARAMETERS                TRN 0106
      RADMIN=BIGNUM                                                     TRN 0107
      RADMAX=0.                                                         TRN 0108
C                                                                       TRN 0109
C     INITIALIZE GROUND EMISSIVITY (ONE MINUS GROUND ALBEDO)            TRN 0110
      EMISS=0.                                                          TRN 0111
      IF(SALB.GE.0. .AND. SALB.LE.1.)EMISS=1.-SALB                      TRN 0112
C                                                                       TRN 0113
C     STORE THE NUMBER OF PATH LAYERS IN IKMX                           TRN 0114
      IKMX=IKMAX                                                        TRN 0115
C                                                                       TRN 0116
C     INITIALIZE INTEGRATED ABSORPTION, RADIANCE, SOLAR IRRADIANCE SUMS.TRN 0117
      SUMA=0.                                                           TRN 0118
      RADSUM=0.                                                         TRN 0119
      SSOL=0.                                                           TRN 0120
C                                                                       TRN 0121
C     INITIALIZE RADIANCE/IRRADIANCE TERMS                              TRN 0122
      BBG=0.                                                            TRN 0123
      SUMSSS=0.                                                         TRN 0124
      SUMSSR=0.                                                         TRN 0125
      RFLSS=0.                                                          TRN 0126
      RFSURF=0.                                                         TRN 0127
      S0=0.                                                             TRN 0128
      TS0=0.                                                            TRN 0129
      TSNOBS=0.                                                         TRN 0130
      TSNREF=0.                                                         TRN 0131
C                                                                       TRN 0132
C     INITIALIZE INTEGRATION WEIGHTING FACTOR                           TRN 0133
      FACTOR=.5                                                         TRN 0134
C                                                                       TRN 0135
C     INITIALIZE ICOUNT, USED TO DETERMINE WHEN HEADER MUST BE PRINTED  TRN 0136
      ICOUNT=IPRINT                                                     TRN 0137
C                                                                       TRN 0138
C     DO NOT PERFORM A MODTRAN CALCULATION IF ALL SOURCES ARE CONTINUUM TRN 0139
      IF(MODTRN.AND.IV1.LT.22680)THEN                                   TRN 0140
C                                                                       TRN 0141
C         WHEN THE BAND MODEL OR LINE-BY-LINE OPTION IS USED, CALL      TRN 0142
C         ROUTINE "BMDATA" TO INITIALIZE PARAMETERS AND TO SET THE      TRN 0143
C         FREQUENCY STEP SIZE "IDVX" TO THE BAND WIDTH (1 CM-1).        TRN 0144
          IDV5=5                                                        TRN 0145
          CALL BMDATA(IV1,IFWHM,IDVX,IKMX,MXFREQ,IMSMX)                 TRN 0146
          IWIDM1=IFWHM/IDVX-1                                           TRN 0147
          IV=5*((IV1-IWIDM1)/5)                                         TRN 0148
          IF(IV.LT.0)IV=0                                               TRN 0149
          IVX=IV-IDVX                                                   TRN 0150
          IV=IV-5                                                       TRN 0151
      ELSE                                                              TRN 0152
          IDV5=IDV                                                      TRN 0153
          IDVX=IDV5                                                     TRN 0154
          IWIDM1=0                                                      TRN 0155
          IV=IV1-IDV5                                                   TRN 0156
          IVX=IV                                                        TRN 0157
          IF(IDV.LT.5)IDV=5                                             TRN 0158
      ENDIF                                                             TRN 0159
      IVXMAX=IV2+IWIDM1                                                 TRN 0160
      IF(IVXMAX.GT.MAXV)IVXMAX=MAXV                                     TRN 0161
      IWRITE=IV1+IWIDM1                                                 TRN 0162
      IWIDTH=IWIDM1+1                                                   TRN 0163
C                                                                       TRN 0164
C     PERFORM TRIANGULAR SLIT INITIALIZATION.  TRANSMITTANCES AT A      TRN 0165
C     GIVEN FREQUENCY CONTRIBUTE TO 2*IWIDTH-1 TRIANGULAR SLITS.        TRN 0166
C     THESE CONTRIBUTIONS ARE STORED IN ARRAY SLIT.  WGT IS THE         TRN 0167
C     NORMALIZED WEIGHT USED TO DEFINE THE TRIANGLE.                    TRN 0168
      NWGT=2*IWIDTH                                                     TRN 0169
      WNORM=1./(IWIDTH*IWIDTH)                                          TRN 0170
      DO 20 I=1,IWIDTH                                                  TRN 0171
          WGT(I)=I*WNORM                                                TRN 0172
   20 WGT(NWGT-I)=WGT(I)                                                TRN 0173
      NWGT=NWGT-1                                                       TRN 0174
      NWGTM1=NWGT-1                                                     TRN 0175
C                                                                       TRN 0176
C     INITIALIZE ICALL (= 0 FOR INITIAL CALL TO ROUTINE SOURCE)         TRN 0177
      ICALL=0                                                           TRN 0178
C                                                                       TRN 0179
C     INITIALIZE TRANSM (.TRUE. FOR TRANSMITTANCE ONLY CALCULATIONS)    TRN 0180
      TRANSM=.TRUE.                                                     TRN 0181
      IF(IEMSCT.EQ.1 .OR. IEMSCT.EQ.2)TRANSM=.FALSE.                    TRN 0182
C                                                                       TRN 0183
C     PRINT HEADERS                                                     TRN 0184
      IF(IEMSCT.EQ.0)THEN                                               TRN 0185
          WRITE(IPU,'((3A))')                                           TRN 0186
     1      '  FREQ TOTAL   H2O  CO2+    O3 TRACE    N2   H2O',         TRN 0187
     2            ' MOLEC   AER  HNO3 AERab  -LOG   CO2    CO',         TRN 0188
     3            '   CH4   N2O    O2   NH3    NO   NO2   SO2',         TRN 0189
     4      '  CM-1 TRANS TRANS TRANS TRANS TRANS  CONT  CONT',         TRN 0190
     5            '  SCAT TRANS TRANS TRANS TOTAL TRANS TRANS',         TRN 0191
     6            ' TRANS TRANS TRANS TRANS TRANS TRANS TRANS'          TRN 0192
      ELSEIF(IEMSCT.EQ.3)THEN                                           TRN 0193
          WRITE(IPU,'(32H  FREQ   TRANS     SOL TR  SOLAR)')            TRN 0194
      ELSE                                                              TRN 0195
          WRITE(IPU,'(3A)')'  FREQ   TOT TRANS  PTH THRML  THRML SCT',  TRN 0196
     1      '  SURF EMIS   SOL SCAT  SING SCAT  GRND RFLT  DRCT RFLT',  TRN 0197
     2      '  TOTAL RAD  REF SOL  SOL@OBS   DEPTH'                     TRN 0198
      ENDIF                                                             TRN 0199
      IF(NOPRNT.EQ.-1)THEN                                              TRN 0200
          IF(IMULT.NE.0)THEN                                            TRN 0201
              WRITE(IPR1,'((3A))')                                      TRN 0202
     1          '  FREQ     ALT        TOTAL    DELTA    THRML UP',     TRN 0203
     2          '    THRML DN   THRML SRC   THRML SUM    SOLAR UP',     TRN 0204
     3          '    SOLAR DN   SOLAR SRC   SOLAR SUM',                 TRN 0205
     4          ' (CM-1)   (KM)  BIN   TRANS    TRANS (          ',     TRN 0206
     5          '       W CM-2 / CM-1               ) (          ',     TRN 0207
     6          '       W CM-2 / CM-1               )'                  TRN 0208
          ELSEIF(IEMSCT.GT.0)THEN                                       TRN 0209
              WRITE(IPR1,'((1X,2A))')                                   TRN 0210
     1          '             ALTITUDES                 B(V,T)        ',TRN 0211
     2          '     TRANSMISSION            RADIANCE',                TRN 0212
     3          '  FREQ  BEGINNING  ENDING  INT   LAYER     BOUNDARY  ',TRN 0213
     4          ' TO BEGIN   IN LAYER    LAYER      TOTAL',             TRN 0214
     5          '(CM-1)    (KM)     (KM)         (W SR-1 CM-2 / CM-1) ',TRN 0215
     6          '                      (W SR-1 CM-2 / CM-1)'            TRN 0216
          ENDIF                                                         TRN 0217
      ENDIF                                                             TRN 0218
C                                                                       TRN 0219
C     INITIALIZE LAYER LOOP VARIABLES                                   TRN 0220
      LOOP0=.TRUE.                                                      TRN 0221
      CALL LOOP(LOOP0,IV,IVX,IKMX,MXFREQ,SUMTMS,SUMMS,TRANSM,           TRN 0222
     1  IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,SUMV,S0,FRAC,UANG,          TRN 0223
     2  DIS,NSTR,GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL)             TRN 0224
      LOOP0=.FALSE.                                                     TRN 0225
C                                                                       TRN 0226
C         DETERMINE LAYER LOOP MAXIMUM                                  TRN 0227
          IF(TRANSM)THEN                                                TRN 0228
C                                                                       TRN 0229
C             FOR TRANSMISSION CALCULATIONS, SKIP OVER LAYER LOOP IN TRATRN 0230
              IKMAX=1                                                   TRN 0231
          ELSEIF(IMULT.NE.0 .AND. .NOT. LSAME)THEN                      TRN 0232
C                                                                       TRN 0233
C             FOR MULTIPLE SCATTERING SET IKMAX TO IMSMX                TRN 0234
              IKMAX=IMSMX                                               TRN 0235
          ELSE                                                          TRN 0236
C                                                                       TRN 0237
C             IF NOT MULTIPLE SCATTERING, RESET IKMAX TO ORIGINAL VALUE TRN 0238
              IKMAX=IKMX                                                TRN 0239
          ENDIF                                                         TRN 0240
C                                                                       TRN 0241
C     END INITIALIZATION, BEGIN OF FREQUENCY LOOP                       TRN 0242
C                                                                       TRN 0243
C     "IVX" IS THE FREQUENCY AT WHICH TRANSMITTANCE WILL BE CALCULATED. TRN 0244
C     DURING THE FIRST PASS, "IVX" AND "IV" MUST BE EQUAL.              TRN 0245
   30 IVX=IVX+IDVX                                                      TRN 0246
          IF(IV.LT.IVX)THEN                                             TRN 0247
              IV=IV+IDV5                                                TRN 0248
              IVTEST=.TRUE.                                             TRN 0249
          ELSE                                                          TRN 0250
              IVTEST=.FALSE.                                            TRN 0251
          ENDIF                                                         TRN 0252
C                                                                       TRN 0253
C         SET INTERPOLATION FRACTION.                                   TRN 0254
          FRAC=FLOAT(IV-IVX)/IDV5                                       TRN 0255
          IF(ICOUNT.EQ.IPRINT)THEN                                      TRN 0256
C                                                                       TRN 0257
C             REINITIALIZE COUNTER AND PRINT HEADER                     TRN 0258
              ICOUNT=0                                                  TRN 0259
              IF(IEMSCT.EQ.0)THEN                                       TRN 0260
                  WRITE(IPR,'(1H1,/33H   FREQ WAVELENGTH  TOTAL     H2O,TRN 0261
     1              47H     CO2+     OZONE    TRACE  N2 CONT  H2O CONT, TRN 0262
     2              47H MOL SCAT  AER-HYD  HNO3    AER-HYD  INTEGRATED, TRN 0263
     3              /43H   1/CM  MICRONS    TRANS    TRANS    TRANS,    TRN 0264
     4              44H    TRANS    TRANS   TRANS    TRANS    TRANS,    TRN 0265
     5              40H     TRANS   TRANS   TAU-ABS  ABSORPTION,/)')    TRN 0266
C                 WRITE(IPR,111)'ALL MINOR SPECIES',                    TRN 0267
C    $                 (CNAMEX(KX),KX=1,11),                            TRN 0268
C    $                 (CNAMEX(KX),KX=12,MIN(22,NSPECX)),               TRN 0269
C    $                 (' ',KX=MIN(22,NSPECX)+1,22),                    TRN 0270
C    $                 ('TRANS',KX=1,12)                                TRN 0271
 111              FORMAT(8X,A,11(1X,A),/,25X,11(1X,A),                  TRN 0272
     $                 /,20X,A,8(4X,A),2(3X,A),4X,A/)                   TRN 0273
              ELSEIF(IEMSCT.EQ.1)THEN                                   TRN 0274
                  WRITE(IPR,'(A,52X,A,3(/3A))')                         TRN 0275
     1              '1','RADIANCE(WATTS/CM2-STER-XXX)',                 TRN 0276
     2              '0  FREQ   WAVLEN   SURF      PATH THERMAL   ',     TRN 0277
     3              'SCAT PART    SURFACE EMISSION   SURFACE REFLECTED',TRN 0278
     4              '     TOTAL RADIANCE   INTEGRAL    TOTAL',          TRN 0279
     5              '  (CM-1) (MICRN) ALBEDO    (CM-1)   (MICRN) ',     TRN 0280
     6              '   (CM-1)    (CM-1)   (MICRN)    (CM-1)   (MICRN)',TRN 0281
     7              '    (CM-1)   (MICRN)    (CM-1)    TRANS'           TRN 0282
              ELSEIF(IEMSCT.EQ.3)THEN                                   TRN 0283
                  WRITE(IPR,'(1H1,22X,27HIRRADIANCE (WATTS/CM2-XXXX),   TRN 0284
     1              /7H0  FREQ,T11,6HWAVLEN,T23,11HTRANSMITTED,         TRN 0285
     2              T45,5HSOLAR,T61,10HINTEGRATED,T80,5HTOTAL,          TRN 0286
     3              /2X,6H(CM-1),T10,7H(MICRN),T20,6H(CM-1),            TRN 0287
     4              T30,7H(MICRN),T40,6H(CM-1),T50,7H(MICRN),           TRN 0288
     5              T60,6HTRANS.,T70,5HSOLAR,T80,5HTRANS)')             TRN 0289
              ELSEIF(IMULT.EQ.0)THEN                                    TRN 0290
                  WRITE(IPR,'(A,45X,A,3(/3A))')                         TRN 0291
     1              '1','RADIANCE(WATTS/CM2-STER-XXX)',                 TRN 0292
     2              '0  FREQ   WAVLEN    PATH THERMAL     SURFACE ',    TRN 0293
     3              'EMISSION    SINGLE SCATTER    GROUND REFLECTED ',  TRN 0294
     4              '   TOTAL RADIANCE   INTEGRAL     TOTAL',           TRN 0295
     5              '  (CM-1) (MICRN)    (CM-1)  (MICRN)    (CM-1)',    TRN 0296
     6              '  (MICRN)    (CM-1)  (MICRN)    (CM-1)  (MICRN)',  TRN 0297
     7              '    (CM-1)  (MICRN)  (CM-1)      TRANS'            TRN 0298
              ELSEIF(DIS)THEN                                           TRN 0299
                  WRITE(IPR,'(A,45X,A,/4(/3A))')                        TRN 0300
     1              '1','RADIANCE(WATTS/CM2-STER-XXX)',                 TRN 0301
     2              '0 FREQ   WAVLEN   PATH THERMAL    SURFACE   ',     TRN 0302
     3              '   PATH SCATTERED SOLAR   GROUND REFLECTED ',      TRN 0303
     4              'RADIANCE   TOTAL RADIANCE   INTEGRAL    TOTAL',    TRN 0304
     5              '                                  EMISSION  ',     TRN 0305
     6              '   TOTAL RAD      SINGLE        TOTAL      ',      TRN 0306
     7              '  DIRECT                                TRANS',    TRN 0307
     8              ' (CM-1) (MICRN)   (CM-1)  (MICRN)   (CM-1)  ',     TRN 0308
     9              ' (CM-1)  (MICRN)   (CM-1)   (CM-1)  (MICRN)',      TRN 0309
     &              '   (CM-1)   (CM-1)  (MICRN)   (CM-1)'              TRN 0310
              ELSE                                                      TRN 0311
                  WRITE(IPR,'(A,52X,A,/4(/3A))')                        TRN 0312
     1              '1','RADIANCE(WATTS/CM2-STER-XXX)',                 TRN 0313
     2              '0 FREQ   WAVLEN   SURF          PATH THERMAL     ',TRN 0314
     3              ' SURFACE   PATH SCAT SOLAR   GROUND REFLECTED',    TRN 0315
     4              '     TOTAL RADIANCE   INTEGRAL   TOTAL',           TRN 0316
     5              '                ALBEDO       TOTAL      SCATTERED',TRN 0317
     6              ' EMISSION    TOTAL   SINGLE    TOTAL   DIRECT',    TRN 0318
     7              '                                 TRANS',           TRN 0319
     8              ' (CM-1) (MICRN)          (CM-1)  (MICRN)   (CM-1)',TRN 0320
     9              '   (CM-1)   (CM-1)   (CM-1)   (CM-1)   (CM-1)',    TRN 0321
     &              '    (CM-1)   (MICRN)   (CM-1)'                     TRN 0322
              ENDIF                                                     TRN 0323
          ENDIF                                                         TRN 0324
C                                                                       TRN 0325
C         DETERMINE LAYER LOOP MAXIMUM                                  TRN 0326
          IF(TRANSM)THEN                                                TRN 0327
C                                                                       TRN 0328
C             LAYER LOOP IS SKIPPED FOR TRANSMISSION CALCULATIONS.      TRN 0329
              IKMAX=1                                                   TRN 0330
          ELSE                                                          TRN 0331
              IF(SALBS.LT.0)CALL RHOEPS(IVX,LSALB,SALB,EMISS)           TRN 0332
              IF(IMULT.NE.0 .AND. .NOT.LSAME)THEN                       TRN 0333
C                                                                       TRN 0334
C                 FOR MULTIPLE SCATTERING SET IKMAX TO IMSMX.           TRN 0335
                  IKMAX=IMSMX                                           TRN 0336
              ELSE                                                      TRN 0337
C                                                                       TRN 0338
C                 IF NOT MULTIPLE SCATTERING, RESET IKMAX TO IKMX.      TRN 0339
                  IKMAX=IKMX                                            TRN 0340
              ENDIF                                                     TRN 0341
          ENDIF                                                         TRN 0342
          SUMV=0.                                                       TRN 0343
C                                                                       TRN 0344
C         INITIALIZE TRANSMISSION ARRAY                                 TRN 0345
          TX(1)=1.                                                      TRN 0346
          TX(2)=1.                                                      TRN 0347
          TX(3)=1.                                                      TRN 0348
          DO 40 K=4,KMAX                                                TRN 0349
   40     TX(K)=0.                                                      TRN 0350
C                                                                       TRN 0351
C         EXTRA-TERRESTRIAL SOLAR IRRADIANCE                            TRN 0352
          IF(IEMSCT.GE.2)S0=SOURCE(IVX,ISOURC,IDAY,ANGLEM)              TRN 0353
C                                                                       TRN 0354
C         CALL LAYER LOOP ROUTINE                                       TRN 0355
          CALL LOOP(LOOP0,IV,IVX,IKMX,MXFREQ,SUMTMS,SUMMS,TRANSM,       TRN 0356
     1      IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,SUMV,S0,FRAC,UANG,      TRN 0357
     2      DIS,NSTR,GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL)         TRN 0358
          IF(IEMSCT.NE.0)THEN                                           TRN 0359
              IF(IEMSCT.EQ.3)THEN                                       TRN 0360
C                                                                       TRN 0361
C                 TRANSMITTED SOLAR IRRADIANCE                          TRN 0362
                  TS0=TX(9)*S0                                          TRN 0363
              ELSE                                                      TRN 0364
                  V=IVX                                                 TRN 0365
                  IF(IVX.EQ.0)V=.5*IDVX                                 TRN 0366
C                                                                       TRN 0367
C                 THERMAL BOUNDARY EMISSION                             TRN 0368
                  IF(TBOUND.GT.0.)BBG=BBFN(TBOUND,V)*TX(9)*(1.-SALB)    TRN 0369
C                                                                       TRN 0370
C                 SURFACE REFLECTED THERMAL SCATTERED RADIANCE          TRN 0371
                  RFSURF=0.                                             TRN 0372
                  IF(IMULT.NE.0 .AND. GROUND)RFSURF=FDNTRT              TRN 0373
                  IF(IEMSCT.EQ.2)THEN                                   TRN 0374
C                                                                       TRN 0375
C                     SINGLE SCATTER SOLAR RADIANCE                     TRN 0376
                      SUMSSS=SUMSSS                                     TRN 0377
C                                                                       TRN 0378
C                     SINGLE + MULTIPLE SOLAR SCATTERED RADIANCE        TRN 0379
                      SUMSSR=SUMSSS+SUMMS                               TRN 0380
C                                                                       TRN 0381
C                     SURFACE REFLECTED SOLAR SCATTER RADIANCES         TRN 0382
                      RFLSS=0.                                          TRN 0383
                      IF(GROUND)THEN                                    TRN 0384
C                                                                       TRN 0385
C                         DIRECT TERM                                   TRN 0386
                          IF(TSNREF.GT.0.)                              TRN 0387
     1                      RFLSS=TSNREF*SALB*COS(ANGSUN/DEG)/PI        TRN 0388
C                                                                       TRN 0389
C                         SOLAR + THERMAL SURFACE REFLECTED RADIANCE    TRN 0390
                          RFSURF=RFSURF+RFLSS                           TRN 0391
                          IF(IMULT.NE.0)RFSURF=RFSURF+FDNSRT            TRN 0392
                      ENDIF                                             TRN 0393
                  ENDIF                                                 TRN 0394
              ENDIF                                                     TRN 0395
          ENDIF                                                         TRN 0396
C                                                                       TRN 0397
C         TRANSMITTANCES, IRRADIANCES [W CM-2 / CM-1], AND RADIANCES    TRN 0398
C         [W SR-1 CM-2 / CM-1] ARE TEMPORARILY STORED IN "TX" SO THAT   TRN 0399
C         THEIR CONVOLUTION OVER THE TRIANGULAR SLIT CAN BE CALCULATED. TRN 0400
          TX(1)=BBG                                                     TRN 0401
          TX(2)=UNIF                                                    TRN 0402
          TX(3)=TRACE                                                   TRN 0403
          TX(8)=SUMV                                                    TRN 0404
          TX(12)=SUMSSR                                                 TRN 0405
          TX(13)=SUMSSS                                                 TRN 0406
          TX(14)=TSNREF                                                 TRN 0407
          TX(15)=RFSURF                                                 TRN 0408
          TX(18)=TSNOBS                                                 TRN 0409
          TX(19)=RFLSS                                                  TRN 0410
          TX(20)=TS0                                                    TRN 0411
          TX(21)=S0                                                     TRN 0412
          TX(22)=SUMTMS                                                 TRN 0413
          TX(23)=SALB                                                   TRN 0414
          DO 60 K=1,KMAX                                                TRN 0415
              IP1=NWGT                                                  TRN 0416
              DO 50 I=NWGTM1,1,-1                                       TRN 0417
                  SLIT(K,IP1)=SLIT(K,I)+WGT(IP1)*TX(K)                  TRN 0418
   50         IP1=I                                                     TRN 0419
              SLIT(K,1)=WGT(1)*TX(K)                                    TRN 0420
   60     TX(K)=SLIT(K,NWGT)                                            TRN 0421
C                                                                       TRN 0422
C         CHECK IF VALUES ARE TO BE PRINTED                             TRN 0423
          IF(IVX.LT.IWRITE)GOTO30                                       TRN 0424
          IWRITE=IWRITE+IDV                                             TRN 0425
          IF(IWRITE.GT.IVXMAX)FACTOR=.5                                 TRN 0426
          ICOUNT=ICOUNT+1                                               TRN 0427
C                                                                       TRN 0428
C         RENORMALIZE IF TRIANGULAR SLIT EXTENDS TO NEGATIVE FREQUENCIESTRN 0429
          IF(IVX.LT.NWGTM1)THEN                                         TRN 0430
              STORE=1.-.5*(NWGTM1-IVX)*(NWGTM1-IVX+1)*WNORM             TRN 0431
              DO 70 K=1,KMAX                                            TRN 0432
   70         TX(K)=TX(K)/STORE                                         TRN 0433
          ENDIF                                                         TRN 0434
          BBG=TX(1)                                                     TRN 0435
          UNIF=TX(2)                                                    TRN 0436
          TRACE=TX(3)                                                   TRN 0437
          SUMV=TX(8)                                                    TRN 0438
          SUMSSR=TX(12)                                                 TRN 0439
          SUMSSS=TX(13)                                                 TRN 0440
          TSNREF=TX(14)                                                 TRN 0441
          RFSURF=TX(15)                                                 TRN 0442
          TSNOBS=TX(18)                                                 TRN 0443
          RFLSS=TX(19)                                                  TRN 0444
          TS0=TX(20)                                                    TRN 0445
          S0=TX(21)                                                     TRN 0446
          SUMTMS=TX(22)                                                 TRN 0447
          SALB=TX(23)                                                   TRN 0448
          V=FLOAT(IVX-IWIDM1)                                           TRN 0449
          ALAM=10000./(V+.000001)                                       TRN 0450
          SUMA=SUMA+FACTOR*IDV*(1.0-TX(9))                              TRN 0451
C                                                                       TRN 0452
C         ALTX9 IS NOW OUTPUT USING AN F FORMAT WHEN IEMSCT = 1 OR 2,   TRN 0453
C         SO THE MAXIMUM IS REDUCED TO 999.999 (STILL ABSURDLY LARGE)   TRN 0454
          ALTX9=999.9                                                   TRN 0455
          IF(TX(9).GT.0.)ALTX9=-LOG(TX(9))                              TRN 0456
          GOTO(80,90,90,100),IEMSCT+1                                   TRN 0457
C                                                                       TRN 0458
C         TRANSMITTANCE ONLY                                            TRN 0459
   80     CONTINUE                                                      TRN 0460
          TX(7)=TX(7)*TX(16)                                            TRN 0461
          WRITE(IPR,'(F8.0,F8.3,11F9.4,F12.3)')V,ALAM,TX(9),TX(17),     TRN 0462
     1      UNIF,TX(31),TRACE,TX(4),TX(5),TX(6),TX(7),TX(11),TX(10),SUMATRN 0463
C                                                                       TRN 0464
C         KEEP 22 ENTRIES WITHIN 132 CHARACTERS                         TRN 0465
          FRMT(1:47)=  '(1X,I5,1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,'TRN 0466
          FRMT(48:95)='1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,'TRN 0467
          FRMT(96:135)=       '1X,F6.3,1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,'TRN 0468
          FRMT(136:175)=      '1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3,1X,F5.3)'TRN 0469
          IF(TX(9) .LT..99995)FRMT( 14: 14)='4'                         TRN 0470
          IF(TX(17).LT..99995)FRMT( 22: 22)='4'                         TRN 0471
          IF(UNIF  .LT..99995)FRMT( 30: 30)='4'                         TRN 0472
          IF(TX(31).LT..99995)FRMT( 38: 38)='4'                         TRN 0473
          IF(TRACE .LT..99995)FRMT( 46: 46)='4'                         TRN 0474
          IF(TX( 4).LT..99995)FRMT( 54: 54)='4'                         TRN 0475
          IF(TX( 5).LT..99995)FRMT( 62: 62)='4'                         TRN 0476
          IF(TX( 6).LT..99995)FRMT( 70: 70)='4'                         TRN 0477
          IF(TX( 7).LT..99995)FRMT( 78: 78)='4'                         TRN 0478
          IF(TX(11).LT..99995)FRMT( 86: 86)='4'                         TRN 0479
          IF(TX(10).LT..99995)FRMT( 94: 94)='4'                         TRN 0480
          IF(ALTX9 .LT..99995)THEN                                      TRN 0481
              FRMT(102:102)='4'                                         TRN 0482
          ELSEIF(ALTX9 .LT.9.9995)THEN                                  TRN 0483
              FRMT(102:102)='3'                                         TRN 0484
          ELSEIF(ALTX9 .LT.99.995)THEN                                  TRN 0485
              FRMT(102:102)='2'                                         TRN 0486
          ELSE                                                          TRN 0487
              FRMT(102:102)='1'                                         TRN 0488
          ENDIF                                                         TRN 0489
          IF(TX(36).LT..99995)FRMT(110:110)='4'                         TRN 0490
          IF(TX(44).LT..99995)FRMT(118:118)='4'                         TRN 0491
          IF(TX(46).LT..99995)FRMT(126:126)='4'                         TRN 0492
          IF(TX(47).LT..99995)FRMT(134:134)='4'                         TRN 0493
          IF(TX(50).LT..99995)FRMT(142:142)='4'                         TRN 0494
          IF(TX(52).LT..99995)FRMT(150:150)='4'                         TRN 0495
          IF(TX(54).LT..99995)FRMT(158:158)='4'                         TRN 0496
          IF(TX(64).LT..99995)FRMT(166:166)='4'                         TRN 0497
          IF(TX(65).LT..99995)FRMT(174:174)='4'                         TRN 0498
          WRITE(IPU,FMT=FRMT)INT(V+.5),TX(9),TX(17),UNIF,TX(31),        TRN 0499
     1      TRACE,TX(4),TX(5),TX(6),TX(7),TX(11),TX(10),ALTX9,TX(36),   TRN 0500
     2      TX(44),TX(46),TX(47),TX(50),TX(52),TX(54),TX(64),TX(65)     TRN 0501
          GOTO110                                                       TRN 0502
C                                                                       TRN 0503
C         RADIANCE PATHS                                                TRN 0504
   90     CONTINUE                                                      TRN 0505
C                                                                       TRN 0506
C         CONVRT IS THE CONVERSION FROM (W SR-1 CM-2 / CM-1)            TRN 0507
C         TO (W SR-1 CM-2 / MICRON).                                    TRN 0508
          CONVRT=1.E-4*V**2                                             TRN 0509
          IF(V.EQ.0.)CONVRT=1.E-4*(.5*IDVX)**2                          TRN 0510
          IF(IEMSCT.EQ.1)THEN                                           TRN 0511
C                                                                       TRN 0512
C             SUMT IS THE TOTAL SPECTRAL RADIANCE, I.E. SUMT            TRN 0513
C             EQUALS THE SUM OF THE DIRECT + MULTIPLY SCATTERED         TRN 0514
C             THERMAL PATH RADIANCE (SUMV), THE SURFACE EMISSION        TRN 0515
C             (BBG), AND THE REFLECTED SURFACE TERM (RFSURF).  SUMTMS   TRN 0516
C             IS THE MULTIPLE SCATTERING CONTRIBUTION TO SUMV.  EACH OF TRN 0517
C             THESE TERMS HAS UNITS (W SR-1 CM-2 / CM-1).  RADSUM IS    TRN 0518
C             THE SPECTRALLY INTEGRATED TOTAL RADIANCE (W SR-1 CM-2).   TRN 0519
              SUMT=SUMV+BBG+RFSURF                                      TRN 0520
              RADSUM=RADSUM+IDV*FACTOR*SUMT                             TRN 0521
              IF(DIS)THEN                                               TRN 0522
                  WRITE(IPR,'(F8.0,F8.3,F7.3,1P,2E10.2,10X,7E10.2,0P,   TRN 0523
     1              F9.5)')V,ALAM,SALB,SUMV,CONVRT*SUMV,BBG,CONVRT*BBG, TRN 0524
     2              RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TX(9)  TRN 0525
              ELSE                                                      TRN 0526
                  WRITE(IPR,'(F8.0,F8.3,F7.3,1P,10E10.2,0P,F9.5)')      TRN 0527
     1              V,ALAM,SALB,SUMV,CONVRT*SUMV,SUMTMS,BBG,CONVRT*BBG, TRN 0528
     2              RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TX(9)  TRN 0529
              ENDIF                                                     TRN 0530
              WRITE(IPU,'(0PF7.0,F11.8,1P3E11.4,11X,2(11X,E11.4),18X,   TRN 0531
     1          0PF8.3)')V,TX(9),SUMV,SUMTMS,BBG,RFSURF,SUMT,ALTX9      TRN 0532
          ELSE                                                          TRN 0533
C                                                                       TRN 0534
C             SUMT IS THE TOTAL SPECTRAL RADIANCE, I.E. SUMT            TRN 0535
C             EQUALS THE SUM OF THE DIRECT + MULTIPLY SCATTERED         TRN 0536
C             THERMAL PATH RADIANCE (SUMV), THE SURFACE EMISSION        TRN 0537
C             (BBG), THE SINGLE + MULTIPLE SOLAR SCATTERED RADIANCE     TRN 0538
C             (SUMSSR) AND THE REFLECTED SURFACE TERM (RFSURF).  EACH   TRN 0539
C             OF THESE TERMS HAS UNITS (W SR-1 CM-2 / CM-1).  RADSUM IS TRN 0540
C             THE SPECTRALLY INTEGRATED TOTAL RADIANCE (W SR-1 CM-2).   TRN 0541
              SUMT=SUMV+BBG+SUMSSR+RFSURF                               TRN 0542
              RADSUM=RADSUM+IDV*FACTOR*SUMT                             TRN 0543
              IF(IMULT.EQ.0)THEN                                        TRN 0544
                  WRITE(IPR,'(F8.0,F8.3,5(1PE10.2,E9.2),E10.2,0PF9.5)') TRN 0545
     1              V,ALAM,SUMV,CONVRT*SUMV,BBG,CONVRT*BBG,             TRN 0546
     2              SUMSSR,CONVRT*SUMSSR,RFSURF,CONVRT*RFSURF,          TRN 0547
     3              SUMT,CONVRT*SUMT,RADSUM,TX(9)                       TRN 0548
              ELSEIF(DIS)THEN                                           TRN 0549
                  WRITE(IPR,'(F7.0,F8.3,1P11E9.2,E10.2,0PF8.5)')        TRN 0550
     1              V,ALAM,SUMV,CONVRT*SUMV,BBG,SUMSSR,CONVRT*SUMSSR,   TRN 0551
     2              SUMSSS,RFSURF,CONVRT*RFSURF,RFLSS,                  TRN 0552
     3              SUMT,CONVRT*SUMT,RADSUM,TX(9)                       TRN 0553
              ELSE                                                      TRN 0554
                  WRITE(IPR,'(F7.0,F8.3,F7.3,1P8E9.2,3E10.3,0PF8.5)')   TRN 0555
     1              V,ALAM,SALB,SUMV,CONVRT*SUMV,SUMTMS,BBG,SUMSSR,     TRN 0556
     2              SUMSSS,RFSURF,RFLSS,SUMT,CONVRT*SUMT,RADSUM,TX(9)   TRN 0557
              ENDIF                                                     TRN 0558
              WRITE(IPU,'(0PF7.0,F11.8,1P8E11.4,2E9.2,0PF8.3)')         TRN 0559
     1          V,TX(9),SUMV,SUMTMS,BBG,SUMSSR,SUMSSS,                  TRN 0560
     2          RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9                   TRN 0561
          ENDIF                                                         TRN 0562
          GOTO110                                                       TRN 0563
C                                                                       TRN 0564
C         DIRECTLY TRANSMITTED SOLAR IRRADIANCE [WATTS/(CM2 MICROMETER)]TRN 0565
  100     CONTINUE                                                      TRN 0566
C                                                                       TRN 0567
C         CONVRT IS THE CONVERSION FROM W/CM2-(CM-1) TO W/CM2-UM        TRN 0568
          CONVRT=1.E-4*V**2                                             TRN 0569
          IF(V.EQ.0)CONVRT=1.E-4*(.5*IDVX)**2                           TRN 0570
C                                                                       TRN 0571
C         RADSUM IS THE INTEGRATED TRANSMITTED SOLAR IRRADIANCE AND SSOLTRN 0572
C         IS THE INTEGRATED EXTRA-TERRESTRIAL SOLAR IRRADIANCE (W/CM2). TRN 0573
          RADSUM=RADSUM+TS0*IDV*FACTOR                                  TRN 0574
          SSOL=SSOL+S0*IDV*FACTOR                                       TRN 0575
          WRITE(IPR,'(F8.0,F8.3,1P6E10.2,0PF9.4)')                      TRN 0576
     1      V,ALAM,TS0,CONVRT*TS0,S0,CONVRT*S0,RADSUM,SSOL,TX(9)        TRN 0577
          WRITE(IPU,'(F7.0,F8.4,1P2E9.2,T96,E10.3)')V,TX(9),TS0,S0,ALTX9TRN 0578
          SUMT=TS0                                                      TRN 0579
  110     CONTINUE                                                      TRN 0580
C                                                                       TRN 0581
C         WRITE OUT plot.dat FILE DATA                                  TRN 0582
          IF(YFLAG.EQ.'T')THEN                                          TRN 0583
              IF(TX(9).GT.YPLTMX)YPLTMX=TX(9)                           TRN 0584
              IF(XFLAG.EQ.'N')THEN                                      TRN 0585
                  WRITE(IPLOT,'(0PF15.3,F15.8)')1000*ALAM,TX(9)         TRN 0586
              ELSEIF(XFLAG.EQ.'M')THEN                                  TRN 0587
                  WRITE(IPLOT,'(0PF15.6,F15.8)')ALAM,TX(9)              TRN 0588
              ELSE                                                      TRN 0589
                  WRITE(IPLOT,'(0PF15.0,F15.8)')V,TX(9)                 TRN 0590
              ENDIF                                                     TRN 0591
          ELSEIF(YFLAG.EQ.'R')THEN                                      TRN 0592
              IF(XFLAG.EQ.'N')THEN                                      TRN 0593
                  CONVRT=1000.*CONVRT                                   TRN 0594
                  IF(CONVRT*SUMT.GT.YPLTMX)YPLTMX=CONVRT*SUMT           TRN 0595
                  WRITE(IPLOT,'(0PF15.3,1PE15.5)')1000.*ALAM,CONVRT*SUMTTRN 0596
              ELSEIF(XFLAG.EQ.'M')THEN                                  TRN 0597
                  IF(CONVRT*SUMT.GT.YPLTMX)YPLTMX=CONVRT*SUMT           TRN 0598
                  WRITE(IPLOT,'(0PF15.6,1PE15.5)')ALAM,CONVRT*SUMT      TRN 0599
              ELSE                                                      TRN 0600
                  IF(SUMT.GT.YPLTMX)YPLTMX=SUMT                         TRN 0601
                  WRITE(IPLOT,'(0PF15.0,1PE15.5)')V,SUMT                TRN 0602
              ENDIF                                                     TRN 0603
          ENDIF                                                         TRN 0604
          IF(IEMSCT.NE.0)THEN                                           TRN 0605
              IF(SUMT.GE.RADMAX)THEN                                    TRN 0606
                  VRMAX=V                                               TRN 0607
                  RADMAX=SUMT                                           TRN 0608
              ENDIF                                                     TRN 0609
              IF(SUMT.LE.RADMIN)THEN                                    TRN 0610
                  VRMIN=V                                               TRN 0611
                  RADMIN=SUMT                                           TRN 0612
              ENDIF                                                     TRN 0613
          ENDIF                                                         TRN 0614
          FACTOR=1.                                                     TRN 0615
      IF(IWRITE.LE.IVXMAX)GOTO30                                        TRN 0616
C                                                                       TRN 0617
C     END OF FREQUENCY LOOP                                             TRN 0618
      IVX=INT(V+.5)                                                     TRN 0619
C                                                                       TRN 0620
C     PRINT plot.dat FILE DELIMITER.                                    TRN 0621
      IF(YFLAG.EQ.'T')THEN                                              TRN 0622
          IF (DLIMIT.EQ.'        ')THEN                                 TRN 0623
                  WRITE(IPLOT,*)                                        TRN 0624
          ELSEIF(XFLAG.EQ.'N')THEN                                      TRN 0625
              WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                      TRN 0626
     1          '   THE MAXIMUM TRANSMITTANCE IS',                      TRN 0627
     2          YPLTMX,' (WAVELENGTHS IN NANOMETERS)'                   TRN 0628
          ELSEIF(XFLAG.EQ.'M')THEN                                      TRN 0629
              WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                      TRN 0630
     1          '   THE MAXIMUM TRANSMITTANCE IS',                      TRN 0631
     2          YPLTMX,' (WAVELENGTHS IN MICRONS)'                      TRN 0632
          ELSE                                                          TRN 0633
              WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                      TRN 0634
     1          '   THE MAXIMUM TRANSMITTANCE IS',                      TRN 0635
     2          YPLTMX,' (FREQUENCIES IN CM-1)'                         TRN 0636
          ENDIF                                                         TRN 0637
      ELSEIF(YFLAG.EQ.'R')THEN                                          TRN 0638
          IF(IEMSCT.EQ.3)THEN                                           TRN 0639
              IF (DLIMIT.EQ.'        ')THEN                             TRN 0640
                  WRITE(IPLOT,*)                                        TRN 0641
              ELSEIF(XFLAG.EQ.'N')THEN                                  TRN 0642
                  WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                  TRN 0643
     1              '   THE MAXIMUM TRANSMITTED SOLAR IRRADIANCE IS',   TRN 0644
     2              YPLTMX,' MICRO-WATTS CM-2 / NANOMETER'              TRN 0645
              ELSEIF(XFLAG.EQ.'M')THEN                                  TRN 0646
                  WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                  TRN 0647
     1              '   THE MAXIMUM TRANSMITTED SOLAR IRRADIANCE IS',   TRN 0648
     2              YPLTMX,' W CM-2 / MICRON'                           TRN 0649
              ELSE                                                      TRN 0650
                  WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                  TRN 0651
     1              '   THE MAXIMUM TRANSMITTED SOLAR IRRADIANCE IS',   TRN 0652
     2              YPLTMX,' W CM-2 / CM-1'                             TRN 0653
              ENDIF                                                     TRN 0654
          ELSE                                                          TRN 0655
              IF (DLIMIT.EQ.'        ')THEN                             TRN 0656
                  WRITE(IPLOT,*)                                        TRN 0657
              ELSEIF(XFLAG.EQ.'N')THEN                                  TRN 0658
                  WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                  TRN 0659
     1              '   THE MAXIMUM SPECTRAL RADIANCE IS',              TRN 0660
     2              YPLTMX,' MICRO-WATTS SR-1 CM-2 / NANOMETER'         TRN 0661
              ELSEIF(XFLAG.EQ.'M')THEN                                  TRN 0662
                  WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                  TRN 0663
     1              '   THE MAXIMUM SPECTRAL RADIANCE IS',              TRN 0664
     2              YPLTMX,' W SR-1 CM-2 / MICRON'                      TRN 0665
              ELSE                                                      TRN 0666
                  WRITE(IPLOT,'(2A,1PE15.5,A)')DLIMIT,                  TRN 0667
     1              '   THE MAXIMUM SPECTRAL RADIANCE IS',              TRN 0668
     2              YPLTMX,' W SR-1 CM-2 / CM-1'                        TRN 0669
              ENDIF                                                     TRN 0670
          ENDIF                                                         TRN 0671
      ENDIF                                                             TRN 0672
CORK  IF(IMULT.NE.0)THEN                                                TRN 0673
C                                                                       TRN 0674
C         INTEGRATED COOLING RATES.                                     TRN 0675
CORK      WRITE(IPR,'(/A)')'1 MULTIPLE SCATTERING CALCULATION RESULTS:' TRN 0676
CORK      CALL WTCOOL                                                   TRN 0677
CORK      WRITE(IPR,'(/A)')'1 NO SCATTERING ATMOSPHERE/SURFACE RESULTS:'TRN 0678
CORK      CALL TCOOL                                                    TRN 0679
CORK      CALL WTCOOL                                                   TRN 0680
CORK  ENDIF                                                             TRN 0681
      WRITE(IPR,'(27H0INTEGRATED ABSORPTION FROM,I6,3H TO,I6,7H CM-1 =, TRN 0682
     1  F10.2,5H CM-1,/24H AVERAGE TRANSMITTANCE =,F6.4,/)')            TRN 0683
     2  IV1,IVX,SUMA,1.-SUMA/(IVX-IV1)                                  TRN 0684
      IF(IEMSCT.EQ.3)THEN                                               TRN 0685
          WRITE(IPR,'(24h0INTEGRATED IRRADIANCE =,1PE11.3,              TRN 0686
     1      11h WATTS CM-2,/24h MINIMUM IRRADIANCE    =,E11.3,          TRN 0687
     2      14h WATTS CM-2 AT,0PF11.1,5h CM-1,/10h  MAXIMUM ,           TRN 0688
     3      14h RADIANCE    =,1PE11.3,24h WATTS CM-2  (CM-1)-1 AT,      TRN 0689
     4      0PF11.1,5H CM-1)')RADSUM,RADMIN,VRMIN,RADMAX,VRMAX          TRN 0690
      ELSEIF(IEMSCT.NE.0)THEN                                           TRN 0691
          WRITE(IPR,'(/A,1PE14.6,A,/(A,1PE14.6,A,0PF10.0,A))')          TRN 0692
     1      ' INTEGRATED TOTAL RADIANCE =',RADSUM,' WATTS CM-2 STER-1', TRN 0693
     2      ' MINIMUM SPECTRAL RADIANCE =',RADMIN,                      TRN 0694
     3                   ' WATTS CM-2 STER-1 / CM-1  AT',VRMIN,' CM-1', TRN 0695
     4      ' MAXIMUM SPECTRAL RADIANCE =',RADMAX,                      TRN 0696
     5                   ' WATTS CM-2 STER-1 / CM-1  AT',VRMAX,' CM-1'  TRN 0697
          IF(SALBS.GE.0)THEN                                            TRN 0698
              WRITE(IPR,'(23H BOUNDARY TEMPERATURE =,F11.2,             TRN 0699
     1          2H K,/22H BOUNDARY EMISSIVITY =,F12.3)')TBOUND,EMISS    TRN 0700
          ELSE                                                          TRN 0701
              WRITE(IPR,'(23H BOUNDARY TEMPERATURE =,F11.2)')TBOUND     TRN 0702
          ENDIF                                                         TRN 0703
      ENDIF                                                             TRN 0704
      RETURN                                                            TRN 0705
      END                                                               TRN 0706
