      SUBROUTINE  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,     DIS 0001
     $                    WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,     DIS 0002
     $                    NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,     DIS 0003
     $                    PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP,       DIS 0004
     $                    TTEMP, TEMIS, DELTAM, PLANK, ONLYFL,          DIS 0005
     $                    ACCUR, PRNT, HEADER, MAXCLY, MAXULV,          DIS 0006
     $                    MAXUMU, MAXCMU, MAXPHI, RFLDIR, RFLDN,        DIS 0007
     $              FLDN, FLUP, DFDT, UAVG, UU, U0U, ALBMED, TRNMED,    DIS 0008
     $                    MSFLAG, WN0, S0CMS, T0CMS,SFDSRT,SFDTRT )     DIS 0009
C=====>        LOWER CASE VARIABLES ADDED                               DIS 0010
              IMPLICIT DOUBLE PRECISION (A-H, O-Z)                      DIS 0011
                                                                        DIS 0012
                                                                        DIS 0013
C **********************************************************************DIS 0014
C       PLANE-PARALLEL DISCRETE ORDINATES RADIATIVE TRANSFER PROGRAM    DIS 0015
C             ( SEE DISORT.DOC FOR COMPLETE DOCUMENTATION )             DIS 0016
C **********************************************************************DIS 0017
                                                                        DIS 0018
C+---------------------------------------------------------------------+DIS 0019
C------------------    I/O VARIABLE SPECIFICATIONS     -----------------DIS 0020
C+---------------------------------------------------------------------+DIS 0021
C+---------------------------------------------------------------------+DIS 0022
C  LOCAL SYMBOLIC DIMENSIONS (HAVE BIG EFFECT ON STORAGE REQUIREMENTS): DIS 0023
                                                                        DIS 0024
C       MXCLY  = MAX NO. OF COMPUTATIONAL LAYERS                        DIS 0025
C       MXULV  = MAX NO. OF OUTPUT LEVELS                               DIS 0026
C       MXCMU  = MAX NO. OF COMPUTATION POLAR ANGLES                    DIS 0027
C       MXUMU  = MAX NO. OF OUTPUT POLAR ANGLES                         DIS 0028
C       MXPHI  = MAX NO. OF OUTPUT AZIMUTHAL ANGLES                     DIS 0029
                                                                        DIS 0030
CJ 7/26 NEED TO MAKE MXCLY AND MXULV CONSISTENT WITH THE                DIS 0031
CJ      PARAMETERIZATION IN MODTRAN3 BY INCLUDING PARAM.LST.            DIS 0032
C                                                                       DIS 0033
CJ      PARAMETER ( MXCLY = 34, MXULV = 34, MXCMU = 16, MXUMU = 1,      DIS 0034
CJ     $            MXPHI = 1, MI = MXCMU/2, MI9M2 = 9*MI-2,            DIS 0035
CJ     $            NNLYRI = MXCMU*MXCLY )                              DIS 0036
C                                                                       DIS 0037
         INCLUDE 'PARAM.LST'                                            DIS 0038
C                                                                       DIS 0039
      PARAMETER (MXCLY =LAYDIM, MXULV =LAYDIM+1, MXCMU =16, MXUMU=1,    DIS 0040
     $            MXPHI = 1, MI = MXCMU/2, MI9M2 = 9*MI-2,              DIS 0041
     $            NNLYRI = MXCMU*MXCLY )                                DIS 0042
C                                                                       DIS 0043
CJ ^                                                                    DIS 0044
                                                                        DIS 0045
C+---------------------------------------------------------------------+DIS 0046
                                                                        DIS 0047
      CHARACTER  HEADER*127                                             DIS 0048
      LOGICAL  DELTAM, LAMBER, PLANK, ONLYFL, PRNT(7), USRANG, USRTAU   DIS 0049
      INTEGER  IBCND, MAXCLY, MAXUMU, MAXULV, MAXCMU, MAXPHI, NLYR,     DIS 0050
     $         NUMU, NSTR, NPHI, NTAU                                   DIS 0051
      REAL*8     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,   DIS 0052
     $         HL( 0:MAXCMU ), PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ), DIS 0053
     $         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP, DIS 0054
     $         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )      DIS 0055
                                                                        DIS 0056
      REAL*8     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),     DIS 0057
     $         UAVG( MAXULV ), DFDT( MAXULV ), U0U( MAXUMU, MAXULV ),   DIS 0058
     $         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),          DIS 0059
     $         TRNMED( MAXUMU ),                                        DIS 0060
     $         S0CMS( MAXUMU,MAXULV ), T0CMS( MAXUMU,MAXULV ),          DIS 0061
     $         FDNSRT, FDNTRT, WN0                                      DIS 0062
                                                                        DIS 0063
          REAL SFDSRT,SFDTRT                                            DIS 0064
                                                                        DIS 0065
C=====>        LOWER CASE VARIABLES ADDED                               DIS 0066
                                                                        DIS 0067
C+---------------------------------------------------------------------+DIS 0068
C      ROUTINES CALLED (IN ORDER):  SLFTST, ZEROAL, CHEKIN, SETDIS,     DIS 0069
C                                   PRTINP, ALBTRN, LEPOLY, SURFAC,     DIS 0070
C                                   SOLEIG, UPBEAM, UPISOT, TERPEV,     DIS 0071
C                                   TERPSO, SETMTX, SOLVE0, FLUXES,     DIS 0072
C                                   USRINT, PRAVIN, PRTINT              DIS 0073
C+---------------------------------------------------------------------+DIS 0074
                                                                        DIS 0075
C  INDEX CONVENTIONS (FOR ALL DO-LOOPS AND ALL VARIABLE DESCRIPTIONS):  DIS 0076
                                                                        DIS 0077
C     IU     :  FOR USER POLAR ANGLES                                   DIS 0078
                                                                        DIS 0079
C  IQ,JQ,KQ  :  FOR COMPUTATIONAL POLAR ANGLES ('QUADRATURE ANGLES')    DIS 0080
                                                                        DIS 0081
C   IQ/2     :  FOR HALF THE COMPUTATIONAL POLAR ANGLES (JUST THE ONES  DIS 0082
C               IN EITHER 0-90 DEGREES, OR 90-180 DEGREES)              DIS 0083
                                                                        DIS 0084
C     J      :  FOR USER AZIMUTHAL ANGLES                               DIS 0085
                                                                        DIS 0086
C     K,L    :  FOR LEGENDRE EXPANSION COEFFICIENTS OR, ALTERNATIVELY,  DIS 0087
C               SUBSCRIPTS OF ASSOCIATED LEGENDRE POLYNOMIALS           DIS 0088
                                                                        DIS 0089
C     LU     :  FOR USER LEVELS                                         DIS 0090
                                                                        DIS 0091
C     LC     :  FOR COMPUTATIONAL LAYERS (EACH HAVING A DIFFERENT       DIS 0092
C               SINGLE-SCATTER ALBEDO AND/OR PHASE FUNCTION)            DIS 0093
                                                                        DIS 0094
C    LEV     :  FOR COMPUTATIONAL LEVELS                                DIS 0095
                                                                        DIS 0096
C    MAZIM   :  FOR AZIMUTHAL COMPONENTS IN FOURIER COSINE EXPANSION    DIS 0097
C               OF INTENSITY AND PHASE FUNCTION                         DIS 0098
                                                                        DIS 0099
C+---------------------------------------------------------------------+DIS 0100
C               I N T E R N A L    V A R I A B L E S                    DIS 0101
                                                                        DIS 0102
C   AMB(IQ/2,IQ/2)    FIRST MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM DIS 0103
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN SOLEIG)   DIS 0104
                                                                        DIS 0105
C   APB(IQ/2,IQ/2)    SECOND MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEMDIS 0106
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN SOLEIG)   DIS 0107
                                                                        DIS 0108
C   ARRAY(IQ,IQ)      SCRATCH MATRIX FOR SOLEIG, UPBEAM AND UPISOT      DIS 0109
C                     (SEE EACH SUBROUTINE FOR DEFINITION)              DIS 0110
                                                                        DIS 0111
C   B()               RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO    DIS 0112
C                     SOLVE0,1;  RETURNS AS SOLUTION VECTOR             DIS 0113
C                     VECTOR  L, THE CONSTANTS OF INTEGRATION           DIS 0114
                                                                        DIS 0115
C   BDR(IQ/2,0:IQ/2)  BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A  DIS 0116
C                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS    DIS 0117
C                     REFERS TO A COMPUTATIONAL ANGLE.  SECOND INDEX:   DIS 0118
C                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE UMU0;      DIS 0119
C                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.     DIS 0120
                                                                        DIS 0121
C   BEM(IQ/2)         BOTTOM-BOUNDARY DIRECTIONAL EMISSIVITY AT COMPU-  DIS 0122
C                     TATIONAL ANGLES.                                  DIS 0123
                                                                        DIS 0124
C   BPLANK            INTENSITY EMITTED FROM BOTTOM BOUNDARY            DIS 0125
                                                                        DIS 0126
C   CBAND()           MATRIX OF LEFT-HAND SIDE OF THE LINEAR SYSTEM     DIS 0127
C                     EQ. SC(5), SCALED BY EQ. SC(12);  IN BANDED       DIS 0128
C                     FORM REQUIRED BY LINPACK SOLUTION ROUTINES        DIS 0129
                                                                        DIS 0130
C   CC(IQ,IQ)         C-SUB-IJ IN EQ. SS(5)                             DIS 0131
                                                                        DIS 0132
C   CMU(IQ)           COMPUTATIONAL POLAR ANGLES (GAUSSIAN)             DIS 0133
                                                                        DIS 0134
C   CWT(IQ)           QUADRATURE WEIGHTS CORRESPONDING TO CMU           DIS 0135
                                                                        DIS 0136
C   DELM0             KRONECKER DELTA, DELTA-SUB-M0, WHERE M = MAZIM    DIS 0137
C                     IS THE NUMBER OF THE FOURIER COMPONENT IN THE     DIS 0138
C                     AZIMUTH COSINE EXPANSION                          DIS 0139
                                                                        DIS 0140
C   EMU(IU)           BOTTOM-BOUNDARY DIRECTIONAL EMISSIVITY AT USER    DIS 0141
C                     ANGLES.                                           DIS 0142
                                                                        DIS 0143
C   EVAL(IQ)          TEMPORARY STORAGE FOR EIGENVALUES OF EQ. SS(12)   DIS 0144
                                                                        DIS 0145
C   EVECC(IQ,IQ)      COMPLETE EIGENVECTORS OF SS(7) ON RETURN FROM     DIS 0146
C                     SOLEIG; STORED PERMANENTLY IN  GC                 DIS 0147
                                                                        DIS 0148
C   EXPBEA(LC)        TRANSMISSION OF DIRECT BEAM IN DELTA-M OPTICAL    DIS 0149
C                     DEPTH COORDINATES                                 DIS 0150
                                                                        DIS 0151
C   FLYR(LC)          TRUNCATED FRACTION IN DELTA-M METHOD              DIS 0152
                                                                        DIS 0153
C   GL(K,LC)          PHASE FUNCTION LEGENDRE POLYNOMIAL EXPANSION      DIS 0154
C                     COEFFICIENTS, CALCULATED FROM PMOM BY             DIS 0155
C                     INCLUDING SINGLE-SCATTERING ALBEDO, FACTOR        DIS 0156
C                     2K+1, AND (IF DELTAM=TRUE) THE DELTA-M            DIS 0157
C                     SCALING                                           DIS 0158
                                                                        DIS 0159
C   GC(IQ,IQ,LC)      EIGENVECTORS AT POLAR QUADRATURE ANGLES,          DIS 0160
C                     G  IN EQ. SC(1)                                   DIS 0161
                                                                        DIS 0162
C   GU(IU,IQ,LC)      EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES    DIS 0163
C                     ( G  IN EQS. SC(3) AND S1(8-9), I.E.              DIS 0164
C                       G WITHOUT THE L FACTOR )                        DIS 0165
                                                                        DIS 0166
C   HLPR()            LEGENDRE COEFFICIENTS OF BOTTOM BIDIRECTIONAL     DIS 0167
C                     REFLECTIVITY (AFTER INCLUSION OF 2K+1 FACTOR)     DIS 0168
                                                                        DIS 0169
C   IPVT(LC*IQ)       INTEGER VECTOR OF PIVOT INDICES FOR LINPACK       DIS 0170
C                     ROUTINES                                          DIS 0171
                                                                        DIS 0172
C   KK(IQ,LC)         EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)         DIS 0173
                                                                        DIS 0174
C   KCONV             COUNTER IN AZIMUTH CONVERGENCE TEST               DIS 0175
                                                                        DIS 0176
C   LAYRU(LU)         COMPUTATIONAL LAYER IN WHICH USER OUTPUT LEVEL    DIS 0177
C                     UTAU(LU) IS LOCATED                               DIS 0178
                                                                        DIS 0179
C   LL(IQ,LC)         CONSTANTS OF INTEGRATION L IN EQ. SC(1),          DIS 0180
C                     OBTAINED BY SOLVING SCALED VERSION OF EQ. SC(5)   DIS 0181
                                                                        DIS 0182
C   LYRCUT            TRUE, RADIATION IS ASSUMED ZERO BELOW LAYER       DIS 0183
C                     NCUT BECAUSE OF ALMOST COMPLETE ABSORPTION        DIS 0184
                                                                        DIS 0185
C   NAZ               NUMBER OF AZIMUTHAL COMPONENTS CONSIDERED         DIS 0186
                                                                        DIS 0187
C   NCUT              COMPUTATIONAL LAYER NUMBER IN WHICH ABSORPTION    DIS 0188
C                     OPTICAL DEPTH FIRST EXCEEDS ABSCUT                DIS 0189
                                                                        DIS 0190
C   OPRIM(LC)         SINGLE SCATTERING ALBEDO AFTER DELTA-M SCALING    DIS 0191
                                                                        DIS 0192
C   PASS1             TRUE ON FIRST ENTRY, FALSE THEREAFTER             DIS 0193
                                                                        DIS 0194
C   PKAG(0:LC)        INTEGRATED PLANCK FUNCTION FOR INTERNAL EMISSION  DIS 0195
                                                                        DIS 0196
C   PSIX(IQ)          SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)       DIS 0197
                                                                        DIS 0198
C   RMU(IU,0:IQ)      BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A  DIS 0199
C                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS    DIS 0200
C                     REFERS TO A USER ANGLE.  SECOND INDEX:            DIS 0201
C                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE UMU0;      DIS 0202
C                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.     DIS 0203
                                                                        DIS 0204
C   TAUC(0:LC)        CUMULATIVE OPTICAL DEPTH (UN-DELTA-M-SCALED)      DIS 0205
                                                                        DIS 0206
C   TAUCPR(0:LC)      CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED IF       DIS 0207
C                     DELTAM = TRUE, OTHERWISE EQUAL TO TAUC)           DIS 0208
                                                                        DIS 0209
C   TPLANK            INTENSITY EMITTED FROM TOP BOUNDARY               DIS 0210
                                                                        DIS 0211
C   UUM(IU,LU,MAZIM)  COMPONENTS OF THE INTENSITY (U-SUPER-M) WHEN      DIS 0212
C                     EXPANDED IN FOURIER COSINE SERIES IN AZIMUTH ANGLEDIS 0213
                                                                        DIS 0214
C   U0C(IQ,LU)        AZIMUTHALLY-AVERAGED INTENSITY                    DIS 0215
                                                                        DIS 0216
C   UTAUPR(LU)        OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M   DIS 0217
C                     COORDINATES;  EQUAL TO  UTAU(LU) IF NO DELTA-M    DIS 0218
                                                                        DIS 0219
C   WK()              SCRATCH ARRAY                                     DIS 0220
                                                                        DIS 0221
C   XR0(LC)           X-SUB-ZERO IN EXPANSION OF THERMAL SOURCE FUNC-   DIS 0222
C                     TION PRECEDING EQ. SS(14) (HAS NO MU-DEPENDENCE)  DIS 0223
                                                                        DIS 0224
C   XR1(LC)           X-SUB-ONE IN EXPANSION OF THERMAL SOURCE FUNC-    DIS 0225
C                     TION;  SEE  EQS. SS(14-16)                        DIS 0226
                                                                        DIS 0227
C   YLM0(L)           NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL         DIS 0228
C                     OF SUBSCRIPT L AT THE BEAM ANGLE (NOT SAVED       DIS 0229
C                     AS FUNCTION OF SUPERSCIPT M)                      DIS 0230
                                                                        DIS 0231
C   YLMC(L,IQ)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL         DIS 0232
C                     OF SUBSCRIPT L AT THE COMPUTATIONAL ANGLES        DIS 0233
C                     (NOT SAVED AS FUNCTION OF SUPERSCIPT M)           DIS 0234
                                                                        DIS 0235
C   YLMU(L,IU)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL         DIS 0236
C                     OF SUBSCRIPT L AT THE USER ANGLES                 DIS 0237
C                     (NOT SAVED AS FUNCTION OF SUPERSCIPT M)           DIS 0238
                                                                        DIS 0239
C   Z()               SCRATCH ARRAY USED IN  SOLVE0,1  TO SOLVE A       DIS 0240
C                     LINEAR SYSTEM FOR THE CONSTANTS OF INTEGRATION    DIS 0241
                                                                        DIS 0242
C   Z0(IQ)            SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)         DIS 0243
                                                                        DIS 0244
C   Z0U(IU,LC)        Z-SUB-ZERO IN EQ. SS(16) INTERPOLATED TO USER     DIS 0245
C                     ANGLES FROM AN EQUATION DERIVED FROM SS(16)       DIS 0246
                                                                        DIS 0247
C   Z1(IQ)            SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)         DIS 0248
                                                                        DIS 0249
C   Z1U(IU,LC)        Z-SUB-ONE IN EQ. SS(16) INTERPOLATED TO USER      DIS 0250
C                     ANGLES FROM AN EQUATION DERIVED FROM SS(16)       DIS 0251
                                                                        DIS 0252
C   ZBEAM(IU,LC)      PARTICULAR SOLUTION FOR BEAM SOURCE               DIS 0253
                                                                        DIS 0254
C   ZJ(IQ)            RIGHT-HAND SIDE VECTOR  X-SUB-ZERO IN             DIS 0255
C                     EQ. SS(19), ALSO THE SOLUTION VECTOR              DIS 0256
C                     Z-SUB-ZERO AFTER SOLVING THAT SYSTEM              DIS 0257
                                                                        DIS 0258
C   ZZ(IQ,LC)         PERMANENT STORAGE FOR THE BEAM SOURCE VECTORS ZJ  DIS 0259
                                                                        DIS 0260
C   ZPLK0(IQ,LC)      PERMANENT STORAGE FOR THE THERMAL SOURCE          DIS 0261
C                     VECTORS  Z0  OBTAINED BY SOLVING  EQ. SS(16)      DIS 0262
                                                                        DIS 0263
C   ZPLK1(IQ,LC)      PERMANENT STORAGE FOR THE THERMAL SOURCE          DIS 0264
C                     VECTORS  Z1  OBTAINED BY SOLVING  EQ. SS(16)      DIS 0265
                                                                        DIS 0266
                                                                        DIS 0267
      LOGICAL LYRCUT, PASS1                                             DIS 0268
      INTEGER IPVT( NNLYRI ), LAYRU( MXULV )                            DIS 0269
      REAL*8    AMB( MI,MI ), APB( MI,MI ),                             DIS 0270
     $        ARRAY( MXCMU,MXCMU ), B( NNLYRI ), BDR( MI,0:MI ),        DIS 0271
     $        BEM( MI ), CBAND( MI9M2,NNLYRI ), CC( MXCMU,MXCMU ),      DIS 0272
     $        CMU( MXCMU ), CWT( MXCMU ), EMU( MXUMU ), EVAL( MI ),     DIS 0273
     $        EVECC( MXCMU, MXCMU ), EXPBEA( 0:MXCLY ), FLYR( MXCLY ),  DIS 0274
     $        FLDN( MXULV ), FLDIR( MXULV ), GL( 0:MXCMU,MXCLY ),       DIS 0275
     $        GC( MXCMU,MXCMU,MXCLY ), GU( MXUMU,MXCMU,MXCLY ),         DIS 0276
     $        HLPR( 0:MXCMU ), KK( MXCMU,MXCLY ), LL( MXCMU,MXCLY ),    DIS 0277
     $        OPRIM( MXCLY ), PHIRAD( MXPHI ), PKAG( 0:MXCLY ),         DIS 0278
     $        PSIX( MXCMU ), RMU( MXUMU,0:MI ), TAUC( 0:MXCLY ),        DIS 0279
     $        TAUCPR( 0:MXCLY ), U0C( MXCMU,MXULV ), UTAUPR( MXULV ),   DIS 0280
     $        UUM( MXUMU,MXULV,0:MXCMU ), WK( MXCMU ), XR0( MXCLY ),    DIS 0281
     $        XR1( MXCLY ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU,MXCMU ),     DIS 0282
     $        YLMU( 0:MXCMU,MXUMU ), Z( NNLYRI ), Z0( MXCMU ),          DIS 0283
     $        Z0U( MXUMU,MXCLY ), Z1( MXCMU ),                          DIS 0284
     $        Z1U( MXUMU,MXCLY ), ZJ( MXCMU ), ZZ( MXCMU,MXCLY ),       DIS 0285
     $        ZPLK0( MXCMU,MXCLY ), ZPLK1( MXCMU,MXCLY ),               DIS 0286
     $        ZBEAM( MXUMU,MXCLY )                                      DIS 0287
C---- ------------------------------------------------------------------DIS 0288
C==== >       LOWER CASE VARIABLES ADDED                                DIS 0289
      REAL*8  CBANDT( MI9M2,NNLYRI ), CBANDS( MI9M2,NNLYRI ),           DIS 0290
     $        DUMMY0( MXCMU,MXCLY ), DUMMY1( MXCMU,MXCLY ),             DIS 0291
     $        DUMBEM( MI ), USAVE( MXUMU ), Z0UMS( MXUMU,MXCLY ),       DIS 0292
     $        Z1UMS( MXUMU,MXCLY ), BEAMMS( MXUMU,MXCLY )               DIS 0293
      LOGICAL MSFLAG                                                    DIS 0294
C-----------------------------------------------------------------------DIS 0295
                                                                        DIS 0296
      DOUBLE PRECISION   D1MACH                                         DIS 0297
      DOUBLE PRECISION   AAD( MI,MI ), EVALD( MI ) , EVECCD( MI,MI ),   DIS 0298
     $                   WKD( MXCMU )                                   DIS 0299
                                                                        DIS 0300
      SAVE  PASS1, PI, EPSIL, RPD                                       DIS 0301
      DATA  PASS1 / .TRUE. /                                            DIS 0302
                                                                        DIS 0303
           PASS1=.FALSE.                                                DIS 0304
         PI = 2. * DASIN(1.0D0)                                         DIS 0305
         EPSIL = 10.*D1MACH(3)                                          DIS 0306
         RPD = PI / 180.0                                               DIS 0307
      IF ( PASS1 )  THEN                                                DIS 0308
C                                ** INSERT INPUT VALUES FOR SELF-TEST   DIS 0309
C                                   NOTE: SELF-TEST MUST NOT USE IBCND=1DIS 0310
                                                                        DIS 0311
         CALL  SLFTST( ACCUR, ALBEDO, BTEMP, DELTAM, DTAUC( 1 ), FBEAM, DIS 0312
     $                 FISOT, IBCND, LAMBER, NLYR, PLANK, NPHI,         DIS 0313
     $                 NUMU, NSTR, NTAU, ONLYFL, PHI( 1 ), PHI0, PMOM,  DIS 0314
     $                 PRNT, SSALB( 1 ), TEMIS, TEMPER, TTEMP, UMU( 1 ),DIS 0315
     $                 USRANG, USRTAU, UTAU( 1 ), UMU0, WVNMHI, WVNMLO, DIS 0316
     $                 .FALSE., DUM, DUM, DUM, DUM )                    DIS 0317
      END IF                                                            DIS 0318
                                                                        DIS 0319
   1  CONTINUE                                                          DIS 0320
                                                                        DIS 0321
C-----------------------------------------------------------------------DIS 0322
C=====>      LOWER CASE VARIABLES ADDED                                 DIS 0323
C            **  INITIALIZE FOR MULTIPLE SCATTERING SOURCE FUNCTIONS USEDIS 0324
       IF ( (.NOT.PASS1).AND. MSFLAG)  THEN                             DIS 0325
          MSFLAG=.TRUE.                                                 DIS 0326
          DO 2 IU = 1, NUMU                                             DIS 0327
2         USAVE(IU) = UMU(IU)                                           DIS 0328
          CALL  ZEROIT( BEAMMS,MXUMU*MXCLY )                            DIS 0329
          CALL  ZEROIT( Z0UMS,MXUMU*MXCLY )                             DIS 0330
          CALL  ZEROIT( Z1UMS,MXUMU*MXCLY )                             DIS 0331
       END IF                                                           DIS 0332
C-----------------------------------------------------------------------DIS 0333
                                                                        DIS 0334
      IF ( PRNT(1) )  WRITE( *,1010 )  HEADER                           DIS 0335
                                                                        DIS 0336
C                         ** ZERO SOME ARRAYS (NOT STRICTLY NECESSARY,  DIS 0337
C                            BUT OTHERWISE UNUSED PARTS OF ARRAYS       DIS 0338
C                            COLLECT GARBAGE)                           DIS 0339
      DO 10 I = 1, NNLYRI                                               DIS 0340
         IPVT(I) = 0                                                    DIS 0341
10    CONTINUE                                                          DIS 0342
      CALL ZEROAL( MXCLY, XR0, XR1, TAUC(1),                            DIS 0343
     $             MXCMU, CMU, CWT, PSIX, EVAL, WK, Z0, Z1, ZJ,         DIS 0344
     $             MXCMU+1, HLPR, YLM0,                                 DIS 0345
     $             MXCMU**2, ARRAY, CC, EVECC, YLMC,                    DIS 0346
     $             (MXCMU+1)*MXUMU, YLMU,                               DIS 0347
     $             MI**2, AMB, APB,                                     DIS 0348
     $             MXCMU*MXCLY, KK, LL, ZZ, ZPLK0, ZPLK1,               DIS 0349
     $             MXUMU*MXCLY, Z0U, Z1U, ZBEAM,                        DIS 0350
     $             MXCMU**2*MXCLY, GC,                                  DIS 0351
     $             MXUMU*MXCMU*MXCLY, GU,                               DIS 0352
     $             NNLYRI, Z )                                          DIS 0353
                                                                        DIS 0354
C                                  ** CALCULATE CUMULATIVE OPTICAL DEPTHDIS 0355
C                                     AND DITHER SINGLE-SCATTER ALBEDO  DIS 0356
C                                     TO IMPROVE NUMERICAL BEHAVIOR OF  DIS 0357
C                                     EIGENVALUE/VECTOR COMPUTATION     DIS 0358
      TAUC( 0 ) = 0.                                                    DIS 0359
      DO 20  LC = 1, NLYR                                               DIS 0360
         IF( SSALB(LC).EQ.1.0 )  SSALB(LC) = 1.0 - EPSIL                DIS 0361
         TAUC(LC) = TAUC(LC-1) + DTAUC(LC)                              DIS 0362
20    CONTINUE                                                          DIS 0363
C                                ** CHECK INPUT DIMENSIONS AND VARIABLESDIS 0364
                                                                        DIS 0365
      CALL  CHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,           DIS 0366
     $              WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,           DIS 0367
     $              NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,           DIS 0368
     $              PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP,             DIS 0369
     $              TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, MAXCLY,         DIS 0370
     $              MAXULV, MAXUMU, MAXCMU, MAXPHI, MXCLY,              DIS 0371
     $              MXULV,  MXUMU,  MXCMU,  MXPHI, TAUC )               DIS 0372
                                                                        DIS 0373
C                                 ** PERFORM VARIOUS SETUP OPERATIONS   DIS 0374
                                                                        DIS 0375
      CALL  SETDIS( CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR,       DIS 0376
     $              GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT,         DIS 0377
     $              MAXUMU, MAXCMU, MXCMU, NCUT, NLYR, NTAU, NN,        DIS 0378
     $              NSTR, PLANK, NUMU, ONLYFL, OPRIM, PMOM, SSALB,      DIS 0379
     $              TAUC, TAUCPR, UTAU, UTAUPR, UMU, UMU0, USRTAU,      DIS 0380
     $              USRANG )                                            DIS 0381
                                                                        DIS 0382
C                                             ** PRINT INPUT INFORMATIONDIS 0383
      IF ( PRNT(1) )                                                    DIS 0384
     $     CALL PRTINP( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,       DIS 0385
     $                  WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, NPHI, PHI, DIS 0386
     $                  IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER,        DIS 0387
     $                  ALBEDO, HL, BTEMP, TTEMP, TEMIS, DELTAM, PLANK, DIS 0388
     $                  ONLYFL, ACCUR, FLYR, LYRCUT, OPRIM, TAUC,       DIS 0389
     $                  TAUCPR, MAXCMU, PRNT(7) )                       DIS 0390
                                                                        DIS 0391
                                                                        DIS 0392
      IF ( IBCND.EQ.1 )  THEN                                           DIS 0393
C                              ** HANDLE SPECIAL CASE FOR GETTING ALBEDODIS 0394
C                                 AND TRANSMISSIVITY OF MEDIUM FOR MANY DIS 0395
C                                 BEAM ANGLES AT ONCE                   DIS 0396
                                                                        DIS 0397
         CALL  ALBTRN( ALBEDO, AMB, APB, ARRAY, B, BDR, CBAND, CC,      DIS 0398
     $                 CMU, CWT, EVAL, EVECC, GL, GC, GU, IPVT,         DIS 0399
     $                 KK, LL, NLYR, NN, NSTR, NUMU, PRNT, TAUCPR,      DIS 0400
     $                 UMU, U0U, WK, YLMC, YLMU, Z, AAD, EVALD,         DIS 0401
     $                 EVECCD, WKD, MI, MI9M2, MAXULV, MAXUMU,          DIS 0402
     $                 MXCMU, MXUMU, NNLYRI, ALBMED, TRNMED )           DIS 0403
         RETURN                                                         DIS 0404
      ENDIF                                                             DIS 0405
C                                   ** CALCULATE PLANCK FUNCTIONS       DIS 0406
      IF ( .NOT.PLANK )  THEN                                           DIS 0407
         BPLANK = 0.0                                                   DIS 0408
         TPLANK = 0.0                                                   DIS 0409
         CALL  ZEROIT( PKAG, MXCLY+1 )                                  DIS 0410
      ELSE                                                              DIS 0411
C                      **  USE DIFFERENT PLANCK FUNCTIONS PLKAVG OR BBFNDIS 0412
C=====>    LOWER CASE VARIABLES ADDED                                   DIS 0413
           IF ( PASS1 .OR. (.NOT.MSFLAG) )  THEN                        DIS 0414
            TPLANK = TEMIS * PLKAVG( WVNMLO, WVNMHI, TTEMP )            DIS 0415
            BPLANK =         PLKAVG( WVNMLO, WVNMHI, BTEMP )            DIS 0416
            DO 40  LEV = 0, NLYR                                        DIS 0417
40            PKAG( LEV ) = PLKAVG( WVNMLO, WVNMHI, TEMPER(LEV) )       DIS 0418
           ELSE                                                         DIS 0419
              TPLANK = TEMIS * DISBBF( TTEMP, WN0 )                     DIS 0420
              BPLANK =         DISBBF( BTEMP, WN0 )                     DIS 0421
              DO 41  LEV = 0, NLYR                                      DIS 0422
            PKAG( LEV ) = DISBBF( TEMPER(LEV), WN0 )                    DIS 0423
41           CONTINUE                                                   DIS 0424
           END IF                                                       DIS 0425
        END IF                                                          DIS 0426
C ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======DIS 0427
C ========  (EQ STWJ 5)                                                 DIS 0428
                                                                        DIS 0429
      KCONV = 0                                                         DIS 0430
      NAZ = NSTR-1                                                      DIS 0431
C                                            ** AZIMUTH-INDEPENDENT CASEDIS 0432
                                                                        DIS 0433
      IF ( FBEAM.EQ.0.0 .OR. (1.-UMU0).LT.1.E-5 .OR. ONLYFL .OR.        DIS 0434
     $     (NUMU.EQ.1.AND.(1.-UMU(1)).LT.1.E-5 ) )                      DIS 0435
     $   NAZ = 0                                                        DIS 0436
                                                                        DIS 0437
      CALL  ZEROIT( UU, MAXUMU*MAXULV*MAXPHI )                          DIS 0438
      DO  200  MAZIM = 0, NAZ                                           DIS 0439
                                                                        DIS 0440
      IF ( MAZIM.EQ.0 )  DELM0 = 1.0                                    DIS 0441
      IF ( MAZIM.GT.0 )  DELM0 = 0.0                                    DIS 0442
C                                  ** GET NORMALIZED ASSOCIATED LEGENDREDIS 0443
C                                     POLYNOMIALS FOR INCIDENT BEAM     DIS 0444
C                                     ANGLE COSINE                      DIS 0445
      IF ( FBEAM.GT.0.0 )                                               DIS 0446
     $     CALL  LEPOLY( 1, MAZIM, MXCMU, NSTR-1, -UMU0, YLM0 )         DIS 0447
                                                                        DIS 0448
C                                  ** GET NORMALIZED ASSOCIATED LEGENDREDIS 0449
C                                         POLYNOMIALS FOR COMPUTATIONAL DIS 0450
C                                         AND USER POLAR ANGLE COSINES  DIS 0451
      IF ( (.NOT.ONLYFL .AND. USRANG) )                                 DIS 0452
     $     CALL  LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, YLMU )        DIS 0453
         CALL  LEPOLY( NN,   MAZIM, MXCMU, NSTR-1, CMU, YLMC )          DIS 0454
                                                                        DIS 0455
C=====> LOWER CASE VARIABLES ADDED                                      DIS 0456
       IF ( (.NOT.PASS1) .AND. MSFLAG )  THEN                           DIS 0457
           DO 42 IU = 1, NUMU                                           DIS 0458
42         UMU(IU) = USAVE(IU)                                          DIS 0459
           CALL  LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, YLMU )        DIS 0460
        END IF                                                          DIS 0461
C                       ** EVALUATE NORMALIZED ASSOCIATED LEGENDRE      DIS 0462
C                          POLYNOMIALS WITH NEGATIVE  CMU  FROM THOSE   DIS 0463
C                          WITH POSITIVE  CMU; DAVE/ARMSTRONG EQ. (15)  DIS 0464
      SGN  = - 1.0                                                      DIS 0465
      DO  50  L = MAZIM, NSTR-1                                         DIS 0466
         SGN = - SGN                                                    DIS 0467
         DO  50  IQ = NN+1, NSTR                                        DIS 0468
            YLMC( L,IQ ) = SGN * YLMC( L,IQ-NN )                        DIS 0469
 50   CONTINUE                                                          DIS 0470
C                                 ** SPECIFY USERS BOTTOM REFLECTIVITY  DIS 0471
C                                    AND EMISSIVITY PROPERTIES          DIS 0472
      IF ( .NOT.LYRCUT )                                                DIS 0473
     $   CALL  SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER,              DIS 0474
     $                 MI, MAZIM, MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL, DIS 0475
     $                 UMU, USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM,    DIS 0476
     $                 RMU )                                            DIS 0477
C ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============DIS 0478
                                                                        DIS 0479
      DO 100  LC = 1, NCUT                                              DIS 0480
                                                                        DIS 0481
C                        ** SOLVE EIGENFUNCTION PROBLEM IN EQ. STWJ(8B);DIS 0482
C                           RETURN EIGENVALUES AND EIGENVECTORS         DIS 0483
                                                                        DIS 0484
         CALL  SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL(0,LC), MI, MAZIM,  DIS 0485
     $                 MXCMU, NN, NSTR, WK, YLMC, CC, EVECC, EVAL,      DIS 0486
     $                 KK(1,LC), GC(1,1,LC), AAD, WKD, EVECCD, EVALD )  DIS 0487
C                                  ** CALCULATE PARTICULAR SOLUTIONS OF DIS 0488
C                                     EQ.SS(18) FOR INCIDENT BEAM SOURCEDIS 0489
         IF ( FBEAM.GT.0.0 )                                            DIS 0490
     $        CALL  UPBEAM( ARRAY, CC, CMU, DELM0, FBEAM, GL(0,LC),     DIS 0491
     $                      IPVT, MAZIM, MXCMU, NN, NSTR, PI, UMU0, WK, DIS 0492
     $                      YLM0, YLMC, ZJ, ZZ(1,LC) )                  DIS 0493
         IF ( PLANK .AND. MAZIM.EQ.0 ) THEN                             DIS 0494
                                                                        DIS 0495
C                              ** CALCULATE PARTICULAR SOLUTIONS OF     DIS 0496
C                                 EQ. SS(15) FOR THERMAL EMISSION SOURCEDIS 0497
                                                                        DIS 0498
            DELTAT = TAUCPR(LC) - TAUCPR(LC-1)                          DIS 0499
            XR1( LC ) = 0.0                                             DIS 0500
                                                                        DIS 0501
            IF ( DELTAT.GT.EPSIL ) XR1( LC ) = ( PKAG(LC) - PKAG(LC-1) )DIS 0502
     $                                       / DELTAT                   DIS 0503
            XR0( LC ) = PKAG(LC-1) - XR1(LC) * TAUCPR(LC-1)             DIS 0504
            CALL UPISOT( ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR,         DIS 0505
     $                     OPRIM(LC), WK, XR0(LC), XR1(LC), Z0, Z1,     DIS 0506
     $                     ZPLK0(1,LC), ZPLK1(1,LC) )                   DIS 0507
         END IF                                                         DIS 0508
C=====>    LOWER CASE VARIABLES ADDED                                   DIS 0509
         IF ( (.NOT.ONLYFL .AND. USRANG) .OR. MSFLAG ) THEN             DIS 0510
                                                                        DIS 0511
C                                            ** INTERPOLATE EIGENVECTORSDIS 0512
C                                               TO USER ANGLES          DIS 0513
                                                                        DIS 0514
            CALL  TERPEV( CWT, EVECC, GL(0,LC), GU(1,1,LC), MAZIM,      DIS 0515
     $                    MXCMU, MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU) DIS 0516
                                                                        DIS 0517
C                                            ** INTERPOLATE SOURCE TERMSDIS 0518
C                                               TO USER ANGLES          DIS 0519
                                                                        DIS 0520
            CALL  TERPSO( CWT, DELM0, FBEAM, GL(0,LC), MAZIM,           DIS 0521
     $                    MXCMU, PLANK, NUMU, NSTR, OPRIM(LC),          DIS 0522
     $                    PI, YLM0, YLMC, YLMU, PSIX, XR0(LC), XR1(LC), DIS 0523
     $                    Z0, ZJ, ZBEAM(1,LC), Z0U(1,LC), Z1U(1,LC),    DIS 0524
     $                    Z0UMS(1,LC), Z1UMS(1,LC), BEAMMS(1,LC) )      DIS 0525
C=====>    LOWER CASE VARIABLES ADDED                                   DIS 0526
         END IF                                                         DIS 0527
100   CONTINUE                                                          DIS 0528
                                                                        DIS 0529
C ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============DIS 0530
                                                                        DIS 0531
C                      ** SET COEFFICIENT MATRIX OF EQUATIONS COMBINING DIS 0532
C                         BOUNDARY AND LAYER INTERFACE CONDITIONS       DIS 0533
      CALL  SETMTX( BDR, CBAND, CMU, CWT, DELM0, GC, KK, LAMBER,        DIS 0534
     $              LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,       DIS 0535
     $              NN, NSTR, TAUCPR, WK )                              DIS 0536
C=====>        LOWER CASE VARIABLES ADDED                               DIS 0537
        DO 110 JK = 1, MI9M2                                            DIS 0538
           DO 110 IK = 1, NNLYRI                                        DIS 0539
           CBANDS(JK,IK)=CBAND(JK,IK)                                   DIS 0540
           CBANDT(JK,IK)=CBAND(JK,IK)                                   DIS 0541
110     CONTINUE                                                        DIS 0542
                                                                        DIS 0543
C                      ** SOLVE FOR CONSTANTS OF INTEGRATION IN HOMO-   DIS 0544
C                         GENEOUS SOLUTION (GENERAL BOUNDARY CONDITIONS)DIS 0545
C               **  COMBINED SOLAR AND THERMAL SOURCES FOR SELFTEST CASEDIS 0546
      CALL  SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,       DIS 0547
     $              FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,             DIS 0548
     $              MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR,      DIS 0549
     $              NNLYRI, PI, TPLANK, TAUCPR, UMU0, Z, ZZ,            DIS 0550
     $              ZPLK0, ZPLK1 )                                      DIS 0551
                                                                        DIS 0552
C                                  ** COMPUTE UPWARD AND DOWNWARD FLUXESDIS 0553
      IF ( MAZIM.EQ.0 )                                                 DIS 0554
     $     CALL FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,     DIS 0555
     $                  MXCMU, MXULV, NCUT, NN, NSTR, NTAU, PI,         DIS 0556
     $                  PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,        DIS 0557
     $                  XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,         DIS 0558
     $                  FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C, MAXULV ) DIS 0559
C-----------------------------------------------------------------------DIS 0560
C            **  COMPUTE MULTIPLE SCATTERING SOURCE FUNCTIONS SEPARATELYDIS 0561
C=====>      LOWER CASE VARIABLES ADDED                                 DIS 0562
       IF ( (.NOT.PASS1) .AND. MSFLAG )  THEN                           DIS 0563
                                                                        DIS 0564
           CALL ZEROIT( DUMMY0, MXCMU*MXCLY )                           DIS 0565
C                            **  FIRST CALL FOR SOLAR MS SOURCE FUNCTIONDIS 0566
                                                                        DIS 0567
           IF ( FBEAM.GT.0. )  THEN                                     DIS 0568
           CALL ZEROIT( DUMMY1, MXCMU*MXCLY )                           DIS 0569
           CALL ZEROIT( DUMBEM, MI )                                    DIS 0570
           CALL  SOLVE0( B,BDR,DUMBEM, 0.D0, CBANDS, CMU, CWT, EXPBEA,  DIS 0571
     $                   FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,        DIS 0572
     $                   MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, DIS 0573
     $                   NNLYRI, PI, 0.D0, TAUCPR, UMU0, Z, ZZ,         DIS 0574
     $                   DUMMY0, DUMMY1 )                               DIS 0575
           CALL  MSSOLR( CMU, CWT, FBEAM, GC, GU, KK, LAYRU, LL,        DIS 0576
     $                   LYRCUT, MAXUMU, MXCMU, MXUMU, NCUT, NN,        DIS 0577
     $                   NSTR, NTAU, NUMU, PI, TAUCPR, UMU0, UTAU,      DIS 0578
     $                   UTAUPR, BEAMMS, ZZ, FDNSRT, S0CMS )            DIS 0579
           END IF                                                       DIS 0580
                                                                        DIS 0581
C                         **  SECOND CALL FOR THERMAL MS SOURCE FUNCTIONDIS 0582
                                                                        DIS 0583
           IF ( PLANK )  THEN                                           DIS 0584
           CALL ZEROIT(DUMMY0,MXCMU*MXCLY)                              DIS 0585
           DO 111 JK = 0, MXCLY                                         DIS 0586
111        EXPBEA(JK) = 0.0                                             DIS 0587
           CALL  SOLVE0( B, BDR, BEM, BPLANK, CBANDT, CMU, CWT, EXPBEA, DIS 0588
     $                   FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,        DIS 0589
     $                   MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, DIS 0590
     $                   NNLYRI, PI, TPLANK, TAUCPR, UMU0, Z, DUMMY0,   DIS 0591
     $                   ZPLK0, ZPLK1 )                                 DIS 0592
                                                                        DIS 0593
           CALL  MSTHML( CMU, CWT, GC, GU, KK, LAYRU, LL, LYRCUT,       DIS 0594
     $                   MAXUMU, MXCMU, MXUMU, NCUT, NN, NSTR, NTAU,    DIS 0595
     $                   NUMU, OPRIM, PI, TAUCPR, UMU0, UTAUPR,         DIS 0596
     $                   XR0, XR1, Z0UMS, Z1UMS, ZPLK0, ZPLK1,          DIS 0597
     $                   FDNTRT, T0CMS )                                DIS 0598
                                                                        DIS 0599
           END IF                                                       DIS 0600
                                                                        DIS 0601
        END IF                                                          DIS 0602
C-----------------------------------------------------------------------DIS 0603
                                                                        DIS 0604
           IF ( ONLYFL )  THEN                                          DIS 0605
         IF( MAXUMU.GE.NSTR )  THEN                                     DIS 0606
C                                         ** SAVE AZIM-AVGD INTENSITIES DIS 0607
C                                            AT QUADRATURE ANGLES       DIS 0608
            DO 120 LU = 1, NTAU                                         DIS 0609
               DO 120 IQ = 1, NSTR                                      DIS 0610
                  U0U( IQ,LU ) = U0C( IQ,LU )                           DIS 0611
120         CONTINUE                                                    DIS 0612
         ELSE                                                           DIS 0613
               CALL  ZEROIT( U0U, MAXUMU*MAXULV )                       DIS 0614
         ENDIF                                                          DIS 0615
         GO TO 210                                                      DIS 0616
      ENDIF                                                             DIS 0617
                                                                        DIS 0618
      IF ( USRANG ) THEN                                                DIS 0619
C                          ** COMPUTE AZIMUTHAL INTENSITY               DIS 0620
C                                        COMPONENTS AT USER ANGLES      DIS 0621
         CALL  USRINT( BPLANK, CMU, CWT, DELM0, EMU, EXPBEA,            DIS 0622
     $                 FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL,     DIS 0623
     $                 LYRCUT, MAZIM, MXCMU, MXULV, MXUMU, NCUT,        DIS 0624
     $                 NLYR, NN, NSTR, PLANK, NUMU, NTAU, PI, RMU,      DIS 0625
     $                 TAUCPR, TPLANK, UMU, UMU0, UTAUPR, WK,           DIS 0626
     $                 ZBEAM, Z0U, Z1U, ZZ, ZPLK0, ZPLK1, UUM )         DIS 0627
                                                                        DIS 0628
      ELSE                                                              DIS 0629
C                                     ** COMPUTE AZIMUTHAL INTENSITY    DIS 0630
C                                        COMPONENTS AT QUADRATURE ANGLESDIS 0631
                                                                        DIS 0632
         CALL  CMPINT( FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM,         DIS 0633
     $                 MXCMU, MXULV, MXUMU, NCUT, NN, NSTR,             DIS 0634
     $                 PLANK, NTAU, TAUCPR, UMU0, UTAUPR,               DIS 0635
     $                 ZZ, ZPLK0, ZPLK1, UUM )                          DIS 0636
      END IF                                                            DIS 0637
                                                                        DIS 0638
      IF( MAZIM.EQ.0 ) THEN                                             DIS 0639
                                                                        DIS 0640
         DO  140  J = 1, NPHI                                           DIS 0641
            PHIRAD( J ) = RPD * ( PHI(J) - PHI0 )                       DIS 0642
 140     CONTINUE                                                       DIS 0643
C                               ** SAVE AZIMUTHALLY AVERAGED INTENSITIESDIS 0644
         DO 160  LU = 1, NTAU                                           DIS 0645
            DO 160  IU = 1, NUMU                                        DIS 0646
               U0U( IU,LU ) = UUM( IU,LU,0 )                            DIS 0647
 160     CONTINUE                                                       DIS 0648
C                              ** PRINT AZIMUTHALLY AVERAGED INTENSITIESDIS 0649
C                                 AT USER ANGLES                        DIS 0650
         IF ( PRNT(4) )                                                 DIS 0651
     $        CALL PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU, U0U )         DIS 0652
                                                                        DIS 0653
      END IF                                                            DIS 0654
C                                ** INCREMENT INTENSITY BY CURRENT      DIS 0655
C                                   AZIMUTHAL COMPONENT (FOURIER        DIS 0656
C                                   COSINE SERIES);  EQ SD(2)           DIS 0657
      AZERR = 0.0                                                       DIS 0658
      DO 180  J = 1, NPHI                                               DIS 0659
         COSPHI = COS( MAZIM * PHIRAD(J) )                              DIS 0660
         DO 180  LU = 1, NTAU                                           DIS 0661
            DO 180  IU = 1, NUMU                                        DIS 0662
               AZTERM = UUM( IU,LU,MAZIM ) * COSPHI                     DIS 0663
               UU( IU,LU,J ) = UU( IU,LU,J ) + AZTERM                   DIS 0664
               AZERR = DMAX1( RATIO( DABS(AZTERM), DABS(UU(IU,LU,J)) ), DIS 0665
     $                        AZERR )                                   DIS 0666
180   CONTINUE                                                          DIS 0667
      IF ( AZERR.LE.ACCUR )  KCONV = KCONV + 1                          DIS 0668
      IF ( KCONV.GE.2 )      GOTO 210                                   DIS 0669
                                                                        DIS 0670
200   CONTINUE                                                          DIS 0671
                                                                        DIS 0672
C ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============DIS 0673
                                                                        DIS 0674
                                                                        DIS 0675
C                                                 ** PRINT INTENSITIES  DIS 0676
                                                                        DIS 0677
 210  IF ( PRNT(5) .AND. .NOT.ONLYFL )                                  DIS 0678
     $     CALL  PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI,          DIS 0679
     $                   MAXULV, MAXUMU )                               DIS 0680
                                                                        DIS 0681
                                                                        DIS 0682
      IF ( PASS1 )  THEN                                                DIS 0683
C                                    ** COMPARE TEST CASE RESULTS WITH  DIS 0684
C                                       CORRECT ANSWERS AND ABORT IF BADDIS 0685
                                                                        DIS 0686
         CALL SLFTST( ACCUR, ALBEDO, BTEMP, DELTAM, DTAUC( 1 ), FBEAM,  DIS 0687
     $                FISOT, IBCND, LAMBER, NLYR, PLANK, NPHI,          DIS 0688
     $                NUMU, NSTR, NTAU, ONLYFL, PHI( 1 ), PHI0, PMOM,   DIS 0689
     $                PRNT, SSALB( 1 ), TEMIS, TEMPER, TTEMP, UMU( 1 ), DIS 0690
     $                USRANG, USRTAU, UTAU( 1 ), UMU0, WVNMHI, WVNMLO,  DIS 0691
     $                .TRUE., FLUP( 1 ), RFLDIR( 1 ), RFLDN( 1 ),       DIS 0692
     $                UU( 1,1,1 ) )                                     DIS 0693
         PASS1 = .FALSE.                                                DIS 0694
         GO TO 1                                                        DIS 0695
      END IF                                                            DIS 0696
                                                                        DIS 0697
                                                                        DIS 0698
                SFDTRT=FDNTRT                                           DIS 0699
                SFDSRT=FDNSRT                                           DIS 0700
                                                                        DIS 0701
                                                                        DIS 0702
                                                                        DIS 0703
      RETURN                                                            DIS 0704
                                                                        DIS 0705
1010  FORMAT ( ////, 1X, 120('*'), /, 25X,                              DIS 0706
     $  'DISCRETE ORDINATES RADIATIVE TRANSFER PROGRAM, VERSION 1.0',   DIS 0707
     $  /, 1X, A, /, 1X, 120('*') )                                     DIS 0708
      END                                                               DIS 0709
