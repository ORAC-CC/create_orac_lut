      SUBROUTINE FLXADD(ML,IEMSCT,SALB)                                 FLX 0001
C                                                                       FLX 0002
C     CALCULATES UPWARD (UPF) AND DOWNWARD (DNF) FLUX PROFILES USING    FLX 0003
C     ADDING METHOD.  LAYER SOURCE FUNCTIONS CALCULATED FROM STREAM     FLX 0004
C     APPROXIMATION [FOR THERMAL SEE ISAACS ET AL, AFGL-TR-86-0073      FLX 0005
C     AND FOR SOLAR SEE MEADOR AND WEAVER, J ATMOS SCI V37, P630].      FLX 0006
C     ***  A.E.R. 1986;  S.S.I. 1995.  ***                              FLX 0007
C                                                                       FLX 0008
C     DECLARE INPUTS                                                    FLX 0009
C       ML       NUMBER OF ATMOSPHERIC LAYER BOUNDARIES                 FLX 0010
C       IEMSCT   FLAG = 1 FOR THERMAL SCATTER ONLY                      FLX 0011
C                     = 2 FOR THERMAL AND SOLAR SCATTER                 FLX 0012
C       SALB     SURFACE SCATTERING ALBEDO                              FLX 0013
      INTEGER ML,IEMSCT                                                 FLX 0014
      REAL SALB                                                         FLX 0015
C                                                                       FLX 0016
C     LIST PARAMETERS                                                   FLX 0017
      INCLUDE 'PARAM.LST'                                               FLX 0018
      INTEGER IBND                                                      FLX 0019
      REAL AA,BB,CC,QA,CPS                                              FLX 0020
      COMMON/AABBCC/AA(11),BB(11),CC(11),IBND(11),QA(11),CPS(11)        FLX 0021
      INTEGER KPOINT                                                    FLX 0022
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     FLX 0023
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   FLX 0024
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   FLX 0025
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     FLX 0026
C                                                                       FLX 0027
C     COMMON /MSRD/                                                     FLX 0028
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            FLX 0029
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             FLX 0030
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               FLX 0031
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    FLX 0032
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  FLX 0033
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. FLX 0034
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               FLX 0035
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         FLX 0036
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          FLX 0037
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              FLX 0038
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  FLX 0039
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               FLX 0040
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           FLX 0041
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     FLX 0042
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   FLX 0043
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       FLX 0044
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     FLX 0045
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              FLX 0046
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 FLX 0047
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          FLX 0048
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             FLX 0049
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                FLX 0050
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      FLX 0051
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         FLX 0052
      REAL EDN,EUP,EUPC,TDF,RUPC,REF,EDNS,EUPS,EUPCS,TDFS,RUPCS,REFS    FLX 0053
      COMMON/FLUX/EDN(LAYDIM),EUP(LAYDIM),EUPC(LAYDIM),TDF(LAYDIM),     FLX 0054
     1  RUPC(LAYDIM),REF(LAYDIM),EDNS(LAYDIM),EUPS(LAYDIM),             FLX 0055
     2  EUPCS(LAYDIM),TDFS(LAYDIM),RUPCS(LAYDIM),REFS(LAYDIM)           FLX 0056
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               FLX 0057
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           FLX 0058
C                                                                       FLX 0059
C     INTERNAL VARIABLES:                                               FLX 0060
C       TAUM(K,N) MOLECULAR OPTICAL THICKNESS OF LAYER N FOR K          FLX 0061
C       EUP(N)    INTRINSIC UPWARD THERMAL FLUX OF LAYER N              FLX 0062
C       EDN(N)    INTRINSIC DOWNWARD THERMAL FLUX OF LAYER N            FLX 0063
C       TDF(N)    INTRINSIC THERMAL TRANSMISSION OF LAYER N             FLX 0064
C       REF(N)    INTRINSIC THERMAL REFLECTANCE OF LAYER N              FLX 0065
C       EUPS(N)   INTRINSIC UPWARD SOLAR FLUX OF LAYER N                FLX 0066
C       EDNS(N)   INTRINSIC DOWNWARD SOLAR FLUX OF LAYER N              FLX 0067
C       TDFS(N)   INTRINSIC SOLAR TRANSMISSION OF LAYER N               FLX 0068
C       REFS(N)   INTRINSIC SOLAR REFLECTANCE OF LAYER N                FLX 0069
C                                                                       FLX 0070
C     DECLARE LOCAL VARIABLES                                           FLX 0071
      INTEGER NLAYRS,IKP1,N,NM1,IK,MOL,K,IB                             FLX 0072
      REAL RDNCN,RUPCN,EUPCN,RDNCNS,RUPCNS,EUPCNS,TAU,A1ME0,A2C,        FLX 0073
     1  EDNCN,EDNCNS,AC,DENOM,RT3,COEF,A2C2,C,HALFE1,E0,EXPAN,          FLX 0074
     2  ONEME0,E1,COEF1,COEF2,ACM1,ONEPE0,ONEME1,E2,EX,SOURC,TWGP,SMDPJ FLX 0075
C                                                                       FLX 0076
C     DECLARE LOCAL ARRAYS                                              FLX 0077
      REAL TAUM(3,LAYDIM),DPJ(3,LAYDIM),GKWJ(3,11),DPWJ(3,11)           FLX 0078
C                                                                       FLX 0079
C     LIST DATA                                                         FLX 0080
      DATA RT3/1.7320508/                                               FLX 0081
C                                                                       FLX 0082
C     NUMBER OF LAYERS                                                  FLX 0083
      NLAYRS=ML-1                                                       FLX 0084
C                                                                       FLX 0085
C     11 MOLECULES BY JOSEPH H PIERLUISSI                               FLX 0086
C       DPWJ   PROBABILITY FOR EACH MOLECULE FIT DOUBLE EXPONENTIAL.    FLX 0087
C       GKWJ   BAND DEPENDENT SCALING OF DENSITIES TO GET K AMOUNT.     FLX 0088
C       CPS    THE STORED VALUES OF PIERLUISSI BAND MODEL COEFFICIENTS. FLX 0089
C       IBND   MAPPING FROM MOLECULE NUMBER TO BAND NUMBER.             FLX 0090
      DO 10 MOL=1,11                                                    FLX 0091
          IF(IBND(MOL).LE.0)GOTO10                                      FLX 0092
          GKWJ(1,MOL)=CC(MOL)*10.**CPS(MOL)                             FLX 0093
          GKWJ(2,MOL)=.09*GKWJ(1,MOL)                                   FLX 0094
          GKWJ(3,MOL)=.015*GKWJ(1,MOL)                                  FLX 0095
          DPWJ(1,MOL)=AA(MOL)                                           FLX 0096
          DPWJ(2,MOL)=BB(MOL)                                           FLX 0097
          DPWJ(3,MOL)=1.-AA(MOL)-BB(MOL)                                FLX 0098
   10 CONTINUE                                                          FLX 0099
C                                                                       FLX 0100
C     TAUM IS DEFINED AS THE SUM OF THE OPTICAL DEPTHS BY MOLECULE      FLX 0101
      IK=LAYTWO+ML                                                      FLX 0102
      DO 40 N=1,NLAYRS                                                  FLX 0103
          IK=IK-1                                                       FLX 0104
          TAUT(1,N)=TCONT(N)                                            FLX 0105
          UPF(1,N)=0.                                                   FLX 0106
          DNF(1,N)=0.                                                   FLX 0107
          UPFS(1,N)=0.                                                  FLX 0108
          DNFS(1,N)=0.                                                  FLX 0109
          SMDPJ=0.                                                      FLX 0110
          DO 30 K=1,3                                                   FLX 0111
              TAUM(K,N)=0.                                              FLX 0112
              TWGP=0.                                                   FLX 0113
              DO 20 MOL=1,11                                            FLX 0114
                  IB=IBND(MOL)                                          FLX 0115
                  IF(IB.LE.0)GOTO20                                     FLX 0116
                  TAU=WPATH(IK,IB)*GKWJ(K,MOL)                          FLX 0117
                  TAUM(K,N)=TAUM(K,N)+TAU                               FLX 0118
                  TWGP=TWGP+TAU*DPWJ(K,MOL)                             FLX 0119
   20         CONTINUE                                                  FLX 0120
C                                                                       FLX 0121
C             EFFECTIVE PROBABILITY BY LAYER DPJ IS BASED ON            FLX 0122
C             MOLECULAR PROBABILITY WEIGHTED BY OPTICAL DEPTH.          FLX 0123
              DPJ(K,N)=.33333333                                        FLX 0124
              IF(TAUM(K,N).NE.0.)DPJ(K,N)=TWGP/TAUM(K,N)                FLX 0125
              SMDPJ=SMDPJ+DPJ(K,N)                                      FLX 0126
   30     CONTINUE                                                      FLX 0127
          DPJ(1,N)=DPJ(1,N)/SMDPJ                                       FLX 0128
          DPJ(2,N)=DPJ(2,N)/SMDPJ                                       FLX 0129
          DPJ(3,N)=DPJ(3,N)/SMDPJ                                       FLX 0130
   40 CONTINUE                                                          FLX 0131
      UPF(1,ML)=0.                                                      FLX 0132
      DNF(1,ML)=0.                                                      FLX 0133
      UPFS(1,ML)=0.                                                     FLX 0134
      DNFS(1,ML)=0.                                                     FLX 0135
C                                                                       FLX 0136
C     PROBABILITY INTEGRATION LOOP                                      FLX 0137
      DO 70 K=1,3                                                       FLX 0138
C                                                                       FLX 0139
C         COMPOSITE DOWNWARD REFLECTION                                 FLX 0140
          RDNCN=0.                                                      FLX 0141
C                                                                       FLX 0142
C         DEFINE INITIAL UPWARD COMPOSITE SURFACE REFLECTANCE           FLX 0143
          RUPCN=SALB                                                    FLX 0144
C                                                                       FLX 0145
C         SURFACE EMISSION                                              FLX 0146
          EUPCN=(1.-RUPCN)*BBGRND                                       FLX 0147
          EUP(ML)=EUPCN                                                 FLX 0148
          EUPC(ML)=EUPCN                                                FLX 0149
          RUPC(ML)=0.                                                   FLX 0150
          IF(IEMSCT.EQ.2)THEN                                           FLX 0151
              RDNCNS=0.                                                 FLX 0152
              RUPCNS=SALB                                               FLX 0153
              EUPCNS=RUPCNS*CSZEN0(1)*S0TRN(1,1)                        FLX 0154
              EUPS(ML)=EUPCNS                                           FLX 0155
              EUPCS(ML)=EUPCNS                                          FLX 0156
              RUPCS(ML)=0.                                              FLX 0157
          ENDIF                                                         FLX 0158
C                                                                       FLX 0159
C         UPWARD ADDING LOOP STARTS FROM BOTTOM OF ATMOSPHERE.          FLX 0160
          IKP1=1                                                        FLX 0161
          DO 50 N=NLAYRS,1,-1                                           FLX 0162
C                                                                       FLX 0163
C             LAYER INDICES IN OPPOSITE DIRECTION IN ROUTINE LOOP.      FLX 0164
              IK=IKP1                                                   FLX 0165
              IKP1=IKP1+1                                               FLX 0166
C                                                                       FLX 0167
C             STORE MOLECULAR OPTICAL THICKNESS FOR CALCULATION         FLX 0168
C             OF SCATTERING ALBEDO IN ROUTINE LOOP.                     FLX 0169
              TAUT(1,IK)=TAUT(1,IK)+TAUM(K,N)*DPJ(K,N)                  FLX 0170
C                                                                       FLX 0171
C             TAU IS THE LAYER OPTICAL THICKNESS                        FLX 0172
              TAU=TAUM(K,N)+TCONT(IK)                                   FLX 0173
C                                                                       FLX 0174
C             USE TWO STREAM APPROXIMATION FOR THERMAL                  FLX 0175
C               EUP AND EDN ARE UPWARD AND DOWNWARD THERMAL FLUX        FLX 0176
C               FOR AN ISOLATED LAYER.  TDF AND REF ARE THE THERMAL     FLX 0177
C               TRANSMISSION AND REFLECTANCE FOR AN ISOLATED LAYER      FLX 0178
              C=RT3*(TAU-TSCAT(IK)*COSBAR(IK))                          FLX 0179
              A2C=RT3*(TAU-TSCAT(IK))                                   FLX 0180
              A2C2=A2C*C                                                FLX 0181
              AC=SQRT(A2C2)                                             FLX 0182
              IF(AC.GT.20.)THEN                                         FLX 0183
                  DENOM=C+AC                                            FLX 0184
                  REF(N)=(C-AC)/DENOM                                   FLX 0185
                  DENOM=.5*DENOM                                        FLX 0186
                  TDF(N)=0.                                             FLX 0187
                  IF(AC.LT.35.)TDF(N)=C*AC*EXP(-AC)/DENOM**2            FLX 0188
                  ACM1=AC-1.                                            FLX 0189
                  EUP(N)=(BBNDRY(IK)+ACM1*BBNDRY(IKP1))/DENOM           FLX 0190
                  EDN(N)=(BBNDRY(IKP1)+ACM1*BBNDRY(IK))/DENOM           FLX 0191
              ELSE                                                      FLX 0192
                  IF(AC.LT..08)THEN                                     FLX 0193
                      EXPAN=1.-AC*(3.-AC)/6.                            FLX 0194
                      ONEME0=AC*EXPAN                                   FLX 0195
                      A1ME0=A2C*EXPAN                                   FLX 0196
                      E0=1.-ONEME0                                      FLX 0197
                      ONEPE0=2.-ONEME0                                  FLX 0198
                      ONEME1=AC*(.5-AC*(4.-AC)/24.)                     FLX 0199
                      E1=1.-ONEME1                                      FLX 0200
                      E2=A2C2*(1.-AC+.55*A2C2)/3.                       FLX 0201
                  ELSE                                                  FLX 0202
                      E0=EXP(-AC)                                       FLX 0203
                      ONEME0=1.-E0                                      FLX 0204
                      ONEPE0=1.+E0                                      FLX 0205
                      E1=ONEME0/AC                                      FLX 0206
                      A1ME0=A2C*E1                                      FLX 0207
                      ONEME1=1.-E1                                      FLX 0208
                      E2=ONEME0-ONEPE0*ONEME1                           FLX 0209
                  ENDIF                                                 FLX 0210
                  HALFE1=.5*E1                                          FLX 0211
                  DENOM=(1.+(C-AC)*HALFE1)*(ONEPE0+A1ME0)               FLX 0212
                  TDF(N)=2*E0/DENOM                                     FLX 0213
                  REF(N)=(C-A2C)*HALFE1*ONEPE0/DENOM                    FLX 0214
                  COEF1=A1ME0*E1+E2                                     FLX 0215
                  COEF2=A1ME0*(ONEME1+E0)+ONEME0**2-E2                  FLX 0216
                  EUP(N)=(COEF1*BBNDRY(IK)+COEF2*BBNDRY(IKP1))/DENOM    FLX 0217
                  EDN(N)=(COEF1*BBNDRY(IKP1)+COEF2*BBNDRY(IK))/DENOM    FLX 0218
              ENDIF                                                     FLX 0219
C                                                                       FLX 0220
C             CALCULATE COMPOSITE FLUXES AND REFLECTANCES               FLX 0221
              COEF=TDF(N)/(1.-RUPCN*REF(N))                             FLX 0222
              EUPCN=EUP(N)+COEF*(EUPCN+EDN(N)*RUPCN)                    FLX 0223
              RUPCN=REF(N)+COEF*TDF(N)*RUPCN                            FLX 0224
              EUPC(N)=EUPCN                                             FLX 0225
              RUPC(N)=RUPCN                                             FLX 0226
              IF(IEMSCT.EQ.2)THEN                                       FLX 0227
C                                                                       FLX 0228
C                 CALCULATE VARIABLES FOR SOLAR HYBRID MODIFIED         FLX 0229
C                 DELTA EDDINGTON 2-STREAM APPROXIMATION.               FLX 0230
                  CALL TRLAY(TAU,TSCAT(IK),COSBAR(IK),CSZEN(IK),        FLX 0231
     1              S0DEP(1,IK),DEPRAT(IK),EX,TDFS(N),REFS(N))          FLX 0232
                  SOURC=CSZENX(IK)*S0TRN(1,IKP1)                        FLX 0233
                  EUPS(N)=SOURC*REFS(N)                                 FLX 0234
                  EDNS(N)=SOURC*TDFS(N)                                 FLX 0235
                  TDFS(N)=EX+TDFS(N)                                    FLX 0236
                  COEF=TDFS(N)/(1.-RUPCNS*REFS(N))                      FLX 0237
                  EUPCNS=EUPS(N)+COEF*(EUPCNS+EDNS(N)*RUPCNS)           FLX 0238
                  RUPCNS=REFS(N)+COEF*TDFS(N)*RUPCNS                    FLX 0239
                  EUPCS(N)=EUPCNS                                       FLX 0240
                  RUPCS(N)=RUPCNS                                       FLX 0241
              ENDIF                                                     FLX 0242
   50     CONTINUE                                                      FLX 0243
C                                                                       FLX 0244
C         NOW ADD DOWNWARD FROM TOP LAYER (N=1)                         FLX 0245
          EDNCN=0.                                                      FLX 0246
          IK=ML                                                         FLX 0247
          DPJ(K,ML)=DPJ(K,NLAYRS)                                       FLX 0248
          DNF(1,ML)=DNF(1,ML)+DPJ(K,ML)*EDNCN                           FLX 0249
          UPF(1,ML)=UPF(1,ML)+DPJ(K,ML)*EUPC(1)                         FLX 0250
          RDNCN=0.                                                      FLX 0251
          IF(IEMSCT.EQ.2)THEN                                           FLX 0252
              EDNCNS=0.                                                 FLX 0253
              DNFS(K,1)=EDNCNS                                          FLX 0254
              UPFS(K,1)=EUPCS(1)                                        FLX 0255
              RDNCNS=0.                                                 FLX 0256
          ENDIF                                                         FLX 0257
          NM1=1                                                         FLX 0258
          DO 60 N=2,ML                                                  FLX 0259
              IK=IK-1                                                   FLX 0260
              DENOM=1.-RDNCN*REF(NM1)                                   FLX 0261
              COEF=TDF(NM1)/DENOM                                       FLX 0262
              EDNCN=EDN(NM1)+COEF*(EDNCN+EUP(NM1)*RDNCN)                FLX 0263
              RDNCN=REF(NM1)+COEF*TDF(NM1)*RDNCN                        FLX 0264
              COEF=DPJ(K,N)/DENOM                                       FLX 0265
              UPF(1,IK)=UPF(1,IK)+COEF*(EUPC(N)+EDNCN*RUPC(N))          FLX 0266
              DNF(1,IK)=DNF(1,IK)+COEF*(EDNCN+EUPC(N)*RDNCN)            FLX 0267
              IF(IEMSCT.EQ.2)THEN                                       FLX 0268
                  DENOM=1.-RDNCNS*REFS(NM1)                             FLX 0269
                  COEF=TDFS(NM1)/DENOM                                  FLX 0270
                  EDNCNS=EDNS(NM1)+COEF*(EDNCNS+EUPS(NM1)*RDNCNS)       FLX 0271
                  RDNCNS=REFS(NM1)+COEF*TDFS(NM1)*RDNCNS                FLX 0272
                  COEF=DPJ(K,N)/DENOM                                   FLX 0273
                  UPFS(1,IK)=UPFS(1,IK)+COEF*(EUPCS(N)+EDNCNS*RUPCS(N)) FLX 0274
                  DNFS(1,IK)=DNFS(1,IK)+COEF*(EDNCNS+EUPCS(N)*RDNCNS)   FLX 0275
              ENDIF                                                     FLX 0276
   60     NM1=N                                                         FLX 0277
   70 CONTINUE                                                          FLX 0278
      WRITE(ISCRCH)(COSBAR(IK),TSCAT(IK),TAUT(1,IK),IK=1,NLAYRS),       FLX 0279
     1  (UPF(1,IK),DNF(1,IK),UPFS(1,IK),DNFS(1,IK),IK=1,ML)             FLX 0280
      RETURN                                                            FLX 0281
      END                                                               FLX 0282
