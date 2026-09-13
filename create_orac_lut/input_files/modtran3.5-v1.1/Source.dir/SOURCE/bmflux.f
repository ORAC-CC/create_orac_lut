      SUBROUTINE BMFLUX(ML,KNTRVL,IEMSCT,SALB,S0)                       BMF 0001
C                                                                       BMF 0002
C     FLUX ADDING ROUTINE                                               BMF 0003
C                                                                       BMF 0004
C     DECLARE INPUTS                                                    BMF 0005
C       ML       NUMBER OF ATMOSPHERIC LAYER BOUNDARIES                 BMF 0006
C       KNTRVL   NUMBER OF INTERVALS IN CORRELATED-K APPROACH           BMF 0007
C       IEMSCT   FLAG = 1 FOR THERMAL SCATTER ONLY                      BMF 0008
C                     = 2 FOR THERMAL AND SOLAR SCATTER                 BMF 0009
C       SALB     SURFACE SCATTERING ALBEDO                              BMF 0010
C       S0       SOURCE IRRADIANCE [W CM-2 / CM-1]                      BMF 0011
      INTEGER ML,KNTRVL,IEMSCT                                          BMF 0012
      REAL SALB,S0                                                      BMF 0013
C                                                                       BMF 0014
C     LIST PARAMETERS                                                   BMF 0015
      INCLUDE 'PARAM.LST'                                               BMF 0016
C                                                                       BMF 0017
C     LIST COMMONS:                                                     BMF 0018
C     COMMON /MSRD/                                                     BMF 0019
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            BMF 0020
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             BMF 0021
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               BMF 0022
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    BMF 0023
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  BMF 0024
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. BMF 0025
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               BMF 0026
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         BMF 0027
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          BMF 0028
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              BMF 0029
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  BMF 0030
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               BMF 0031
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           BMF 0032
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     BMF 0033
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   BMF 0034
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       BMF 0035
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     BMF 0036
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              BMF 0037
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 BMF 0038
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          BMF 0039
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             BMF 0040
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                BMF 0041
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      BMF 0042
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         BMF 0043
      REAL EDN,EUP,EUPC,TDF,RUPC,REF,EDNS,EUPS,EUPCS,TDFS,RUPCS,REFS    BMF 0044
      COMMON/FLUX/EDN(LAYDIM),EUP(LAYDIM),EUPC(LAYDIM),TDF(LAYDIM),     BMF 0045
     1  RUPC(LAYDIM),REF(LAYDIM),EDNS(LAYDIM),EUPS(LAYDIM),             BMF 0046
     2  EUPCS(LAYDIM),TDFS(LAYDIM),RUPCS(LAYDIM),REFS(LAYDIM)           BMF 0047
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               BMF 0048
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           BMF 0049
C                                                                       BMF 0050
C       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.       BMF 0051
C       UPFLX    LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].    BMF 0052
C       DNFLX    LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].  BMF 0053
C       UPFLXS   BOUNDARY UPWARD SCATTERED SOLAR FLUX [W CM-2 / CM-1].  BMF 0054
C       DNFLXS   BOUNDARY DOWNWARD SCATTERED SOLAR FLUX [W CM-2 / CM-1].BMF 0055
C       NTFLX    LAYER BOUNDARY NET (THERMAL PLUS SCATTERED SOLAR       BMF 0056
C                PLUS DIRECT SOLAR) UPWARD FLUX [W CM-2 / CM-1].        BMF 0057
      REAL SUBINT,UPFLX,DNFLX,UPFLXS,DNFLXS,NTFLX                       BMF 0058
      COMMON/NETFLX/SUBINT(NKSUB),UPFLX(LAYDIM),DNFLX(LAYDIM),          BMF 0059
     1  UPFLXS(LAYDIM),DNFLXS(LAYDIM),NTFLX(LAYDIM)                     BMF 0060
C                                                                       BMF 0061
C     INTERNAL VARIABLES:                                               BMF 0062
C       EUP(N)    INTRINSIC UPWARD THERMAL FLUX OF LAYER N              BMF 0063
C       EDN(N)    INTRINSIC DOWNWARD THERMAL FLUX OF LAYER N            BMF 0064
C       TDF(N)    INTRINSIC THERMAL TRANSMISSION OF LAYER N             BMF 0065
C       REF(N)    INTRINSIC THERMAL REFLECTANCE OF LAYER N              BMF 0066
C       EUPS(N)   INTRINSIC UPWARD SOLAR FLUX OF LAYER N                BMF 0067
C       EDNS(N)   INTRINSIC DOWNWARD SOLAR FLUX OF LAYER N              BMF 0068
C       TDFS(N)   INTRINSIC SOLAR TRANSMISSION OF LAYER N               BMF 0069
C       REFS(N)   INTRINSIC SOLAR REFLECTANCE OF LAYER N                BMF 0070
C                                                                       BMF 0071
C     DECLARE LOCAL VARIABLES                                           BMF 0072
      INTEGER NLAYRS,IKP1,N,NM1,IK,INTRVL                               BMF 0073
      REAL RDNCN,RUPCN,EUPCN,RDNCNS,RUPCNS,EUPCNS,TAU,A1ME0,A2C,        BMF 0074
     1  EDNCN,EDNCNS,AC,DENOM,RT3,COEF,A2C2,C,HALFE1,E0,EXPAN,          BMF 0075
     2  ONEME0,E1,COEF1,COEF2,ACM1,ONEPE0,ONEME1,E2,EX,SOURC            BMF 0076
C                                                                       BMF 0077
C     LIST DATA                                                         BMF 0078
      DATA RT3/1.7320508/                                               BMF 0079
C                                                                       BMF 0080
C     NUMBER OF LAYERS                                                  BMF 0081
      NLAYRS=ML-1                                                       BMF 0082
C                                                                       BMF 0083
C     LOOP OVER K BINS                                                  BMF 0084
      DO 30 INTRVL=1,KNTRVL                                             BMF 0085
C                                                                       BMF 0086
C         COMPOSITE DOWNWARD REFLECTION                                 BMF 0087
          RDNCN=0.                                                      BMF 0088
C                                                                       BMF 0089
C         DEFINE INITIAL UPWARD COMPOSITE SURFACE REFLECTANCE           BMF 0090
          RUPCN=SALB                                                    BMF 0091
C                                                                       BMF 0092
C         SURFACE EMISSION                                              BMF 0093
          EUPCN=(1.-RUPCN)*BBGRND                                       BMF 0094
          EUP(ML)=EUPCN                                                 BMF 0095
          EUPC(ML)=EUPCN                                                BMF 0096
          RUPC(ML)=0.                                                   BMF 0097
          IF(IEMSCT.EQ.2)THEN                                           BMF 0098
              RDNCNS=0.                                                 BMF 0099
              RUPCNS=SALB                                               BMF 0100
              EUPCNS=RUPCNS*CSZEN0(1)*S0TRN(INTRVL,1)                   BMF 0101
              EUPS(ML)=EUPCNS                                           BMF 0102
              EUPCS(ML)=EUPCNS                                          BMF 0103
              RUPCS(ML)=0.                                              BMF 0104
          ENDIF                                                         BMF 0105
C                                                                       BMF 0106
C         UPWARD ADDING LOOP STARTS FROM BOTTOM OF ATMOSPHERE.          BMF 0107
          IKP1=1                                                        BMF 0108
          DO 10 N=NLAYRS,1,-1                                           BMF 0109
C                                                                       BMF 0110
C             LAYER INDICES IN OPPOSITE DIRECTION IN ROUTINE LOOP.      BMF 0111
              IK=IKP1                                                   BMF 0112
              IKP1=IKP1+1                                               BMF 0113
              TAU=TAUT(INTRVL,IK)                                       BMF 0114
C                                                                       BMF 0115
C             USE TWO STREAM APPROXIMATION FOR THERMAL                  BMF 0116
              C=RT3*(TAU-TSCAT(IK)*COSBAR(IK))                          BMF 0117
              A2C=RT3*(TAU-TSCAT(IK))                                   BMF 0118
              A2C2=A2C*C                                                BMF 0119
              AC=SQRT(A2C2)                                             BMF 0120
              IF(AC.GT.20.)THEN                                         BMF 0121
                  DENOM=C+AC                                            BMF 0122
                  REF(N)=(C-AC)/DENOM                                   BMF 0123
                  DENOM=.5*DENOM                                        BMF 0124
                  TDF(N)=0.                                             BMF 0125
                  IF(AC.LT.35.)TDF(N)=C*AC*EXP(-AC)/DENOM**2            BMF 0126
                  ACM1=AC-1.                                            BMF 0127
                  EUP(N)=(BBNDRY(IK)+ACM1*BBNDRY(IKP1))/DENOM           BMF 0128
                  EDN(N)=(BBNDRY(IKP1)+ACM1*BBNDRY(IK))/DENOM           BMF 0129
              ELSE                                                      BMF 0130
                  IF(AC.LT..08)THEN                                     BMF 0131
                      EXPAN=1.-AC*(3.-AC)/6.                            BMF 0132
                      ONEME0=AC*EXPAN                                   BMF 0133
                      A1ME0=A2C*EXPAN                                   BMF 0134
                      E0=1.-ONEME0                                      BMF 0135
                      ONEPE0=2.-ONEME0                                  BMF 0136
                      ONEME1=AC*(.5-AC*(4.-AC)/24.)                     BMF 0137
                      E1=1.-ONEME1                                      BMF 0138
                      E2=A2C2*(1.-AC+.55*A2C2)/3.                       BMF 0139
                  ELSE                                                  BMF 0140
                      E0=EXP(-AC)                                       BMF 0141
                      ONEME0=1.-E0                                      BMF 0142
                      ONEPE0=1.+E0                                      BMF 0143
                      E1=ONEME0/AC                                      BMF 0144
                      A1ME0=A2C*E1                                      BMF 0145
                      ONEME1=1.-E1                                      BMF 0146
                      E2=ONEME0-ONEPE0*ONEME1                           BMF 0147
                  ENDIF                                                 BMF 0148
                  HALFE1=.5*E1                                          BMF 0149
                  DENOM=(1.+(C-AC)*HALFE1)*(ONEPE0+A1ME0)               BMF 0150
                  TDF(N)=2*E0/DENOM                                     BMF 0151
                  REF(N)=(C-A2C)*HALFE1*ONEPE0/DENOM                    BMF 0152
                  COEF1=A1ME0*E1+E2                                     BMF 0153
                  COEF2=A1ME0*(ONEME1+E0)+ONEME0**2-E2                  BMF 0154
                  EUP(N)=(COEF1*BBNDRY(IK)+COEF2*BBNDRY(IKP1))/DENOM    BMF 0155
                  EDN(N)=(COEF1*BBNDRY(IKP1)+COEF2*BBNDRY(IK))/DENOM    BMF 0156
              ENDIF                                                     BMF 0157
C                                                                       BMF 0158
C             CALCULATE COMPOSITE FLUXES AND REFLECTANCES               BMF 0159
              COEF=TDF(N)/(1.-RUPCN*REF(N))                             BMF 0160
              EUPCN=EUP(N)+COEF*(EUPCN+EDN(N)*RUPCN)                    BMF 0161
              RUPCN=REF(N)+COEF*TDF(N)*RUPCN                            BMF 0162
              EUPC(N)=EUPCN                                             BMF 0163
              RUPC(N)=RUPCN                                             BMF 0164
              IF(IEMSCT.EQ.2)THEN                                       BMF 0165
C                                                                       BMF 0166
C                 CALCULATE VARIABLES FOR SOLAR HYBRID MODIFIED         BMF 0167
C                 DELTA EDDINGTON 2-STREAM APPROXIMATION.               BMF 0168
                  CALL TRLAY(TAU,TSCAT(IK),COSBAR(IK),CSZEN(IK),        BMF 0169
     1              S0DEP(INTRVL,IK),DEPRAT(IK),EX,TDFS(N),REFS(N))     BMF 0170
                  IF(EX.GT.1.E-30)THEN                                  BMF 0171
                      SOURC=S0TRN(INTRVL,IK)/EX                         BMF 0172
                  ELSE                                                  BMF 0173
                      SOURC=S0*EXP((DEPRAT(IK)-1.)*S0DEP(INTRVL,IK))    BMF 0174
                  ENDIF                                                 BMF 0175
                  SOURC=CSZENX(IK)*SOURC                                BMF 0176
                  EUPS(N)=SOURC*REFS(N)                                 BMF 0177
                  EDNS(N)=SOURC*TDFS(N)                                 BMF 0178
                  TDFS(N)=EX+TDFS(N)                                    BMF 0179
                  COEF=TDFS(N)/(1.-RUPCNS*REFS(N))                      BMF 0180
                  EUPCNS=EUPS(N)+COEF*(EUPCNS+EDNS(N)*RUPCNS)           BMF 0181
                  RUPCNS=REFS(N)+COEF*TDFS(N)*RUPCNS                    BMF 0182
                  EUPCS(N)=EUPCNS                                       BMF 0183
                  RUPCS(N)=RUPCNS                                       BMF 0184
              ENDIF                                                     BMF 0185
   10     CONTINUE                                                      BMF 0186
C                                                                       BMF 0187
C         NOW ADD DOWNWARD FROM TOP LAYER (N=1)                         BMF 0188
          EDNCN=0.                                                      BMF 0189
          IK=ML                                                         BMF 0190
          UPF(INTRVL,IK)=EUPC(1)                                        BMF 0191
          DNF(INTRVL,IK)=0.                                             BMF 0192
          UPFLX(IK)=UPFLX(IK)+SUBINT(INTRVL)*EUPC(1)                    BMF 0193
          RDNCN=0.                                                      BMF 0194
          IF(IEMSCT.EQ.2)THEN                                           BMF 0195
              EDNCNS=0.                                                 BMF 0196
              DNFS(INTRVL,IK)=0.                                        BMF 0197
              UPFS(INTRVL,IK)=EUPCS(1)                                  BMF 0198
              UPFLXS(IK)=UPFLXS(IK)+SUBINT(INTRVL)*EUPCS(1)             BMF 0199
              NTFLX(IK)=NTFLX(IK)+                                      BMF 0200
     1          SUBINT(INTRVL)*(EUPCS(1)-CSZEN0(IK)*S0TRN(INTRVL,IK))   BMF 0201
              RDNCNS=0.                                                 BMF 0202
          ENDIF                                                         BMF 0203
          NM1=1                                                         BMF 0204
          DO 20 N=2,ML                                                  BMF 0205
              IK=IK-1                                                   BMF 0206
              DENOM=1.-RDNCN*REF(NM1)                                   BMF 0207
              COEF=TDF(NM1)/DENOM                                       BMF 0208
              EDNCN=EDN(NM1)+COEF*(EDNCN+EUP(NM1)*RDNCN)                BMF 0209
              RDNCN=REF(NM1)+COEF*TDF(NM1)*RDNCN                        BMF 0210
              UPF(INTRVL,IK)=(EUPC(N)+EDNCN*RUPC(N))/DENOM              BMF 0211
              DNF(INTRVL,IK)=(EDNCN+EUPC(N)*RDNCN)/DENOM                BMF 0212
              UPFLX(IK)=UPFLX(IK)+SUBINT(INTRVL)*UPF(INTRVL,IK)         BMF 0213
              DNFLX(IK)=DNFLX(IK)+SUBINT(INTRVL)*DNF(INTRVL,IK)         BMF 0214
              NTFLX(IK)=NTFLX(IK)+SUBINT(INTRVL)*                       BMF 0215
     1          (EUPC(N)*(1.-RDNCN)-EDNCN*(1.-RUPC(N)))/DENOM           BMF 0216
              IF(IEMSCT.EQ.2)THEN                                       BMF 0217
                  DENOM=1.-RDNCNS*REFS(NM1)                             BMF 0218
                  COEF=TDFS(NM1)/DENOM                                  BMF 0219
                  EDNCNS=EDNS(NM1)+COEF*(EDNCNS+EUPS(NM1)*RDNCNS)       BMF 0220
                  RDNCNS=REFS(NM1)+COEF*TDFS(NM1)*RDNCNS                BMF 0221
                  UPFS(INTRVL,IK)=(EUPCS(N)+EDNCNS*RUPCS(N))/DENOM      BMF 0222
                  DNFS(INTRVL,IK)=(EDNCNS+EUPCS(N)*RDNCNS)/DENOM        BMF 0223
                  UPFLXS(IK)=UPFLXS(IK)+SUBINT(INTRVL)*UPFS(INTRVL,IK)  BMF 0224
                  DNFLXS(IK)=DNFLXS(IK)+SUBINT(INTRVL)*DNFS(INTRVL,IK)  BMF 0225
                  NTFLX(IK)=NTFLX(IK)-                                  BMF 0226
     1              SUBINT(INTRVL)*(CSZEN0(IK)*S0TRN(INTRVL,IK)-        BMF 0227
     2              (EUPCS(N)*(1.-RDNCNS)-EDNCNS*(1.-RUPCS(N)))/DENOM)  BMF 0228
              ENDIF                                                     BMF 0229
   20     NM1=N                                                         BMF 0230
   30 CONTINUE                                                          BMF 0231
      NTFLX(ML)=NTFLX(ML)+UPFLX(ML)                                     BMF 0232
C                                                                       BMF 0233
C     PRESENTLY, FLUXES ARE NOT SAVED IF THE CORRELATED-K METHOD        BMF 0234
C     IS USED SINCE THE SCRATCH FILE CAN GROW TOO LARGE.                BMF 0235
      IF(KNTRVL.EQ.1)WRITE(ISCRCH)(COSBAR(IK),TSCAT(IK),                BMF 0236
     1  (TAUT(INTRVL,IK),INTRVL=1,KNTRVL),IK=1,NLAYRS),                 BMF 0237
     2  ((UPF(INTRVL,IK),DNF(INTRVL,IK),                                BMF 0238
     3  UPFS(INTRVL,IK),DNFS(INTRVL,IK),INTRVL=1,KNTRVL),IK=1,ML)       BMF 0239
      RETURN                                                            BMF 0240
      END                                                               BMF 0241
