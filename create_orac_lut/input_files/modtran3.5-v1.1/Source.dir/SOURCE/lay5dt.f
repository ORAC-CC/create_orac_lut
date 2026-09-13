      SUBROUTINE LAY5DT(V,IKOFF,IPATH)                                  L5D 0001
C                                                                       L5D 0002
C     THIS ROUTINE DEFINES THE LAYER DEPENDENT 5 CM-1 DATA              L5D 0003
      INTEGER IKOFF,IPATH                                               L5D 0004
      REAL V                                                            L5D 0005
C                                                                       L5D 0006
C     THE TXNEW AND TXOLD ARRAYS CONTAIN:                               L5D 0007
C        1   ASYMMETRY PARAMETER WEIGHTED BY SCATTERING DEPTH           L5D 0008
C        2   INCREMENTAL AEROSOL SCATTERING OPTICAL DEPTH               L5D 0009
C        3   TOTAL O2 CONTINUUM TRANSMITTANCE                           L5D 0010
C        4   N2 CONTINUUM TRANSMITTANCE                                 L5D 0011
C        5   TOTAL H2O CONTINUUM TRANSMITTANCE                          L5D 0012
C        6   RAYLEIGH MOLECULAR SCATTERED TRANSMITTANCE                 L5D 0013
C        7   AEROSOL EXTINCTION                                         L5D 0014
C        8   TOTAL OZONE CONTINUUM TRANSMITTANCE                        L5D 0015
C        9   PRODUCT OF ALL CONTINUUM TRANSMITTANCES EXCEPT O2 AND HNO3 L5D 0016
C       10   AEROSOL ABSORPTION                                         L5D 0017
C       11   HNO3 TRANSMITTANCE                                         L5D 0018
C       12   MOLECULAR CONTINUUM OPTICAL DEPTH                          L5D 0019
C       13   INCREMENTAL AEROSOL + CLOUD EXTINCTION OPTICAL DEPTH       L5D 0020
C       14   TOTAL CONTINUUM OPTICAL DEPTH                              L5D 0021
C       15   RAYLEIGH MOLECULAR SCATTERING OPTICAL DEPTH                L5D 0022
C       16   CIRRUS CLOUD TRANSMITTANCE (ICLD = 20 ONLY)                L5D 0023
C       17   UV/VIS NO2 TRANSMISSION                                    L5D 0024
C       18   UV/VIS SO2 TRANSMISSION                                    L5D 0025
C       19   INCREMENTAL WATER DROPLET SCATTERING OPTICAL DEPTH         L5D 0026
C       20   INCREMENTAL ICE PARTICLE SCATTERING OPTICAL DEPTH          L5D 0027
C                                                                       L5D 0028
C     FOR LOWTRAN RUNS, TX CONTAINS MOLECULAR LINE CENTER TRANSMITTANCESL5D 0029
C       17=H2O  36=CO2  31=O3   47=N2O  44=CO   46=CH4                  L5D 0030
C       50=O2   54=NO   56=SO2  55=NO2  52=NH3                          L5D 0031
C                                                                       L5D 0032
      INCLUDE 'PARAM.LST'                                               L5D 0033
      INTEGER KPOINT                                                    L5D 0034
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     L5D 0035
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   L5D 0036
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   L5D 0037
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     L5D 0038
C                                                                       L5D 0039
C       PI       THE CONSTANT PI                                        L5D 0040
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       L5D 0041
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       L5D 0042
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         L5D 0043
      REAL PI,DEG,BIGNUM,BIGEXP                                         L5D 0044
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                L5D 0045
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     L5D 0046
      REAL TBOUND,SALB                                                  L5D 0047
      LOGICAL MODTRN                                                    L5D 0048
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   L5D 0049
     1  MODTRN                                                          L5D 0050
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA                       L5D 0051
      REAL VIS,WSS,WHH,RAINRT                                           L5D 0052
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,     L5D 0053
     1  RAINRT                                                          L5D 0054
      REAL EXTV,ABSV,ASYV                                               L5D 0055
      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       L5D 0056
      INTEGER IBND                                                      L5D 0057
      REAL A1,B1,C1,QA,CPS                                              L5D 0058
      COMMON/AABBCC/A1(11),B1(11),C1(11),IBND(11),QA(11),CPS(11)        L5D 0059
      REAL SIGO20,SIGO2A,SIGO2B,C0,CT1,CT2,ABB,CRSNO2,CRSSO2            L5D 0060
      COMMON/FRQ5/SIGO20,SIGO2A,SIGO2B,C0,CT1,CT2,ABB(19),CRSNO2,CRSSO2 L5D 0061
      REAL TXNEW,TXOLD                                                  L5D 0062
      COMMON/LAY5/TXNEW(20,LAYTHR,3),TXOLD(20,LAYTHR,3)                 L5D 0063
      SAVE /LAY5/                                                       L5D 0064
C                                                                       L5D 0065
C     DECLARE LOCAL VARIABLES                                           L5D 0066
      INTEGER IAER,K                                                    L5D 0067
      REAL ABSSM,SCTSM,ASYMSM,EXTSM,ASYMR,SCTTRM,                       L5D 0068
     1  SH2OT0,SH2OT1,SH2O,FH2O,TAVE                                    L5D 0069
C                                                                       L5D 0070
C     DECLARE FUNCTIONS                                                 L5D 0071
      REAL DBLTX                                                        L5D 0072
C                                                                       L5D 0073
C     LIST DATA                                                         L5D 0074
      INTEGER MAP(NAER)                                                 L5D 0075
      DATA MAP/7,12,13,14,16,66,67/                                     L5D 0076
C                                                                       L5D 0077
C     CONTRIBUTION OF RAIN TO EXTINCTION, SCATTERING AND ABSORPTION     L5D 0078
      EXTSM=0.                                                          L5D 0079
      ABSSM=0.                                                          L5D 0080
      SCTSM=0.                                                          L5D 0081
      ASYMR=0.                                                          L5D 0082
      IF(W(3).NE.0.)CALL RAIN(V,EXTSM,ABSSM,SCTSM,ASYMR)                L5D 0083
C                                                                       L5D 0084
C     CALCULATE INCREMENTAL OPTICAL DEPTH VARIABLES                     L5D 0085
      IF(IPATH.EQ.3)THEN                                                L5D 0086
          TXOLD( 1,IKOFF,IPATH)=TXNEW( 1,IKOFF,IPATH)                   L5D 0087
          TXOLD( 2,IKOFF,IPATH)=TXNEW( 2,IKOFF,IPATH)                   L5D 0088
          TXOLD(13,IKOFF,IPATH)=TXNEW(13,IKOFF,IPATH)                   L5D 0089
          TXOLD(19,IKOFF,IPATH)=TXNEW(19,IKOFF,IPATH)                   L5D 0090
          TXOLD(20,IKOFF,IPATH)=TXNEW(20,IKOFF,IPATH)                   L5D 0091
          ASYMSM=SCTSM*ASYMR                                            L5D 0092
          DO 10 IAER=1,NAER                                             L5D 0093
              SCTTRM=W(MAP(IAER))*(EXTV(IAER)-ABSV(IAER))               L5D 0094
              ASYMSM=ASYMSM+SCTTRM*ASYV(IAER)                           L5D 0095
C                                                                       L5D 0096
C             FOR THE SINGLE SCATTER SOLAR RADIANCE CALCULATIONS        L5D 0097
C             (ROUTINE SSRAD), THE AEROSOL AND CLOUD INCREMENTAL        L5D 0098
C             OPTICAL DEPTH CONTRIBUTIONS ARE KEPT SEPARATE.            L5D 0099
              IF(IAER.EQ.6)THEN                                         L5D 0100
                  TXNEW(19,IKOFF,IPATH)=SCTTRM                          L5D 0101
              ELSEIF(IAER.EQ.7)THEN                                     L5D 0102
                  TXNEW(20,IKOFF,IPATH)=SCTTRM                          L5D 0103
              ELSE                                                      L5D 0104
                  SCTSM=SCTSM+SCTTRM                                    L5D 0105
              ENDIF                                                     L5D 0106
              EXTSM=EXTSM+W(MAP(IAER))*EXTV(IAER)                       L5D 0107
              ABSSM=ABSSM+W(MAP(IAER))*ABSV(IAER)                       L5D 0108
   10     CONTINUE                                                      L5D 0109
          TXNEW( 1,IKOFF,IPATH)=ASYMSM                                  L5D 0110
          TXNEW( 2,IKOFF,IPATH)=SCTSM                                   L5D 0111
          TXNEW(13,IKOFF,IPATH)=EXTSM                                   L5D 0112
      ELSE                                                              L5D 0113
C                                                                       L5D 0114
C         ADD AEROSOL EXTINCTION AND ABSORPTION TO RAIN CONTRIBUTION    L5D 0115
          DO 20 IAER=1,NAER                                             L5D 0116
              EXTSM=EXTSM+W(MAP(IAER))*EXTV(IAER)                       L5D 0117
              ABSSM=ABSSM+W(MAP(IAER))*ABSV(IAER)                       L5D 0118
   20     CONTINUE                                                      L5D 0119
      ENDIF                                                             L5D 0120
C                                                                       L5D 0121
C     DETERMINE OPTICAL DEPTHS                                          L5D 0122
C       TX(4)   N2 CONTINUUM                                            L5D 0123
C       TX(5)   SUM OF H2O CONTINUUM CONTRIBUTIONS                      L5D 0124
C       TX(6)   RAYLEIGH MOLECULAR SCATTERING                           L5D 0125
C       TX(8)   VISIBLE AND ULTRAVIOLET O3 (= 0 IF ABBUV > 0)           L5D 0126
C       TX(11)  HNO3 ABSORPTION (LOWTRAN CALCULATIONS ONLY)             L5D 0127
C       TX(16)  CIRRUS CLOUDS                                           L5D 0128
C       TX(9)   SUM OF CONTINUUM CONTRIBUTIONS EXCLUDING O2 AND HNO3    L5D 0129
C       TX(3)   O2 CONTINUUM CONTRIBUTIONS                              L5D 0130
      TX(4)=ABB(4)*W(4)                                                 L5D 0131
      TX(5)=0.                                                          L5D 0132
      IF(W(5).GT.0.)THEN                                                L5D 0133
          SH2OT0=ABB(5)                                                 L5D 0134
          SH2OT1=ABB(9)                                                 L5D 0135
          FH2O=ABB(10)                                                  L5D 0136
          TAVE=296.-36.*W(9)/W(5)                                       L5D 0137
          CALL CKD(SH2OT0,SH2OT1,TAVE,SH2O,FH2O,V)                      L5D 0138
          TX(5)=1.E-20*(SH2O*W(5)+FH2O*W(10))                           L5D 0139
      ENDIF                                                             L5D 0140
      TX(6)=ABB(6)*W(6)                                                 L5D 0141
      TX(8)=ABB(8)*W(8)+C0*(.269*W(8)+CT1*W(59)+CT2*W(60))              L5D 0142
      TX(11)=ABB(11)*W(11)                                              L5D 0143
      TX(16)=0.                                                         L5D 0144
      IF(ICLD.EQ.20)TX(16)=2*W(16)                                      L5D 0145
      TX(3)=W(58)*ABB(17)+                                              L5D 0146
     1  SIGO20*(W(63)+SIGO2A*(W(1)-220.*W(63))+SIGO2B*W(2))             L5D 0147
      TX(64)=CRSNO2*W(64)                                               L5D 0148
      TX(65)=CRSSO2*W(65)                                               L5D 0149
      TX(9)=TX(4)+TX(5)+TX(6)+TX(8)+EXTSM+TX(16)                        L5D 0150
C                                                                       L5D 0151
C     STORE CUMULATIVE AEROSOL PARAMETERS FOR DIFFERENT VERTICAL REGIONSL5D 0152
      TX(10)=ABSSM                                                      L5D 0153
      TX(7)=EXTSM                                                       L5D 0154
C                                                                       L5D 0155
C     STORE OPTICAL THICKNESS PARAMETERS                                L5D 0156
      TXOLD(12,IKOFF,IPATH)=TXNEW(12,IKOFF,IPATH)                       L5D 0157
      TXNEW(12,IKOFF,IPATH)=TX(4)+TX(5)+TX(8)+TX(3)                     L5D 0158
      TXOLD(14,IKOFF,IPATH)=TXNEW(14,IKOFF,IPATH)                       L5D 0159
      TXNEW(14,IKOFF,IPATH)=TX(9)+TX(3)+TX(64)+TX(65)                   L5D 0160
      TXOLD(15,IKOFF,IPATH)=TXNEW(15,IKOFF,IPATH)                       L5D 0161
      TXNEW(15,IKOFF,IPATH)=TX(6)                                       L5D 0162
      DO 30 K=3,11                                                      L5D 0163
          TXOLD(K,IKOFF,IPATH)=TXNEW(K,IKOFF,IPATH)                     L5D 0164
          IF(TX(K).LE.BIGEXP)THEN                                       L5D 0165
              TXNEW(K,IKOFF,IPATH)=EXP(-TX(K))                          L5D 0166
          ELSE                                                          L5D 0167
              TXNEW(K,IKOFF,IPATH)=1./BIGNUM                            L5D 0168
          ENDIF                                                         L5D 0169
   30 CONTINUE                                                          L5D 0170
      TXOLD(16,IKOFF,IPATH)=TXNEW(16,IKOFF,IPATH)                       L5D 0171
      IF(TX(16).LE.BIGEXP)THEN                                          L5D 0172
          TXNEW(16,IKOFF,IPATH)=EXP(-TX(16))                            L5D 0173
      ELSE                                                              L5D 0174
          TXNEW(16,IKOFF,IPATH)=1./BIGNUM                               L5D 0175
      ENDIF                                                             L5D 0176
      TXOLD(17,IKOFF,IPATH)=TXNEW(17,IKOFF,IPATH)                       L5D 0177
      IF(TX(64).LE.BIGEXP)THEN                                          L5D 0178
          TXNEW(17,IKOFF,IPATH)=EXP(-TX(64))                            L5D 0179
      ELSE                                                              L5D 0180
          TXNEW(17,IKOFF,IPATH)=1./BIGNUM                               L5D 0181
      ENDIF                                                             L5D 0182
C                                                                       L5D 0183
      TXOLD(18,IKOFF,IPATH)=TXNEW(18,IKOFF,IPATH)                       L5D 0184
      IF(TX(65).LE.BIGEXP)THEN                                          L5D 0185
          TXNEW(18,IKOFF,IPATH)=EXP(-TX(65))                            L5D 0186
      ELSE                                                              L5D 0187
          TXNEW(18,IKOFF,IPATH)=1./BIGNUM                               L5D 0188
      ENDIF                                                             L5D 0189
C                                                                       L5D 0190
      IF(.NOT.MODTRN)THEN                                               L5D 0191
C                                                                       L5D 0192
C         LOWTRAN7 DOUBLE EXPONENTIAL MOLECULAR TRANSMITTANCES          L5D 0193
          DO 40 K=1,11                                                  L5D 0194
              TX(KPOINT(K))=1.                                          L5D 0195
              IF(CPS(K).GT.-20.)                                        L5D 0196
     1          TX(KPOINT(K))=DBLTX(W(IBND(K)),CPS(K),QA(K))            L5D 0197
   40     CONTINUE                                                      L5D 0198
      ENDIF                                                             L5D 0199
      RETURN                                                            L5D 0200
      END                                                               L5D 0201
