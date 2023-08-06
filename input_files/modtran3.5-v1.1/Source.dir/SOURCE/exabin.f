      SUBROUTINE EXABIN(ICH)                                            EXA 0001
      INTEGER ICH(4)                                                    EXA 0002
      INCLUDE 'PARAM.LST'                                               EXA 0003
C                                                                       EXA 0004
C      LOADS EXTINCTION, ABSORPTION AND ASYMMETRY COEFFICIENTS          EXA 0005
C      FOR THE FOUR AEROSOL ALTITUDE REGIONS                            EXA 0006
C                                                                       EXA 0007
C      MODIFIED FOR ASYMMETRY - JAN 1986 (A.E.R. INC.)                  EXA 0008
C                                                                       EXA 0009
      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   EXA 0010
     1    RAINRT                                                        EXA 0011
      COMMON /CARD2D/ IREG(4),ALTB(4),IREGC(4)                          EXA 0012
C                                                                       EXA 0013
      INTEGER KPOINT                                                    EXA 0014
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     EXA 0015
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   EXA 0016
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   EXA 0017
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     EXA 0018
      COMMON/EXTD/VX0(NWAVLN),                                          EXA 0019
     1  RUREXT(NWAVLN,4),RURABS(NWAVLN,4),RURSYM(NWAVLN,4),             EXA 0020
     2  URBEXT(NWAVLN,4),URBABS(NWAVLN,4),URBSYM(NWAVLN,4),             EXA 0021
     3  OCNEXT(NWAVLN,4),OCNABS(NWAVLN,4),OCNSYM(NWAVLN,4),             EXA 0022
     4  TROEXT(NWAVLN,4),TROABS(NWAVLN,4),TROSYM(NWAVLN,4),             EXA 0023
     5  FG1EXT(NWAVLN),FG1ABS(NWAVLN),FG1SYM(NWAVLN),                   EXA 0024
     6  FG2EXT(NWAVLN),FG2ABS(NWAVLN),FG2SYM(NWAVLN),                   EXA 0025
     7  BSTEXT(NWAVLN),BSTABS(NWAVLN),BSTSYM(NWAVLN),                   EXA 0026
     8  AVOEXT(NWAVLN),AVOABS(NWAVLN),AVOSYM(NWAVLN),                   EXA 0027
     9  FVOEXT(NWAVLN),FVOABS(NWAVLN),FVOSYM(NWAVLN),                   EXA 0028
     &  DMEEXT(NWAVLN),DMEABS(NWAVLN),DMESYM(NWAVLN),                   EXA 0029
     1  CCUEXT(NWAVLN),CCUABS(NWAVLN),CCUSYM(NWAVLN),                   EXA 0030
     2  CALEXT(NWAVLN),CALABS(NWAVLN),CALSYM(NWAVLN),                   EXA 0031
     3  CSTEXT(NWAVLN),CSTABS(NWAVLN),CSTSYM(NWAVLN),                   EXA 0032
     4  CSCEXT(NWAVLN),CSCABS(NWAVLN),CSCSYM(NWAVLN),                   EXA 0033
     5  CNIEXT(NWAVLN),CNIABS(NWAVLN),CNISYM(NWAVLN)                    EXA 0034
      COMMON/CIRR/CI64XT(NWAVLN),CI64AB(NWAVLN),CI64G(NWAVLN),          EXA 0035
     1            CIR4XT(NWAVLN),CIR4AB(NWAVLN),CIR4G(NWAVLN)           EXA 0036
      DIMENSION RHZONE(4)                                               EXA 0037
      DIMENSION ELWCR(4),ELWCU(4),ELWCM(4),ELWCT(4)                     EXA 0038
      REAL MDLWC                                                        EXA 0039
      DATA RHZONE/0.,70.,80.,99./                                       EXA 0040
      DATA ELWCR/3.517E-04,3.740E-04,4.439E-04,9.529E-04/               EXA 0041
      DATA ELWCM/4.675E-04,6.543E-04,1.166E-03,3.154E-03/               EXA 0042
      DATA ELWCU/3.102E-04,3.802E-04,4.463E-04,9.745E-04/               EXA 0043
      DATA ELWCT/1.735E-04,1.820E-04,2.020E-04,2.408E-04/               EXA 0044
      DATA AFLWC/1.295E-02/,RFLWC/1.804E-03/,CULWC/7.683E-03/           EXA 0045
      DATA ASLWC/4.509E-03/,STLWC/5.272E-03/,SCLWC/4.177E-03/           EXA 0046
      DATA SNLWC/7.518E-03/,BSLWC/1.567E-04/,FVLWC/5.922E-04/           EXA 0047
      DATA AVLWC/1.675E-04/,MDLWC/4.775E-04/                            EXA 0048
      DO 2 I = 1,4                                                      EXA 0049
2     AWCCON(I) = 0                                                     EXA 0050
CC    I1=1                                                              EXA 0051
C                                                                       EXA 0052
C     "NWAVLN" VALUES CALCULATED IN AEREXT USING ROUTINE GAMFOG         EXA 0053
      NE=NWAVLN-1                                                       EXA 0054
C                                                                       EXA 0055
CC    IF (IHAZE.EQ.7) I1=2                                              EXA 0056
CC    IF(IHAZE.EQ.3)  I1 = 2                                            EXA 0057
      DO 185 M=1 ,4                                                     EXA 0058
      IF(IREG(M) .NE. 0)        GO TO 185                               EXA 0059
      ITA=ICH(M)                                                        EXA 0060
      ITC=ICH(M)-7                                                      EXA 0061
      ITAS = ITA                                                        EXA 0062
      IF(IREGC(M) .NE. 0) GO TO 100                                     EXA 0063
12    CONTINUE                                                          EXA 0064
      WRH=W(15)                                                         EXA 0065
      IF (ICH(M).EQ.6.AND.M.NE.1) WRH=70.                               EXA 0066
C     THIS CODING  DOES NOT ALLOW TROP RH DEPENDENT  ABOVE EH(7,I)      EXA 0067
C     DEFAULTS TO TROPOSPHERIC AT 70. PERCENT                           EXA 0068
      DO 10 I=2,4                                                       EXA 0069
      IF (WRH.LT.RHZONE(I)) GO TO 15                                    EXA 0070
   10 CONTINUE                                                          EXA 0071
      I=4                                                               EXA 0072
   15 II=I-1                                                            EXA 0073
      IF(WRH.GT.0.0.AND.WRH.LT.99.)X=ALOG(100.0-WRH)                    EXA 0074
      X1=ALOG(100.0-RHZONE(II))                                         EXA 0075
      X2=ALOG(100.0-RHZONE(I))                                          EXA 0076
      IF (WRH.GE.99.0) X=X2                                             EXA 0077
      IF (WRH.LE.0.0) X=X1                                              EXA 0078
17    DO 80 N=1,NE                                                      EXA 0079
      ITA = ITAS                                                        EXA 0080
      IF(ITA.EQ.3. AND. M.EQ.1) GO TO 18                                EXA 0081
      ABSC(M,N)=0.                                                      EXA 0082
      EXTC(M,N)=0.                                                      EXA 0083
      ASYM(M,N)=0.                                                      EXA 0084
      IF(ITA.GT.6) GO TO 45                                             EXA 0085
      IF(ITA.LE.0) GO TO 80                                             EXA 0086
18    IF(N.GE.41. AND. ITA.EQ.3) ITA = 4                                EXA 0087
C     RH DEPENDENT AEROSOLS                                             EXA 0088
      GO TO (20,20,22,25,30,35), ITA                                    EXA 0089
   20 Y2=ALOG(RUREXT(N,I))                                              EXA 0090
      Y1=ALOG(RUREXT(N,II))                                             EXA 0091
      Z2=ALOG(RURABS(N,I))                                              EXA 0092
      Z1=ALOG(RURABS(N,II))                                             EXA 0093
      A2=ALOG(RURSYM(N,I))                                              EXA 0094
      A1=ALOG(RURSYM(N,II))                                             EXA 0095
      E2=ALOG(ELWCR(I))                                                 EXA 0096
      E1=ALOG(ELWCR(II))                                                EXA 0097
      GO TO 40                                                          EXA 0098
  22  IF(M.GT.1) GO TO 25                                               EXA 0099
      A2=ALOG(OCNSYM(N,I))                                              EXA 0100
      A1=ALOG(OCNSYM(N,II))                                             EXA 0101
      A=A1+(A2-A1)*(X-X1)/(X2-X1)                                       EXA 0102
      ASYM(M,N)=EXP(A)                                                  EXA 0103
      E2=ALOG(ELWCM(I))                                                 EXA 0104
      E1=ALOG(ELWCM(II))                                                EXA 0105
C                                                                       EXA 0106
C     NAVY MARITIME AEROSOL CHANGES TO MARINE IN MICROWAVE              EXA 0107
C     NO NEED TO DEFINE EQUIVALENT WATER                                EXA 0108
C                                                                       EXA 0109
      GO TO 80                                                          EXA 0110
   25 Y2=ALOG(OCNEXT(N,I))                                              EXA 0111
      Y1=ALOG(OCNEXT(N,II))                                             EXA 0112
      Z2=ALOG(OCNABS(N,I))                                              EXA 0113
      Z1=ALOG(OCNABS(N,II))                                             EXA 0114
      A2=ALOG(OCNSYM(N,I))                                              EXA 0115
      A1=ALOG(OCNSYM(N,II))                                             EXA 0116
      E2=ALOG(ELWCM(I))                                                 EXA 0117
      E1=ALOG(ELWCM(II))                                                EXA 0118
      GO TO 40                                                          EXA 0119
   30 Y2=ALOG(URBEXT(N,I))                                              EXA 0120
      Y1=ALOG(URBEXT(N,II))                                             EXA 0121
      Z2=ALOG(URBABS(N,I))                                              EXA 0122
      Z1=ALOG(URBABS(N,II))                                             EXA 0123
      A2=ALOG(URBSYM(N,I))                                              EXA 0124
      A1=ALOG(URBSYM(N,II))                                             EXA 0125
      E2=ALOG(ELWCU(I))                                                 EXA 0126
      E1=ALOG(ELWCU(II))                                                EXA 0127
      GO TO 40                                                          EXA 0128
   35 Y2=ALOG(TROEXT(N,I))                                              EXA 0129
      Y1=ALOG(TROEXT(N,II))                                             EXA 0130
      Z2=ALOG(TROABS(N,I))                                              EXA 0131
      Z1=ALOG(TROABS(N,II))                                             EXA 0132
      A2=ALOG(TROSYM(N,I))                                              EXA 0133
      A1=ALOG(TROSYM(N,II))                                             EXA 0134
      E2=ALOG(ELWCT(I))                                                 EXA 0135
      E1=ALOG(ELWCT(II))                                                EXA 0136
   40 Y=Y1+(Y2-Y1)*(X-X1)/(X2-X1)                                       EXA 0137
      ZK=Z1+(Z2-Z1)*(X-X1)/(X2-X1)                                      EXA 0138
      A=A1+(A2-A1)*(X-X1)/(X2-X1)                                       EXA 0139
      ABSC(M,N)=EXP(ZK)                                                 EXA 0140
      EXTC(M,N)=EXP(Y)                                                  EXA 0141
      ASYM(M,N)=EXP(A)                                                  EXA 0142
      IF(N.EQ.1) EC=E1+(E2-E1)*(X-X1)/(X2-X1)                           EXA 0143
      IF(N.EQ.1) AWCCON(M)=EXP(EC)                                      EXA 0144
      GO TO 80                                                          EXA 0145
   45 IF (ITA.GT.19) GO TO 75                                           EXA 0146
      IF (ITC.LT.1) GO TO 80                                            EXA 0147
      GO TO (50,55,80,60,65,70,65,70,60,60,70,75), ITC                  EXA 0148
   50 ABSC(M,N)=FG1ABS(N)                                               EXA 0149
      EXTC(M,N)=FG1EXT(N)                                               EXA 0150
      ASYM(M,N)=FG1SYM(N)                                               EXA 0151
      IF(N.EQ.1) AWCCON(M)=AFLWC                                        EXA 0152
      GO TO 80                                                          EXA 0153
   55 ABSC(M,N)=FG2ABS(N)                                               EXA 0154
      EXTC(M,N)=FG2EXT(N)                                               EXA 0155
      ASYM(M,N)=FG2SYM(N)                                               EXA 0156
      IF(N.EQ.1) AWCCON(M)=RFLWC                                        EXA 0157
      GO TO 80                                                          EXA 0158
   60 ABSC(M,N)=BSTABS(N)                                               EXA 0159
      EXTC(M,N)=BSTEXT(N)                                               EXA 0160
      ASYM(M,N)=BSTSYM(N)                                               EXA 0161
      IF(N.EQ.1) AWCCON(M)=BSLWC                                        EXA 0162
      GO TO 80                                                          EXA 0163
   65 ABSC(M,N)=AVOABS(N)                                               EXA 0164
      EXTC(M,N)=AVOEXT(N)                                               EXA 0165
      ASYM(M,N)=AVOSYM(N)                                               EXA 0166
      IF(N.EQ.1) AWCCON(M)=AVLWC                                        EXA 0167
      GO TO 80                                                          EXA 0168
   70 ABSC(M,N)=FVOABS(N)                                               EXA 0169
      EXTC(M,N)=FVOEXT(N)                                               EXA 0170
      ASYM(M,N)=FVOSYM(N)                                               EXA 0171
      IF(N.EQ.1) AWCCON(M)=FVLWC                                        EXA 0172
      GO TO 80                                                          EXA 0173
   75 ABSC(M,N)=DMEABS(N)                                               EXA 0174
      EXTC(M,N)=DMEEXT(N)                                               EXA 0175
      ASYM(M,N)=DMESYM(N)                                               EXA 0176
      IF(N.EQ.1) AWCCON(M)=MDLWC                                        EXA 0177
   80 CONTINUE                                                          EXA 0178
      GO TO 185                                                         EXA 0179
100   CONTINUE                                                          EXA 0180
CCC                                                                     EXA 0181
CCC       SECTION TO LOAD EXTINCTION AND ABSORPTION COEFFICIENTS        EXA 0182
CCC       FOR CLOUD AND OR RAIN MODELS                                  EXA 0183
CCC                                                                     EXA 0184
      DO 150 N=1,NE                                                     EXA 0185
      ABSC(M,N)=0.0                                                     EXA 0186
      EXTC(M,N)=0.0                                                     EXA 0187
      ASYM(M,N)=0.0                                                     EXA 0188
      IC=ICLD                                                           EXA 0189
      GO TO (125,130,135,140,145,135,145,145,125,125,125), IC           EXA 0190
125   ABSC(M,N)=CCUABS(N)                                               EXA 0191
      EXTC(M,N)=CCUEXT(N)                                               EXA 0192
      ASYM(M,N)=CCUSYM(N)                                               EXA 0193
      IF(N.EQ.1) AWCCON(M)=CULWC                                        EXA 0194
      GO TO 150                                                         EXA 0195
130   ABSC(M,N)=CALABS(N)                                               EXA 0196
      EXTC(M,N)=CALEXT(N)                                               EXA 0197
      ASYM(M,N)=CALSYM(N)                                               EXA 0198
      IF(N.EQ.1) AWCCON(M)=ASLWC                                        EXA 0199
      GO TO 150                                                         EXA 0200
135   ABSC(M,N)=CSTABS(N)                                               EXA 0201
      EXTC(M,N)=CSTEXT(N)                                               EXA 0202
      ASYM(M,N)=CSTSYM(N)                                               EXA 0203
      IF(N.EQ.1) AWCCON(M)=STLWC                                        EXA 0204
      GO TO 150                                                         EXA 0205
140   ABSC(M,N)=CSCABS(N)                                               EXA 0206
      EXTC(M,N)=CSCEXT(N)                                               EXA 0207
      ASYM(M,N)=CSCSYM(N)                                               EXA 0208
      IF(N.EQ.1) AWCCON(M)=SCLWC                                        EXA 0209
      GO TO 150                                                         EXA 0210
145   ABSC(M,N)=CNIABS(N)                                               EXA 0211
      EXTC(M,N)=CNIEXT(N)                                               EXA 0212
      ASYM(M,N)=CNISYM(N)                                               EXA 0213
      IF(N.EQ.1) AWCCON(M)=SNLWC                                        EXA 0214
150   CONTINUE                                                          EXA 0215
185   CONTINUE                                                          EXA 0216
      DO 200 N=1,NWAVLN                                                 EXA 0217
      ABSC(5,N)=0.                                                      EXA 0218
      EXTC(5,N)=0.                                                      EXA 0219
      ASYM(5,N)=0.                                                      EXA 0220
      AWCCON(5)=0.                                                      EXA 0221
      IF(ICLD .EQ. 18) THEN                                             EXA 0222
           EXTC(5,N)= CI64XT(N)                                         EXA 0223
           ABSC(5,N)= CI64AB(N)                                         EXA 0224
           ASYM(5,N)= CI64G(N)                                          EXA 0225
           AWCCON(5)=3.446E-3                                           EXA 0226
      ENDIF                                                             EXA 0227
      IF(ICLD .EQ. 19) THEN                                             EXA 0228
           EXTC(5,N)= CIR4XT(N)                                         EXA 0229
           ABSC(5,N)= CIR4AB(N)                                         EXA 0230
           ASYM(5,N)= CIR4G(N)                                          EXA 0231
           AWCCON(5)=5.811E-2                                           EXA 0232
      ENDIF                                                             EXA 0233
200   CONTINUE                                                          EXA 0234
      RETURN                                                            EXA 0235
C                                                                       EXA 0236
      END                                                               EXA 0237
