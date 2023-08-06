      REAL FUNCTION PHASEF(V,ALT,SANGLE,RH)                             PHF 0001
      INCLUDE 'PARAM.LST'                                               PHF 0002
C                                                                       PHF 0003
C     THIS ROUTINE IS A BIT DIFFERENT FROM AND REPLACES                 PHF 0004
C     THE OLD PHASEF                                                    PHF 0005
C                                                                       PHF 0006
C     RETURNS THE AEROSOL PHASE FUNCTION FROM THE STORED DATA BASE      PHF 0007
C                                                                       PHF 0008
C     THE TRUTH TABLE MNUM(27,26) STORED IN COMMON/MNMPHS/              PHF 0009
C     IN SUBROUTINE PHSDTA IS QUERIED TO DETERMINE THE PROPER PHASE     PHF 0010
C     FUNCTION NEEDED.                                                  PHF 0011
C     THE 27 POSITIONS REPRESENT THE 27 SPECIFIC FREQUENCIES SHOWN IN   PHF 0012
C     DATA STATEMENT WAVE  .2-40 MICRONS.                               PHF 0013
C     THE NUMBERS STORED IN THESE 27 POSITIONS REPRESENT THE CORRECT    PHF 0014
C     PHASE FUNCTIONS CHOSEN FROM THE DATA STATEMENT PHSFNC'S 1-70      PHF 0015
C     POSSIBLE CHOICES.                                                 PHF 0016
C     THE 26 DATA STATEMENTS EACH HAVING 27 FREQUENCIES REPRESENT THE   PHF 0017
C     FOLLOWING 26 MODELS;                                              PHF 0018
C      1=RURAL     0%RH   2=RURAL    70%RH   3=RURAL    80%RH           PHF 0019
C      4=RURAL    99%RH   5=MARITIME  0%RH   6=MARITIME 70%RH           PHF 0020
C      7=MARITIME 80%RH   8=MARITIME 99%RH   9=URBAN     0%RH           PHF 0021
C     10=URBAN    70%RH  11=URBAN    80%RH  12=URBAN    99%RH           PHF 0022
C     13=OCEANIC   0%RH  14=OCEANIC  70%RH  15=OCEANIC  80%RH           PHF 0023
C     16=OCEANIC  99%RH  17=TROPOSPH  0%RH  18=TROPOSPH 70%RH           PHF 0024
C     19=TROPOSPH 80%RH  20=TROPOSPH 99%RH  21=STRATOSPHERIC            PHF 0025
C     22=AGED VOLCANIC   23=FRESH VOLCANIC  24=RADIATION FOG            PHF 0026
C     25=ADVECTIVE FOG   26=METEORIC DUST                               PHF 0027
C                                                                       PHF 0028
C     IN THE PRESENT VERSION THE 4 OCEANIC MODELS 13-16                 PHF 0029
C     ARE NOT UTILIZED.                                                 PHF 0030
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     PHF 0031
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                PHF 0032
C                                                                       PHF 0033
C     FILE UNIT NUMBERS                                                 PHF 0034
C       IRD      MODTRAN INPUT FILE, tape5, UNIT NUMBER (1).            PHF 0035
C       IPR      MODTRAN STANDARD OUTPUT FILE, tape6, UNIT NUMBER (2).  PHF 0036
C       IPU      MODTRAN SPECTRAL DATA FILE, tape7, UNIT NUMBER (7).    PHF 0037
C       NPR      PRINTOUT LEVEL SWITCH (1=small,0=normal,-1=large).     PHF 0038
C       IPR1     MODTRAN FLUX OUTPUT FILE, tape8, UNIT NUMBER (8).      PHF 0039
C       ISCRCH   MULTIPLE SCATTERING SCRATCH FILE UNIT NUMBER (10).     PHF 0040
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               PHF 0041
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           PHF 0042
C                                                                       PHF 0043
C       PI       THE CONSTANT PI                                        PHF 0044
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       PHF 0045
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       PHF 0046
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         PHF 0047
      REAL PI,DEG,BIGNUM,BIGEXP                                         PHF 0048
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                PHF 0049
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,     PHF 0050
     1  RAINRT                                                          PHF 0051
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)                            PHF 0052
      REAL EXTV,ABSV,ASYV                                               PHF 0053
      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       PHF 0054
      DIMENSION  RHPTS(4),WAVE(27),ANG(34)                              PHF 0055
      DATA ANG/0.,2.,4.,6.,8.,10.,12.,16.,20.,24.,28.,32.,36.,40.,      PHF 0056
     1  50.,60.,70.,80.,90.,100.,110.,120.,125.,130.,135.,140.,145.,    PHF 0057
     2  150.,155.,160.,165.,170.,175.,180./                             PHF 0058
      DATA WAVE/.2,.3,.55,.6943,1.06,1.536,2.0,2.5,2.7,3.,3.2,3.39,5.,  PHF 0059
     1  6.,7.2,7.9,8.7,9.2,10.0,10.59,12.5,15.0,17.2,18.5,21.3,30.,40./ PHF 0060
      DATA RHPTS /0.,70.,80.,99./                                       PHF 0061
      PHASEF=0.                                                         PHF 0062
      ALAM=1.E4/(V+1.E-5)                                               PHF 0063
      IF(SANGLE.LT.0. .OR. SANGLE.GT.180.)GOTO900                       PHF 0064
      IF(ALAM.GT.WAVE(27))THEN                                          PHF 0065
          COSANG=COS(SANGLE/DEG)                                        PHF 0066
          PHASEF=HENGNS(ASYV(1),COSANG)                                 PHF 0067
          RETURN                                                        PHF 0068
      ENDIF                                                             PHF 0069
C                                                                       PHF 0070
C     DETERMINE THE AEROSOL MODEL NUMBER                                PHF 0071
      IF(ALT.GT.ALTB(1))GOTO95                                          PHF 0072
      IF(IHAZE.EQ.0)GOTO400                                             PHF 0073
      IF(IHAZE .EQ.10) THEN                                             PHF 0074
           COSANG=COS(SANGLE/DEG)                                       PHF 0075
           PHASEF=HENGNS(ASYV(1),COSANG)                                PHF 0076
           RETURN                                                       PHF 0077
      ENDIF                                                             PHF 0078
C                                                                       PHF 0079
C     CHECK IF CLOUD,RAIN OR DESERT MODEL IS REQUESTED                  PHF 0080
      IF(IHAZE.EQ.7)THEN                                                PHF 0081
           COSANG=COS(SANGLE/DEG)                                       PHF 0082
           PHASEF=HENGNS(ASYV(1),COSANG)                                PHF 0083
           RETURN                                                       PHF 0084
      ENDIF                                                             PHF 0085
      IF(IHAZE.GE.8)GOTO90                                              PHF 0086
C                                                                       PHF 0087
C     0-2KM BOUNDARY LAYER MODELS, RH DEPENDENT                         PHF 0088
      DO 50 I1=1,4                                                      PHF 0089
          I=I1                                                          PHF 0090
          IF(RHPTS(I).EQ.RH) GOTO70                                     PHF 0091
          IF(RHPTS(I).GT.RH) GOTO60                                     PHF 0092
50    CONTINUE                                                          PHF 0093
60    IRHLO=I-1                                                         PHF 0094
      IRHHI=I                                                           PHF 0095
      GOTO80                                                            PHF 0096
70    IRHLO=I                                                           PHF 0097
      IRHHI=I                                                           PHF 0098
80    CONTINUE                                                          PHF 0099
C                                                                       PHF 0100
C     RURAL MODEL                                                       PHF 0101
      IF(IHAZE.EQ.1 .OR. IHAZE.EQ.2)NN0=0                               PHF 0102
C                                                                       PHF 0103
C     MARITIME MODEL                                                    PHF 0104
      IF(IHAZE.EQ.3 .OR. IHAZE.EQ.4)NN0=4                               PHF 0105
C                                                                       PHF 0106
C     URBAN MODEL                                                       PHF 0107
      IF(IHAZE.EQ.5)NN0=8                                               PHF 0108
C                                                                       PHF 0109
C     TROPOSPHERIC MODEL                                                PHF 0110
      IF(IHAZE.EQ.6)NN0=16                                              PHF 0111
      NN=NN0+IRHLO                                                      PHF 0112
      GOTO130                                                           PHF 0113
C     0-2KM FOG MODELS, NO RH DEPENDENCE                                PHF 0114
90    IF(IHAZE.EQ.8)NN=24                                               PHF 0115
      IF(IHAZE.EQ.9)NN=25                                               PHF 0116
      GOTO130                                                           PHF 0117
95    IF(ALT.GT.ALTB(2)) GOTO110                                        PHF 0118
C     2-10KM TROPOSPHERIC MODEL                                         PHF 0119
      IF(IHAZE.EQ.7)THEN                                                PHF 0120
           COSANG=COS(SANGLE/DEG)                                       PHF 0121
           PHASEF=HENGNS(ASYV(2),COSANG)                                PHF 0122
              RETURN                                                    PHF 0123
        ENDIF                                                           PHF 0124
      NN=18                                                             PHF 0125
      GOTO130                                                           PHF 0126
110   IF(ALT.GT.ALTB(3)) GOTO120                                        PHF 0127
C     10-30KM STRATOSPHERIC MODELS                                      PHF 0128
C     BACKGROUND MODEL                                                  PHF 0129
      IF(IVULCN.EQ.0.OR.IVULCN.EQ.1) NN=21                              PHF 0130
C     AGED VOLCANIC MODEL                                               PHF 0131
      IF(IVULCN.EQ.2.OR.IVULCN.EQ.4) NN=22                              PHF 0132
C     FRESH VOLCANIC                                                    PHF 0133
      IF(IVULCN.EQ.3.OR.IVULCN.EQ.5 .OR.IVULCN.EQ.8) NN=23              PHF 0134
C     BACKGROUND STRATO                                                 PHF 0135
      IF(IVULCN.EQ.6.OR.IVULCN.EQ.7) NN=21                              PHF 0136
      GOTO130                                                           PHF 0137
C     30-100KM METEORIC MODEL                                           PHF 0138
120   NN=26                                                             PHF 0139
130   IRH=0                                                             PHF 0140
C                                                                       PHF 0141
C     DETERMINE THE BOUNDING ANGLE INDICES                              PHF 0142
140   DO 210 I1=1,ML                                                    PHF 0143
      I=I1                                                              PHF 0144
      IF(ANG(I).EQ.SANGLE) GOTO230                                      PHF 0145
      IF(ANG(I).GT.SANGLE) GOTO220                                      PHF 0146
210   CONTINUE                                                          PHF 0147
220   IANG1=I-1                                                         PHF 0148
      IANG2=I                                                           PHF 0149
      GOTO240                                                           PHF 0150
230   IANG1=I                                                           PHF 0151
      IANG2=I                                                           PHF 0152
240   CONTINUE                                                          PHF 0153
C                                                                       PHF 0154
C     DETERMINE THE BOUNDING WAVELENGTH INDICES                         PHF 0155
      DO 250 I1=1,27                                                    PHF 0156
      I=I1                                                              PHF 0157
      IF(WAVE(I).EQ.ALAM) GOTO270                                       PHF 0158
      IF(WAVE(I).GT.ALAM) GOTO260                                       PHF 0159
250   CONTINUE                                                          PHF 0160
260   IWAV1=I-1                                                         PHF 0161
      IWAV2=I                                                           PHF 0162
      IF(IWAV1.LT.1) THEN                                               PHF 0163
          IWAV1=1                                                       PHF 0164
      ENDIF                                                             PHF 0165
      GOTO280                                                           PHF 0166
270   IWAV1=I                                                           PHF 0167
      IWAV2=I                                                           PHF 0168
280   CONTINUE                                                          PHF 0169
C                                                                       PHF 0170
C     FUNCTION PF CHOOSES DESIRED PHASE FUNCTION FROM LOOK UP TABLE     PHF 0171
C     MNUM(IWAV,NN)  WHERE IWAV IS FREQ. AND NN IS MODEL NO.            PHF 0172
C                                                                       PHF 0173
C     WAVELENGTH INTERPOLATION ONLY USES PF11 AND PF21                  PHF 0174
C     ANGLE INTERPOLATION ONLY USES PF11 AND PF12                       PHF 0175
C     WAVELENGTH AND ANGLE INTERPOLATION USES PF11,PF21 AND PF12,PF22.  PHF 0176
C                                                                       PHF 0177
      PF11=PF(NN,IWAV1,IANG1)                                           PHF 0178
      PF21=PF(NN,IWAV2,IANG1)                                           PHF 0179
      PF12=PF(NN,IWAV1,IANG2)                                           PHF 0180
      PF22=PF(NN,IWAV2,IANG2)                                           PHF 0181
C     INTERPOLATE IN WAVELENGTH THEN ANGLE                              PHF 0182
      IF(IWAV1.EQ.IWAV2) GOTO310                                        PHF 0183
      IF(IANG1.EQ.IANG2) GOTO290                                        PHF 0184
C     BOTH INTERPOLATIONS ARE NECESSARY                                 PHF 0185
      CALL INTERP(2,ALAM,WAVE(IWAV1),WAVE(IWAV2),YANG1,                 PHF 0186
     1PF11,PF21)                                                        PHF 0187
      CALL INTERP(2,ALAM,WAVE(IWAV1),WAVE(IWAV2),YANG2,                 PHF 0188
     1PF12,PF22)                                                        PHF 0189
      CALL INTERP(2,SANGLE,ANG(IANG1),ANG(IANG2),Y,YANG1,YANG2)         PHF 0190
      GOTO330                                                           PHF 0191
C     ONLY WAVELENGTH INTERPOLATION IS NECESSARY                        PHF 0192
290   CALL INTERP(2,ALAM,WAVE(IWAV1),WAVE(IWAV2),Y,PF11,                PHF 0193
     1PF21)                                                             PHF 0194
      GOTO330                                                           PHF 0195
310   IF(IANG1.EQ.IANG2) GOTO320                                        PHF 0196
C     ONLY ANGLE INTERPOLATION IS NECESSARY                             PHF 0197
      CALL INTERP(2,SANGLE,ANG(IANG1),ANG(IANG2),Y,PF11,PF12)           PHF 0198
      GOTO330                                                           PHF 0199
C     NO INTERPOLATION IS NECESSARY                                     PHF 0200
320   Y=PF(NN,IWAV1,IANG1)                                              PHF 0201
330   CONTINUE                                                          PHF 0202
      PHASEF=Y                                                          PHF 0203
C                                                                       PHF 0204
C     HUMIDITY DEPENDENCE                                               PHF 0205
      IF(ALT.GT.ALTB(1).OR.NN.GE.17.OR.IRHLO.EQ.IRHHI) GOTO400          PHF 0206
      IF(IRH.EQ.1) GOTO340                                              PHF 0207
      NN=NN0+IRHHI                                                      PHF 0208
      PHFA1=PHASEF                                                      PHF 0209
      IRH=1                                                             PHF 0210
      GOTO280                                                           PHF 0211
340   CONTINUE                                                          PHF 0212
      PHFA2=PHASEF                                                      PHF 0213
      CALL INTERP(1,RH,RHPTS(IRHLO),RHPTS(IRHHI),PHFA,PHFA1,PHFA2)      PHF 0214
      PHASEF=PHFA                                                       PHF 0215
400   CONTINUE                                                          PHF 0216
      RETURN                                                            PHF 0217
  900 WRITE(IPR,901) SANGLE                                             PHF 0218
  901 FORMAT('0FROM PHASEF- SCATTERING ANGLE IS OUT OF RANGE, '         PHF 0219
     1    ,'ANGLE = ',E12.5)                                            PHF 0220
      STOP                                                              PHF 0221
      END                                                               PHF 0222
