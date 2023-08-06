      SUBROUTINE ABCDTA(IV)                                             ABC 0001
C                                                                       ABC 0002
      COMMON /ABC/ FACTOR(3),ANH3(2),ACO2(10),ACO(3),                   ABC 0003
     X             ACH4(4),ANO2(3),AN2O(11),AO2(6),AO3(5),              ABC 0004
     X             ASO2(4),AH2O(14),ANO,                                ABC 0005
     X             AANH3(2),BBNH3(2),CCNH3(2),                          ABC 0006
     X             AACO2(10),BBCO2(10),CCCO2(10),                       ABC 0007
     X             AACO(3),BBCO(3),CCCO(3),                             ABC 0008
     X             AACH4(4),BBCH4(4),CCCH4(4),                          ABC 0009
     X             AANO2(3),BBNO2(3),CCNO2(3),                          ABC 0010
     X             AAN2O(11),BBN2O(11),CCN2O(11),                       ABC 0011
     X             AAO2(6),BBO2(6),CCO2(6),                             ABC 0012
     X             AAO3(5),BBO3(5),CCO3(5),                             ABC 0013
     X             AASO2(4),BBSO2(4),CCSO2(4),                          ABC 0014
     X             AAH2O(14),BBH2O(14),CCH2O(14),                       ABC 0015
     X             AANO     ,BBNO     ,CCNO                             ABC 0016
C                                                                       ABC 0017
      COMMON/AABBCC/AA(11),BB(11),CC(11),IBND(11),QA(11),CPS(11)        ABC 0018
C                                                                       ABC 0019
C    MOL                                                                ABC 0020
C     1    H2O (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0021
C     2    CO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0022
C     3    O3  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0023
C     4    N2O (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0024
C     5    CO  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0025
C     6    CH4 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0026
C     7    O2  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0027
C     8    NO  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0028
C     9    SO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0029
C    10    NO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0030
C    11    NH3 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)                ABC 0031
C                                                                       ABC 0032
C  ---H2O                                                               ABC 0033
      IMOL=1                                                            ABC 0034
      IW=-1                                                             ABC 0035
      IF(IV.GE.     0.AND.IV.LE.   345) IW=17                           ABC 0036
      IF(IV.GE.   350.AND.IV.LE.  1000) IW=18                           ABC 0037
      IF(IV.GE.  1005.AND.IV.LE.  1640) IW=19                           ABC 0038
      IF(IV.GE.  1645.AND.IV.LE.  2530) IW=20                           ABC 0039
      IF(IV.GE.  2535.AND.IV.LE.  3420) IW=21                           ABC 0040
      IF(IV.GE.  3425.AND.IV.LE.  4310) IW=22                           ABC 0041
      IF(IV.GE.  4315.AND.IV.LE.  6150) IW=23                           ABC 0042
      IF(IV.GE.  6155.AND.IV.LE.  8000) IW=24                           ABC 0043
      IF(IV.GE.  8005.AND.IV.LE.  9615) IW=25                           ABC 0044
      IF(IV.GE.  9620.AND.IV.LE. 11540) IW=26                           ABC 0045
      IF(IV.GE. 11545.AND.IV.LE. 13070) IW=27                           ABC 0046
      IF(IV.GE. 13075.AND.IV.LE. 14860) IW=28                           ABC 0047
      IF(IV.GE. 14865.AND.IV.LE. 16045) IW=29                           ABC 0048
      IF(IV.GE. 16340.AND.IV.LE. 17860) IW=30                           ABC 0049
      IBAND=IW - 16                                                     ABC 0050
      IBND(IMOL)=IW                                                     ABC 0051
      IF(IW .GT.  0) THEN                                               ABC 0052
           QA(IMOL)= AH2O(IBAND)                                        ABC 0053
           AA(IMOL) =AAH2O(IBAND)                                       ABC 0054
           BB(IMOL) =BBH2O(IBAND)                                       ABC 0055
           CC(IMOL) =CCH2O(IBAND)                                       ABC 0056
      ENDIF                                                             ABC 0057
C  ---O3                                                                ABC 0058
      IMOL=3                                                            ABC 0059
      IW=-1                                                             ABC 0060
      IF (IV .GE.     0 .AND. IV .LE.   200)  IW=31                     ABC 0061
      IF (IV .GE.   515 .AND. IV .LE.  1275)  IW=32                     ABC 0062
      IF (IV .GE.  1630 .AND. IV .LE.  2295)  IW=33                     ABC 0063
      IF (IV .GE.  2670 .AND. IV .LE.  2845)  IW=34                     ABC 0064
      IF (IV .GE.  2850 .AND. IV .LE.  3260)  IW=35                     ABC 0065
      IBAND     =IW - 30                                                ABC 0066
      IBND(IMOL)=IW                                                     ABC 0067
      IF(IW .GT.  0) THEN                                               ABC 0068
           QA(IMOL)=AO3(IBAND)                                          ABC 0069
           AA(IMOL)=AAO3(IBAND)                                         ABC 0070
           BB(IMOL)=BBO3(IBAND)                                         ABC 0071
           CC(IMOL)=CCO3(IBAND)                                         ABC 0072
      ENDIF                                                             ABC 0073
C  ---CO2                                                               ABC 0074
      IMOL=2                                                            ABC 0075
      IW=-1                                                             ABC 0076
      IF (IV .GE.   425 .AND. IV .LE.   835)  IW=36                     ABC 0077
      IF (IV .GE.   840 .AND. IV .LE.  1440)  IW=37                     ABC 0078
      IF (IV .GE.  1805 .AND. IV .LE.  2855)  IW=38                     ABC 0079
      IF (IV .GE.  3070 .AND. IV .LE.  3755)  IW=39                     ABC 0080
      IF (IV .GE.  3760 .AND. IV .LE.  4065)  IW=40                     ABC 0081
      IF (IV .GE.  4530 .AND. IV .LE.  5380)  IW=41                     ABC 0082
      IF (IV .GE.  5905 .AND. IV .LE.  7025)  IW=42                     ABC 0083
      IF((IV .GE.  7395 .AND. IV .LE.  7785) .OR.                       ABC 0084
     *   (IV .GE.  8030 .AND. IV .LE.  8335) .OR.                       ABC 0085
     *   (IV .GE.  9340 .AND. IV .LE.  9670)) IW=43                     ABC 0086
      IBAND=IW - 35                                                     ABC 0087
      IBND(IMOL)=IW                                                     ABC 0088
      IF(IW .GT.  0) THEN                                               ABC 0089
           QA(IMOL)=ACO2(IBAND)                                         ABC 0090
           AA(IMOL)=AACO2(IBAND)                                        ABC 0091
           BB(IMOL)=BBCO2(IBAND)                                        ABC 0092
           CC(IMOL)=CCCO2(IBAND)                                        ABC 0093
      ENDIF                                                             ABC 0094
C  ---CO                                                                ABC 0095
      IMOL=5                                                            ABC 0096
      IW=-1                                                             ABC 0097
      IF (IV .GE.     0 .AND. IV .LE.   175) IW=44                      ABC 0098
      IF((IV .GE.  1940 .AND. IV .LE.  2285) .OR.                       ABC 0099
     *   (IV .GE.  4040 .AND. IV .LE.  4370)) IW=45                     ABC 0100
      IBAND=IW - 43                                                     ABC 0101
      IBND(IMOL)=IW                                                     ABC 0102
      IF(IW .GT.  0) THEN                                               ABC 0103
           QA(IMOL)=ACO(IBAND)                                          ABC 0104
           AA(IMOL)=AACO(IBAND)                                         ABC 0105
           BB(IMOL)=BBCO(IBAND)                                         ABC 0106
           CC(IMOL)=CCCO(IBAND)                                         ABC 0107
      ENDIF                                                             ABC 0108
C  ---CH4                                                               ABC 0109
      IMOL=6                                                            ABC 0110
      IW=-1                                                             ABC 0111
      IF((IV .GE.  1065 .AND. IV .LE.  1775) .OR.                       ABC 0112
     *   (IV .GE.  2345 .AND. IV .LE.  3230) .OR.                       ABC 0113
     *   (IV .GE.  4110 .AND. IV .LE.  4690) .OR.                       ABC 0114
     *   (IV .GE.  5865 .AND. IV .LE.  6135))IW=46                      ABC 0115
      IBAND=IW - 45                                                     ABC 0116
      IBND(IMOL)=IW                                                     ABC 0117
      IF(IW .GT.  0) THEN                                               ABC 0118
           QA(IMOL)=ACH4(IBAND)                                         ABC 0119
           AA(IMOL)=AACH4(IBAND)                                        ABC 0120
           BB(IMOL)=BBCH4(IBAND)                                        ABC 0121
           CC(IMOL)=CCCH4(IBAND)                                        ABC 0122
      ENDIF                                                             ABC 0123
C  ---N2O                                                               ABC 0124
      IMOL=4                                                            ABC 0125
      IW=-1                                                             ABC 0126
      IF (IV .GE.     0 .AND. IV .LE.   120)  IW=47                     ABC 0127
      IF((IV .GE.   490 .AND. IV .LE.   775) .OR.                       ABC 0128
     *   (IV .GE.   865 .AND. IV .LE.   995) .OR.                       ABC 0129
     *   (IV .GE.  1065 .AND. IV .LE.  1385) .OR.                       ABC 0130
     *   (IV .GE.  1545 .AND. IV .LE.  2040) .OR.                       ABC 0131
     *   (IV .GE.  2090 .AND. IV .LE.  2655)) IW=48                     ABC 0132
      IF((IV .GE.  2705 .AND. IV .LE.  2865) .OR.                       ABC 0133
     *   (IV .GE.  3245 .AND. IV .LE.  3925) .OR.                       ABC 0134
     *   (IV .GE.  4260 .AND. IV .LE.  4470) .OR.                       ABC 0135
     *   (IV .GE.  4540 .AND. IV .LE.  4785) .OR.                       ABC 0136
     *   (IV .GE.  4910 .AND. IV .LE.  5165)) IW=49                     ABC 0137
      IBAND=IW - 46                                                     ABC 0138
      IBND(IMOL)=IW                                                     ABC 0139
C                                                                       ABC 0140
      IF(IW .EQ. 49)IBAND=7                                             ABC 0141
C                                                                       ABC 0142
C     THIS CORRECTION IS ONLY FOR N2O AS CURRENTLY WRITTEN              ABC 0143
C                                                                       ABC 0144
      IF(IW .GT.  0) THEN                                               ABC 0145
           QA(IMOL)=AN2O(IBAND)                                         ABC 0146
           AA(IMOL)=AAN2O(IBAND)                                        ABC 0147
           BB(IMOL)=BBN2O(IBAND)                                        ABC 0148
           CC(IMOL)=CCN2O(IBAND)                                        ABC 0149
      ENDIF                                                             ABC 0150
C  ---O2                                                                ABC 0151
      IMOL=7                                                            ABC 0152
      IW=-1                                                             ABC 0153
      IF (IV .GE.     0 .AND. IV .LE.   265)  IW=50                     ABC 0154
      IF((IV .GE.  7650 .AND. IV .LE.  8080) .OR.                       ABC 0155
     *   (IV .GE.  9235 .AND. IV .LE.  9490) .OR.                       ABC 0156
     *   (IV .GE. 12850 .AND. IV .LE. 13220) .OR.                       ABC 0157
     *   (IV .GE. 14300 .AND. IV .LE. 14600) .OR.                       ABC 0158
     *   (IV .GE. 15695 .AND. IV .LE. 15955)) IW=51                     ABC 0159
       IF(IV .GE. 49600 .AND. IV. LE. 52710)  IW=51                     ABC 0160
      IBAND=IW - 49                                                     ABC 0161
      IBND(IMOL)=IW                                                     ABC 0162
      IF(IW .GT.  0) THEN                                               ABC 0163
           QA(IMOL)=AO2(IBAND)                                          ABC 0164
           IF(IV .GE. 49600 .AND. IV. LE. 52710)  QA(IMOL) =.4704       ABC 0165
           AA(IMOL)=AAO2(IBAND)                                         ABC 0166
           BB(IMOL)=BBO2(IBAND)                                         ABC 0167
           CC(IMOL)=CCO2(IBAND)                                         ABC 0168
      ENDIF                                                             ABC 0169
C  ---NH3                                                               ABC 0170
      IMOL=11                                                           ABC 0171
      IW=-1                                                             ABC 0172
      IF (IV .GE.     0 .AND. IV .LE.   385)  IW=52                     ABC 0173
      IF (IV .GE.   390 .AND. IV .LE.  2150)  IW=53                     ABC 0174
      IBAND=IW - 51                                                     ABC 0175
      IBND(IMOL)=IW                                                     ABC 0176
      IF(IW .GT.  0) THEN                                               ABC 0177
           QA(IMOL)=ANH3(IBAND)                                         ABC 0178
           AA(IMOL)=AANH3(IBAND)                                        ABC 0179
           BB(IMOL)=BBNH3(IBAND)                                        ABC 0180
           CC(IMOL)=CCNH3(IBAND)                                        ABC 0181
      ENDIF                                                             ABC 0182
C  ---NO                                                                ABC 0183
      IMOL=8                                                            ABC 0184
      IW=-1                                                             ABC 0185
      IF (IV .GE.  1700 .AND. IV .LE.  2005) IW =54                     ABC 0186
      IBAND=IW - 53                                                     ABC 0187
      IBND(IMOL)=IW                                                     ABC 0188
      IF(IW .GT.  0) THEN                                               ABC 0189
           QA(IMOL)=ANO                                                 ABC 0190
           AA(IMOL)=AANO                                                ABC 0191
           BB(IMOL)=BBNO                                                ABC 0192
           CC(IMOL)=CCNO                                                ABC 0193
      ENDIF                                                             ABC 0194
C  ---NO2                                                               ABC 0195
      IW=-1                                                             ABC 0196
      IMOL=10                                                           ABC 0197
      IF((IV .GE.   580 .AND. IV .LE.   925) .OR.                       ABC 0198
     *   (IV .GE.  1515 .AND. IV .LE.  1695) .OR.                       ABC 0199
     *   (IV .GE.  2800 .AND. IV .LE.  2970)) IW=55                     ABC 0200
      IBAND=IW - 54                                                     ABC 0201
      IBND(IMOL)=IW                                                     ABC 0202
      IF(IW .GT.  0) THEN                                               ABC 0203
           QA(IMOL)=ANO2(IBAND)                                         ABC 0204
           AA(IMOL)=AANO2(IBAND)                                        ABC 0205
           BB(IMOL)=BBNO2(IBAND)                                        ABC 0206
           CC(IMOL)=CCNO2(IBAND)                                        ABC 0207
      ENDIF                                                             ABC 0208
C  ---SO2                                                               ABC 0209
      IMOL=9                                                            ABC 0210
      IW=-1                                                             ABC 0211
      IF (IV .GE.     0 .AND. IV .LE.   185)  IW=56                     ABC 0212
      IF((IV .GE.   400 .AND. IV .LE.   650) .OR.                       ABC 0213
     *   (IV .GE.   950 .AND. IV .LE.  1460) .OR.                       ABC 0214
     *   (IV .GE.  2415 .AND. IV .LE.  2580)) IW=57                     ABC 0215
      IBAND=IW - 55                                                     ABC 0216
      IBND(IMOL)=IW                                                     ABC 0217
      IF(IW .GT.  0) THEN                                               ABC 0218
           QA(IMOL)=ASO2(IBAND)                                         ABC 0219
           AA(IMOL)=AASO2(IBAND)                                        ABC 0220
           BB(IMOL)=BBSO2(IBAND)                                        ABC 0221
           CC(IMOL)=CCSO2(IBAND)                                        ABC 0222
      ENDIF                                                             ABC 0223
      RETURN                                                            ABC 0224
      END                                                               ABC 0225
