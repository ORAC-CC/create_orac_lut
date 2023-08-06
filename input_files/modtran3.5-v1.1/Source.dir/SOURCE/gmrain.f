      FUNCTION   GMRAIN(FREQ,T,RATE)                                    GMR 0001
C                                                                       GMR 0002
C        COMPUTES ATTENUATION OF CONDENSED WATER IN FORM OF RAIN        GMR 0003
C                                                                       GMR 0004
C        FREQ = WAVENUMBER (CM-1)                                       GMR 0005
C        T    = TEMPERATURE (DEGREES KELVIN)                            GMR 0006
C        RATE = PRECIPITATION RATE (MM/HR)                              GMR 0007
C       WVLTH = WAVELENGTH IN CM                                        GMR 0008
C                                                                       GMR 0009
C     TABLES ATTAB AND FACTOR CALCULATED FROM FULL MIE THEORY           GMR 0010
C     UTILIZING MARSHALL-PALMER SIZE DISTRIBUTION WITH RAYS INDEX       GMR 0011
C     OF REFRACTION                                                     GMR 0012
C                                                                       GMR 0013
C     ATTAB IS ATTENUATION DATA TABLE IN NEPERS FOR 20 DEG CELSIUS      GMR 0014
C     WITH RADIATION FIELD REMOVED                                      GMR 0015
C                                                                       GMR 0016
C     WVNTBL IS WAVENUMBER TABLE FOR WAVENUMBERS USED IN TABLE ATTAB    GMR 0017
C     TMPTAB IS INTERPOLATION DATA TABLE FOR TEMPERATURES IN DEG KELVIN GMR 0018
C                                                                       GMR 0019
C     TLMDA IS INTERPOLATION DATA TABLE FOR WAVELENGTH IN CM            GMR 0020
C     TFREQ IS INTERPOLATION DATA TABLE FOR WAVENUMBER IN CM-1          GMR 0021
C                                                                       GMR 0022
C     RATTAB IS RAIN RATE TABLE IN MM/HR                                GMR 0023
C                                                                       GMR 0024
C     FACTOR IS TABLE OF TEMPERATURE CORRECTION FACTORS FOR             GMR 0025
C     TABLE ATTAB FOR REPRESENTATIVE RAINS WITHOUT RADIATION FIELD      GMR 0026
C                                                                       GMR 0027
C                                                                       GMR 0028
C     AITKEN INTERPOLATION SCHEME WRITTEN BY                            GMR 0029
C           E.T. FLORANCE O.N.R. PASADENA CA.                           GMR 0030
C                                                                       GMR 0031
C                                                                       GMR 0032
      DIMENSION ATTAB1(35),ATTAB2(35),ATTAB3(35),ATTAB4(35),ATTAB5(35)  GMR 0033
      DIMENSION ATTAB6(35),ATTAB7(35),ATTAB8(35),ATTAB9(35)             GMR 0034
      DIMENSION ATTAB(35,9),WVLTAB(27),RATTAB(9),FACTOR(5,8,5)          GMR 0035
      DIMENSION X(4),Y(4),ATTN(4),RATES(4)                              GMR 0036
CCC   DIMENSION X(3),Y(3),ATTN(3),RATES(3)                              GMR 0037
      DIMENSION TMPTAB(5),TLMDA(6),FACIT(5),TFACT(5)                    GMR 0038
      DIMENSION TFREQ(8),WVNTBL(35)                                     GMR 0039
      DIMENSION FACEQ1(5,8),FACEQ2(5,8),FACEQ3(5,8),FACEQ4(5,8)         GMR 0040
      DIMENSION FACEQ5(5,8)                                             GMR 0041
      EQUIVALENCE (ATTAB1(1),ATTAB(1,1)),(ATTAB2(1),ATTAB(1,2))         GMR 0042
      EQUIVALENCE (ATTAB3(1),ATTAB(1,3)),(ATTAB4(1),ATTAB(1,4))         GMR 0043
      EQUIVALENCE (ATTAB5(1),ATTAB(1,5)),(ATTAB6(1),ATTAB(1,6))         GMR 0044
      EQUIVALENCE (ATTAB7(1),ATTAB(1,7)),(ATTAB8(1),ATTAB(1,8))         GMR 0045
      EQUIVALENCE (ATTAB9(1),ATTAB(1,9))                                GMR 0046
      EQUIVALENCE (FACEQ1(1,1),FACTOR(1,1,1))                           GMR 0047
      EQUIVALENCE (FACEQ2(1,1),FACTOR(1,1,2))                           GMR 0048
      EQUIVALENCE (FACEQ3(1,1),FACTOR(1,1,3))                           GMR 0049
      EQUIVALENCE (FACEQ4(1,1),FACTOR(1,1,4))                           GMR 0050
      EQUIVALENCE (FACEQ5(1,1),FACTOR(1,1,5))                           GMR 0051
      DATA WVLTAB/.03,.033,.0375,.043,.05,.06,.075,.1,.15,.2,.25,.3,.5, GMR 0052
     1.8,1.,2.,3.,4.,5.,5.5,6.,6.5,7.,8.,9.,10.,15./                    GMR 0053
      DATA WVNTBL/ 0.0000,                                              GMR 0054
     1    .0667,.1000,.1111,.1250,.1429,.1538,                          GMR 0055
     2  .1667,.1818,.2000,.2500,.3333,.5000,1.0000,                     GMR 0056
     3 1.2500,2.0000,3.3333,4.0000,5.0000,6.6667,10.0000,               GMR 0057
     4 13.3333,16.6667,20.0000,23.2558,26.6667,30.3030,33.3333,         GMR 0058
     5 50.0,80.0,120.0,180.0,250.0,300.0,350.0/                         GMR 0059
      DATA RATTAB /.25,1.25,2.5,5.,12.5,25.,50.,100.,150./              GMR 0060
      DATA TLMDA/.03,.1,.5,1.25,3.2,10./                                GMR 0061
      DATA TFREQ/0.0,0.1,0.3125,0.8,2.0,10.0,33.3333,350.0/             GMR 0062
      DATA TMPTAB/273.15,283.15,293.15,303.15,313.15/                   GMR 0063
      DATA ATTAB1/                                                      GMR 0064
     1 1.272E+00,1.332E+00,1.361E+00,1.368E+00,1.393E+00,1.421E+00,     GMR 0065
     2 1.439E+00,1.466E+00,1.499E+00,1.541E+00,1.682E+00,1.951E+00,     GMR 0066
     3 2.571E+00,3.575E+00,3.808E+00,4.199E+00,3.665E+00,3.161E+00,     GMR 0067
     4 2.462E+00,1.632E+00,8.203E-01,4.747E-01,3.052E-01,2.113E-01,     GMR 0068
     5 1.551E-01,1.168E-01,8.958E-02,7.338E-02,3.174E-02,1.178E-02,     GMR 0069
     6 5.016E-03,2.116E-03,1.123E-03,8.113E-04,6.260E-04/               GMR 0070
      DATA ATTAB2/                                                      GMR 0071
     1 4.915E+00,5.257E+00,5.518E+00,5.632E+00,5.807E+00,6.069E+00,     GMR 0072
     2 6.224E+00,6.452E+00,6.756E+00,7.132E+00,8.453E+00,1.132E+01,     GMR 0073
     3 1.685E+01,2.177E+01,2.246E+01,2.156E+01,1.470E+01,1.167E+01,     GMR 0074
     4 8.333E+00,5.089E+00,2.356E+00,1.320E+00,8.315E-01,5.705E-01,     GMR 0075
     5 4.151E-01,3.119E-01,2.385E-01,1.955E-01,8.373E-02,3.138E-02,     GMR 0076
     6 1.351E-02,5.789E-03,3.090E-03,2.236E-03,1.725E-03/               GMR 0077
      DATA ATTAB3/                                                      GMR 0078
     1 8.798E+00,9.586E+00,1.023E+01,1.049E+01,1.093E+01,1.159E+01,     GMR 0079
     2 1.205E+01,1.263E+01,1.343E+01,1.450E+01,1.832E+01,2.627E+01,     GMR 0080
     3 3.904E+01,4.664E+01,4.702E+01,4.152E+01,2.542E+01,1.959E+01,     GMR 0081
     4 1.363E+01,8.087E+00,3.660E+00,2.028E+00,1.274E+00,8.710E-01,     GMR 0082
     5 6.340E-01,4.757E-01,3.634E-01,2.971E-01,1.275E-01,4.795E-02,     GMR 0083
     6 2.072E-02,8.936E-03,4.780E-03,3.460E-03,2.670E-03/               GMR 0084
      DATA ATTAB4/                                                      GMR 0085
     1 1.575E+01,1.750E+01,1.914E+01,1.991E+01,2.108E+01,2.276E+01,     GMR 0086
     2 2.399E+01,2.561E+01,2.785E+01,3.097E+01,4.204E+01,6.334E+01,     GMR 0087
     3 8.971E+01,9.853E+01,9.609E+01,7.718E+01,4.290E+01,3.220E+01,     GMR 0088
     4 2.188E+01,1.271E+01,5.641E+00,3.110E+00,1.947E+00,1.327E+00,     GMR 0089
     5 9.657E-01,7.242E-01,5.539E-01,4.528E-01,1.942E-01,7.335E-02,     GMR 0090
     6 3.181E-02,1.380E-02,7.394E-03,5.354E-03,4.132E-03/               GMR 0091
      DATA ATTAB5/                                                      GMR 0092
     1 3.400E+01,3.927E+01,4.523E+01,4.796E+01,5.207E+01,5.886E+01,     GMR 0093
     2 6.383E+01,7.060E+01,8.005E+01,9.360E+01,1.381E+02,2.069E+02,     GMR 0094
     3 2.620E+02,2.534E+02,2.366E+02,1.673E+02,8.285E+01,6.059E+01,     GMR 0095
     4 4.013E+01,2.280E+01,9.939E+00,5.439E+00,3.400E+00,2.315E+00,     GMR 0096
     5 1.685E+00,1.263E+00,9.664E-01,7.914E-01,3.397E-01,1.288E-01,     GMR 0097
     6 5.611E-02,2.450E-02,1.316E-02,9.536E-03,7.360E-03/               GMR 0098
      DATA ATTAB6/                                                      GMR 0099
     1 6.087E+01,7.347E+01,8.886E+01,9.653E+01,1.081E+02,1.283E+02,     GMR 0100
     2 1.435E+02,1.649E+02,1.947E+02,2.346E+02,3.543E+02,4.991E+02,     GMR 0101
     3 5.705E+02,5.048E+02,4.510E+02,2.900E+02,1.335E+02,9.607E+01,     GMR 0102
     4 6.269E+01,3.520E+01,1.519E+01,8.295E+00,5.182E+00,3.529E+00,     GMR 0103
     5 2.569E+00,1.927E+00,1.474E+00,1.208E+00,5.191E-01,1.975E-01,     GMR 0104
     6 8.627E-02,3.784E-02,2.037E-02,1.476E-02,1.139E-02/               GMR 0105
      DATA ATTAB7/                                                      GMR 0106
     1 1.090E+02,1.396E+02,1.811E+02,2.029E+02,2.396E+02,3.039E+02,     GMR 0107
     2 3.536E+02,4.189E+02,5.081E+02,6.217E+02,9.038E+02,1.165E+03,     GMR 0108
     3 1.212E+03,9.731E+02,8.330E+02,4.901E+02,2.123E+02,1.507E+02,     GMR 0109
     4 9.718E+01,5.408E+01,2.316E+01,1.264E+01,7.896E+00,5.377E+00,     GMR 0110
     5 3.915E+00,2.939E+00,2.249E+00,1.844E+00,7.940E-01,3.029E-01,     GMR 0111
     6 1.327E-01,5.846E-02,3.151E-02,2.284E-02,1.763E-02/               GMR 0112
      DATA ATTAB8/                                                      GMR 0113
     1 1.950E+02,2.703E+02,3.904E+02,4.614E+02,5.825E+02,7.909E+02,     GMR 0114
     2 9.475E+02,1.142E+03,1.380E+03,1.656E+03,2.237E+03,2.610E+03,     GMR 0115
     3 2.500E+03,1.820E+03,1.491E+03,8.103E+02,3.336E+02,2.344E+02,     GMR 0116
     4 1.495E+02,8.273E+01,3.524E+01,1.922E+01,1.203E+01,8.182E+00,     GMR 0117
     5 5.961E+00,4.477E+00,3.429E+00,2.812E+00,1.216E+00,4.651E-01,     GMR 0118
     6 2.043E-01,9.033E-02,4.874E-02,3.534E-02,2.728E-02/               GMR 0119
      DATA ATTAB9/                                                      GMR 0120
     1 2.742E+02,4.012E+02,6.353E+02,7.829E+02,1.027E+03,1.439E+03,     GMR 0121
     2 1.725E+03,2.071E+03,2.475E+03,2.909E+03,3.738E+03,4.104E+03,     GMR 0122
     3 3.776E+03,2.589E+03,2.070E+03,1.078E+03,4.326E+02,3.023E+02,     GMR 0123
     4 1.918E+02,1.059E+02,4.499E+01,2.454E+01,1.539E+01,1.045E+01,     GMR 0124
     5 7.615E+00,5.722E+00,4.384E+00,3.596E+00,1.561E+00,5.978E-01,     GMR 0125
     6 2.630E-01,1.165E-01,6.292E-02,4.562E-02,3.522E-02/               GMR 0126
      DATA FACEQ1/                                                      GMR 0127
     1 1.606,1.252,1.000, .816, .680,1.603,1.246,1.000, .817, .684,     GMR 0128
     2 1.444,1.207,1.000, .838, .694,1.016, .985,1.000,1.034,1.058,     GMR 0129
     3  .950, .976,1.000,1.034,1.068, .922, .956,1.000,1.044,1.090,     GMR 0130
     4  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/     GMR 0131
      DATA FACEQ2/                                                      GMR 0132
     1 1.606,1.252,1.000, .816, .680,1.612,1.256,1.000, .817, .684,     GMR 0133
     2 1.193,1.101,1.000, .889, .769, .885, .927,1.000,1.086,1.175,     GMR 0134
     3  .941, .976,1.000,1.024,1.047, .932, .966,1.000,1.034,1.079,     GMR 0135
     4  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/     GMR 0136
      DATA FACEQ3/                                                      GMR 0137
     1 1.606,1.252,1.000, .816, .680,1.621,1.256,1.000, .817, .673,     GMR 0138
     2  .969, .995,1.000, .982, .940, .895, .937,1.000,1.075,1.143,     GMR 0139
     3  .950, .976,1.000,1.024,1.036, .932, .966,1.000,1.034,1.079,     GMR 0140
     4  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/     GMR 0141
      DATA FACEQ4/                                                      GMR 0142
     1 1.606,1.252,1.000, .816, .680,1.631,1.265,1.000, .807, .662,     GMR 0143
     2  .848, .927,1.000,1.044,1.079, .922, .956,1.000,1.055,1.111,     GMR 0144
     3  .950, .976,1.000,1.013,1.036, .932, .966,1.000,1.034,1.079,     GMR 0145
     4  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/     GMR 0146
      DATA FACEQ5/                                                      GMR 0147
     1 1.606,1.252,1.000, .816, .680,1.603,1.265,1.000, .807, .662,     GMR 0148
     2  .820, .918,1.000,1.075,1.132, .941, .966,1.000,1.034,1.079,     GMR 0149
     3  .960, .976,1.000,1.013,1.036, .932, .966,1.000,1.034,1.079,     GMR 0150
     4  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/     GMR 0151
      DATA RATLIM /.05/                                                 GMR 0152
C         GIVE ZERO ATTN IF RATE FALLS BELOW LIMIT                      GMR 0153
      IF(RATE.GT.RATLIM) GO TO 12                                       GMR 0154
      GMRAIN = 0.                                                       GMR 0155
      RETURN                                                            GMR 0156
12    CONTINUE                                                          GMR 0157
CC 12 WVLTH =  1.0  /FREQ                                               GMR 0158
CCC   JMAX=3                                                            GMR 0159
      JMAX=4                                                            GMR 0160
CCC   IF(WVLTH.GT.WVLTAB(1)) GO TO      14                              GMR 0161
CCC   ILOW=0                                                            GMR 0162
CCC   JMAX=2                                                            GMR 0163
CCC   GO TO 18                                                          GMR 0164
CCC   THIS DO LOOP IS 2 LESS THAN NO. OF WVLTAB INPUT                   GMR 0165
CCC14 DO 15 I=2,25                                                      GMR 0166
14    DO 15 I=3,33                                                      GMR 0167
CCC   IF(WVLTH.LT.(.5*(WVLTAB(I)+WVLTAB(I+1)))) GO TO 16                GMR 0168
      IF(FREQ.LT.WVNTBL(I)) GO TO 16                                    GMR 0169
   15 CONTINUE                                                          GMR 0170
CCC   SET ILOW EQUAL TO 1 LESS THAN DO MAX                              GMR 0171
CCC   ILOW=24                                                           GMR 0172
      I=34                                                              GMR 0173
CCC   GO TO 18                                                          GMR 0174
CCC16 ILOW = I-2                                                        GMR 0175
16    ILOW=I-3                                                          GMR 0176
   18 CONTINUE                                                          GMR 0177
CCC   DO 190 I=2,7                                                      GMR 0178
      DO 190 K=3,7                                                      GMR 0179
CCC   IF (RATE. LT.(.5*(RATTAB(I)+RATTAB(I+1))))GO TO 195               GMR 0180
      IF(RATE.LT.RATTAB(K)) GO TO 195                                   GMR 0181
  190 CONTINUE                                                          GMR 0182
CCC   KMIN=6                                                            GMR 0183
      K=8                                                               GMR 0184
CCC   GO TO 198                                                         GMR 0185
CC195 KMIN=I-2                                                          GMR 0186
195   KMIN=K-3                                                          GMR 0187
  198 CONTINUE                                                          GMR 0188
      DO 20 J=1,JMAX                                                    GMR 0189
      IJ = ILOW + J                                                     GMR 0190
      X(J) =       WVNTBL(IJ)                                           GMR 0191
   20 CONTINUE                                                          GMR 0192
C        INTERPOLATE                                                    GMR 0193
CCC   Z = -ALOG(FREQ)                                                   GMR 0194
CCC   DO 25 K=1,3                                                       GMR 0195
      DO 25 K=1,4                                                       GMR 0196
      KJ=KMIN+K                                                         GMR 0197
      RATES(K)=RATTAB(KJ)                                               GMR 0198
      DO 24 J=1,JMAX                                                    GMR 0199
      IJ = ILOW + J                                                     GMR 0200
      Y(J)= ALOG(ATTAB(IJ,KJ))                                          GMR 0201
   24 CONTINUE                                                          GMR 0202
      ATTN(K)=EXP(AITK(X,Y,FREQ,JMAX) )                                 GMR 0203
   25 CONTINUE                                                          GMR 0204
C        APPLY TEMPERATURE CORRECTION                                   GMR 0205
      DO 31 I=2,5                                                       GMR 0206
      IF(T.LT.TMPTAB(I)) GO TO 33                                       GMR 0207
   31 CONTINUE                                                          GMR 0208
      ILOW = 4                                                          GMR 0209
      GO TO 35                                                          GMR 0210
   33 ILOW = I-1                                                        GMR 0211
   35 CONTINUE                                                          GMR 0212
      DO 41 J=2,8                                                       GMR 0213
      IF(FREQ.LT.TFREQ(J)) GO TO 43                                     GMR 0214
   41 CONTINUE                                                          GMR 0215
CCC   JLOW IS 2 LESS THAN DO MAX                                        GMR 0216
      JLOW=6                                                            GMR 0217
      GO TO 45                                                          GMR 0218
   43 JLOW = J-2                                                        GMR 0219
   45 CONTINUE                                                          GMR 0220
      DO 50 K=1,2                                                       GMR 0221
      DO 49 J=1,2                                                       GMR 0222
C        INTERPOLATE IN TEMPERATURE                                     GMR 0223
CCC   KJ=(KMIN/2)+K                                                     GMR 0224
      KJ=K+(KMIN+1)/2                                                   GMR 0225
      JI = JLOW + J                                                     GMR 0226
      FAC = ((TMPTAB(ILOW)-T)*FACTOR(ILOW+1,JI,KJ)+(T-TMPTAB(ILOW+1))*  GMR 0227
     1 FACTOR(ILOW,JI,KJ))/(TMPTAB(ILOW)-TMPTAB(ILOW+1))                GMR 0228
      JI = JLOW +3-J                                                    GMR 0229
      FACIT(J) = (TFREQ(JI)-FREQ )*FAC                                  GMR 0230
   49 CONTINUE                                                          GMR 0231
      TFACT(K) = (FACIT(2)-FACIT(1))/(TFREQ(JLOW+1)-TFREQ(JLOW+2))      GMR 0232
   50 CONTINUE                                                          GMR 0233
C        COMPUTE ATTENUATION (DB/KM)                                    GMR 0234
CCC   KJ=2*KMIN/2+1                                                     GMR 0235
      KJ=2*((KMIN+1)/2)+1                                               GMR 0236
CCC   GMRAIN=AITK(RATES,ATTN,RATE,3)*                                   GMR 0237
      GMRAIN=AITK(RATES,ATTN,RATE,4)*                                   GMR 0238
     1((RATE-RATTAB(KJ))*TFACT(2)+(RATTAB(KJ+2)-RATE)*                  GMR 0239
     2TFACT(1))/(RATTAB(KJ+2)-RATTAB(KJ))                               GMR 0240
CCC                                                                     GMR 0241
CCC    APPLY CONVERSION TO NEPERS                                       GMR 0242
CCC                                                                     GMR 0243
      RETURN                                                            GMR 0244
      END                                                               GMR 0245
