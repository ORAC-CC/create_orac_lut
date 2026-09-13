      BLOCK DATA XMLATM                                                 XML 0001
C                                                                       XML 0002
C***********************************************************************XML 0003
C     THIS BLOCK DATA SUBROUTINE INITIALIZES THE STANDARD PROFILES      XML 0004
C     FOR THE "CROSS-SECTION" MOLECULES, THAT IS, THE MOLECULES FOR     XML 0005
C     WHICH THE SPECTRAL DATA IS IN THE FORM OF CROSS-SECTIONS          XML 0006
C     (ABSORPTION COEFFICIENTS) INSTEAD OF LINE PARAMETERS.             XML 0007
C     THE PROFILES OF VOLUME MIXING RATIOS GIVEN HERE ARE FROM:         XML 0008
C                                                                       XML 0009
C     M. ALLEN, SPRING EQUINOX, DIUNRNAL AVERAGE, 1990.                 XML 0010
C     (PRIVATE COMMUNICATION)                                           XML 0011
C***********************************************************************XML 0012
C                                                                       XML 0013
C**   COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES       XML 0014
C**   FOR THE CROSS-SECTION MOLECULES.                                  XML 0015
C                                                                       XML 0016
C**   MOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH   XML 0017
C**   LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL, LAYXMX LEVELS MAX     XML 0018
C                                                                       XML 0019
C                                                                       XML 0020
C                                                                       XML 0021
C     CONVENTION                                                        XML 0022
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         XML 0023
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            XML 0024
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             XML 0025
C                                                                       XML 0026
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              XML 0027
C     NSPEC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL      XML 0028
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     XML 0029
C                                                                       XML 0030
      INCLUDE 'PARAM.LST'                                               XML 0031
C     PARAMETER (MMOLX=18)                                              XML 0032
C     PARAMETER (MMOL=12)                                               XML 0033
C     PARAMETER (MMOLT=MMOL+MMOLX, MMOLT2=2*MMOLT)                      XML 0034
C                                                                       XML 0035
C     PARAMETER(NSPEC=12, NSPECX=13, NSPECT=NSPEC+NSPECX)               XML 0036
C                                                                       XML 0037
      CHARACTER*8 CNAMEX(MMOLX)                                         XML 0038
C                                                                       XML 0039
C                                                                       XML 0040
C                                                                       XML 0041
      REAL AMWTX(MMOLX)                                                 XML 0042
      REAL ALTX(50),                                                    XML 0043
     +     MOLX1(50), MOLX2(50), MOLX3(50), MOLX4(50), MOLX5(50),       XML 0044
     +     MOLX6(50), MOLX7(50), MOLX8(50), MOLX9(50), MOLX10(50),      XML 0045
     +     MOLX11(50),MOLX12(50),MOLX13(50),MOLX14(50),MOLX15(50),      XML 0046
     +     MOLX16(50),MOLX17(50),MOLX18(50),MOLX19(50),MOLX20(50),      XML 0047
     +     MOLX21(50),MOLX22(50),MOLX23(50),MOLX24(50),MOLX25(50),      XML 0048
     +     MOLX26(50),MOLX27(50),MOLX28(50),MOLX29(50),MOLX30(50),      XML 0049
     +     MOLX31(50),MOLX32(50),MOLX33(50),MOLX34(50),MOLX35(50)       XML 0050
C                                                                       XML 0051
      COMMON /NAMEX/CNAMEX                                              XML 0052
      COMMON /ATMWTX/AMWTX                                              XML 0053
      COMMON /MLATMX/ LAYXMX,ALTX,                                      XML 0054
     +     MOLX1,  MOLX2,  MOLX3,  MOLX4,  MOLX5,                       XML 0055
     +     MOLX6,  MOLX7,  MOLX8,  MOLX9,  MOLX10,                      XML 0056
     +     MOLX11, MOLX12, MOLX13, MOLX14, MOLX15,                      XML 0057
     +     MOLX16, MOLX17, MOLX18, MOLX19, MOLX20,                      XML 0058
     +     MOLX21, MOLX22, MOLX23, MOLX24, MOLX25,                      XML 0059
     +     MOLX26, MOLX27, MOLX28, MOLX29, MOLX30,                      XML 0060
     +     MOLX31, MOLX32, MOLX33, MOLX34, MOLX35                       XML 0061
C                                                                       XML 0062
      DATA (CNAMEX(I),I=1,NSPECX)/' CFC-11 ',                           XML 0063
     $     ' CFC-12 ',                                                  XML 0064
     $     ' CFC-13 ',                                                  XML 0065
     $     ' CFC-14 ',                                                  XML 0066
     $     ' CFC-22 ',                                                  XML 0067
     $     ' CFC-113',                                                  XML 0068
     $     ' CFC-114',                                                  XML 0069
     $     ' CFC-115',                                                  XML 0070
     $     ' CLONO2 ',                                                  XML 0071
     $     '  HNO4  ',                                                  XML 0072
     $     ' CHCL2F ',                                                  XML 0073
     $     '  CCL4  ',                                                  XML 0074
     $     '  N2O5  '/                                                  XML 0075
C                                                                       XML 0076
C                                                                       XML 0077
C     DATA AMWTX/CFC-11,CFC-12,CFC-13,CFC-14,CFC-22,CFC-113,CFC-114,CFC-XML 0078
C     CLONO2,HNO4,CHCL2F,CCL4,N2O5/                                     XML 0079
      DATA (AMWTX(I),I=1,NSPECX)                                        XML 0080
     $     /137.3683, 120.9138, 104.4592, 88.00461, 86.4687,            XML 0081
     $     187.3762, 170.9216, 154.467,                                 XML 0082
     $     97.45790, 79.01227, 102.9239, 153.8235,108.0104/             XML 0083
C                                                                       XML 0084
      DATA LAYXMX/50/                                                   XML 0085
C                                                                       XML 0086
      DATA ALTX/                                                        XML 0087
     1     0.0,    1.0,   2.0,   3.0,   4.0,  5.0,   6.0,   7.0,   8.0, XML 0088
     2     9.0,   10.0,  11.0,  12.0,  13.0, 14.0,  15.0,  16.0,  17.0, XML 0089
     3     18.0,  19.0,  20.0,  21.0,  22.0, 23.0,  24.0,  25.0,  27.5, XML 0090
     4     30.0,  32.5,  35.0,  37.5,  40.0, 42.5,  45.0,  47.5,  50.0, XML 0091
     5     55.0,  60.0,  65.0,  70.0,  75.0, 80.0,  85.0,  90.0,  95.0, XML 0092
     6     100.0,105.0, 110.0, 115.0, 120.0/                            XML 0093
C                                                                       XML 0094
C     DATA MOLX1 / CCL3F, AKA CFC-11 /                                  XML 0095
C                                                                       XML 0096
      DATA MOLX1 /                                                      XML 0097
     X     1.400E-04, 1.400E-04, 1.400E-04, 1.400E-04, 1.400E-04,       XML 0098
     X     1.400E-04, 1.400E-04, 1.400E-04, 1.400E-04, 1.398E-04,       XML 0099
     X     1.395E-04, 1.393E-04, 1.390E-04, 1.350E-04, 1.310E-04,       XML 0100
     X     1.270E-04, 1.230E-04, 1.123E-04, 1.016E-04, 9.090E-05,       XML 0101
     X     8.020E-05, 6.865E-05, 5.710E-05, 4.555E-05, 3.400E-05,       XML 0102
     X     2.728E-05, 1.048E-05, 3.902E-06, 6.031E-07, 1.984E-07,       XML 0103
     X     2.351E-08, 1.850E-09, 7.681E-10, 9.182E-11, 2.389E-11,       XML 0104
     X     5.965E-12, 6.108E-13, 4.240E-14, 4.972E-15, 5.090E-16,       XML 0105
     X     4.235E-17, 1.980E-18, 1.862E-19, 1.540E-20, 1.015E-21,       XML 0106
     X     3.110E-23, 2.073E-24, 1.356E-25, 7.233E-27, 5.570E-28/       XML 0107
C                                                                       XML 0108
C     DATA MOLX2 / CCL2F2, AKA CFC-12 /                                 XML 0109
C                                                                       XML 0110
      DATA MOLX2 /                                                      XML 0111
     X     2.400E-04, 2.400E-04, 2.400E-04, 2.400E-04, 2.400E-04,       XML 0112
     X     2.400E-04, 2.400E-04, 2.400E-04, 2.400E-04, 2.398E-04,       XML 0113
     X     2.395E-04, 2.393E-04, 2.390E-04, 2.355E-04, 2.320E-04,       XML 0114
     X     2.285E-04, 2.250E-04, 2.148E-04, 2.045E-04, 1.943E-04,       XML 0115
     X     1.840E-04, 1.710E-04, 1.580E-04, 1.450E-04, 1.320E-04,       XML 0116
     X     1.188E-04, 8.580E-05, 5.885E-05, 3.564E-05, 2.133E-05,       XML 0117
     X     1.203E-05, 6.090E-06, 3.834E-06, 2.133E-06, 1.264E-06,       XML 0118
     X     8.520E-07, 4.160E-07, 1.940E-07, 9.038E-08, 3.905E-08,       XML 0119
     X     1.525E-08, 5.190E-09, 1.829E-09, 5.670E-10, 1.450E-10,       XML 0120
     X     2.480E-11, 4.303E-12, 6.140E-13, 6.613E-14, 9.590E-15/       XML 0121
C                                                                       XML 0122
C     DATA MOLX3 / CCLF3, AKA CFC-13  /                                 XML 0123
C                                                                       XML 0124
      DATA MOLX3 /                                                      XML 0125
     X     50*1.E-12                                           /        XML 0126
C                                                                       XML 0127
C     DATA MOLX4 / CF4, AKA CFC-14 /                                    XML 0128
C                                                                       XML 0129
      DATA MOLX4 /                                                      XML 0130
     X     50*1.E-12                                            /       XML 0131
C                                                                       XML 0132
C     DATA MOLX5 / CHF2CL, AKA CFC-22 /  ????                           XML 0133
C                                                                       XML 0134
      DATA MOLX5 /                                                      XML 0135
     X     6.000E-05, 5.995E-05, 5.990E-05, 5.985E-05, 5.980E-05,       XML 0136
     X     5.978E-05, 5.975E-05, 5.973E-05, 5.970E-05, 5.965E-05,       XML 0137
     X     5.960E-05, 5.955E-05, 5.950E-05, 5.893E-05, 5.835E-05,       XML 0138
     X     5.778E-05, 5.720E-05, 5.563E-05, 5.405E-05, 5.248E-05,       XML 0139
     X     5.090E-05, 4.893E-05, 4.695E-05, 4.498E-05, 4.300E-05,       XML 0140
     X     4.085E-05, 3.548E-05, 3.025E-05, 2.520E-05, 2.070E-05,       XML 0141
     X     1.703E-05, 1.390E-05, 1.196E-05, 1.038E-05, 9.330E-06,       XML 0142
     X     8.780E-06, 8.118E-06, 7.700E-06, 7.423E-06, 7.185E-06,       XML 0143
     X     6.910E-06, 6.520E-06, 5.838E-06, 4.790E-06, 3.360E-06,       XML 0144
     X     1.790E-06, 7.138E-07, 2.235E-07, 5.860E-08, 1.920E-08/       XML 0145
C                                                                       XML 0146
C     DATA MOLX5 / CHCLF2, AKA CFC-22 /  ????                           XML 0147
C     DATA MOLX5 /                                                      XML 0148
C     X     6.000E-05, 5.994E-05, 5.987E-05, 5.982E-05, 5.977E-05,      XML 0149
C     X     5.974E-05, 5.970E-05, 5.968E-05, 5.966E-05, 5.963E-05,      XML 0150
C     X     5.960E-05, 5.955E-05, 5.949E-05, 5.921E-05, 5.893E-05,      XML 0151
C     X     5.808E-05, 5.723E-05, 5.582E-05, 5.441E-05, 5.265E-05,      XML 0152
C     X     5.089E-05, 4.897E-05, 4.705E-05, 4.502E-05, 4.298E-05,      XML 0153
C     X     4.084E-05, 3.548E-05, 3.021E-05, 2.514E-05, 2.062E-05,      XML 0154
C     X     1.686E-05, 1.392E-05, 1.184E-05, 1.036E-05, 9.356E-06,      XML 0155
C     X     8.784E-06, 8.163E-06, 7.741E-06, 7.449E-06, 7.201E-06,      XML 0156
C     X     6.919E-06, 6.524E-06, 5.872E-06, 4.867E-06, 3.396E-06,      XML 0157
C     X     1.808E-06, 6.935E-07, 2.066E-07, 5.485E-08, 1.930E-08/      XML 0158
C                                                                       XML 0159
C     DATA MOLX6 / C2CL3F3, AKA CFC-113 /                               XML 0160
C                                                                       XML 0161
      DATA MOLX6 /                                                      XML 0162
     X     1.900E-05, 1.900E-05, 1.900E-05, 1.900E-05, 1.900E-05,       XML 0163
     X     1.900E-05, 1.900E-05, 1.900E-05, 1.900E-05, 1.898E-05,       XML 0164
     X     1.895E-05, 1.893E-05, 1.890E-05, 1.855E-05, 1.820E-05,       XML 0165
     X     1.785E-05, 1.750E-05, 1.653E-05, 1.555E-05, 1.458E-05,       XML 0166
     X     1.360E-05, 1.238E-05, 1.116E-05, 9.932E-06, 8.710E-06,       XML 0167
     X     7.603E-06, 4.834E-06, 2.915E-06, 1.411E-06, 7.168E-07,       XML 0168
     X     3.205E-07, 1.230E-07, 7.013E-08, 3.215E-08, 1.653E-08,       XML 0169
     X     9.905E-09, 3.920E-09, 1.420E-09, 5.268E-10, 1.757E-10,       XML 0170
     X     5.143E-11, 1.230E-11, 3.263E-12, 7.440E-13, 1.344E-13,       XML 0171
     X     1.310E-14, 1.238E-15, 7.866E-17, 2.326E-18, 7.910E-20/       XML 0172
C                                                                       XML 0173
C     DATA MOLX7 / C2CL2F4, AKA CFC-114 /                               XML 0174
C                                                                       XML 0175
      DATA MOLX7 /                                                      XML 0176
     X     1.200E-05, 1.200E-05, 1.200E-05, 1.200E-05, 1.200E-05,       XML 0177
     X     1.200E-05, 1.200E-05, 1.200E-05, 1.200E-05, 1.200E-05,       XML 0178
     X     1.200E-05, 1.200E-05, 1.200E-05, 1.188E-05, 1.175E-05,       XML 0179
     X     1.163E-05, 1.150E-05, 1.118E-05, 1.085E-05, 1.053E-05,       XML 0180
     X     1.020E-05, 9.755E-06, 9.310E-06, 8.865E-06, 8.420E-06,       XML 0181
     X     7.948E-06, 6.766E-06, 5.635E-06, 4.559E-06, 3.653E-06,       XML 0182
     X     2.926E-06, 2.320E-06, 1.939E-06, 1.618E-06, 1.386E-06,       XML 0183
     X     1.250E-06, 1.040E-06, 8.690E-07, 7.193E-07, 5.855E-07,       XML 0184
     X     4.645E-07, 3.550E-07, 2.535E-07, 1.590E-07, 7.688E-08,       XML 0185
     X     2.130E-08, 3.514E-09, 3.470E-10, 1.612E-11, 8.970E-13/       XML 0186
C                                                                       XML 0187
C     DATA MOLX8 / C2CLF5, AKA CFC-115 /                                XML 0188
C                                                                       XML 0189
      DATA MOLX8 /                                                      XML 0190
     X     50*1.E-12                                            /       XML 0191
C                                                                       XML 0192
C     DATA MOLX9 / CLONO2 /                                             XML 0193
C                                                                       XML 0194
      DATA MOLX9 /                                                      XML 0195
     X     5.750E-06, 5.128E-06, 4.505E-06, 3.883E-06, 3.260E-06,       XML 0196
     X     2.755E-06, 2.250E-06, 1.745E-06, 1.240E-06, 1.208E-06,       XML 0197
     X     1.175E-06, 1.143E-06, 1.110E-06, 8.007E-06, 1.491E-05,       XML 0198
     X     2.180E-05, 2.870E-05, 1.015E-04, 1.744E-04, 2.472E-04,       XML 0199
     X     3.200E-04, 4.600E-04, 6.000E-04, 7.400E-04, 8.800E-04,       XML 0200
     X     9.725E-04, 1.204E-03, 1.140E-03, 9.456E-04, 5.238E-04,       XML 0201
     X     2.379E-04, 4.270E-05, 1.739E-05, 1.680E-06, 3.544E-07,       XML 0202
     X     4.668E-08, 1.194E-09, 1.710E-11, 8.568E-13, 7.235E-14,       XML 0203
     X     5.415E-15, 4.040E-17, 3.793E-20, 1.172E-22, 8.674E-25,       XML 0204
     X     8.080E-28, 8.438E-30, 7.600E-32, 0.       , 0.       /       XML 0205
C                                                                       XML 0206
C     DATA MOLX10 / HNO4 /                                              XML 0207
C                                                                       XML 0208
      DATA MOLX10 /                                                     XML 0209
     X     4.340E-07, 2.136E-06, 3.837E-06, 5.539E-06, 7.240E-06,       XML 0210
     X     1.816E-05, 2.907E-05, 3.999E-05, 5.090E-05, 4.910E-05,       XML 0211
     X     4.730E-05, 4.550E-05, 4.370E-05, 5.655E-05, 6.940E-05,       XML 0212
     X     8.225E-05, 9.510E-05, 1.426E-04, 1.901E-04, 2.375E-04,       XML 0213
     X     2.850E-04, 3.518E-04, 4.185E-04, 4.853E-04, 5.520E-04,       XML 0214
     X     5.808E-04, 6.526E-04, 5.100E-04, 3.187E-04, 1.469E-04,       XML 0215
     X     5.261E-05, 9.960E-06, 4.435E-06, 8.735E-07, 2.573E-07,       XML 0216
     X     7.630E-08, 6.803E-09, 4.790E-10, 7.730E-11, 3.210E-11,       XML 0217
     X     1.266E-11, 8.430E-13, 1.883E-15, 1.734E-18, 4.466E-21,       XML 0218
     X     3.330E-24, 6.306E-26, 1.522E-27, 3.538E-29, 2.240E-30/       XML 0219
C                                                                       XML 0220
C     DATA MOLX11 / CHCL2F, AKA CFC-21/                                 XML 0221
C                                                                       XML 0222
      DATA MOLX11 /                                                     XML 0223
     X     50*1.E-12/                                                   XML 0224
C                                                                       XML 0225
C     DATA MOLX12 / CCL4 /                                              XML 0226
C                                                                       XML 0227
      DATA MOLX12 /                                                     XML 0228
     X     1.300E-04, 1.300E-04, 1.300E-04, 1.300E-04, 1.300E-04,       XML 0229
     X     1.298E-04, 1.295E-04, 1.293E-04, 1.290E-04, 1.290E-04,       XML 0230
     X     1.290E-04, 1.290E-04, 1.290E-04, 1.250E-04, 1.210E-04,       XML 0231
     X     1.170E-04, 1.130E-04, 1.018E-04, 9.060E-05, 7.940E-05,       XML 0232
     X     6.820E-05, 5.715E-05, 4.610E-05, 3.505E-05, 2.400E-05,       XML 0233
     X     1.884E-05, 5.940E-06, 1.755E-06, 1.316E-07, 3.954E-08,       XML 0234
     X     1.715E-09, 3.960E-11, 1.531E-11, 5.572E-13, 1.100E-13,       XML 0235
     X     1.093E-14, 3.775E-16, 5.250E-18, 2.454E-19, 9.875E-21,       XML 0236
     X     3.051E-22, 2.930E-24, 1.112E-25, 3.611E-27, 8.790E-29,       XML 0237
     X     5.100E-31, 1.358E-32, 0.       , 0.       , 0.       /       XML 0238
C                                                                       XML 0239
C     DATA MOLX13 / N2O5 /                                              XML 0240
C                                                                       XML 0241
      DATA MOLX13 /                                                     XML 0242
C                                                                       XML 0243
     X     2.420E-10, 9.540E-10, 1.666E-09, 2.378E-09, 3.090E-09,       XML 0244
     X     3.608E-09, 4.125E-09, 4.643E-09, 5.160E-09, 5.763E-09,       XML 0245
     X     6.365E-09, 6.968E-09, 7.570E-09, 4.407E-07, 8.738E-07,       XML 0246
     X     1.307E-06, 1.740E-06, 8.780E-06, 1.582E-05, 2.286E-05,       XML 0247
     X     2.990E-05, 4.145E-05, 5.300E-05, 6.455E-05, 7.610E-05,       XML 0248
     X     7.703E-05, 7.934E-05, 5.520E-05, 2.730E-05, 1.082E-05,       XML 0249
     X     2.724E-06, 2.130E-07, 8.313E-08, 3.938E-09, 7.580E-10,       XML 0250
     X     6.289E-11, 1.044E-12, 5.120E-15, 1.355E-16, 8.200E-18,       XML 0251
     X     7.258E-19, 5.090E-21, 3.961E-24, 2.392E-27, 8.388E-30,       XML 0252
     X     0.       , 0.       , 0.       , 0.       , 0.       /       XML 0253
C                                                                       XML 0254
C     DATA MOLX14 / HNO3 /                                              XML 0255
C                                                                       XML 0256
      DATA MOLX14 /                                                     XML 0257
     X     5.550E-05, 6.403E-05, 7.255E-05, 8.108E-05, 8.960E-05,       XML 0258
     X     9.320E-05, 9.680E-05, 1.004E-04, 1.040E-04, 1.205E-04,       XML 0259
     X     1.370E-04, 1.535E-04, 1.700E-04, 4.800E-04, 7.900E-04,       XML 0260
     X     1.100E-03, 1.410E-03, 2.163E-03, 2.915E-03, 3.668E-03,       XML 0261
     X     4.420E-03, 5.010E-03, 5.600E-03, 6.190E-03, 6.780E-03,       XML 0262
     X     6.790E-03, 6.815E-03, 5.875E-03, 4.643E-03, 3.205E-03,       XML 0263
     X     1.981E-03, 9.000E-04, 4.356E-04, 1.229E-04, 3.756E-05,       XML 0264
     X     1.170E-05, 1.048E-06, 7.040E-08, 1.099E-08, 4.315E-09,       XML 0265
     X     9.327E-10, 1.320E-10, 1.112E-12, 4.992E-15, 2.180E-17,       XML 0266
     X     1.840E-20, 3.313E-22, 9.205E-24, 5.043E-25, 1.260E-25/       XML 0267
C                                                                       XML 0268
C     DATA MOLX15 / ?????? /                                            XML 0269
C                                                                       XML 0270
      DATA MOLX15 /                                                     XML 0271
     X     50*-99.                                              /       XML 0272
C                                                                       XML 0273
C     DATA MOLX16 / ?????? /                                            XML 0274
C                                                                       XML 0275
      DATA MOLX16 /                                                     XML 0276
     X     50*-99.                                              /       XML 0277
C                                                                       XML 0278
C     DATA MOLX17 / ?????? /                                            XML 0279
C                                                                       XML 0280
      DATA MOLX17 /                                                     XML 0281
     X     50*-99.                                              /       XML 0282
C                                                                       XML 0283
C     DATA MOLX18 / ?????? /                                            XML 0284
C                                                                       XML 0285
      DATA MOLX18 /                                                     XML 0286
     X     50*-99.                                              /       XML 0287
C                                                                       XML 0288
C     DATA MOLX19 / ?????? /                                            XML 0289
C                                                                       XML 0290
      DATA MOLX19 /                                                     XML 0291
     X     50*-99.                                              /       XML 0292
C                                                                       XML 0293
C     DATA MOLX20 / ?????? /                                            XML 0294
C                                                                       XML 0295
      DATA MOLX20 /                                                     XML 0296
     X     50*-99.                                              /       XML 0297
C                                                                       XML 0298
C     DATA MOLX21 / ?????? /                                            XML 0299
C                                                                       XML 0300
      DATA MOLX21 /                                                     XML 0301
     X     50*-99.                                              /       XML 0302
C                                                                       XML 0303
C     DATA MOLX22 / ?????? /                                            XML 0304
C                                                                       XML 0305
      DATA MOLX22 /                                                     XML 0306
     X     50*-99.                                              /       XML 0307
C                                                                       XML 0308
C     DATA MOLX23 / ?????? /                                            XML 0309
C                                                                       XML 0310
      DATA MOLX23 /                                                     XML 0311
     X     50*-99.                                              /       XML 0312
C                                                                       XML 0313
C     DATA MOLX24 / ?????? /                                            XML 0314
C                                                                       XML 0315
      DATA MOLX24 /                                                     XML 0316
     X     50*-99.                                              /       XML 0317
C                                                                       XML 0318
C     DATA MOLX25 / ?????? /                                            XML 0319
C                                                                       XML 0320
      DATA MOLX25 /                                                     XML 0321
     X     50*-99.                                              /       XML 0322
C                                                                       XML 0323
C     DATA MOLX26 / ?????? /                                            XML 0324
C                                                                       XML 0325
      DATA MOLX26 /                                                     XML 0326
     X     50*-99.                                              /       XML 0327
C                                                                       XML 0328
C     DATA MOLX27 / ?????? /                                            XML 0329
C                                                                       XML 0330
      DATA MOLX27 /                                                     XML 0331
     X     50*-99.                                              /       XML 0332
C                                                                       XML 0333
C     DATA MOLX28 / ?????? /                                            XML 0334
C                                                                       XML 0335
      DATA MOLX28 /                                                     XML 0336
     X     50*-99.                                              /       XML 0337
C                                                                       XML 0338
C     DATA MOLX29 / ?????? /                                            XML 0339
C                                                                       XML 0340
      DATA MOLX29 /                                                     XML 0341
     X     50*-99.                                              /       XML 0342
C                                                                       XML 0343
C     DATA MOLX30 / ?????? /                                            XML 0344
C                                                                       XML 0345
      DATA MOLX30 /                                                     XML 0346
     X     50*-99.                                              /       XML 0347
C                                                                       XML 0348
C     DATA MOLX31 / ?????? /                                            XML 0349
C                                                                       XML 0350
      DATA MOLX31 /                                                     XML 0351
     X     50*-99.                                              /       XML 0352
C                                                                       XML 0353
C     DATA MOLX32 / ?????? /                                            XML 0354
C                                                                       XML 0355
      DATA MOLX32 /                                                     XML 0356
     X     50*-99.                                              /       XML 0357
C                                                                       XML 0358
C     DATA MOLX33 / ?????? /                                            XML 0359
C                                                                       XML 0360
      DATA MOLX33 /                                                     XML 0361
     X     50*-99.                                              /       XML 0362
C                                                                       XML 0363
C     DATA MOLX34 / ?????? /                                            XML 0364
C                                                                       XML 0365
      DATA MOLX34 /                                                     XML 0366
     X     50*-99.                                              /       XML 0367
C                                                                       XML 0368
      DATA MOLX35 /                                                     XML 0369
     X     50*-99.                                              /       XML 0370
C                                                                       XML 0371
      END                                                               XML 0372
