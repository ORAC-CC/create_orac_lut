       PROGRAM MDTRN4                                                   MOD 0001
C                                                                       MOD 0002
C     THIS VERSION WAS UPDATED  SEPT 26 1994                            MOD 0003
C     NOW CONTAINS NEW VSA AND WATER VAPOR CONTINUM CORRECTION          MOD 0004
C     CORRECTION TRANSMISSION FOR RADIANCE WITH SOLAR                   MOD 0005
C     IMPROVED RADIANCE ALGORITHIM                                      MOD 0006
C     SPECTRALLY  DEPENDENT ALBEDO                                      MOD 0007
C     USE NEG VALUES OF ALBEDO TO TRIGGER                               MOD 0008
C     COMMON AND DIMENSIONS VIA PARAMETER LIST                          MOD 0009
C     ADDED X_SECTIONS                                                  MOD 0010
C     CORRECT BLOCK DATA                                                MOD 0011
C     RATCSZ     JULY 11                                                MOD 0012
C     MORE X_SECTIONS AUG 1  (13 IN ALL)                                MOD 0013
C     DISORT          AUG19                                             MOD 0014
C     NEW SOLAR       SEPT 21 1994                                      MOD 0015
C     CO2 CONTINUM  NOW PART OF 1995 VERSION OF DIRAC                   MOD 0016
C                                                                       MOD 0017
C--------------------------------------------------------------         MOD 0018
C     MODTRAN IS A COMPUTER CODE DEVELOPED BY SPECTRAL SCIENCES, INC.   MOD 0019
C     AND DESIGNED TO DETERMINE ATMOSPHERIC TRANSMISSION AND RADIANCE   MOD 0020
C     AT MODERATE RESOLUTION (FWHM = 1 CM-1) FROM 0 TO 50,000 CM-1.     MOD 0021
C     MODTRAN IS BASED ON AFGL'S LOWTRAN 7 CODE.  UNLESS THE CODE IS    MOD 0022
C     USED TO RUN LOWTRAN 7, AN UNFORMATTED BAND MODEL TAPE (UNIT=9),   MOD 0023
C     UFTAPE, MUST EXIST.  THE PROGRAM "MAKEUF" CREATES THIS TAPE FROM  MOD 0024
C     THE FORMATTED BAND MODEL TAPE, BMTAPE.                            MOD 0025
C                                                                       MOD 0026
C     MOST OF THE ROUTINES FROM LOWTRAN 7 REMAIN UNCHANGED.  SOME       MOD 0027
C     ROUTINES HAVE BEEN CHANGED ONLY IN THAT A FEW OF THEIR COMMONS    MOD 0028
C     STATEMENTS HAVE BEEN MODIFIED AND THE LOGICAL VARIABLE "MODTRN"   MOD 0029
C     HAS BEEN DECLARED:                                                MOD 0030
C                                                                       MOD 0031
C        ROUTINE    COMMONS                                             MOD 0032
C        -------    -------                                             MOD 0033
C        ABCDTA     "BLANK"                                             MOD 0034
C        AEREXT     "BLANK", CARD1                                      MOD 0035
C        AERNSM     "BLANK", CARD1                                      MOD 0036
C        AERTMP     "BLANK"                                             MOD 0037
C        CIRR18     "BLANK", CARD1, CARD4                               MOD 0038
C        CIRRUS     "BLANK", CARD1, CARD4                               MOD 0039
C        CLDPRF     "BLANK"                                             MOD 0040
C        DESATT     "BLANK"                                             MOD 0041
C        EQULWC     "BLANK"                                             MOD 0042
C        EXABIN     "BLANK", CARD1, CARD4                               MOD 0043
C        FLXADD     "BLANK", CARD1, TRAN                                MOD 0044
C        LAYVSA     CARD1                                               MOD 0045
C        PHASEF     CARD1                                               MOD 0046
C        RDEXA      "BLANK"                                             MOD 0047
C        RDNSM      CARD1                                               MOD 0048
C        DPRFPA (FORMERLY RFPATH)                                       MOD 0049
C                   CARD1, SOLS                                         MOD 0050
C        SSRAD      "BLANK", CARD1, SOLS                                MOD 0051
C        VSANSM     "BLANK", CARD1                                      MOD 0052
C                                                                       MOD 0053
C     THE OTHER ROUTINES WITH MINOR CHANGES ARE:                        MOD 0054
C                                                                       MOD 0055
C        ROUTINE    NUMBER OF ADDED AND SUBTRACTED LINES                MOD 0056
C        -------    ------------------------------------                MOD 0057
C        GEO        108                                                 MOD 0058
C        MSRAD      34                                                  MOD 0059
C        SSGEO      53                                                  MOD 0060
C        STDMDL     39                                                  MOD 0061
C                                                                       MOD 0062
C     THE MOST SIGNIFICANT CHANGES WERE MADE TO THE MAIN ROUTINE AND    MOD 0063
C     THE SUBROUTINE TRANS.  THE MAIN HAS BEEN SPLIT INTO 2 ROUTINES.   MOD 0064
C     THE NEW MAIN ROUTINE JUST CONTAINS COMMENT CARDS AND A CALL TO    MOD 0065
C     ROUTINE DRIVER.  SUBROUTINE DRIVER IS THE DRIVER FOR MODTRAN.     MOD 0066
C     SUBROUTINE TRANS HAS BEEN COMPLETELY REVAMPED WITH ITS NEW FORM   MOD 0067
C     CONTAINING 1012 LINES.  TRANS NOW STEPS IN 1 CM-1 INCREMENTS,     MOD 0068
C     CALLING THE BAND MODEL ROUTINES, INTERPOLATING BETWEEN 5 CM-1     MOD 0069
C     CONTINUUM TRANSMITTANCES, AND PERFORMING THE TRIANGULAR SLIT      MOD 0070
C     CALCULATION.                                                      MOD 0071
C                                                                       MOD 0072
C     A NUMBER OF NEW ROUTINES HAS BEEN ADDED TO ACCOMMODATE THE HIGHER MOD 0073
C     RESOLUTION OPTION:                                                MOD 0074
C                                                                       MOD 0075
C        ROUTINE (# LINES)    PURPOSE                                   MOD 0076
C        -----------------    --------------------------                MOD 0077
C        BMDATA    (135)      MAKE INITIAL CALLS TO BAND MODEL TAPE     MOD 0078
C        BMERFU     (14)      RETURNS EXP(-Y*Y) + Y*SQRT(PI)*ERF(Y) - 1 MOD 0079
C        BMFLUX    (277)      PERFORMS ADDING METHOD CALCULATION        MOD 0080
C        BMLOAD     (91)      LOADS THE BAND MODEL PARAMETERS           MOD 0081
C        BMOD      (287)      PERFORMS CURTIS-GODSON SUMS               MOD 0082
C        BMTRAN     (57)      CALCULATES MOLECULAR TRANSMITTANCE        MOD 0083
C                                                                       MOD 0084
C     TWO CARDS FROM THE LOWTRAN 7 INPUT STREAM HAVE BEEN CHANGED.      MOD 0085
C     LOGICAL VARIABLE, MODTRN, IS NOW READ IN ON CARD 1, AS SHOWN:     MOD 0086
C                                                                       MOD 0087
C     READ(IRD,'(L1,I4,12I5,F8.3,F7.2)')MODTRN,MODEL,ITYPE,IEMSCT,      MOD 0088
C    1  IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB              MOD 0089
C                                                                       MOD 0090
C     IF MODTRN IS FALSE, THE CODE DEFAULTS TO LOWTRAN 7.  WHEN MODTRN  MOD 0091
C     IS TRUE, THE MODERATE RESOLUTION OPTION IS INVOKED.  CARD 4 HAS   MOD 0092
C     ALSO BEEN CHANGED.  INSTEAD OF READING THE REAL VARIABLES, V1, V2 MOD 0093
C     AND DV, THE INTEGERS IV1, IV2, IDV AND IFWHM ARE READ, AS SHOWN:  MOD 0094
C                                                                       MOD 0095
C     READ(IRD,'(4I10)')IV1,IV2,IDV,IFWHM                               MOD 0096
C                                                                       MOD 0097
C     IV1 AND IV2 ARE THE LOWER AND UPPER BOUNDS ON THE FREQUENCY       MOD 0098
C     (CM-1) RANGE, RESPECTIVELY.  IV1 IS A NON-NEGATIVE INTEGER AND    MOD 0099
C     IV2 MUST EXCEED IV1.  IDV, A POSITIVE INTEGER, IS THE STEP SIZE   MOD 0100
C     USED FOR PRINTING OF TRANSMITTANCES AND RADIANCES.  FINALLY,      MOD 0101
C     IFWHM IS THE FULL WIDTH AT HALF MAXIMUM USED BY THE TRIANGULAR    MOD 0102
C     SLIT FUNCTION.  IFWHM MUST BE A POSITIVE INTEGER NOT EXCEEDING 50 MOD 0103
C                                                                       MOD 0104
C***********************************************************************MOD 0105
C     LOWTRAN7  (LAST REVISED JAN  30 1989) REVISION 3.6                MOD 0106
C                                                                       MOD 0107
C               AUTHORS                                                 MOD 0108
C                                                                       MOD 0109
C               F.X.KNEIZYS                                             MOD 0110
C               E. P. SHETTLE                                           MOD 0111
C               G.P. ANDERSON                                           MOD 0112
C               L. W. ABREU                                             MOD 0113
C               J. H. CHETWYND                                          MOD 0114
C               J. E. A. SELBY    (GRUMMAN AEROSPACE)                   MOD 0115
C               S. A. CLOUGH      (AER INC)                             MOD 0116
C               W. O. GALLERY     (OPTIMETRICS)                         MOD 0117
C                                                                       MOD 0118
C   PROGRAM LOWTRAN  CALCULATES THE TRANSMITTANCE AND/OR RADIANCE       MOD 0119
C   OF THE ATMOSPHERE  FROM   0 CM-1 TO 50000 CM-1 (0.20 TO INFINITY    MOD 0120
C   MICRONS) AT 20 CM-1 SPECTRAL RESOLUTION ON A LINEAR                 MOD 0121
C   WAVENUMBER SCALE WITH 5CM-1 SAMPLING                                MOD 0122
C                                                                       MOD 0123
C   LOWTRAN 7 IS A LOW-RESOLUTION PROPAGATION MODEL FOR CALCULATING     MOD 0124
C   ATMOSPHERIC TRANSMITTANCE AND BACKGROUND RADIANCE FROM 0 TO         MOD 0125
C   50,000 CM-1 AT A RESOLUTION OF 20 CM-1 WITH A MINIMUM OF 5 CM-1     MOD 0126
C   SAMPLING.  THE MODEL IS BASED ON THE LOWTRAN 6 (1983) MODEL.        MOD 0127
C   THE PROGRAM CALCULATES SINGLE SCATTERED SOLAR (OR LUNAR)            MOD 0128
C   RADIATION.  MULTIPLE SCATTERED RADIATION HAS BEEN ADDED TO THE      MOD 0129
C   MODEL AS WELL AS NEW MOLECULAR BAND MODEL PARAMETERS AND NEW OR     MOD 0130
C   UPDATED OZONE AND MOLECULAR OXYGEN ABSORPTION PARAMETERS FOR THE    MOD 0131
C   UV.  OTHER MODIFICATIONS INCLUDE A WIND-DEPENDENT DESERT MODEL, NEW MOD 0132
C   CIRRUS CLOUD MODELS, AND NEW CLOUD AND RAIN MODELS.  THE MODEL ALSO MOD 0133
C   INCLUDES NEW REPRESENTATIVE (GEOGRAPHICAL AND SEASONAL) ATMOSPHERIC MOD 0134
C   MODELS AND UPDATED AEROSOL MODELS WITH OPTIONS TO REPLACE THEM WITH MOD 0135
C   USER-DERIVED VALUES.  SIX MODES OF PROGRAM EXECUTION ARE ALLOWED    MOD 0136
C   WITH THE NEW MODEL AND COMPUTER CODE FOR A GIVEN SLANT PATH         MOD 0137
C   UTILIZING SPHERICAL-REFRACTIVE GEOMETRY.  THE ARMY VERTICAL         MOD 0138
C   STRUCTURE ALGORITHM HAS BEEN MODIFIED TO INCLUDE THE NEW PEDESTAL   MOD 0139
C   MODEL BELOW THE CLOUD BASE.  A NEW OPTION HAS BEEN  ADDED TO        MOD 0140
C   MODIFY THE AEROSOL PROFILE, IF THE GROUND IS NOT AT SEA LEVEL.      MOD 0141
C                                                                       MOD 0142
C***********************************************************************MOD 0143
C                                                                       MOD 0144
C     THE FOLLOWING INFORMATION SHOULD BE PROVIDED BY THE USER          MOD 0145
C     AND MAILED TO   L.W ABREU  ,AFGL/OPI,HANSCOM AFB,MASS 01731       MOD 0146
C     THIS WILL BE USED TO UPDATE THE AFGL MAILING LIST                 MOD 0147
C     AND FOR NOTIFICATION TO THE USER OF ERRORS IN THE CODE            MOD 0148
C                                                                       MOD 0149
C                                                                       MOD 0150
C           MY NAME IS                                                  MOD 0151
C           COMPANY                                                     MOD 0152
C           ADDRESS                                                     MOD 0153
C           MY COMPUTER IS                                              MOD 0154
C                                                                       MOD 0155
C                                                                       MOD 0156
C***********************************************************************MOD 0157
C   THE USE OF THE WORD 'CARD' IS EQUIVALENT TO EDITING WITH 80 COLUMNS MOD 0158
C                                                                       MOD 0159
C     PROGRAM ACTIVATED BY SUBMISSION OF A FIVE  (OR MORE)              MOD 0160
C      CARD SEQUENCE AS FOLLOWS                                         MOD 0161
C                                                                       MOD 0162
C     CARD 1    MODEL,ITYPE,IEMSCT,IMULT,M1,M2,M3,                      MOD 0163
C               M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB                     MOD 0164
C                          FORMAT(13I5,F8.3,F7.2)                       MOD 0165
C                                                                       MOD 0166
C     CARD 2    IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT, MOD 0167
C               GNDALT                                                  MOD 0168
C                          FORMAT(6I5,5F10.3)                           MOD 0169
C                                                                       MOD 0170
C               CARD 2A    CTHIK,CALT,CEXT,ISEED       (ICLD=18,19,20)  MOD 0171
C                          FORMAT(3F10.3,I10)                           MOD 0172
C                                                                       MOD 0173
C               CARD 2B    ZCVSA,ZTVSA,ZINVSA    (IVSA=1)               MOD 0174
C                          FORMAT(3F10.3)                               MOD 0175
C                                                                       MOD 0176
C               CARD 2C    ML,IRD1,IRD2,TITLE (MODEL=0 / 7,IM=1)        MOD 0177
C                          FORMAT(3I5,18A4)                             MOD 0178
C                                                                       MOD 0179
CC------------------------ BEGIN ML LOOP                                MOD 0180
CC-                                                                     MOD 0181
CC-             CARD 2C1   ZMDL,P,T,WMOL(1),WMOL(2),WMOL(3),JCHAR       MOD 0182
CC-                        FORMAT(F10.3,5E10.3,15A1)                    MOD 0183
CC-                                                                     MOD 0184
CC-             CARD 2C2   (WMOL(J),J=4,11)                             MOD 0185
CC-                        FORMAT(8E10.3)                               MOD 0186
CC-                                                                     MOD 0187
CC-             CARD 2C2   WMOL(12)                                     MOD 0188
CC-                        FORMAT(8E10.3)                               MOD 0189
CC-                                                                     MOD 0190
CC-             CARD 2C3   AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,               MOD 0191
CC-                        IVUL1,ISEA1,ICHR                             MOD 0192
CC-                        FORMAT(10X,3F10.3,5I5)                       MOD 0193
CC-                                                                     MOD 0194
CC------------------------ END ML LOOP                                  MOD 0195
C                                                                       MOD 0196
C               CARD 2D    IREG(1 TO 4) (IHAZE=7 OR ICLD = 11)          MOD 0197
C                          FORMAT(4I5)                                  MOD 0198
C                                                                       MOD 0199
C               CARD 2D1   AWCCON,TITLE                                 MOD 0200
C                          FORMAT(E10.3,18A4)                           MOD 0201
C                                                                       MOD 0202
C               CARD 2D2   (VX(I),EXTC(N,I),ABSC(N,I),ASYM(N,I),I=1,47) MOD 0203
C                          (IHAZE=7 OR ICLD=11)                         MOD 0204
C                          FORMAT(3(F6.2,2F7.5,F6.4))                   MOD 0205
C                                                                       MOD 0206
C     CARD 3    H1,H2,ANGLE,RANGE,BETA,RO,LEN                           MOD 0207
C                          FORMAT(6F10.3,I5)                            MOD 0208
C                                                                       MOD 0209
C               ALTERNATE  CARD 3 (IEMSCT=3)                            MOD 0210
C                          H1,H2,ANGLE,IDAY,RO,ISOURC,ANGLEM            MOD 0211
C                          FORMAT(3F10.3,I5,5X,F10.3,I5,F10.3)          MOD 0212
C                                                                       MOD 0213
C               CARD 3A1   IPARM,IPH,IDAY,ISOURC           (IEMSCT=2)   MOD 0214
C                          FORMAT(4I5)                                  MOD 0215
C                                                                       MOD 0216
C               CARD 3A2    PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G MOD 0217
C                           FORMAT(8F10.3)               (IEMSCT=2)     MOD 0218
C                                                                       MOD 0219
C               CARD 3B1    NANGLS       (IPH=1)                        MOD 0220
C                           FORMAT(I5)                                  MOD 0221
C                                                                       MOD 0222
C               CARD 3B2(1 TO NANGLS)    (IPH=1)                        MOD 0223
C                      (ANGF(I),F(1,I),F(2,I),F(3,I),F(4,I),I=1,NANGLS) MOD 0224
C                           FORMAT(5E10.3)                              MOD 0225
C                                                                       MOD 0226
C     CARD 4    V1, V2, DV                                              MOD 0227
C                           FORMAT(3F10.3)                              MOD 0228
C                                                                       MOD 0229
C     CARD 5    IRPT                                                    MOD 0230
C                           FORMAT(I5)                                  MOD 0231
C                                                                       MOD 0232
C***********************************************************************MOD 0233
C   ** FOLLOWING IS A FULL DESCRIPTION OF EACH CARD                     MOD 0234
C                                                                       MOD 0235
C     CARD 1    MODEL,ITYPE,IEMSCT,IMULT,M1,M2,M3,                      MOD 0236
C               M4,M5,M6,MDEF,IM,NOPRNT,TBOUND,SALB                     MOD 0237
C                         FORMAT(13I5,F8.3,F7.2)                        MOD 0238
C                                                                       MOD 0239
C             'MODEL' SELECTS ONE OF SIX GEOGRAPHICAL MODEL ATMOSPHERES MOD 0240
C              OR SPECIFIES THAT USER-DEFINED METEOROLOGICAL            MOD 0241
C              DATA ARE TO BE USED.                                     MOD 0242
C                                                                       MOD 0243
C                                                                       MOD 0244
C     MODEL=0 IF METEOROLOGICAL DATA ARE SPECIFIED(HORIZONTAL PATH ONLY)MOD 0245
C           1 TROPICAL ATMOSPHERE                                       MOD 0246
C           2 MIDLATITUDE SUMMER                                        MOD 0247
C           3 MIDLATITUDE WINTER                                        MOD 0248
C           4 SUBARCTIC   SUMMER                                        MOD 0249
C           5 SUBARCTIC   WINTER                                        MOD 0250
C           6 1976 U.S. STANDARD ATMOSPHERE                             MOD 0251
C           7 IF A NEW MODEL ATMOSPHERE( OR RADIOSONDE DATA) IS TO BE   MOD 0252
C             READ IN.                                                  MOD 0253
C                                                                       MOD 0254
C     [NOTE: MODEL=0  USED FOR HORIZONTAL PATH ONLY]                    MOD 0255
C                                                                       MOD 0256
C                                                                       MOD 0257
C           'ITYPE' INDICATES THE TYPE OF ATMOSPHERIC PATH              MOD 0258
C                                                                       MOD 0259
C     ITYPE=1 FOR A HORIZONTAL (CONSTANT-PRESSURE) PATH                 MOD 0260
C           2 VERTICAL OR SLANT PATH BETWEEN TWO ALTITUDES              MOD 0261
C           3 FOR A VERTICAL OR SLANT PATH TO SPACE                     MOD 0262
C                                                                       MOD 0263
C                                                                       MOD 0264
C           'IEMSCT' DETERMINES THE MODE OF EXECUTION OF THE PROGRAM    MOD 0265
C                                                                       MOD 0266
C     IEMSCT=0    PROGRAM EXECUTION IN TRANSMITTANCE MODE.              MOD 0267
C            1    PROGRAM EXECUTION IN RADIANCE MODE.                   MOD 0268
C            2    PROGRAM EXECUTION IN RADIANCE MODE WITH SOLAR/LUNAR   MOD 0269
C                  SCATTERED RADIANCE INCLUDED.                         MOD 0270
C            3    DIRECT SOLAR IRRADIANCE                               MOD 0271
C                                                                       MOD 0272
C           'IMULT' DETERMINES EXECUTION WITH MULTIPLE SCATTERING       MOD 0273
C                                                                       MOD 0274
C     IMULT = 0 PROGRAM EXECUTED WITHOUT MULTIPLE SCATTERING            MOD 0275
C             1 PROGRAM EXECUTED WITH MULTIPLE SCATTERING               MOD 0276
C              [NOTE: IEMSCT MUST EQUAL 1 OR 2 FOR MULTIPLE SCATTERING] MOD 0277
C                                                                       MOD 0278
C                                                                       MOD 0279
C           'M1,M2,M3' ARE USED TO MODIFY OR SUPPLEMENT THE ALTITUDE    MOD 0280
C            PROFILES OF TEMPERATURE AND PRESSURE,WATER VAPOR,AND OZONE MOD 0281
C                                                                       MOD 0282
C           'M4,M5,M6'  SEASONAL DEPENDENCE CH4,N2O,CO                  MOD 0283
C           'MDEF'       USE DEFAULT FOR OTHER GASES                    MOD 0284
C                                                                       MOD 0285
C     FOR NORMAL OPERATION OF PROGRAM   (MODEL 1 TO 6)                  MOD 0286
C     SET M1=M2=M3=0 , M4=M5=M6=MDEF = 0                                MOD 0287
C                                                                       MOD 0288
C     THESE PARAMETERS ARE RESET TO DEFAULT VALUES BY MODEL             MOD 0289
C     WHEN THEY ARE EQUAL TO ZERO                                       MOD 0290
C                                                                       MOD 0291
C      EXCEPT FOR MODEL 0 AND 7                                         MOD 0292
C      WHEN M1 = 0 M1 RESET TO 'MODEL'                                  MOD 0293
C      WHEN M2 = 0 M2 RESET TO 'MODEL'                                  MOD 0294
C      WHEN M3 = 0 M3 RESET TO 'MODEL'                                  MOD 0295
C      WHEN M4 = 0 M4 RESET TO 'MODEL'                                  MOD 0296
C      WHEN M5 = 0 M5 RESET TO 'MODEL'                                  MOD 0297
C      WHEN M6 = 0 M6 RESET TO 'MODEL'                                  MOD 0298
C      WHEN MDEF=0 MDEF RESET TO 1  FOR ALL REMAINING                   MOD 0299
C                                                                       MOD 0300
C     M1=1-6 DEFAULT TEMP. AND PRESSURE TO SPECIFIED MODEL ATM.         MOD 0301
C                                                                       MOD 0302
C     M2=1-6 DEFAULT H2O   TO SPECIFIED MODEL ATM.                      MOD 0303
C                                                                       MOD 0304
C     M3=1-6 DEFAULT OZONE TO SPECIFIED MODEL ATM.                      MOD 0305
C                                                                       MOD 0306
C     M4=1-6 DEFAULT CH4   TO SPECIFIED MODEL ATM.                      MOD 0307
C                                                                       MOD 0308
C     M5=1-6 DEFAULT N2O   TO SPECIFIED MODEL ATM.                      MOD 0309
C                                                                       MOD 0310
C     M6=1-6 DEFAULT CO    TO SPECIFIED MODEL ATM.                      MOD 0311
C                                                                       MOD 0312
C     MDEF=1     USE DEFAULT   PROFILE  FOR CO2,O2,NO,SO2,NO2,NH3,HNO3  MOD 0313
C                NOT NEEDED WITH MODEL 1 TO 6                           MOD 0314
C                                                                       MOD 0315
C                                                                       MOD 0316
C     IF 'MODEL' 0 OR 'MODEL' 7  THE PROGRAM EXPECTS TO READ            MOD 0317
C     "USER SUPPLIED" ATMOSPHERIC PROFILES. SET:'IM' = 1 FOR            MOD 0318
C      FIRST RUN. TO RERUN THE SAME "USER-ATMOSPHERE" FOR A SERIES      MOD 0319
C      OF CASES SET:'IM' = 0 TO REUSE THE PREVIOUSLY READ DATA.         MOD 0320
C                                                                       MOD 0321
C     IM=0    FOR  NORMAL OPERATION OF PROGRAM OR WHEN SUBSEQUENT       MOD 0322
C                  CALCULATIONS ARE TO BE RUN WITH MODEL =7             MOD 0323
C        1    WHEN RADIOSONDE DATA ARE TO BE READ INITIALLY.            MOD 0324
C                                                                       MOD 0325
C     NOPRNT=0 FOR NORMAL OPERATION OF PROGRAM.                         MOD 0326
C                                                                       MOD 0327
C            1 TO MINIMIZE PRINTING OF TRANSMITTANCE /OR RADIANCE TABLE MOD 0328
C                   AND ATMOSPHERIC PROFILES                            MOD 0329
C                                                                       MOD 0330
C                                                                       MOD 0331
C     TBOUND =BOUNDARY TEMPERATURE ( K),USED IN THE RADIATION MODE      MOD 0332
C             (IEMSCT = 1 OR 2) FOR SLANT PATHS THAT INTERSECT THE      MOD 0333
C             EARTH OR TERMINATE AT A GREY BOUNDARY (FOR EXAMPLE        MOD 0334
C             CLOUD,TARGET).  IF TBOUND IS LEFT BLANK AND THE PATH      MOD 0335
C             INTERSECTS THE EARTH, THE PROGRAM WILL USE THE            MOD 0336
C             TEMPERATURE OF THE FIRST ATMOSPHERIC LEVEL AS THE         MOD 0337
C             BOUNDARY TEMPERATURE.                                     MOD 0338
C                                                                       MOD 0339
C      SALB = SURFACE ALBEDO OF THE EARTH AT THE LOCATION               MOD 0340
C             AND AVERAGE FREQUENCY OF THE CALCULATION (0 TO 1.)        MOD 0341
C             IF SALB IS LEFT BLANK THE PROGRAM ASSUMES                 MOD 0342
C             THE SURFACE IS A BLACKBODY.                               MOD 0343
C             NEGITIVE VALUE USE SPECTRALLY DEPENDENT VALUES FROM       MOD 0344
C             REFBKG  SALB  -1 USES THE 1ST FILE                        MOD 0345
C                                                                       MOD 0346
C***********************************************************************MOD 0347
C                                                                       MOD 0348
C     CARD 2   IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,  MOD 0349
C              GNDALT                                                   MOD 0350
C                          FORMAT(6I5,5F10.3)                           MOD 0351
C                                                                       MOD 0352
C     'IHAZE' SELECTS THE TYPE OF EXTINCTION AND A DEFAULT              MOD 0353
C     METEOROLOGICAL RANGE FOR THE BOUNDARY-LAYER AEROSOL MODEL         MOD 0354
C     (0 TO 2KM ALTITUDE)                                               MOD 0355
C     IF 'VIS' IS ALSO SPECIFIED ON CARD 2 IT WILL OVERRIDE THE         MOD 0356
C     DEFAULT 'IHAZE' VALUE  FOR METEOROLOGICAL RANGE                   MOD 0357
C                                                                       MOD 0358
C     IHAZE=0  NO AEROSOL ATTENUATION INCLUDED IN CALCULATION.          MOD 0359
C          =1  RURAL EXTINCTION, 23-KM VIS.                             MOD 0360
C          =2  RURAL EXTINCTION, 5-KM VIS.                              MOD 0361
C          =3  NAVY MARITIME EXTINCTION,SETS OWN VIS.                   MOD 0362
C          =4  MARITIME EXTINCTION, 23-KM VIS.    (LOWTRAN 5 MODEL)     MOD 0363
C          =5  URBAN EXTINCTION, 5-KM VIS.                              MOD 0364
C          =6  TROPOSPHERIC EXTINCTION, 50-KM VIS.                      MOD 0365
C          =7  USER DEFINED  AEROSOL EXTINCTION COEFFICIENTS            MOD 0366
C              TRIGGERS READING IREG FOR UP TO 4 REGIONS OF             MOD 0367
C              USER DEFINED EXTINCTION ABSORPTION AND ASYMMETRY         MOD 0368
C          =8  FOG1 (ADVECTIVE FOG) EXTINCTION, 0.2-KM VIS.             MOD 0369
C          =9  FOG2 (RADIATIVE FOG) EXTINCTION, 0.5-KM VIS.             MOD 0370
C          =10 DESERT EXTINCTION  SETS OWN VISIBILITY FROM WIND SPEED   MOD 0371
C                                                                       MOD 0372
C     'ISEASN' SELECTS THE SEASONAL DEPENDENCE OF THE PROFILES          MOD 0373
C     FOR BOTH THE TROPOSPHERIC (2 TO 10 KM) AND                        MOD 0374
C     STRATOSPHERIC(10 TO 30 KM) AEROSOLS.                              MOD 0375
C                                                                       MOD 0376
C     ISEASN=0 DEFAULTS TO SEASON OF 'MODEL'                            MOD 0377
C              (MODEL 0,1,2,4,6,7) SUMMER                               MOD 0378
C              (MODEL 3,5)         WINTER                               MOD 0379
C           =1 SPRING-SUMMER                                            MOD 0380
C           =2 FALL - WINTER                                            MOD 0381
C                                                                       MOD 0382
C     'IVULCN' SELECTS BOTH THE PROFILE AND EXTINCTION TYPE             MOD 0383
C     FOR THE STRATOSPHERIC AEROSOLS AND DETERMINES TRANSITION          MOD 0384
C     PROFILES ABOVE THE STRATOSPHERE TO 100 KM.                        MOD 0385
C                                                                       MOD 0386
C     IVULCN=0 DEFAULT TO STRATOSPHERIC BACKGROUND                      MOD 0387
C           =1 STRATOSPHERIC BACKGROUND                                 MOD 0388
C           =2 AGED VOLCANIC TYPE/MODERATE VOLCANIC PROFILE             MOD 0389
C           =3 FRESH VOLCANIC TYPE/HIGH VOLCANIC PROFILE                MOD 0390
C           =4 AGED VOLCANIC TYPE/HIGH VOLCANIC PROFILE                 MOD 0391
C           =5 FRESH VOLCANIC TYPE/MODERATE VOLCANIC PROFILE            MOD 0392
C           =6 BACKGROUND STRATOSPHERIC TYPE/MODERATE VOLCANIC PROFILE  MOD 0393
C           =7 BACKGROUND STRATOSPHERIC TYPE/HIGH VOLCANIC PROFILE      MOD 0394
C           =8 FRESH VOLCANIC TYPE/EXTREME VOLCANIC PROFILE             MOD 0395
C                                                                       MOD 0396
C     'ICSTL' IS THE AIR MASS CHARACTER(1 TO 10) ONLY USED WITH         MOD 0397
C     NAVY MARITIME MODEL(IHAZE=3)                                      MOD 0398
C                                                                       MOD 0399
C     ICSTL = 1 OPEN OCEAN                                              MOD 0400
C            .                                                          MOD 0401
C            .                                                          MOD 0402
C            .                                                          MOD 0403
C           10 STRONG CONTINENTAL INFLUENCE                             MOD 0404
C                                                                       MOD 0405
C                                                                       MOD 0406
C     'ICLD' SPECIFIES WHICH OF THE CLOUD MODELS AND THE RAIN RATES     MOD 0407
C     ARE USED                                                          MOD 0408
C                                                                       MOD 0409
C     ICLD  FOR CLOUD AND OR RAIN                                       MOD 0410
C     ICLD = 0   NO CLOUDS OR RAIN                                      MOD 0411
C          = 1  CUMULUS CLOUD BASE .66KM TOP 2.7KM                      MOD 0412
C          = 2  ALTOSTRATUS CLOUD BASE 2.4KM TOP 3.0KM                  MOD 0413
C          = 3  STRATUS CLOUD BASE .33KM TOP 1.0KM                      MOD 0414
C          = 4  STRATUS/STRATO CU BASE .66KM TOP 2.0KM                  MOD 0415
C          = 5  NIMBOSTRATUS CLOUD BASE .16KM TOP .66KM                 MOD 0416
C          = 6  2.0MM/HR DRIZZLE (MODELED WITH CLOUD  3)                MOD 0417
C               RAIN  2. MM HR AT 0KM TO .22 MM HR AT 1.5KM             MOD 0418
C          = 7  5.0MM/HR LIGHT RAIN (MODELED WITH CLOUD  5)             MOD 0419
C               RAIN  5. MM HR AT 0KM TO .2  MM HR AT 1.5KM             MOD 0420
C          = 8  12.5MM/HR MODERATE RAIN (MODELED WITH CLOUD  5)         MOD 0421
C               RAIN 12.5MM HR AT 0KM TO .2  MM HR AT 2.0KM             MOD 0422
C          = 9  25.0MM/HR HEAVY RAIN (MODELED WITH CLOUD  1)            MOD 0423
C               RAIN 25. MM HR AT 0KM TO .2  MM HR AT 3.0KM             MOD 0424
C          =10  75.0MM/HR EXTREME RAIN (MODELED WITH CLOUD  1)          MOD 0425
C               RAIN 75. MM HR AT 0KM TO .2  MM HR AT 3.5KM             MOD 0426
C          =11  READ IN USER DEFINED CLOUD EXTINCTION AND ABSORPTION    MOD 0427
C              USER DEFINED  AEROSOL EXTINCTION COEFFICIENTS            MOD 0428
C              TRIGGERS READING IREG FOR UP TO 4 REGIONS OF             MOD 0429
C              USER DEFINED EXTINCTION ABSORPTION AND ASYMMETRY         MOD 0430
C          =18  STANDARD   CIRRUS MODEL                                 MOD 0431
C          =19  SUB VISUAL CIRRUS MODEL                                 MOD 0432
C          =20  NOAA       CIRRUS MODEL  (LOWTRAN 6 MODEL)              MOD 0433
C                                                                       MOD 0434
C                                                                       MOD 0435
C     IVSA DETERMINES THE USE OF THE ARMY VERTICAL STRUCTURE            MOD 0436
C     ALGORITHM FOR AEROSOLS IN THE BOUNDARY LAYER.                     MOD 0437
C                                                                       MOD 0438
C     IVSA=0   NOT USED                                                 MOD 0439
C         =1   VERTICAL STRUCTURE ALGORITHM                             MOD 0440
C                                                                       MOD 0441
C     'VIS'   SPECIFIES THE METEOROLIGICAL RANGE                        MOD 0442
C     VIS =    METEOROLOGICAL RANGE (KM) (WHEN SPECIFIED,SUPERSEDES     MOD 0443
C              DEFAULT VALUE SET BY IHAZE)                              MOD 0444
C                                                                       MOD 0445
C     'WSS'     SPECIFIES THE CURRENT WIND SPEED                        MOD 0446
C     WSS =    CURRENT WIND SPEED (M/S).    WITH (IHAZE=3/IHAZE=10)     MOD 0447
C                                                                       MOD 0448
C     'WHH'    SPECIFIES THE 24 HOUR AVERAGE WIND SPEED                 MOD 0449
C     WHH =    24 HOUR AVERAGE WIND SPEED (M/S).  ONLY WITH (IHAZE=3)   MOD 0450
C                                                                       MOD 0451
C     'RAINRT' SPECIFIES THE RAIN RATE                                  MOD 0452
C     RAINRT = RAIN RATE (MM/HR).             DEFAULT VALUE IS ZERO.    MOD 0453
C     USED  TO TOP OF CLOUD WHEN CLOUD IS PRESENT                       MOD 0454
C     WHEN NO CLOUDS RAIN RATE USED TO 6KM                              MOD 0455
C                                                                       MOD 0456
C     'GNDALT' SPECIFIES THE ALTITUDE OF SURFACE RELATIVE TO SEA LEVEL  MOD 0457
C     GNDALT = ALTITUDE OF SURFACE RELATIVE TO SEA LEVEL (KM)           MOD 0458
C              USED TO MODIFY  AEROSOL PROFILES BELOW 6 KM ALTITUDE     MOD 0459
C                                                                       MOD 0460
C***********************************************************************MOD 0461
C                                                                       MOD 0462
C     OPTIONAL INPUT CARDS AFTER CARD 2                                 MOD 0463
C     SELECTED BY PARAMETERS ICLD,IVSA,MODEL,AND IHAZE ON CARDS 2       MOD 0464
C                                                                       MOD 0465
C                                                                       MOD 0466
C     CARD 2A   CTHIK,CALT,CEXT,ISEED     (ICLD=18,19,20)               MOD 0467
C                          FORMAT(3F10.3,I10)                           MOD 0468
C                   INPUT CARD FOR CIRRUS ALTITUDE PROFILE              MOD 0469
C                   SUBROUTINE WHEN ICLD = 18,19,20                     MOD 0470
C                                                                       MOD 0471
C     CHTIK    = CIRRUS THICKNESS (KM)                                  MOD 0472
C                0  USE THICKNESS STATISTICS                            MOD 0473
C                                                                       MOD 0474
C     CALT     = CIRRUS BASE ALTITUDE(KM)                               MOD 0475
C                0 USE DEFAULT DETERMINED BY 'MODEL'                    MOD 0476
C                                                                       MOD 0477
C     CEXT     = EXTINCTION COEFFIENT(KM-1) AT 0.55                     MOD 0478
C                0 USE 0.14 * CTHIK                                     MOD 0479
C                                                                       MOD 0480
C     ISEED    = RANDOM NUMBER INITIALIZATION FLAG.                     MOD 0481
C                0 USE DEFAULT MEAN VALUES FOR CIRRUS                   MOD 0482
C                .NE. 0 INITIAL VALUE OF SEED FOR                       MOD 0483
C                RANDOM NUMBER GENERATOR FUNCTION                       MOD 0484
C                CHANGE SEED VALUE EACH RUN FOR DIFFERENT               MOD 0485
C                RANDOM NUMBER SEQUENCES.  THIS PROVIDES FOR            MOD 0486
C                STATISTICAL DETERMINATION OF CIRRUS BASE               MOD 0487
C                ALTITUDE AND THICKNESS.                                MOD 0488
C                                                                       MOD 0489
C   NOTE: RANDOM NUMBERS GENERATION IS SYSTEM DEPENDENT                 MOD 0490
C                                                                       MOD 0491
C***********************************************************************MOD 0492
C                                                                       MOD 0493
C     CARD 2B             ZCVSA,ZTVSA,ZINVSA     (IVSA=1)               MOD 0494
C                          FORMAT(3F10.3)                               MOD 0495
C               INPUT CARD FOR ARMY VERTICAL STRUCTURE                  MOD 0496
C               ALGORITHM SUBROUTINE WHEN IVSA=1.                       MOD 0497
C                                                                       MOD 0498
C     ZCVSA = CLOUD CEILING HEIGHT (KM)                                 MOD 0499
C             LT 0 NO CLOUD CEILING                                     MOD 0500
C             GT 0 KNOWN CLOUD CEILING                                  MOD 0501
C                0 UNKNOWN CLOUD CEILING HEIGHT                         MOD 0502
C                  PROGRAM CALCULATES CLOUD HEIGHT                      MOD 0503
C                                                                       MOD 0504
C     ZTVSA = THICKNESS OF CLOUD OR FOG (KM),                           MOD 0505
C               0 DEFAULTS TO 200 METERS                                MOD 0506
C                                                                       MOD 0507
C     ZINVSA= HEIGHT OF THE INVERSION (KM)                              MOD 0508
C                 0 DEFAULTS TO 100 METERS                              MOD 0509
C             LT  0 NO INVERSION LAYER                                  MOD 0510
C                                                                       MOD 0511
C***********************************************************************MOD 0512
C                                                                       MOD 0513
C     CARD 2C  ML,IRD1,IRD2,TITLE   (MODEL=0 / 7,IM=1)                  MOD 0514
C                          FORMAT(3I5,18A4)                             MOD 0515
C              ADDITIONAL ATMOSPHERIC MODEL       (MODEL=7)             MOD 0516
C              NEW MODEL ATMOSPHERE CAN BE INSERTED PROVIDED THE        MOD 0517
C              PARAMETERS 'MODEL' AND 'IM' ARE SET EQUAL TO 7 AND 1     MOD 0518
C              RESPECTIVELY ON CARD 1.                                  MOD 0519
C                                                                       MOD 0520
C     ML=      NUMBER OF ATMOSPHERIC LEVELS TO BE INSERTED              MOD 0521
C                   (MAXIMUM OF LAYDIM)                                 MOD 0522
C                                                                       MOD 0523
C     CARD 2C1 IS READ AUTOMATICALLY FOR MODEL 0 AND 7                  MOD 0524
C                                                                       MOD 0525
C     IRD1 CONTROL READING WN2O,WCO ... AND WNH3,WHNO3      CARD        MOD 0526
C                                                                       MOD 0527
C     IRD1 = 0     NO READ  (MOLECULAR DENSITIES BY LAYER)              MOD 0528
C     IRD1 = 1     READ                                                 MOD 0529
C                                                                       MOD 0530
C     IRD2  CONTROL READING AHAZE,EQLWCZ,... CARD                       MOD 0531
C                                                                       MOD 0532
C     IRD2 = 0     NO READ  (AEROSOL CONTROL BY LAYER)                  MOD 0533
C     IRD2 = 1     READ                                                 MOD 0534
C                                                                       MOD 0535
C     JCHAR  ON CARD 2C1 IS USUALLY USED TO DEFINE MOLECULES 4 TO 12    MOD 0536
C     IHAZE  ON CARD 2   IS USUALLY USED TO DEFINE AEROSOL PROFILES     MOD 0537
C     IRD1 = 1 OR IRD2 = 1   SELDOM USED                                MOD 0538
C                                                                       MOD 0539
C     TITLE=   IDENTIFICATION OF NEW MODEL ATMOSPHERE                   MOD 0540
C                                                                       MOD 0541
C                                                                       MOD 0542
C     THE FOLLOWING CARDS ARE READ IN SUBROUTINE AERNSM                 MOD 0543
C                                                                       MOD 0544
CC------------------------ BEGIN ML LOOP                                MOD 0545
CC-                                                                     MOD 0546
CC-   CARD 2C1 ZMDL,P,T,WMOL(1),WMOL(2),WMOL(3),JCHAR                   MOD 0547
CC-   LAYER VARIABLES        WH,   WCO2,     WO,JCHAR (1 TO 13)         MOD 0548
CC-                        FORMAT(F10.3,5E10.3,15A1)                    MOD 0549
CC-                                                                     MOD 0550
CC-   ZMDL     ALTITUDE OF LAYER (KM)                                   MOD 0551
CC-   P        PRESSURE AT LAYER                                        MOD 0552
CC-   T        TEMPERATURE                                              MOD 0553
CC-   WMOL     READ, INTERPRETED AND MOVED INTO LAYER VARIABLES         MOD 0554
CC-   WH =     WATER VAPOR                                              MOD 0555
CC-   WCO2 =   CO2                                                      MOD 0556
CC-   WO =     OZONE                                                    MOD 0557
CC-                                                                     MOD 0558
CC-   JCHAR    FLAGS TO SPECIFY UNITS OR DEFAULTS FOR                   MOD 0559
CC-   P,T,WH,WCO2,WO,WN2O,WCO,.. AND WNH3,WHNO3                         MOD 0560
CC-   JCHAR BLANK DEFAULT TO M1,M2,M3,M4,M5,M6,MDEF WHEN AMOUNT ZERO    MOD 0561
CC-                                                                     MOD 0562
CC-      PARAMETERS - JCHAR = INPUT KEY                                 MOD 0563
CC-                                                                     MOD 0564
CC-   **  ACCEPTS VARIABLE UNITS ON PRESS AND TEMP                      MOD 0565
CC-                                                                     MOD 0566
CC-     JCHAR(1)                                                        MOD 0567
CC-                                                                     MOD 0568
CC-    " ",A           PRESSURE IN (MB)    OR BLANK                     MOD 0569
CC-        B              "     "  (ATM)                                MOD 0570
CC-        C              "     "  (TORR)                               MOD 0571
CC-       1-6          DEFAULT TO SPECIFIED MODEL ATMOSPHERE            MOD 0572
CC-                                                                     MOD 0573
CC-     JCHAR(2)                                                        MOD 0574
CC-    " ",A           AMBIENT TEMPERATURE IN DEG(K)  OR BLANK          MOD 0575
CC-        B              "         "       "  " (C)                    MOD 0576
CC-       1-6          DEFAULT TO SPECIFIED MODEL ATMOSPHERE            MOD 0577
CC-                                                                     MOD 0578
CC-   ****************************************************************  MOD 0579
CC-   FOR MOLECULAR SPECIES ONLY                                        MOD 0580
CC-                                                                     MOD 0581
CC-     JCHAR             JCHAR(3-13)                                   MOD 0582
CC-                                                                     MOD 0583
CC-   " ",A            VOLUME MIXING RATIO (PPMV)                       MOD 0584
CC-       B            NUMBER DENSITY (CM-3)                            MOD 0585
CC-       C            MASS MIXING RATIO (GM(K)/KG(AIR))                MOD 0586
CC-       D            MASS DENSITY (GM M-3)                            MOD 0587
CC-       E            PARTIAL PRESSURE (MB)                            MOD 0588
CC-       F            DEW POINT TEMP (TD IN T(K)) - H2O ONLY           MOD 0589
CC-       G             "    "     "  (TD IN T(C)) - H2O ONLY           MOD 0590
CC-       H            RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY (3) MOD 0591
CC-       I            AVAILABLE FOR USER DEFINITION                    MOD 0592
CC-      1-6           DEFAULT TO SPECIFIED MODEL ATMOSPHERE            MOD 0593
CC-                                                                     MOD 0594
CC-   ****************************************************************  MOD 0595
CC-                                                                     MOD 0596
CC-   CARD 2C2   (WMOL(J),J=4,11)                                       MOD 0597
CC-   VARIABLES  WN2O,WCO,WCH4,WO2,WNO,WSO2,WNO2,WNH3                   MOD 0598
CC-                        FORMAT(8E10.3)                               MOD 0599
CC-                                                                     MOD 0600
CC-   CARD 2C2   WMOL(12)             (CONT)                            MOD 0601
CC-   VARIABLE   WHNO3                                                  MOD 0602
CC-                        FORMAT(8E10.3)                               MOD 0603
CC-                                                                     MOD 0604
CC-                                                                     MOD 0605
CC-   WMOL     READ, INTERPRETED AND MOVED INTO LAYER VARIABLES         MOD 0606
CC-   WN2O =   N2O                                                      MOD 0607
CC-   WCO  =   CO                                                       MOD 0608
CC-   WCH4 =   CH4                                                      MOD 0609
CC-   WO2  =   O2                                                       MOD 0610
CC-   WNO  =   NO                                                       MOD 0611
CC-   WSO2 =   SO2                                                      MOD 0612
CC-   WNO2 =   NO2                                                      MOD 0613
CC-   WNH3 =   NH3                                                      MOD 0614
CC-   WHNO3 =  HNO3                                                     MOD 0615
CC-                                                                     MOD 0616
CC-   ****************************************************************  MOD 0617
CC-                                                                     MOD 0618
CC- CARD 2C3     AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,ICHR         MOD 0619
CC-                        FORMAT(10X,3F10.3,5I5)                       MOD 0620
CC-                                                                     MOD 0621
CC-  'AHAZE' AEROSOL SCALING FACTOR (EQUAL TO THE VISIBLE [0.55UM]      MOD 0622
CC-          EXTINCTION COEFFICIENT AT ALTZ)                            MOD 0623
CC-                                                                     MOD 0624
CC-          [NOTE ** ONE OF AHAZE OR EQLWCZ IS ALLOWED ]               MOD 0625
CC-                                                                     MOD 0626
CC-  'EQLWCZ' EQUIVALENT LIQUID WATER CONTENT ( GM/M3) AT ALT Z         MOD 0627
CC-           FOR THE AEROSOL, CLOUD OR FOG MODELS                      MOD 0628
CC-                                                                     MOD 0629
CC-   RRATZ=RAIN RATE (MM/HR) AT ALT Z                                  MOD 0630
CC-                                                                     MOD 0631
CC-   'IHA1' AEROSOL EXTINCTION AND METEOROLOGICAL RANGE CONTROL FOR    MOD 0632
CC-          THE ALTITUDE, Z                                            MOD 0633
CC-                                                                     MOD 0634
CC-   'ICLD1' CLOUD EXTINCTION CONTROL FOR THE ALTITUDE, Z              MOD 0635
CC-                                                                     MOD 0636
CC-   WHEN USING 'ICLD1' IT IS NECESSARY TO SET 'ICLD' (NON ZERO)       MOD 0637
CC-                                                                     MOD 0638
CC-   'IVUL1' STRATOSPHERIC AEROSOL PROFILE AND EXTINCTION CONTROL FOR  MOD 0639
CC-   THE ALTITUDE Z                                                    MOD 0640
CC-                                                                     MOD 0641
CC-   ONLY ONE OF 'IHA1','ICLD1' OR 'IVUL1' IS ALLOWED                  MOD 0642
CC-   IF 'IHA1' NE 0  THEN OTHERS IGNORED                               MOD 0643
CC-   IF 'IHA1' EQ 0 AND 'ICLD1' NE 0 THEN USE 'ICLD1'                  MOD 0644
CC-                                                                     MOD 0645
CC-   IF 'AHAZE' AND 'EQLWCZ' ARE BOTH ZERO  DEFAULT PROFILE LOADED     MOD 0646
CC-      FROM 'IHAZ1','ICLD1','IVUL1'                                   MOD 0647
CC-                                                                     MOD 0648
CC-   'ISEA1'  AEROSOL SEASON CONTROL FOR THE ALTITUDE,Z                MOD 0649
CC-                                                                     MOD 0650
CC-   'ICHR '  CHANGE PROFILE REGION  AT ALTITUDE Z                     MOD 0651
CC-            USED WHEN IHA1 IS 7 IN TWO ADJACENT ALTITUDE REGIMES     MOD 0652
CC-                                                                     MOD 0653
CC------------- END ML LOOP                                             MOD 0654
C***********************************************************************MOD 0655
C                                                                       MOD 0656
C   IHAZE = 7 OR ICLD = 11 INPUT                                        MOD 0657
C                                                                       MOD 0658
C   CARD 2D  (IREG(II),II=1,4)                                          MOD 0659
C                        FORMAT(4I5)                                    MOD 0660
C                                                                       MOD 0661
C             'IREG' SPECIFIES WHICH OF THE FOUR AEROSOL REGIONS        MOD 0662
C              A USER DEFINED AEROSOL MODEL IS USED (IHAZE=7/ICLD=11)   MOD 0663
C                                                                       MOD 0664
C             [NOTE   REGIONS DEFAULT TO                                MOD 0665
C             0-2 ,3-10,11-30,35-100 KM                                 MOD 0666
C             AND CAN BE OVERRIDDEN WITH 'IHA1' SETTINGS IN MODEL 7]    MOD 0667
C                                                                       MOD 0668
C             IREG = 0  USE DEFAULT VALUES FOR THIS REGION              MOD 0669
C                                                                       MOD 0670
C             IREG = 1  READ EXTINCTION ABSORPTION ASYMMETRY            MOD 0671
C             FOR A REGION                                              MOD 0672
C                                                                       MOD 0673
C   CARD 2D1 AWCCON,TITLE                                               MOD 0674
C                        FORMAT(E10.3,18A4)                             MOD 0675
C                                                                       MOD 0676
C            'AWWCON' EQUIVALENT LIQUID WATER CONTENT(GM/M3)            MOD 0677
C             FOR A REGION                                              MOD 0678
C                                                                       MOD 0679
C             'TITLE' FOR A AEROSOL REGION                              MOD 0680
C                                                                       MOD 0681
C    CARD 2D2 (VX(I),EXTC(N,I),ABSC(N,I),ASYM(N,I),I=1,47)              MOD 0682
C                   FORMAT(4(F6.2,2F7.5,F6.4))                          MOD 0683
C                                                                       MOD 0684
C             WHERE  N = IREG(II)   FOR UP TO 4 ALTITUDE REGIONS        MOD 0685
C             USER DEFINED AEROSOL  OR CLOUD EXTINCTION AND ABSORPTION  MOD 0686
C              COEFFICIENTS WHEN IHAZE = 7 OR ICLD = 11                 MOD 0687
C                                                                       MOD 0688
C     VX(I)    = WAVELENGTH OF AEROSOL COEFFICIENT                      MOD 0689
C                   (NOT USED BY PROGRAM BUT CORRESPONDING TO           MOD 0690
C                   WAVELENGTHS DEFINED IN ARRAY VX0                    MOD 0691
C                   IN SUBROUTINE EXTDTA)                               MOD 0692
C                                                                       MOD 0693
C     EXTC(N,I) = AEROSOL EXTINCTION COEFFICIENT                        MOD 0694
C     ABSC(N,I) = AEROSOL ABSORPTION COEFFICIENT                        MOD 0695
C     ASYM(N,I) = AEROSOL ASYMMETRY PARAMETER                           MOD 0696
C     WHERE  N = IREG(II)   FOR UP TO 4 ALTITUDE REGIONS                MOD 0697
C                                                                       MOD 0698
C***********************************************************************MOD 0699
C                                                                       MOD 0700
C     CARD 3    H1,H2,ANGLE,RANGE,BETA,RO,LEN    FORMAT(6F10.3,I5)      MOD 0701
C            USED TO DEFINE THE GEOMETRICAL PATH PARAMETERS FOR A GIVEN MOD 0702
C            PROBLEM.                                                   MOD 0703
C                                                                       MOD 0704
C     OR IF IEMSCT=3 ; CARD 3 H1,H2,ANGLE,IDAY,RO,ISOURC,ANGLEM         MOD 0705
C                                                                       MOD 0706
C     H1  =  INITIAL ALTITUDE(KM)                                       MOD 0707
C     H2  =  FINAL ALTITUDE(KM)                                         MOD 0708
C                                                                       MOD 0709
C                   IN THE RADIANCE MODE OF THE PROGRAM EXECUTION       MOD 0710
C            H1, THE INITIAL ALTITUDE,ALWAYS DEFINES THE POSITION OF    MOD 0711
C            THE OBSERVER (OR SENSOR).                                  MOD 0712
C                                                                       MOD 0713
C     ANGLE =INITIAL ZENITH ANGLE (DEGREES) AS MEASURED FROM H1         MOD 0714
C     [NOTE: ANGLE = 0 LOOKS STRAIGHT UP.                               MOD 0715
C            ANGLE IS DEFINED  FROM 0 TO 180 DEGREES ]                  MOD 0716
C                                                                       MOD 0717
C     RANGE =PATH LENGTH (KM)                                           MOD 0718
C     BETA  =EARTH CENTER ANGLE SUBTENDED BY H1 AND H2 (DEGREES)        MOD 0719
C                                                                       MOD 0720
C     RO =   RADIUS OF THE EARTH (KM) AT THE PARTICULAR GEOGRAPHICAL    MOD 0721
C            LOCATION AT WHICH THE CALCULATION IS TO BE PERFORMED.      MOD 0722
C              IF RO BLANK PROGRAM USES RADIUS FOR APPROPRIATE MODEL    MOD 0723
C              ATMOSPHERE. (MODEL 0 OR 7 DEFAULT = 6371.23 KM)          MOD 0724
C                                                                       MOD 0725
C     LEN =0 FOR NORMAL OPERATION OF PROGRAM                            MOD 0726
C         =1 LONG PATH THROUGH THE TANGENT HEIGHT                       MOD 0727
C                                                                       MOD 0728
C            IT IS NOT NECESSARY TO SPECIFY EVERY QUANTITY GIVEN ABOVE  MOD 0729
C            ONLY THOSE THAT ADEQUATELY DESCRIBE THE PROBLEM ACCORDING  MOD 0730
C            TO THE PARAMETER ITYPE                                     MOD 0731
C                                                                       MOD 0732
C            ITYPE=1 READ H1,RANGE                                      MOD 0733
C                 =2 READ H1,H2,ANGLE   OR H1,H2,RANGE   OR H1,H2,BETA  MOD 0734
C                    OR H1,ANGLE,RANGE                                  MOD 0735
C                 =3 READ H1,ANGLE OR H1,H2                             MOD 0736
C                    [NOTE: H2 IS INTERPRETED AS HMIN FOR THIS CASE]    MOD 0737
C                                                                       MOD 0738
C--------------                                                         MOD 0739
C     CARD 3    OPTION  IEMSCT = 3                                      MOD 0740
C    'IDAY'     DAY OF THE YEAR, I.E. FROM 1 TO 365  (IEMSCT = 3)       MOD 0741
C                                                                       MOD 0742
C     ISOURC=0  EXTRATERRESTRIAL SOURCE IS THE SUN                      MOD 0743
C           =1  EXTRATERRESTRIAL SOURCE IS THE MOON                     MOD 0744
C                                                                       MOD 0745
C     ANGLEM=PHASE ANGLE OF THE MOON, I.E. THE ANGLE FORMED             MOD 0746
C            BY THE SUN, MOON AND EARTH (REQUIRED IF ISOURC=1)          MOD 0747
C                                                                       MOD 0748
C***********************************************************************MOD 0749
C                                                                       MOD 0750
C     CARD 3A1   IPARM,IPH,IDAY,ISOURC           (IEMSCT=2)             MOD 0751
C                          FORMAT(4I5)                                  MOD 0752
C              INPUT CARD FOR SOLAR/LUNAR SCATTERED RADIATION WHEN      MOD 0753
C              IEMSCT = 2                                               MOD 0754
C                                                                       MOD 0755
C              IPARM =0,1,2 AND CONTROLS THE METHOD OF SPECIFYING THE   MOD 0756
C              SOLAR/LUNAR GEOMETRY ON CARD 3A2.                        MOD 0757
C                     (SEE DEFINITION BELOW FOR CARD 3A2)               MOD 0758
C                                                                       MOD 0759
C              IPH DETERMINES THE TYPE OF PHASE FUNCTION USED IN THE    MOD 0760
C               CALCULATION                                             MOD 0761
C                                                                       MOD 0762
C     IPH=0     HENYEY-GREENSTEIN AEROSOL PHASE FUNCTION                MOD 0763
C        =1     USER SUPPLIED AEROSOL PHASE FUNCTION (SEE CARD 3B)      MOD 0764
C        =2     MIE GENERATED DATA BASE OF AEROSOL PHASE FUNCTIONS FOR  MOD 0765
C               THE LOWTRAN MODELS.                                     MOD 0766
C                                                                       MOD 0767
C     IDAY=     DAY OF THE YEAR, I.E. FROM 1 TO 365   (REQUIRED)        MOD 0768
C                                                                       MOD 0769
C     ISOURC=0  EXTRATERRESTRIAL SOURCE IS THE SUN                      MOD 0770
C           =1  EXTRATERRESTRIAL SOURCE IS THE MOON                     MOD 0771
C                                                                       MOD 0772
C***********************************************************************MOD 0773
C                                                                       MOD 0774
C     CARD 3A2  PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G             MOD 0775
C                          FORMAT(8F10.3)                (IEMSCT=2)     MOD 0776
C              INPUT CARD FOR SOLAR/LUNAR SCATTERED RADIATION WHEN      MOD 0777
C              IEMSCT = 2                                               MOD 0778
C              DEFINITIONS OF PARM1,PARM2,PARM3,PARM4 DETERMINED BY     MOD 0779
C              VALUE OF IPARM ON CARD 3A1.                              MOD 0780
C                                                                       MOD 0781
C                       FOR IPARM=0                                     MOD 0782
C                                                                       MOD 0783
C     PARM1= OBSERVER LATITUDE (-90 TO +90)                             MOD 0784
C          NOTE- IF ABS(PARM1) IS GREATER THAN 89.5 THE OBSERVER IS     MOD 0785
C          ASSUMED TO BE AT EITHER THE NORTH OR THE SOUTH POLE.  IN     MOD 0786
C          THAT CASE THE PATH AZIMUTH IS UNDEFINED.  THE DIRECTION OF   MOD 0787
C          LINE OF SIGHT MUST BE SPECIFIED AS THE LONGITUDE ALONG WHICH MOD 0788
C          THE PATH LIES. THIS QUANTITY RATHER THAN THE USUAL AZIMUTH   MOD 0789
C          IS READ IN                                                   MOD 0790
C     PARM2= OBSERVER LONGITUDE (0 TO 360)                              MOD 0791
C     PARM3= SOURCE (SUN OR MOON) LATITUDE                              MOD 0792
C     PARM4= SOURCE (SUN OR MOON) LONGITUDE                             MOD 0793
C                                                                       MOD 0794
C                       FOR IPARM=1                                     MOD 0795
C        (IDAY AND TIME MUST BE SPECIFIED,CANNOT BE USED WITH ISOURC=1) MOD 0796
C                                                                       MOD 0797
C                                                                       MOD 0798
C     PARM1= OBSERVER LATITUDE (-90 TO +90)                             MOD 0799
C     PARM2= OBSERVER LONGITUDE (0 TO 360)                              MOD 0800
C             PARM3,PARM4 ARE NOT REQUIRED                              MOD 0801
C                                                                       MOD 0802
C     [NOTE: THAT THE CALCULATED APPARENT SOLAR ZENITH                  MOD 0803
C            ANGLE IS THE ZENITH ANGLE AT H1 OF THE REFRACTED           MOD 0804
C            PATH TO THE SUN AND IS LESS THAN THE ASTRONOMICAL          MOD 0805
C            SOLAR ZENITH ANGLE.  THE DIFFERENCE BETWEEN THE            MOD 0806
C            TWO ANGLES IS NEGLIGIBLE FOR ANGLES LESS THAN 80           MOD 0807
C            DEGREES.]                                                  MOD 0808
C                                                                       MOD 0809
C                        FOR IPARM=2                                    MOD 0810
C                                                                       MOD 0811
C     PARM1= AZIMUTHAL ANGLE BETWEEN THE OBSERVER'S LINE OF SIGHT       MOD 0812
C            AND THE OBSERVER-TO-SUN PATH, MEASURED FROM THE LINE       MOD 0813
C            OF SIGHT, POSITIVE EAST OF NORTH, BETWEEN -180 AND 180     MOD 0814
C     PARM2= THE SUN'S ZENITH ANGLE                                     MOD 0815
C                                                                       MOD 0816
C              PARM3,PARM4 ARE NOT REQUIRED                             MOD 0817
C                                                                       MOD 0818
C                                                                       MOD 0819
C              REMAINING CONTROL PARAMETERS                             MOD 0820
C                                                                       MOD 0821
C     TIME=  GREENWICH TIME IN DECIMAL HOURS, I.E. 8:45 AM IS 8.75,     MOD 0822
C            5:20 PM IS 17.33 ETC. (USED WITH IPARM = 1)                MOD 0823
C                                                                       MOD 0824
C     PSIPO= PATH AZIMUTH (DEGREES EAST OF NORTH, I.E. DUE NORTH IS 0.0 MOD 0825
C            DUE EAST IS 90.0 ETC.  (USED WITH IPARM = 0 OR 1)          MOD 0826
C                                                                       MOD 0827
C     ANGLEM=PHASE ANGLE OF THE MOON, I.E. THE ANGLE FORMED             MOD 0828
C            BY THE SUN, MOON AND EARTH (REQUIRED IF ISOURC=1)          MOD 0829
C                                                                       MOD 0830
C     G=     ASYMMETRY FACTOR FOR USE WITH HENYEY-GREENSTEIN            MOD 0831
C            PHASE FUNCTION (USED WITH IPH = 0)                         MOD 0832
C                                                                       MOD 0833
C***********************************************************************MOD 0834
C     CARD 3B1 NANGLS          (IPH=1)                                  MOD 0835
C                   FORMAT(I5)                                          MOD 0836
C                                                                       MOD 0837
C              INPUT CARD FOR USER DEFINED PHASE FUNCTIONS WHEN IPH=1.  MOD 0838
C                                                                       MOD 0839
C     NANGLS=  NUMBER OF ANGLES FOR THE USER DEFINED PHASE              MOD 0840
C              FUNCTIONS(MAXIMUM OF 50)                                 MOD 0841
C                                                                       MOD 0842
C***********************************************************************MOD 0843
C                                                                       MOD 0844
C     CARD 3B2(1 TO NANGLS)    (IPH=1)                                  MOD 0845
C             (ANGF(I),F(1,I),F(2,I),F(3,I),F(4,I),I=1,NANGLS)          MOD 0846
C                   FORMAT(5E10.3)                                      MOD 0847
C              INPUT CARD FOR USER DEFINED PHASE FUNCTION WHEN IPH=1.   MOD 0848
C              FOR AVERAGE FREQUENCY OF CALCULATION                     MOD 0849
C                                                                       MOD 0850
C     ANGF(I)= PHASE ANGLE IN DECIMAL DEGREES                           MOD 0851
C                   (0.0 TO 180.0)                                      MOD 0852
C                                                                       MOD 0853
C     F(1,I)=  USER DEFINED PHASE FUNCTION AT ANGF(I)                   MOD 0854
C              BOUNDARY LAYER DEFAULTS TO (0 TO 2KM))                   MOD 0855
C                                                                       MOD 0856
C     F(2,I)=  USER DEFINED PHASE FUNCTION AT ANGF(I)                   MOD 0857
C              TROPOSPHERE DEFAULTS TO (2 TO 10 KM)                     MOD 0858
C                                                                       MOD 0859
C     F(3,I)=  USER DEFINED PHASE FUNCTION AT ANGF(I)                   MOD 0860
C              STRATOSPHERE DEFAULTS TO (10 TO 30 KM)                   MOD 0861
C                                                                       MOD 0862
C     F(4,I)=  USER DEFINED PHASE FUNCTION AT ANGF(I)                   MOD 0863
C              MESOSPHERE DEFAULTS TO (30 TO 100 KM)                    MOD 0864
C                                                                       MOD 0865
C***********************************************************************MOD 0866
C                                                                       MOD 0867
C     CARD 4    V1, V2, DV                       FORMAT(3F10.3)         MOD 0868
C                                                                       MOD 0869
C              THE SPECTRAL RANGE OVER WHICH DATA ARE REQUIRED AND      MOD 0870
C              THE SPECTRAL INCREMENTS AT WHICH THE DATA ARE TO BE      MOD 0871
C              CALCULATED  IS DETERMINED BY CARD 4.                     MOD 0872
C                                                                       MOD 0873
C     V1 =     INITIAL FREQUENCY (WAVENUMBER CM-1 )                     MOD 0874
C     V2 =     FINAL FREQUENCY(WAVENUMBER CM-1 )                        MOD 0875
C     DV =     FREQUENCY INCREMENT (OR STEP SIZE) (CM-1)                MOD 0876
C              NOTE: DV MUST BE A MULTIPLE OF 5 CM-1                    MOD 0877
C              ANY STEP SIZE .GT. 5 CM-1 WILL UNDERSAMPLE THE RESULTS   MOD 0878
C                                                                       MOD 0879
C              SCANNING FUNCTION IS AVAILABLE TO PROPERLY HANDLE DATA   MOD 0880
C              WITH LOWER RESOLUTION THAN 20CM-1 LOWTRAN 7              MOD 0881
C                                                                       MOD 0882
C***********************************************************************MOD 0883
C                                                                       MOD 0884
C     CARD 5    IRPT                             FORMAT(I5)             MOD 0885
C     IRPT=0  TO END PROGRAM                                            MOD 0886
C         =1  READ ALL DATA CARDS (1,2,3,4,5)                           MOD 0887
C         =2  NOT USED  (WILL STOP PROGRAM)                             MOD 0888
C         =3  READ CARD 3   THE GEOMETRY CARD AND CARD 5                MOD 0889
C         =4  READ CARD 4 TO CHANGE FREQUENCY AND CARD 5                MOD 0890
C     GT 4 OR IRPT=2 WILL CAUSE PROGRAM TO STOP                         MOD 0891
C                                                                       MOD 0892
C     IRPT GE 1 USED FOR MULTIPLE RUNS OF LOWTRAN                       MOD 0893
C     WARNING IRPT 3 CANNOT BE USED WHEN RUNNING MULTIPLE SCATTERING    MOD 0894
C     CASES WITH SOLAR SCATTERING                                       MOD 0895
C                                                                       MOD 0896
C     REFERENCES                                                        MOD 0897
C                                                                       MOD 0898
C       (1980) ATMOSPHERIC TRANSMITTANCE/RADIANCE- COMPUTER CODE        MOD 0899
C       LOWTRAN 5 AFGL-TR-80-0067                                       MOD 0900
C       KNEIZYS, F. X.,SHETTLE, E. P. ,GALLERY, W. O.,CHETWYND, J. H.,  MOD 0901
C       ABREU, L. W., SELBY, J. E. A., FENN, R. W. ,MCCLATCHEY R. A.    MOD 0902
C                                                                       MOD 0903
C       (1983) ATMOSPHERIC TRANSMITTANCE/RADIANCE- COMPUTER CODE        MOD 0904
C       LOWTRAN 6  AFGL TR 83 0187                                      MOD 0905
C       KNEIZYS, F. X.,SHETTLE, E. P. ,GALLERY, W. O.,CHETWYND, J. H.,  MOD 0906
C       ABREU, L. W., SELBY, J. E. A., CLOUGH, S. A., FENN, R. W.       MOD 0907
C                                                                       MOD 0908
C       (1988) ATMOSPHERIC TRANSMITTANCE/RADIANCE- COMPUTER CODE        MOD 0909
C       LOWTRAN 7 AFGL-TR-88-XXXX                                       MOD 0910
C       KNEIZYS, F. X.,SHETTLE, E. P. ,ANDERSON G. P.,ABREU ,L. W.      MOD 0911
C       CHETWYND, J H,SELBY, J. E. A. ,CLOUGH, S. A.,GALLERY, W. O      MOD 0912
C                                                                       MOD 0913
C       (1988) LOWTRAN 7 COMPUTER CODE : USER'S MANUAL AFGL-TR-88-XXXX  MOD 0914
C       KNEIZYS, F. X.,SHETTLE, E. P. ,ANDERSON G. P.,ABREU ,L. W.      MOD 0915
C       CHETWYND, J H,SELBY, J. E. A. ,CLOUGH, S. A.,GALLERY, W. O      MOD 0916
C                                                                       MOD 0917
C       MOLECULAR TRANSMISSION BAND MODELS FOR LOWTRAN AFGL-TR-86-0272  MOD 0918
C       PIERLUISSI, J. H., MARAGOUDAKIS, C. E.                          MOD 0919
C                                                                       MOD 0920
C       MULTIPLE SCATTERING TREATMENT FOR USE IN                        MOD 0921
C       THE LOWTRAN AND FASCODE MODELS  AFGL-TR-86-0073                 MOD 0922
C       ISAACS, R. G., WANG, W. C., WORSHAM, R. D.,GOLDENBERG S.        MOD 0923
C                                                                       MOD 0924
C       AFGL ATMOSPHERIC CONSTITUENT PROFILES (0 TO 120KM)              MOD 0925
C                                              AFGL-TR-86-0110          MOD 0926
C       ANDERSON, G. P., CLOUGH, S. A., KNEIZYS, F. X.,                 MOD 0927
C       CHETWYND, J. H., SHETTLE, E. P.                                 MOD 0928
C                                                                       MOD 0929
C       PROGRAM FOR ATMOSPHERIC TRANSMITTANCE RADIANCE/CALCULATIONS     MOD 0930
C       FSCATM                                  AFGL-TR-83-0065         MOD 0931
C       GALLERY, W. O., KNEIZYS, F. X., AND CLOUGH, S. A.               MOD 0932
C                                                                       MOD 0933
C       AFGL HANDBOOK OF GEOPHYSICS AND THE SPACE ENVIRONMENT           MOD 0934
C       EDITOR, A. S. JURSA  CHAPTER 18 1985                            MOD 0935
C                                                                       MOD 0936
C       MODELS OF THE AEROSOLS OF THE LOWER ATMOSPHERE AND THE EFFECTS  MOD 0937
C       OF HUMIDITY VARIATIONS ON THEIR OPTICAL PROPERTIES              MOD 0938
C       SHETTLE, E.P. AND FENN, R. W.            AFGL-TR-79-0214        MOD 0939
C                                                                       MOD 0940
C       OPTICAL PROPAGATION IN THE ATMOSPHERE    AGARD-CP-183  1975     MOD 0941
C       SHETTLE, E. P., AND FENN, R. W.          NTIS (NO AD-A028-615)  MOD 0942
C                                                                       MOD 0943
C                                                                       MOD 0944
C       ATMOSPHERIC ATTENUATION OF MILLIMETER AND SUBMILLIMETER         MOD 0945
C       WAVES:  MODEL AND COMPUTER CODE          AFGL-TR-79-0253        MOD 0946
C       FALCONE,V. J. JR.,ABREU,L. W. AND SHETTLE, E. P.                MOD 0947
C                                                                       MOD 0948
C       LOWTRAN  PLUS ULTRAVIOLET O2 ABSORPTION                         MOD 0949
C                                                                       MOD 0950
C       REFERENCES- JOHNSTON, ET AL, J GEOPHYS RES, 89,11661-11665,1984.MOD 0951
C                                                                       MOD 0952
C       FREQUENCY RANGE: 50000-36000CM-1 FOR HERZBERG CALCULATION       MOD 0953
C                                                                       MOD 0954
C       THE SCHUMANN-RUNGE BANDS (PARTICULARLY THE 1,0 AND 0,0) ARE NOT MOD 0955
C       INCLUDED IN THE CALCULATIONS (50000 AND 49400 CM-1).            MOD 0956
C       THE HERZBERG BANDS ARE APPROXIMATED BY AN EXTRAPOLATION OF THE  MOD 0957
C       HERZBERG CONTINUUM (41322-36000 CM-1).                          MOD 0958
C                                                                       MOD 0959
C***********************************************************************MOD 0960
      CALL DRIVER                                                       MOD 0961
      STOP                                                              MOD 0962
C@    END                                                               MOD 0963
C@    THE FOLLOWING TIME AND DATE SUBROUTINES APPLY TO A CDC 6600       MOD 0964
C@    SUBROUTINE FDATE(HDATE)                                           MOD 0965
C@    CALL DATE(GDATE)                                                  MOD 0966
C@    HDATE=SHIFT(GDATE,6)                                              MOD 0967
C@    RETURN                                                            MOD 0968
C@    END                                                               MOD 0969
C@    SUBROUTINE FCLOCK(HTIME)                                          MOD 0970
C@    CALL CLOCK(GTIME)                                                 MOD 0971
C@    HTIME=SHIFT(GTIME,6)                                              MOD 0972
C@    RETURN                                                            MOD 0973
      END                                                               MOD 0974
