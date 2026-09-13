      SUBROUTINE VSA(IHAZE,VIS,CEILHT,DEPTH,ZINVHT,Z,RH,AHAZE,IH)       VSA 0001
C                                                                       VSA 0002
C     VERTICAL STRUCTURE ALGORITHM                                      VSA 0003
C                                                                       VSA 0004
C     FROM U.S. ARMY ATMOSPHERIC SCIENCES LAB                           VSA 0005
C     WHITE SANDS MISSILE RANGE, NM                                     VSA 0006
C                                                                       VSA 0007
C     CREATES A PROFILE OF AEROSOL DENSITY NEAR THE GROUND,INCLUDING    VSA 0008
C     CLOUDS AND FOG                                                    VSA 0009
C                                                                       VSA 0010
C     THESE PROFILES ARE AT 9 HEIGHTS BETWEEN 0 KM AND 2 KM             VSA 0011
C                                                                       VSA 0012
C                                                                       VSA 0013
C  ***VISIBILITY IS ASSUMED TO BE THE SURFACE VISIBILITY***             VSA 0014
C                                                                       VSA 0015
C     IHAZE  = THE TYPE OF AEROSOL                                      VSA 0016
C     VIS    = VISIBILITY IN KM AT THE SURFACE                          VSA 0017
C     CEILHT = THE CLOUD CEILING HEIGHT IN KM                           VSA 0018
C     DEPTH  = THE CLOUD/FOG DEPTH IN KM                                VSA 0019
C     ZINVHT = THE HEIGHT OF INVERSION OR BOUNDARY LAYER IN KM          VSA 0020
C                                                                       VSA 0021
C     VARIABLES USED IN VSA                                             VSA 0022
C                                                                       VSA 0023
C     ZC     = CLOUD CEILING HEIGHT IN M                                VSA 0024
C     ZT     = CLOUD DEPTH IN M                                         VSA 0025
C     ZINV   = INVERSION HEIGHT IN M                                    VSA 0026
C           SEE BELOW FOR MORE INFORMATION ABOUT ZC, ZT, AND ZINV       VSA 0027
C     D      = INITIAL EXTINCTION AT THE SURFACE (D=3.912/VIS-0.012)    VSA 0028
C     ZALGO  = THE DEPTH OF THE LAYER FOR THE ALGORITHM                 VSA 0029
C                                                                       VSA 0030
C     OUTPUT FROM VSA:                                                  VSA 0031
C                                                                       VSA 0032
C     Z      = HEIGHT IN KM                                             VSA 0033
C     RH     = RELATIVE HUMIDITY AT HEIGHT Z IN PERCENT                 VSA 0034
C     AHAZE  = EXTINCTION AT HEIGHT Z IN KM**-1                         VSA 0035
C     IH     = AEROSAL TYPE FOR HEIGHT Z                                VSA 0036
C     HMAX   = MAXIMUM HEIGHT IN KM USED IN VSA, NOT NECESSARILY 2.0 KM VSA 0037
C                                                                       VSA 0038
C                                                                       VSA 0039
C     THE SLANT PATH CALCULATION USES THE FOLLOWING FUNCTION:           VSA 0040
C                                                                       VSA 0041
C                 EXT55=A*EXP(B*EXP(C*Z))                               VSA 0042
C                                                                       VSA 0043
C     WHERE 'Z' IS THE HEIGHT IN KILOMETERS,                            VSA 0044
C           'A' IS A FUNCTION OF EXT55 AT Z=0.0 AND IS ALWAYS POSITIVE, VSA 0045
C           'B' AND 'C' ARE FUNCTIONS OF CLOUD CONDITIONS AND SURFACE   VSA 0046
C               VISIBILITY (EITHER A OR B CAN BE POSITIVE OR NEGATIVE), VSA 0047
C           'EXT55' IS THE VISIBILE EXTINCTION COEFFICIENT IN KM**-1.   VSA 0048
C                                                                       VSA 0049
C     THEREFORE, THERE ARE 4 CASES DEPENDING ON THE SIGNS OF 'B' AND 'C'VSA 0050
C     CEILHT AND ZINVHT ARE USED AS SWITCHES TO DETERMINE WHICH CASE    VSA 0051
C     TO USE.  THE SURFACE EXTINCTION 'D' IS CALCULATED FROM THE        VSA 0052
C     VISIBILITY USING  D=3.912/VIS-0.012 AS FOLLOWS-                   VSA 0053
C                                                                       VSA 0054
C         CASE=1  FOG/CLOUD CONDITIONS                                  VSA 0055
C                 'B' LT 0.0, 'C' LT 0.0                                VSA 0056
C                 'D' GE 7.0   KM**-1                                   VSA 0057
C                 FOR A CLOUD 7.    KM**-1 IS THE BOUNDARY VALUE AT     VSA 0058
C                 THE CLOUD BASE AND 'Z' IS THE VERTICAL DISTANCE       VSA 0059
C                 INTO THE CLOUD.                                       VSA 0060
C                 VARIABLE USED:   DEPTH                                VSA 0061
C                 ** DEFAULT:  DEPTH OF FOG/CLOUD IS 0.2 KM WHEN        VSA 0062
C                              'DEPTH' IS 0.0                           VSA 0063
C                                                                       VSA 0064
C             =2  CLOUD CEILING PRESENT                                 VSA 0065
C                 'B' GT 0.0, 'C' GT 0.0                                VSA 0066
C                 VARIABLE USED:   CEILHT (MUST BE GE 0.0)              VSA 0067
C                 ** DEFAULTS:  CASE 2 - CEILHT IS CALCULATED FROM      VSA 0068
C                               SURFACE EXTINCTION                      VSA 0069
C                                                                       VSA 0070
C             =3  RADIATION FOG OR INVERSION OR BOUNDARY LAYER PRESENT  VSA 0071
C                 'B' LT 0.0, 'C' GT 0.0                                VSA 0072
C                 VIS LE 2.0 KM DEFAULTS TO A RADIATION FOG AT THE      VSA 0073
C                     GROUND AND OVERRIDES INPUT BOUNDARY AEROSOL TYPE  VSA 0074
C                 VIS GT 2.0 KM FOR AN INVERSION OR BOUNDARY LAYER      VSA 0075
C                     WITH INPUT BOUNDARY AEROSOL TYPE                  VSA 0076
C                 ** IHAZE=9 (RADIATION FOG) ALWAYS DEFAULTS TO A       VSA 0077
C                    RADIATION FOG NO MATTER WHAT THE VISIBILITY IS.    VSA 0078
C                 SWITCH VARIABLE: CEILHT (MUST BE LT 0.0)              VSA 0079
C                 VARIABLE USED:   ZINVHT (MUST BE GE 0.0)              VSA 0080
C                 ** CEILHT MUST BE LT 0.0 FOR ZINVHT TO BE USED **     VSA 0081
C                    HOWEVER, IF DEPTH IS GT 0.0 AND ZINVHT IS EQ 0.0,  VSA 0082
C                    THE PROGRAM WILL SUBSTITUTE DEPTH FOR ZINVHT.      VSA 0083
C                 ** DEFAULT:  FOR A RADIATION FOG ZINVHT IS 0.05 K     VSA 0084
C                              FOR AN INVERSION LAYER ZINVHT IS 2.0 KM  VSA 0085
C                                                                       VSA 0086
C           NOTE: IF IHAZE = 9, BUT VIS GT 2.0 KM RECOMEND              VSA 0087
C           THAT IHAZE DEFAULT TO RURAL AEROSOL                         VSA 0088
C                                                                       VSA 0089
C             =4  NO CLOUD CEILING, INVERSION LAYER, OR BOUNDARY        VSA 0090
C                 LAYER PRESENT, I.E. CLEAR SKIES                       VSA 0091
C                 EXTINCTION PROFILE CONSTANT WITH HEIGHT A SHORT       VSA 0092
C                 DISTANCE ABOVE THE SURFACE                            VSA 0093
C                                                                       VSA 0094
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          VSA 0095
      DIMENSION Z(10),RH(10),AHAZE(10),IH(10)                           VSA 0096
      DIMENSION AA(2),CC(3),EE(4),A(2),B(2),C(2),FAC1(9),FAC2(9)        VSA 0097
      REAL KMTOM                                                        VSA 0098
      DATA AA/92.1,0.3981/,CC/-0.014,0.0125,-0.03 /,KMTOM/1000.0/       VSA 0099
C     THE LAST 3 VALUES OF EE BELOW ARE EXTINCTIONS FOR VISIBILITIES    VSA 0100
C     EQUAL TO 5.0, 23.0, AND 50.0 KM, RESPECTIVELY.                    VSA 0101
      DATA EE/7.0  ,0.7824,0.17009,0.012  /                             VSA 0102
      DATA FAC1/0.0,0.03,0.05,0.075,0.1,0.18,0.3,0.45,1.0/              VSA 0103
      DATA FAC2/0.0,0.03,0.1,0.18,0.3,0.45,0.6,0.78,1.0/                VSA 0104
      WRITE(IPR,599)                                                    VSA 0105
C                                                                       VSA 0106
C     UPPER LIMIT ON VERTICAL DISTANCE - 2 KM                           VSA 0107
      ZHIGH=2000.                                                       VSA 0108
      HMAX=ZHIGH                                                        VSA 0109
      IF(VIS.GT.0.0)GO TO 5                                             VSA 0110
C     DEFAULT FOR VISIBILITY DEPENDS ON THE VALUE OF IHAZE.             VSA 0111
      IF(IHAZE.EQ.8)VIS=0.2                                             VSA 0112
      IF(IHAZE.EQ.9)VIS=0.5                                             VSA 0113
      IF(IHAZE.EQ.2.OR.IHAZE.EQ.5)VIS=5.0                               VSA 0114
      IF(IHAZE.EQ.1.OR.IHAZE.EQ.4.OR.IHAZE.EQ.7)VIS=23.0                VSA 0115
      IF(IHAZE.EQ.6)VIS=50.0                                            VSA 0116
C     IF(IHAZE.EQ.3)VIS= OR IHAZE = 10 VIS IS DETERMINED ELSEWHERE      VSA 0117
   5  D=3.912/VIS-0.012                                                 VSA 0118
C                                                                       VSA 0119
      ZC=CEILHT*KMTOM                                                   VSA 0120
      ZT=DEPTH*KMTOM                                                    VSA 0121
      ZINV=ZINVHT*KMTOM                                                 VSA 0122
C     IHAZE=9 (RADIATION FOG) IS ALWAYS CALCULATED AS A RADIATION FOG.  VSA 0123
      IF(IHAZE.EQ.9)ZC=-1.0                                             VSA 0124
C     ALSO, CHECK TO SEE IF THE FOG DEPTH FOR A RADIATION FOG           VSA 0125
C     WAS INPUT TO DEPTH INSTEAD OF THE CORRECT VARIABLE ZINVHT.        VSA 0126
      IF(IHAZE.EQ.9.AND.ZT.GT.0.0.AND.ZINV.EQ.0.0)ZINV=ZT               VSA 0127
C                                                                       VSA 0128
C     'IC' DEFINES WHICH CASE TO USE.                                   VSA 0129
      IC=2                                                              VSA 0130
      IF(D.GE.EE(1).AND.ZC.GE.0.0)IC=1                                  VSA 0131
C                                                                       VSA 0132
      IF(ZC.LT.0.0.AND.IC.EQ.2)IC=3                                     VSA 0133
      IF(ZINV.LT.0.0.AND.IC.EQ.3)IC=4                                   VSA 0134
      K=1                                                               VSA 0135
      GO TO (10,20,40,50),IC                                            VSA 0136
C                                                                       VSA 0137
C     CASE 1:  DEPTH FOG/CLOUD; INCREASING EXTINCTION WITH HEIGHT FROM  VSA 0138
C              CLOUD/FOG BASE TO CLOUD/FOG TOP.                         VSA 0139
 10   CONTINUE                                                          VSA 0140
      IF(ZC.LT.HMAX.AND.IC.EQ.2)K=2                                     VSA 0141
C     IC=-1 WHEN A CLOUD IS PRESENT AND THE PATH GOES INTO IT.          VSA 0142
C     USE CASE 2 OR 2' BELOW CLOUD AND CASE 1 INSIDE IT.                VSA 0143
      IF(K.EQ.2)IC=(-1)                                                 VSA 0144
C     THE BASE OF THE CLOUD HAS AN EXTINCTION COEFFICIENT OF 7.0   KM-1.VSA 0145
      IF(K.EQ.2)D=EE(1)                                                 VSA 0146
      A(K)=AA(1)                                                        VSA 0147
C     IF THE SURFACE EXTINCTION IS GREATER THAN THE UPPER LIMIT OF 92.1 VSA 0148
C     KM**-1, RUN THE ALGORITHM WITH AN UPPER LIMIT OF 'D+10'.          VSA 0149
      IF(D.GE.AA(1))A(K)=D+10.0                                         VSA 0150
      C(K)=CC(1)                                                        VSA 0151
      IF(ZT.LE.0.0)WRITE(IPR  ,603)                                     VSA 0152
      IF(ZT.LE.0.0)WRITE(IPR  ,604)                                     VSA 0153
      IF(ZT.GT.0.0)WRITE(IPR  ,611)ZT                                   VSA 0154
C     IF THE DISTANCE FROM THE GROUND TO THE CLOUD/FOG TOP IS LESS      VSA 0155
C     THAN 2.0 KM, VSA WILL ONLY CALCULATE UP TO THE CLOUD TOP.         VSA 0156
      IF(ZT.LE.0.0)ZT=200.                                              VSA 0157
      HMAX=AMIN1(ZT+ZC,HMAX)                                            VSA 0158
      GO TO 60                                                          VSA 0159
C                                                                       VSA 0160
C     CASE 2:  CLEAR/HAZY/LIGHTLY FOGGY; INCREASING EXTINCTION WITH HEIGVSA 0161
C              UP TO THE CLOUD BASE.                                    VSA 0162
 20   A(K)=AA(2)                                                        VSA 0163
      E=EE(1)                                                           VSA 0164
      IF(ZC.EQ.0.0)WRITE(IPR  ,600)                                     VSA 0165
      IF(ZC.EQ.0.0)THEN                                                 VSA 0166
        EAK =  ALOG(E/A(K) )                                            VSA 0167
        DAK =  ALOG(D/A(K) )                                            VSA 0168
        ANUM = EAK / DAK                                                VSA 0169
        IF(ANUM . GT. 0) THEN                                           VSA 0170
              CEIL = ALOG(ANUM)/CC(2)                                   VSA 0171
         ELSE                                                           VSA 0172
               CEIL = 2000.                                             VSA 0173
          ENDIF                                                         VSA 0174
      ENDIF                                                             VSA 0175
CC    IF(ZC.EQ.0.0)CEIL=ALOG(ALOG(E/A(K))/(ALOG(D/A(K))))/CC(2)         VSA 0176
      IF(ZC.EQ.0.0)WRITE(IPR  ,602)CEIL                                 VSA 0177
      IF(ZC.GT.0.0)WRITE(IPR  ,610)ZC                                   VSA 0178
      IF(ZC.EQ.0.0)ZC=CEIL                                              VSA 0179
      F = (VIS * ZC/350.0)**2                                           VSA 0180
C                                                                       VSA 0181
C     F IS A SCALING FACTOR USED IN CASE 2                              VSA 0182
C                                                                       VSA 0183
      DF = D/F                                                          VSA 0184
      IF(DF .LT. 1.0E-5) THEN                                           VSA 0185
           A(K)=D - (D*D/(2.*F))                                        VSA 0186
      ELSE                                                              VSA 0187
           A(K) = F*(1.0 - EXP(-D/F))                                   VSA 0188
      ENDIF                                                             VSA 0189
C                                                                       VSA 0190
C     THE COEFFICIENT A IS RECALCULATED BASED UPON THE SCALING FACTOR   VSA 0191
C                                                                       VSA 0192
      GO TO 60                                                          VSA 0193
C                                                                       VSA 0194
C                                                                       VSA 0195
C     CASE 3:  NO CLOUD CEILING BUT A RADIATION FOG OR AN INVERSION     VSA 0196
C              OR BOUNDARY LAYER PRESENT; DECREASING EXTINCTION WITH    VSA 0197
C              HEIGHT UP TO THE HEIGHT OF THE FOG OR LAYER.             VSA 0198
 40   A(K)=D*1.1                                                        VSA 0199
      E=EE(3)                                                           VSA 0200
      IF(IHAZE.EQ.2.OR.IHAZE.EQ.5)E=EE(2)                               VSA 0201
      IF(IHAZE.EQ.6.OR.(VIS.GT.2.0.AND.IHAZE.NE.9))E=EE(4)              VSA 0202
      IF(E.GT.D)E=D*0.99999                                             VSA 0203
      IF(ZT.GT.0.0.AND.ZINV.EQ.0.0.AND.VIS.LE.2.0)ZINV=ZT               VSA 0204
      IF(ZINV.EQ.0.0.AND.VIS.GT.2.0.AND.IHAZE.NE.9)WRITE(IPR  ,601)     VSA 0205
      IF(ZINV.EQ.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9))WRITE(IPR  ,605)    VSA 0206
      IF(ZINV.EQ.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9))WRITE(IPR  ,604)    VSA 0207
      IF(ZINV.GT.0.0.AND.VIS.GT.2.0.AND.IHAZE.NE.9)WRITE(IPR  ,612)ZINV VSA 0208
      IF(ZINV.GT.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9))WRITE(IPR,614)ZINV  VSA 0209
      IF(ZINV.EQ.0.0.AND.VIS.GT.2.0.AND.IHAZE.NE.9)ZINV=2000            VSA 0210
      IF(ZINV.EQ.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9))ZINV= 50            VSA 0211
      HMAX=AMIN1(ZINV,HMAX)                                             VSA 0212
      ZC=0.0                                                            VSA 0213
      GO TO 60                                                          VSA 0214
C                                                                       VSA 0215
C     CASE 4:  NO CLOUD CEILING OR INVERSION LAYER;                     VSA 0216
C              CONSTANT EXTINCTION WITH HEIGHT.                         VSA 0217
C                                                                       VSA 0218
50     A(K) = EE(4)                                                     VSA 0219
       C(K) = CC(3)                                                     VSA 0220
C                                                                       VSA 0221
60               B(K)=ALOG(D/A(K))                                      VSA 0222
      IF(IC.EQ.2)C(K)=ALOG(ALOG(E/A(K))/B(K))/ZC                        VSA 0223
      IF(IC.EQ.3)C(K)=ALOG(ALOG(E/A(K))/B(K))/ZINV                      VSA 0224
      IF(ZC.LT.HMAX.AND.K.EQ.1.AND.IC.EQ.2)GO TO 10                     VSA 0225
      IF(IC.EQ.2)HMAX=AMIN1(ZC,HMAX)                                    VSA 0226
      ZALGO=HMAX                                                        VSA 0227
      IF(IC.LT.0)ZALGO=ZC                                               VSA 0228
      WRITE(IPR  ,619)                                                  VSA 0229
      IF(IC.LT.0)K=1                                                    VSA 0230
C                                                                       VSA 0231
      DO 70 I=1,9                                                       VSA 0232
      IF(IC.LT.0.AND.I.EQ.5)K=2                                         VSA 0233
      IF(IC.LT.0.AND.I.EQ.5)ZALGO=HMAX-ZC                               VSA 0234
      Z(I)=ZALGO*(1.0-FAC2(10-I))                                       VSA 0235
      IF(IC.EQ.1)Z(I)=ZALGO*FAC1(I)                                     VSA 0236
      IF(IC.EQ.4)Z(I)=ZALGO*FLOAT(I-1)/8.0                              VSA 0237
      IF(IC.LT.0.AND.I.LT.5)Z(I)=ZALGO*(1.0-FAC2(11-2*I))               VSA 0238
      IF(IC.LT.0.AND.I.GE.5)Z(I)=ZALGO*FAC1(2*I-9)                      VSA 0239
C     IF(IC.LT.0.AND.(I.EQ.7.OR.I.EQ.8))Z(I)=ZALGO*FAC1(2*I-10)         VSA 0240
                 AHAZE(I)=A(K)*EXP(B(K)*EXP(C(K)*Z(I)))                 VSA 0241
      IF(IC.LE.0.AND.I.GE.5)Z(I)=Z(I)+ZC                                VSA 0242
      Z(I)=Z(I)/KMTOM                                                   VSA 0243
      RH(I)=6.953*ALOG(AHAZE(I))+86.407                                 VSA 0244
      IF(AHAZE(I).GE.EE(1))RH(I)=100.0                                  VSA 0245
      VISIB=3.912/(AHAZE(I)+0.012)                                      VSA 0246
      IH(I)=IHAZE                                                       VSA 0247
C     IF A RADIATION FOG IS PRESENT (I.E. VIS<=2.0 KM AND IC=3),        VSA 0248
C     IH IS SET TO 9 FOR ALL LEVELS.                                    VSA 0249
      IF(VISIB.LE.2.0.AND.IC.EQ.3)IH(I)=9                               VSA 0250
C     FOR A DEPTH FOG/CLOUD CASE, IH=8 DENOTING AN ADVECTION FOG.       VSA 0251
      IF(IC.EQ.1.OR.(IC.LT.0.AND.I.GE.5))IH(I)=8                        VSA 0252
      WRITE(IPR  ,620)Z(I),RH(I),AHAZE(I),VISIB,IH(I)                   VSA 0253
   70 CONTINUE                                                          VSA 0254
      HMAX=HMAX/KMTOM                                                   VSA 0255
      RETURN                                                            VSA 0256
C                                                                       VSA 0257
599   FORMAT('0 VERTICAL STRUCTURE ALGORITHM (VSA) USED')               VSA 0258
600   FORMAT(1H ,50X,28HCLOUD CEILING HEIGHT UNKNOWN)                   VSA 0259
601   FORMAT(1H ,50X,42HINVERSION OR BOUNDARY LAYER HEIGHT UNKNOWN,/,   VSA 0260
     1  1H ,50X,39HVSA WILL USE A DEFAULT OF 2000.0 METERS,/)           VSA 0261
605   FORMAT(1H ,50X,27HRADIATION FOG DEPTH UNKNOWN)                    VSA 0262
619   FORMAT(5X,10HHEIGHT(KM),5X,7HR.H.(%),5X,16HEXTINCTION(KM-1),      VSA 0263
     1   5X,15HVIS(3.912/EXTN),5X,5HIHAZE,/)                            VSA 0264
620   FORMAT(7X,F7.4,7X,F5.1,8X,E12.4,11X,F7.4,10X,I2)                  VSA 0265
602   FORMAT(1H ,39X,35HVSA WILL USE A CALCULATED VALUE OF ,F7.1,       VSA 0266
     1       7H METERS,/)                                               VSA 0267
603   FORMAT(1H ,50X,19HCLOUD DEPTH UNKNOWN)                            VSA 0268
604   FORMAT(1H ,50X,38HVSA WILL USE A DEFAULT OF 200.0 METERS,/)       VSA 0269
610   FORMAT(1H ,50X,24HCLOUD CEILING HEIGHT IS ,F9.1,7H METERS,/)      VSA 0270
611   FORMAT(1H ,50X,15HCLOUD DEPTH IS ,F14.1,7H METERS,/)              VSA 0271
612   FORMAT(1H ,50X,38HINVERSION OR BOUNDARY LAYER HEIGHT IS ,F7.1,    VSA 0272
     1 7H METERS,/)                                                     VSA 0273
614   FORMAT(1H ,50X,26HDEPTH OF RADIATION FOG IS ,F7.1,7H METERS,/)    VSA 0274
613   FORMAT(1H ,50X,43HTHERE IS NO INVERSION OR BOUNDARY LAYER OR ,    VSA 0275
     1 13HCLOUD PRESENT,/)                                              VSA 0276
      END                                                               VSA 0277
