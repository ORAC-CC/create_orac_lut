      SUBROUTINE CRPROF(CTHICK,CWDCOL,CIPCOL)                           CRP 0001
C                                                                       CRP 0002
C     THIS ROUTINE DEFINES CLOUD/RAIN MODEL PROFILES.                   CRP 0003
C                                                                       CRP 0004
C     LIST PARAMETERS                                                   CRP 0005
      INTEGER NMDLS                                                     CRP 0006
      PARAMETER(NMDLS=10)                                               CRP 0007
      INCLUDE 'PARAM.LST'                                               CRP 0008
C                                                                       CRP 0009
C     LIST ARGUMENTS                                                    CRP 0010
C     CTHICK   OUTPUT CLOUD THICKNESS [KM]                              CRP 0011
C     CWDCOL   OUTPUT WATER DROPLET VERTICAL COLUMN DENSITY [KM GM/M3]  CRP 0012
C     CIPCOL   OUTPUT ICE PARTICLE VERTICAL COLUMN DENSITY [KM GM/M3]   CRP 0013
      REAL CTHICK,CWDCOL,CIPCOL                                         CRP 0014
C                                                                       CRP 0015
C     LIST COMMONS                                                      CRP 0016
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               CRP 0017
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           CRP 0018
      REAL ZM,PM,TM,RFNDX,DENSTY                                        CRP 0019
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    CRP 0020
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               CRP 0021
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     CRP 0022
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                CRP 0023
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA                       CRP 0024
      REAL VIS,WSS,WHH,RAINRT                                           CRP 0025
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 CRP 0026
     1  VIS,WSS,WHH,RAINRT                                              CRP 0027
      INTEGER NCRALT,NCRSPC                                             CRP 0028
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    CRP 0029
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      CRP 0030
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       CRP 0031
      REAL ZCLD,CLD,CLDICE,RR                                           CRP 0032
      COMMON/CLDRR/ZCLD(1:NZCLD,0:1),CLD(1:NZCLD,0:5),                  CRP 0033
     1  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)                             CRP 0034
      INTEGER NBND                                                      CRP 0035
      REAL ZCLDRN(NZCLD),DRPWAT(NZCLD),PRTICE(NZCLD),RNPROF(NZCLD)      CRP 0036
      COMMON/CLDRN/NBND,ZCLDRN,DRPWAT,PRTICE,RNPROF                     CRP 0037
C                                                                       CRP 0038
C     LIST LOCAL VARIABLES AND ARRAYS                                   CRP 0039
      INTEGER ICLD0,ICE0,IRAIN0,IZ0,NUMALT,N,NBASE,NTOP,NBNDM1          CRP 0040
      REAL CBASE,FACLEN,CMIN,FAC,RNRATE(NZCLD)                          CRP 0041
C                                                                       CRP 0042
C     LIST DATA                                                         CRP 0043
      INTEGER MDLCLD(NMDLS)                                             CRP 0044
      CHARACTER*37 CLDTIT(NMDLS)                                        CRP 0045
      DATA MDLCLD/1,2,3,4,5,3,5,5,1,1/                                  CRP 0046
      DATA CLDTIT/'CUMULUS CLOUD                        ',              CRP 0047
     1            'ALTOSTRATUS CLOUD                    ',              CRP 0048
     2            'STRATUS CLOUD                        ',              CRP 0049
     3            'STRATUS/STRATO-CUMULUS CLOUD         ',              CRP 0050
     4            'NIMBOSTRATUS CLOUD                   ',              CRP 0051
     5            'STRATUS CLOUD WITH DRIZZLE           ',              CRP 0052
     6            'NIMBOSTRATUS CLOUD WITH LIGHT RAIN   ',              CRP 0053
     7            'NIMBOSTRATUS CLOUD WITH MODERATE RAIN',              CRP 0054
     8            'CUMULUS CLOUD WITH HEAVY RAIN        ',              CRP 0055
     9            'CUMULUS CLOUD WITH EXTREME RAIN      '/              CRP 0056
C                                                                       CRP 0057
C     RETURN IF CLOUD/RAIN MODEL 1 THROUGH NMDLS WAS NOT SELECTED       CRP 0058
      IF(ICLD.LT.1 .OR. ICLD.GT.NMDLS)RETURN                            CRP 0059
C                                                                       CRP 0060
C     ICLD0    INDEX OF CLOUD/RAIN MODEL WATER DROPLET DENSITY PROFILE  CRP 0061
C     ICE0     INDEX OF CLOUD/RAIN MODEL ICE PARTICLE DENSITY PROFILE   CRP 0062
C     IRAIN0   INDEX OF CLOUD/RAIN MODEL RAIN RATE PROFILE              CRP 0063
C     IZ0      INDEX OF CLOUD/RAIN MODEL ALTITUDE PROFILE               CRP 0064
C     NUMALT   NUMBER OF CLOUD/RAIN MODEL PROFILE BOUNDARY ALTITUDES    CRP 0065
      IF(NCRALT.LT.2)THEN                                               CRP 0066
C                                                                       CRP 0067
C         USE A BUILT-IN CLOUD/RAIN PROFILE DATA                        CRP 0068
          ICLD0=MDLCLD(ICLD)                                            CRP 0069
          ICE0=1                                                        CRP 0070
          IRAIN0=ICLD-5                                                 CRP 0071
          IZ0=1                                                         CRP 0072
          NUMALT=NZCLD                                                  CRP 0073
      ELSE                                                              CRP 0074
C                                                                       CRP 0075
C         READ USER-DEFINED CLOUD/RAIN PROFILE DATA                     CRP 0076
          IZ0=0                                                         CRP 0077
          ICLD0=0                                                       CRP 0078
          ICE0=0                                                        CRP 0079
          IRAIN0=0                                                      CRP 0080
          NUMALT=NCRALT                                                 CRP 0081
          CALL CRUPRF(NCRALT)                                           CRP 0082
      ENDIF                                                             CRP 0083
C                                                                       CRP 0084
C     NBASE    INDEX OF CLOUD BASE                                      CRP 0085
C     NTOP     INDEX OF CLOUD TOP                                       CRP 0086
      IF(CLD(NUMALT,ICLD0).NE.0. .OR. CLDICE(NUMALT,ICE0).NE.0.)THEN    CRP 0087
          WRITE(IPR,'(/2A,/14X,A)')' FATAL ERROR:  CLOUD',              CRP 0088
     1      ' WATER DROPLET AND ICE PARTICLE PROFILES ARE NOT',         CRP 0089
     2      ' BOTH ZERO AT THE HIGHEST CLOUD PROFILE ALTITUDE.'         CRP 0090
          STOP 'CLOUD PROFILES ARE NOT PROPERLY DEFINED AT TOP.'        CRP 0091
      ENDIF                                                             CRP 0092
      DO 10 NBASE=1,NUMALT                                              CRP 0093
          IF(CLD(NBASE,ICLD0).GT.0. .OR. CLDICE(NBASE,ICE0).GT.0.)GOTO20CRP 0094
   10 CONTINUE                                                          CRP 0095
      WRITE(IPR,'(/A)')                                                 CRP 0096
     1  ' FATAL ERROR:  NO WATER DROPLETS OR ICE PARTICLES IN CLOUD.'   CRP 0097
      STOP ' FATAL ERROR:  NO WATER DROPLETS OR ICE PARTICLES IN CLOUD.'CRP 0098
   20 CONTINUE                                                          CRP 0099
      DO 30 NTOP=NUMALT-1,NBASE+1,-1                                    CRP 0100
          IF(CLD(NTOP,ICLD0).GT.0. .OR. CLDICE(NTOP,ICE0).GT.0.)GOTO40  CRP 0101
   30 CONTINUE                                                          CRP 0102
      WRITE(IPR,'(/A)')                                                 CRP 0103
     1  ' FATAL ERROR:  CLOUD BASE AND TOP ALTITUDES ARE EQUAL.'        CRP 0104
      STOP ' FATAL ERROR:  CLOUD BASE AND TOP ALTITUDES ARE EQUAL.'     CRP 0105
   40 CONTINUE                                                          CRP 0106
C                                                                       CRP 0107
C     PRIORITY FOR INPUT RAIN RATES IS AS FOLLOWS:                      CRP 0108
C       IF "RAINRT" > 0, USE CONSTANT RAIN RATE FROM GROUND             CRP 0109
C       TO CLOUD TOP.  OTHERWISE, USE THE RAIN PROFILE.                 CRP 0110
      IF(RAINRT.GT.0. .OR. (NCRALT.LT.2 .AND. IRAIN0.LE.0))THEN         CRP 0111
          IF(RAINRT.LT.0.)RAINRT=0.                                     CRP 0112
          DO 50 N=1,NTOP                                                CRP 0113
              RNRATE(N)=RAINRT                                          CRP 0114
   50     CONTINUE                                                      CRP 0115
      ELSE                                                              CRP 0116
          DO 60 N=1,NTOP                                                CRP 0117
              RNRATE(N)=RR(N,IRAIN0)                                    CRP 0118
   60     CONTINUE                                                      CRP 0119
      ENDIF                                                             CRP 0120
      RNRATE(NTOP+1)=0.                                                 CRP 0121
C                                                                       CRP 0122
C     CBASE    CLOUD BASE ALTITUDE (RELATIVE TO SEA LEVEL) [KM]         CRP 0123
      CBASE=CALT                                                        CRP 0124
      IF(CBASE.LT.0.)CBASE=ZCLD(NBASE,IZ0)                              CRP 0125
      CBASE=CBASE+ZM(1)                                                 CRP 0126
C                                                                       CRP 0127
C     FACLEN   CLOUD STRETCH/COMPRESSION FACTOR                         CRP 0128
C     CTHICK   CLOUD THICKNESS [KM]                                     CRP 0129
      FACLEN=1.                                                         CRP 0130
      CTHICK=ZCLD(NTOP,IZ0)-ZCLD(NBASE,IZ0)                             CRP 0131
      IF(CTHIK.GT.0.)THEN                                               CRP 0132
          FACLEN=CTHIK/CTHICK                                           CRP 0133
          CTHICK=CTHIK                                                  CRP 0134
      ENDIF                                                             CRP 0135
C                                                                       CRP 0136
C     NBND     NUMBER OF CLOUD/RAIN PROFILE BOUNDARY ALTITUDES          CRP 0137
C     RNPROF   TRUE RAIN RATE PROFILE [MM/HR]                           CRP 0138
C     DRPWAT   CLOUD WATER DROPLET PROFILE [GM/M3]                      CRP 0139
C     PRTICE   CLOUD ICE PARTICLE PROFILE [GM/M3]                       CRP 0140
C     ZCLDRN   CLOUD/RAIN PROFILE ALTITUDES [KM]                        CRP 0141
C     CWDCOL   CLOUD WATER DROPLET VERTICAL COLUMN DENSITY [KM GM/M3]   CRP 0142
C     CIPCOL   CLOUD ICE PARTICLE VERTICAL COLUMN DENSITY [KM GM/M3]    CRP 0143
C     CMIN     ALTITUDE NON-ZERO CLOUD DENSITY BEGINS [KM]              CRP 0144
C                                                                       CRP 0145
C     RAIN/CLOUD PROFILES UP TO AND INCLUDING THE CLOUD BASE            CRP 0146
      RNPROF(1)=RNRATE(1)                                               CRP 0147
      NBND=1                                                            CRP 0148
      CWDCOL=0.                                                         CRP 0149
      IF(CBASE.GT.ZM(1))THEN                                            CRP 0150
C                                                                       CRP 0151
C         CLOUD BASE IS ABOVE THE GROUND                                CRP 0152
          ZCLDRN(1)=ZM(1)                                               CRP 0153
          DRPWAT(1)=0.                                                  CRP 0154
          PRTICE(1)=0.                                                  CRP 0155
          IF(NBASE.EQ.1)THEN                                            CRP 0156
C                                                                       CRP 0157
C             ASSUME CLOUD WATER DROPLET AND ICE PARTICLE               CRP 0158
C             DENSITIES DROP TO ZERO ONE METER BELOW THE CLOUD          CRP 0159
C             BASE SINCE THE FALL OFF WAS NOT EXPLICITLY DICTATED.      CRP 0160
              CMIN=CBASE-.001                                           CRP 0161
          ELSE                                                          CRP 0162
              CMIN=CBASE+FACLEN*(ZCLD(NBASE-1,IZ0)-ZCLD(NBASE,IZ0))     CRP 0163
          ENDIF                                                         CRP 0164
          NBND=2                                                        CRP 0165
          IF(CMIN.LE.ZM(1))THEN                                         CRP 0166
C                                                                       CRP 0167
C             CLOUD BASE IS NEAR THE GROUND.  ASSUME RAIN RATE          CRP 0168
C             VARIES LINEARLY FROM CLOUD BASE TO GROUND.                CRP 0169
              CMIN=ZM(1)                                                CRP 0170
          ELSEIF(RNRATE(1).LE.0.)THEN                                   CRP 0171
C                                                                       CRP 0172
C             CLOUD WITH NO RAIN.  BEGIN PROFILE AT CMIN.               CRP 0173
              ZCLDRN(1)=CMIN                                            CRP 0174
          ELSEIF(NBASE.GT.2)THEN                                        CRP 0175
C                                                                       CRP 0176
C             CLOUD WITH RAIN.  STRETCH/COMPRESS RAIN                   CRP 0177
C             PROFILE BETWEEN GROUND AND CMIN.                          CRP 0178
              FAC=(CMIN-ZM(1))/(ZCLD(NBASE-1,IZ0)-ZCLD(1,IZ0))          CRP 0179
              DO 70 NBND=2,NBASE-1                                      CRP 0180
                  ZCLDRN(NBND)=ZM(1)+FAC*(ZCLD(NBND,IZ0)-ZCLD(1,IZ0))   CRP 0181
                  DRPWAT(NBND)=0.                                       CRP 0182
                  PRTICE(NBND)=0.                                       CRP 0183
                  RNPROF(NBND)=RNRATE(NBND)                             CRP 0184
   70         CONTINUE                                                  CRP 0185
              NBND=NBASE                                                CRP 0186
          ENDIF                                                         CRP 0187
          RNPROF(NBND)=RNRATE(NBASE)                                    CRP 0188
          FAC=.5*(CBASE-ZCLDRN(NBND-1))                                 CRP 0189
          CWDCOL=FAC*CLD(NBASE,ICLD0)                                   CRP 0190
          CIPCOL=FAC*CLDICE(NBASE,ICE0)                                 CRP 0191
      ENDIF                                                             CRP 0192
      ZCLDRN(NBND)=CBASE                                                CRP 0193
      DRPWAT(NBND)=CLD(NBASE,ICLD0)                                     CRP 0194
      PRTICE(NBND)=CLDICE(NBASE,ICE0)                                   CRP 0195
C                                                                       CRP 0196
C     RAIN/CLOUD PROFILES ABOVE THE CLOUD BASE                          CRP 0197
      NBNDM1=NBND                                                       CRP 0198
      DO 80 N=NBASE+1,NTOP+1                                            CRP 0199
          NBND=NBND+1                                                   CRP 0200
          ZCLDRN(NBND)=CBASE+FACLEN*(ZCLD(N,IZ0)-ZCLD(NBASE,IZ0))       CRP 0201
          DRPWAT(NBND)=CLD(N,ICLD0)                                     CRP 0202
          PRTICE(NBND)=CLDICE(N,ICE0)                                   CRP 0203
          RNPROF(NBND)=RNRATE(N)                                        CRP 0204
          FAC=.5*(ZCLDRN(NBND)-ZCLDRN(NBNDM1))                          CRP 0205
          CWDCOL=CWDCOL+FAC*(DRPWAT(NBND)+DRPWAT(NBNDM1))               CRP 0206
          CIPCOL=CIPCOL+FAC*(PRTICE(NBND)+PRTICE(NBNDM1))               CRP 0207
   80 NBNDM1=NBND                                                       CRP 0208
C                                                                       CRP 0209
C     CHECK THAT CLOUD TOP IS NOT TOO HIGH                              CRP 0210
      IF(ZCLDRN(NBND).GT.ZM(ML))THEN                                    CRP 0211
          WRITE(IPR,'(/3(A,F8.3))')                                     CRP 0212
     1      ' FATAL ERROR:  TOP OF CLOUD PROFILE (',ZCLDRN(NBND),       CRP 0213
     2      'KM) IS ABOVE THE TOP OF THE ATMOSPHERE (',ZM(ML),'KM).'    CRP 0214
          STOP                                                          CRP 0215
      ENDIF                                                             CRP 0216
                                                                        CRP 0217
C                                                                       CRP 0218
C     CHECK INPUT CLOUD WATER DROPLET AND ICE PARTICLE COLUMN DENSITIES CRP 0219
      IF(CCOLWD.EQ.0. .AND. CCOLIP.EQ.0.)THEN                           CRP 0220
          WRITE(IPR,'(/2A,/(10X,A))')' WARNING: ',                      CRP 0221
     1      ' BOTH THE CLOUD WATER DROPLET AND ICE PARTICLE COLUMN',    CRP 0222
     2      ' DENSITIES WERE INPUT AS ZERO.  IT IS ASSUMED THAT THE',   CRP 0223
     3      ' INTENT OF THE USER WAS TO NOT SCALE THE CLOUD COLUMN',    CRP 0224
     4      ' DENSITIES (I.E. TO USE DEFAULTS).  THEREFORE, THESE',     CRP 0225
     5      ' INPUTS (CCOLWD AND CCOLIP) HAVE BEEN RESET TO -1.000'     CRP 0226
          CCOLWD=-1.                                                    CRP 0227
          CCOLIP=-1.                                                    CRP 0228
      ENDIF                                                             CRP 0229
C                                                                       CRP 0230
C     SCALE WATER DROPLET PROFILE TO INPUT VERTICAL COLUMN DENSITY      CRP 0231
      IF(CCOLWD.GE.0. .AND. CWDCOL.GT.0.)THEN                           CRP 0232
          FAC=CCOLWD/CWDCOL                                             CRP 0233
          CWDCOL=CCOLWD                                                 CRP 0234
          DO 90 N=1,NBND                                                CRP 0235
              DRPWAT(N)=FAC*DRPWAT(N)                                   CRP 0236
   90     CONTINUE                                                      CRP 0237
      ENDIF                                                             CRP 0238
C                                                                       CRP 0239
C     SCALE ICE PARTICLE PROFILE TO INPUT VERTICAL COLUMN DENSITY       CRP 0240
      IF(CCOLIP.GE.0. .AND. CIPCOL.GT.0.)THEN                           CRP 0241
          FAC=CCOLIP/CIPCOL                                             CRP 0242
          CIPCOL=CCOLIP                                                 CRP 0243
          DO 100 N=1,NBND                                               CRP 0244
              PRTICE(N)=FAC*PRTICE(N)                                   CRP 0245
  100     CONTINUE                                                      CRP 0246
      ENDIF                                                             CRP 0247
C                                                                       CRP 0248
C     RETURN IF CLOUD/RAIN PROFILES ARE NOT BE TO OUTPUT                CRP 0249
      IF(NPR.GE.0)RETURN                                                CRP 0250
C                                                                       CRP 0251
C     WRITE OUT CLOUD/RAIN PROFILES                                     CRP 0252
      IF(NCRALT.GE.2)THEN                                               CRP 0253
          WRITE(IPR,'(//A,4X,A)')'1','USER-DEFINED CLOUD/RAIN PROFILES' CRP 0254
      ELSEIF(RAINRT.LE.0.)THEN                                          CRP 0255
          WRITE(IPR,'(//A,4X,A)')'1',CLDTIT(ICLD)                       CRP 0256
      ELSE                                                              CRP 0257
          N=INDEX(CLDTIT(ICLD0),'CLOUD')+4                              CRP 0258
          WRITE(IPR,'(//A,4X,2A)')                                      CRP 0259
     1      '1',CLDTIT(ICLD0)(1:N),' WITH CONSTANT RAIN RATE'           CRP 0260
      ENDIF                                                             CRP 0261
      WRITE(IPR,'(5X,2A,F12.5,A)')'(CLOUD WATER DROPLET',               CRP 0262
     1  ' VERTICAL COLUMN DENSITY:',CWDCOL,' KM GM/M3)'                 CRP 0263
      WRITE(IPR,'(5X,2A,F12.5,A)')'(CLOUD ICE PARTICLE ',               CRP 0264
     1  ' VERTICAL COLUMN DENSITY:',CIPCOL,' KM GM/M3)'                 CRP 0265
      WRITE(IPR,'(4(/A),/(I5,4F14.5))')                                 CRP 0266
     1  ' BOUNDARY                 WATER         ICE',                  CRP 0267
     2  ' LAYER                    DROPLET       PARTICLE        RAIN', CRP 0268
     3  ' NUMBER    ALTITUDE       DENSITY       DENSITY         RATE', CRP 0269
     4  '             (KM)         (GM/M3)       (GM/M3)       (MM/HR)',CRP 0270
     5  (N,ZCLDRN(N),DRPWAT(N),PRTICE(N),RNPROF(N),N=1,NBND)            CRP 0271
      WRITE(IPR,'(/A,//)')' END OF CLOUD/RAIN PROFILES'                 CRP 0272
      RETURN                                                            CRP 0273
      END                                                               CRP 0274
