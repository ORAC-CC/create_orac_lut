      SUBROUTINE NEWSRC(H1SAV,H2SAV,ANGSAV,RNGSAV,BETASV,LENNSV,        NSC 0001
     1  IPARM,PARM1,PARM2,PSIPO,BETAH2)                                 NSC 0002
C                                                                       NSC 0003
C     IF THE OBSERVER IS ABOVE THE TOP OF THE DEFINED ATMOSPHERE,       NSC 0004
C     THIS ROUTINE DEFINES THE SOLAR/LUNAR ANGLES FOR THE POINT         NSC 0005
C     AT WHICH THE LINE-OF-SIGHT PATH ENTERS THE ATMOSPHERE.            NSC 0006
C                                                                       NSC 0007
C     DECLARE INPUTS/OUTPUTS                                            NSC 0008
C       H1SAV    OBSERVER ALTITUDE [KM]                                 NSC 0009
C       H2SAV    FINAL ALTITUDE [KM]                                    NSC 0010
C       ANGSAV   PATH ZENITH ANGLE AT THE OBSERVER [DEG]                NSC 0011
C       RNGSAV   PATH RANGE [KM]                                        NSC 0012
C       BETASV   EARTH CENTER ANGLE SUBTENDED BY PATH [DEG]             NSC 0013
C       LENNSV   FLAG EQUAL TO 1 FOR PATHS THROUGH A TANGENT HEIGHT     NSC 0014
C       IPARM    SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG                NSC 0015
C       PARM1    OBSERVER LATITUDE [DEG NORTH OF EQUATOR] (IPARM<2)     NSC 0016
C                RELATIVE SOLAR AZIMUTH [DEG EAST OF NORTH] (IPARM=2)   NSC 0017
C       PARM2    OBSERVER LONGITUDE [DEG WEST OF GREENWICH] (IPARM<2)   NSC 0018
C                SOLAR ZENITH ANGLE [DEG] (IPARM=2)                     NSC 0019
C       PSIPO    PATH TRUE AZIMUTH ANGLE [DEG EAST OF NORTH]            NSC 0020
C       BETAH2   EARTH CENTER ANGLE BETWEEN H1 AND H2                   NSC 0021
      REAL H1SAV,H2SAV,ANGSAV,RNGSAV,BETASV,PARM1,PARM2,PSIPO,BETAH2    NSC 0022
      INTEGER LENNSV,IPARM                                              NSC 0023
C                                                                       NSC 0024
C     INCLUDE PARAMETERS                                                NSC 0025
      INCLUDE 'PARAM.LST'                                               NSC 0026
C                                                                       NSC 0027
C     INCLUDE COMMONS                                                   NSC 0028
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     NSC 0029
      REAL TBOUND,SALB                                                  NSC 0030
      LOGICAL MODTRN                                                    NSC 0031
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   NSC 0032
     1  MODTRN                                                          NSC 0033
      INTEGER LENN                                                      NSC 0034
      REAL H1,H2,ANGLE,RANGE,BETA,REE                                   NSC 0035
      COMMON/CARD3/H1,H2,ANGLE,RANGE,BETA,REE,LENN                      NSC 0036
C                                                                       NSC 0037
C       PI       THE CONSTANT PI                                        NSC 0038
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       NSC 0039
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       NSC 0040
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         NSC 0041
      REAL PI,DEG,BIGNUM,BIGEXP                                         NSC 0042
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                NSC 0043
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     NSC 0044
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                NSC 0045
      REAL RE,ZMAX                                                      NSC 0046
      INTEGER IMAX,IMOD,IPATH                                           NSC 0047
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             NSC 0048
      REAL ZM,PM,TM,RFNDX,DENSTY                                        NSC 0049
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    NSC 0050
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               NSC 0051
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               NSC 0052
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           NSC 0053
      LOGICAL LSMALL                                                    NSC 0054
      COMMON/SMALL2/LSMALL                                              NSC 0055
C                                                                       NSC 0056
C     DECLARE LOCAL VARIABLES                                           NSC 0057
      REAL HMIN,PHI,STORE,ECA,RLAT,RLONG,OLAT,OLONG,                    NSC 0058
     1  PSISO,PSINEW,PRM1NW,PRM2NW                                      NSC 0059
      DOUBLE PRECISION DPH1,DPH2,DPANG,DPPHI,DPHMIN,DPBETA,DPRANG,DPBENDNSC 0060
      LOGICAL LOGLOS                                                    NSC 0061
      INTEGER ISLCT,IERROR,NPRSAV,ISSGSV,IAMT,LENGTH                    NSC 0062
C                                                                       NSC 0063
C     DEFINE DATA                                                       NSC 0064
      DATA LOGLOS/.FALSE./,ISLCT/0/                                     NSC 0065
C                                                                       NSC 0066
C     DEFINE THE TOP OF THE ATMOSPHERE AND EARTH RADIUS                 NSC 0067
      ZMAX=ZM(ML)                                                       NSC 0068
C                                                                       NSC 0069
C     RETURN IF H1 IS BELOW THE TOP OF THE ATMOSPHERE                   NSC 0070
C     UNLESS IMULT EQUALS -1 (IF IMULT EQUALS -1, THE                   NSC 0071
C     EARTH CENTER ANGLE TO H2 MUST BE DETERMINED).                     NSC 0072
      IF(H1.LE.ZMAX .AND. IMULT.NE.-1)RETURN                            NSC 0073
C                                                                       NSC 0074
C     CALL GEOINP TO DETERMINE THE SET (H1,H2,HMIN,ANGLE,PHI,LENN).     NSC 0075
      IERROR=0                                                          NSC 0076
      RE=REE                                                            NSC 0077
      IMOD=ML                                                           NSC 0078
      IMAX=ML                                                           NSC 0079
      NPRSAV=NPR                                                        NSC 0080
      NPR=2                                                             NSC 0081
      ISSGSV=ISSGEO                                                     NSC 0082
      ISSGEO=0                                                          NSC 0083
      CALL GEOINP(H1,H2,ANGLE,RANGE,BETA,ITYPE,                         NSC 0084
     1  LENN,HMIN,PHI,IERROR,ISLCT,LOGLOS)                              NSC 0085
      NPR=NPRSAV                                                        NSC 0086
      ISSGEO=ISSGSV                                                     NSC 0087
C                                                                       NSC 0088
C     CALCULATE THE EARTH CENTER ANGLE, ECA, BETWEEN THE ORIGINAL       NSC 0089
C     OBSERVER ALTITUDE, H1SAV, AND THE NEW OBSERVER ALTITUDE, H1.      NSC 0090
C     ASSUME NO REFRACTION ABOVE, H1, THE TOP OF THE ATMOSPHERE.        NSC 0091
      STORE=SIN(ANGLE/DEG)/(RE+H1SAV)                                   NSC 0092
      ECA=180.-(DEG*ASIN((RE+H1)*STORE)+ANGLE)                          NSC 0093
      IF(ECA.LT.0.)ECA=0.                                               NSC 0094
C                                                                       NSC 0095
C     RESET VALUES AND RETURN IF EITHER AN INCONSISTENCY WAS DETECTED   NSC 0096
C     OR THE EARTH CENTER ANGLE IS ZERO (NADIR PATH).                   NSC 0097
      IF(IERROR.NE.0)THEN                                               NSC 0098
          H1=H1SAV                                                      NSC 0099
          H2=H2SAV                                                      NSC 0100
          ANGLE=ANGSAV                                                  NSC 0101
          RANGE=RNGSAV                                                  NSC 0102
          BETA=BETASV                                                   NSC 0103
          LENN=LENNSV                                                   NSC 0104
          BETAH2=0.                                                     NSC 0105
          RETURN                                                        NSC 0106
      ENDIF                                                             NSC 0107
C                                                                       NSC 0108
C     IF IMULT = -1, DETERMINE THE EARTH CENTER ANGLE BETWEEN H1 AND H2 NSC 0109
      IF(IMULT.EQ.-1)THEN                                               NSC 0110
          IAMT=2                                                        NSC 0111
          LSMALL=.FALSE.                                                NSC 0112
          DPHMIN=DBLE(HMIN)                                             NSC 0113
          IF(ITYPE.EQ.3 .AND. H2SAV.NE.0.)THEN                          NSC 0114
C                                                                       NSC 0115
C             H2 IS THE TANGENT POINT (CASE 3B).                        NSC 0116
              DPH1=DPHMIN                                               NSC 0117
              DPH2=DBLE(H1)                                             NSC 0118
              DPANG=DBLE(90.)                                           NSC 0119
              DPPHI=DBLE(ANGLE)                                         NSC 0120
              LENGTH=0                                                  NSC 0121
          ELSE                                                          NSC 0122
C                                                                       NSC 0123
C             H2 IS THE FINAL ALTITUDE (NO MORE THAN ZMAX).             NSC 0124
              DPH1=DBLE(H1)                                             NSC 0125
              DPH2=DBLE(H2)                                             NSC 0126
              DPANG=DBLE(ANGLE)                                         NSC 0127
              DPPHI=DBLE(PHI)                                           NSC 0128
              LENGTH=LENN                                               NSC 0129
          ENDIF                                                         NSC 0130
          CALL DPRFPA(DPH1,DPH2,DPANG,DPPHI,LENGTH,                     NSC 0131
     1      DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)                           NSC 0132
          BETAH2=REAL(DPBETA)                                           NSC 0133
          IF(H1SAV.LE.ZMAX)THEN                                         NSC 0134
              H1=H1SAV                                                  NSC 0135
              H2=H2SAV                                                  NSC 0136
              ANGLE=ANGSAV                                              NSC 0137
              RANGE=RNGSAV                                              NSC 0138
              BETA=BETASV                                               NSC 0139
              LENN=LENNSV                                               NSC 0140
              RETURN                                                    NSC 0141
          ENDIF                                                         NSC 0142
      ENDIF                                                             NSC 0143
C                                                                       NSC 0144
C     SET H2 TO ITS ORIGINAL VALUE                                      NSC 0145
      H2=H2SAV                                                          NSC 0146
C                                                                       NSC 0147
C     SLANT PATH BRANCHING                                              NSC 0148
      IF(ITYPE.EQ.2)THEN                                                NSC 0149
C                                                                       NSC 0150
C         IF PATH GEOMETRY WAS DEFINED BY CASE 2C (H1,H2,RANGE) OR      NSC 0151
C         CASE 2D (H1,H2,BETA), CHANGE TO CASE 2A (H1,H2,ANGLE).        NSC 0152
C         IF PATH WAS DEFINED BY CASE 2B (H1,ANGLE,RANGE) DETERMINE     NSC 0153
C         THE NEW RANGE.                                                NSC 0154
          BETA=0.                                                       NSC 0155
          RANGE=0.                                                      NSC 0156
          IF(BETASV.EQ.0. .AND. RNGSAV.NE.0. .AND. ANGSAV.NE.0.)        NSC 0157
     1      RANGE=RNGSAV-SIN(ECA/DEG)/STORE                             NSC 0158
      ENDIF                                                             NSC 0159
C                                                                       NSC 0160
C     WRITE OUT CHANGES IN GEOMETRY.                                    NSC 0161
      WRITE(IPR,'(/A,//2A,/2A,/A,2F13.5)')                              NSC 0162
     1  ' OBSERVER ALTITUDE WAS LOWERED TO THE TOP OF THE ATMOSPHERE:', NSC 0163
     2  ' NAME                   DESCRIPTION                ',          NSC 0164
     3  '      INPUT      REVISION',                                    NSC 0165
     4  ' -----   ------------------------------------------',          NSC 0166
     5  '   ----------   ----------',                                   NSC 0167
     6  ' H1      OBSERVER ALTITUDE (KM)                    ',          NSC 0168
     7  H1SAV,H1                                                        NSC 0169
      IF(ITYPE.EQ.3)THEN                                                NSC 0170
          WRITE(IPR,'(2A,2F13.5)')' H2      ',                          NSC 0171
     1      'TANGENT ALTITUDE (KM)                     ',H2SAV,H2       NSC 0172
      ELSE                                                              NSC 0173
          WRITE(IPR,'(2A,2F13.5)')' H2      ',                          NSC 0174
     1      'FINAL ALTITUDE (KM)                       ',H2SAV,H2       NSC 0175
      ENDIF                                                             NSC 0176
      WRITE(IPR,'(3(A,2F13.5,/),A,2I13)')                               NSC 0177
     1  ' ANGLE   OBSERVER ZENITH (DEG)                     ',          NSC 0178
     2  ANGSAV,ANGLE,                                                   NSC 0179
     3  ' RANGE   SLANT RANGE (KM)                          ',          NSC 0180
     4  RNGSAV,RANGE,                                                   NSC 0181
     5  ' BETA    EARTH CENTER ANGLE (DEG)                  ',          NSC 0182
     6  BETASV,BETA,                                                    NSC 0183
     7  ' LENN    SHORT/LONG PATH FLAG                      ',          NSC 0184
     8  LENNSV,LENN                                                     NSC 0185
C                                                                       NSC 0186
C     REPLACE PATH GEOMETRY INPUTS                                      NSC 0187
      H1SAV=H1                                                          NSC 0188
      H2SAV=H2                                                          NSC 0189
      ANGSAV=ANGLE                                                      NSC 0190
      RNGSAV=RANGE                                                      NSC 0191
      BETASV=BETA                                                       NSC 0192
      LENNSV=LENN                                                       NSC 0193
C                                                                       NSC 0194
C     BRANCH BASED ON SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG           NSC 0195
      IF(IPARM.LT.2)THEN                                                NSC 0196
C                                                                       NSC 0197
C         OBSERVER LATITUDE, LONGITUDE, AND PATH AZIMUTH WERE           NSC 0198
C         INPUTS.  DETERMINE THEIR VALUES AT THE NEW H1.                NSC 0199
          CALL LOCATE(PARM1,PARM2,PSIPO,ECA,PRM1NW,PRM2NW)              NSC 0200
          CALL PSIECA(PRM1NW,PRM2NW,PARM1,PARM2,PSINEW,ECA)             NSC 0201
          PSINEW=PSINEW+180.                                            NSC 0202
          IF(PSINEW.GE.360.)PSINEW=PSINEW-360.                          NSC 0203
          WRITE(IPR,'((A,2F13.5))')                                     NSC 0204
     1      ' PARM1   OBSERVER LATITUDE (DEG NORTH OF EQUATOR)  ',      NSC 0205
     2      PARM1,PRM1NW,                                               NSC 0206
     3      ' PARM2   OBSERVER LONGITUDE (DEG WEST OF GREENWICH)',      NSC 0207
     4      PARM2,PRM2NW,                                               NSC 0208
     5      ' PSIPO   PATH AZIMUTH (DEG EAST OF NORTH)          ',      NSC 0209
     6      PSIPO,PSINEW                                                NSC 0210
          PARM1=PRM1NW                                                  NSC 0211
          PARM2=PRM2NW                                                  NSC 0212
          PSIPO=PSINEW                                                  NSC 0213
      ELSE                                                              NSC 0214
C                                                                       NSC 0215
C         RELATIVE SOLAR AZIMUTH AND SOLAR ZENITH ANGLES WERE INPUTS.   NSC 0216
C         ASSUME THAT THE OBSERVER AND SUN WERE ORIGINALLY BOTH ON THE  NSC 0217
C         EQUATOR, THAT THE OBSERVER WAS AT 0 DEG LONGITUDE, AND THAT   NSC 0218
C         THE SUN WAS WEST OF THE OBSERVER (I.E. AT PARM2 DEG LONGITUDE)NSC 0219
C         STEP 1:  DETERMINE THE ORIGINAL PATH AZIMUTH, (PSIPO = THE    NSC 0220
C                  TRUE SOLAR AZIMUTH - THE RELATIVE SOLAR AZIMUTH).    NSC 0221
          PSIPO=270.-PARM1                                              NSC 0222
C                                                                       NSC 0223
C         STEP 2:  DETERMINE THE NEW OBSERVER LATITUDE AND LONGITUDE.   NSC 0224
          OLAT=0.                                                       NSC 0225
          OLONG=0.                                                      NSC 0226
          CALL LOCATE(OLAT,OLONG,PSIPO,ECA,RLAT,RLONG)                  NSC 0227
C                                                                       NSC 0228
C         STEP 3:  DETERMINE THE NEW PATH AZIMUTH.                      NSC 0229
          CALL PSIECA(RLAT,RLONG,OLAT,OLONG,PSIPO,ECA)                  NSC 0230
          PSIPO=PSIPO+180.                                              NSC 0231
C                                                                       NSC 0232
C         STEP 4:  DETERMINE THE NEW TRUE SOLAR AZIMUTH AND SOLAR ZENITHNSC 0233
          CALL PSIECA(RLAT,RLONG,OLAT,PARM2,PSISO,PRM2NW)               NSC 0234
C                                                                       NSC 0235
C         STEP 5:  DETERMINE THE NEW RELATIVE SOLAR AZIMUTH             NSC 0236
          PRM1NW=PSISO-PSIPO                                            NSC 0237
          IF(PRM1NW.GT.180.)PRM1NW=PRM1NW-360.                          NSC 0238
          IF(PRM1NW.LE.-180.)PRM1NW=PRM1NW+360.                         NSC 0239
C                                                                       NSC 0240
C         STEP 6:  WRITE OUT OLD AND NEW VALUES.                        NSC 0241
          WRITE(IPR,'((A,2F13.5))')                                     NSC 0242
     1      ' PARM1   RELATIVE SOLAR AZIMUTH (DEG EAST OF NORTH)',      NSC 0243
     2      PARM1,PRM1NW,                                               NSC 0244
     3      ' PARM2   SOLAR ZENITH ANGLE (DEG)                  ',      NSC 0245
     4      PARM2,PRM2NW                                                NSC 0246
C                                                                       NSC 0247
C         STEP 7:  REPLACE OLD VALUES                                   NSC 0248
          PARM1=PRM1NW                                                  NSC 0249
          PARM2=PRM2NW                                                  NSC 0250
      ENDIF                                                             NSC 0251
C                                                                       NSC 0252
C     RETURN TO DRIVER                                                  NSC 0253
      RETURN                                                            NSC 0254
      END                                                               NSC 0255
