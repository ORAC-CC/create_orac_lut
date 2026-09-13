      SUBROUTINE GEOINP(SPH1,SPH2,SPANGL,SPRANG,SPBETA,                 GIN 0001
     1  ITYPE,LENN,SPHMIN,SPPHI,IERROR,ISLCT,LOGLOS)                    GIN 0002
C                                                                       GIN 0003
C     THIS ROUTINE INTERPRETS THE ALLOWABLE COMBINATIONS OF INPUT PATH  GIN 0004
C     PARAMETERS INTO THE STANDARD SET:  H1, H2, ANGLE, PHI, HMIN AND   GIN 0005
C     LENN.  THE ALLOWABLE COMBINATIONS OF INPUT PARAMETERS ARE         GIN 0006
C                                                                       GIN 0007
C          FOR ITYPE = 2  (SLANT PATH H1 TO H2)                         GIN 0008
C              A.  H1, H2 AND ANGLE;                                    GIN 0009
C              B.  H1, ANGLE AND RANGE;                                 GIN 0010
C              C.  H1, H2 AND RANGE; OR                                 GIN 0011
C              D.  H1, H2, AND BETA.                                    GIN 0012
C                                                                       GIN 0013
C          FOR ITYPE = 3  (SLANT PATH H1 TO SPACE/GROUND)               GIN 0014
C              A. H1 AND ANGLE; OR                                      GIN 0015
C              B. H1 AND HMIN (INPUT AS H2).                            GIN 0016
C                                                                       GIN 0017
C     THIS ROUTINE ALSO DETECTS BAD INPUT (IMPOSSIBLE GEOMETRY) AND     GIN 0018
C     ITYPE = 2 CASES WHICH INTERSECT THE EARTH, AND RETURNS THESE      GIN 0019
C     CASES WITH ERROR FLAGS.                                           GIN 0020
C                                                                       GIN 0021
C     DECLARE ARGUMENTS:                                                GIN 0022
C       SPH1     ALTITUDE OF OBSERVER [KM].                             GIN 0023
C       SPH2     FINAL OR TANGENT ALTITUDE [KM].                        GIN 0024
C       SPANGL   PATH ZENITH ANGLE FROM OBSERVER TOWARDS SPH2 [DEG].    GIN 0025
C       SPRANG   SLANT PATH RANGE [KM].                                 GIN 0026
C       SPBETA   SLANT PATH EARTH CENTER ANGLE [DEG].                   GIN 0027
C       ITYPE    FLAG FOR GEOMETRY TYPE.                                GIN 0028
C       LENN     LENGTH FLAG USED FOR CASE 2A WHEN H2<H1 AND ANGLE>90   GIN 0029
C                (=0 FOR SHORT PATH, =1 FOR PATH THROUGH TANGENT POINT).GIN 0030
C       SPHMIN   SLANT PATH MINIMUM ALTITUDE [KM].                      GIN 0031
C       SPPHI    PATH ZENITH ANGLE FROM FINAL ALTITUDE TOWARDS          GIN 0032
C                OBSERVER [DEG].                                        GIN 0033
C       IERROR   ERROR FLAG (=0 FOR ACCEPTABLE GEOMETRY INPUTS).        GIN 0034
C       ISLCT    ????                                                   GIN 0035
C       LOGLOS   LOGICAL FLAG(.TRUE. FOR ITYPE=2 WITH NO SUN).          GIN 0036
      REAL SPH1,SPH2,SPANGL,SPRANG,SPBETA,SPHMIN,SPPHI                  GIN 0037
      INTEGER ITYPE,LENN,IERROR,ISLCT                                   GIN 0038
      LOGICAL LOGLOS                                                    GIN 0039
C                                                                       GIN 0040
C     LIST COMMONS:                                                     GIN 0041
C                                                                       GIN 0042
C     FILE UNIT NUMBERS                                                 GIN 0043
C       IRD      MODTRAN INPUT FILE, tape5, UNIT NUMBER (1).            GIN 0044
C       IPR      MODTRAN STANDARD OUTPUT FILE, tape6, UNIT NUMBER (2).  GIN 0045
C       IPU      MODTRAN SPECTRAL DATA FILE, tape7, UNIT NUMBER (7).    GIN 0046
C       NPR      PRINTOUT LEVEL SWITCH (1=small,0=normal,-1=large).     GIN 0047
C       IPR1     MODTRAN FLUX OUTPUT FILE, tape8, UNIT NUMBER (8).      GIN 0048
C       ISCRCH   MULTIPLE SCATTERING SCRATCH FILE UNIT NUMBER (10).     GIN 0049
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               GIN 0050
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           GIN 0051
      REAL RE,ZMAX                                                      GIN 0052
      INTEGER IMAX,IMOD,IPATH                                           GIN 0053
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             GIN 0054
      REAL GNDALT                                                       GIN 0055
      COMMON/GRAUND/GNDALT                                              GIN 0056
C                                                                       GIN 0057
C     DECLARE LOCAL VARIABLES:                                          GIN 0058
      DOUBLE PRECISION H1,H2,ANGLE,RANGE,BETA,HMIN,PHI,H2ST,STORE       GIN 0059
C                                                                       GIN 0060
C     DEFINE DATA                                                       GIN 0061
C       ITER    COUNTER USED IN CALL TO DPFNMN.                         GIN 0062
      INTEGER ITER                                                      GIN 0063
      DATA ITER/0/                                                      GIN 0064
C                                                                       GIN 0065
C     DEFINE DOUBLE PRECISION GEOMETRY PARAMETERS:                      GIN 0066
      H1=DBLE(SPH1)                                                     GIN 0067
      H2=DBLE(SPH2)                                                     GIN 0068
      ANGLE=DBLE(SPANGL)                                                GIN 0069
      RANGE=DBLE(SPRANG)                                                GIN 0070
      BETA=DBLE(SPBETA)                                                 GIN 0071
      HMIN=DBLE(SPHMIN)                                                 GIN 0072
      PHI=DBLE(SPPHI)                                                   GIN 0073
C                                                                       GIN 0074
C     BRANCH BASED ON PATH TYPE.                                        GIN 0075
      IF(ITYPE.EQ.3)THEN                                                GIN 0076
C                                                                       GIN 0077
C         SLANT PATH TO SPACE.                                          GIN 0078
          IF(H2.EQ.0.)THEN                                              GIN 0079
C                                                                       GIN 0080
C             CASE 3A:  H1 AND ANGLE.                                   GIN 0081
              IF(NPR.LT.1)WRITE(IPR,'(//A)')                            GIN 0082
     1          ' CASE 3A:  GIVEN H1,H2=SPACE,ANGLE'                    GIN 0083
              H2=DBLE(ZMAX)                                             GIN 0084
              IF(ANGLE.GT.90.)LENN=1                                    GIN 0085
              CALL DPFNMN(H1,ANGLE,H2,LENN,ITER,HMIN,PHI,IERROR)        GIN 0086
          ELSE                                                          GIN 0087
C                                                                       GIN 0088
C             CASE 3B:  H1 AND HMIN                                     GIN 0089
              IF(NPR.LT.1)WRITE(IPR,'(//A)')                            GIN 0090
     1          ' CASE 3B:  GIVEN H1, HMIN, H2=SPACE'                   GIN 0091
              HMIN=H2                                                   GIN 0092
              H2=DBLE(ZMAX)                                             GIN 0093
              IF(H1.LT.HMIN)THEN                                        GIN 0094
                  WRITE(IPR,'(/2A,//10X,2(A,F13.6))')' GEOINP,',        GIN 0095
     1              ' CASE 3B (H1,HMIN,SPACE):  ERROR IN INPUT DATA',   GIN 0096
     2              'H1 =',H1,'    IS LESS THAN HMIN =',HMIN            GIN 0097
                  IERROR=1                                              GIN 0098
              ELSE                                                      GIN 0099
                  STORE=90.D0                                           GIN 0100
                  CALL DPFNMN(HMIN,STORE,H1,LENN,ITER,HMIN,ANGLE,IERROR)GIN 0101
                  STORE=90.D0                                           GIN 0102
                  CALL DPFNMN(HMIN,STORE,H2,LENN,ITER,HMIN,PHI,IERROR)  GIN 0103
                  LENN=1                                                GIN 0104
                  IF(HMIN.EQ.H1)LENN=0                                  GIN 0105
              ENDIF                                                     GIN 0106
          ENDIF                                                         GIN 0107
      ELSEIF(ITYPE.EQ.2)THEN                                            GIN 0108
C                                                                       GIN 0109
C         SLANT PATH TO BETWEEN ALTITUDES.                              GIN 0110
          IF((LOGLOS .AND. ISLCT.EQ.21) .OR.                            GIN 0111
     1      (.NOT.LOGLOS .AND. RANGE.LE.0. .AND. BETA.LE.0.))THEN       GIN 0112
C                                                                       GIN 0113
C             CASE 2A:  H1, H2, ANGLE                                   GIN 0114
              IF(NPR.LT.1)WRITE(IPR,'(//A)')                            GIN 0115
     1          ' CASE 2A:  GIVEN H1, H2 AND ANGLE'                     GIN 0116
              IF(H1.GE.H2 .AND. ANGLE.LE.90.)THEN                       GIN 0117
                  WRITE(IPR,'(/2A,//(10X,2(A,F13.6),A))')' GEOINP,',    GIN 0118
     1              ' CASE 2A (H1,H2,ANGLE):  ERROR IN INPUT DATA',     GIN 0119
     2              'H1 (=',H1,'KM) IS GREATER THAN OR EQUAL TO H2 (= ',GIN 0120
     3              H2,'KM)','AND ANGLE (=',ANGLE,                      GIN 0121
     4              'DEG) DOES NOT EXCEED 90DEG.'                       GIN 0122
                  IERROR=1                                              GIN 0123
              ELSEIF(H1.LE.GNDALT .AND. ANGLE.GT.90.)THEN               GIN 0124
                  WRITE(IPR,'(/2A)')                                    GIN 0125
     1              ' GEOINP, ITYPE = 2: SLANT PATH INTERSECTS',        GIN 0126
     1              ' THE EARTH OR GNDALT, AND CANNOT REACH H2.'        GIN 0127
                  IERROR=1                                              GIN 0128
              ELSE                                                      GIN 0129
                  IF(H2.LT.H1 .AND. ANGLE.GT.90. .AND. NPR.LT.1)        GIN 0130
     1              WRITE(IPR,'(//3A,I3)')' EITHER A SHORT PATH',       GIN 0131
     2              ' (LENN=0) OR A LONG PATH THROUGH A TANGENT',       GIN 0132
     3              ' HEIGHT (LENN=1) IS POSSIBLE: LENN = ',LENN        GIN 0133
                  H2ST=H2                                               GIN 0134
                  CALL DPFNMN(H1,ANGLE,H2,LENN,ITER,HMIN,PHI,IERROR)    GIN 0135
                  IF(H2.NE.H2ST)THEN                                    GIN 0136
                      WRITE(IPR,'(/2A)')                                GIN 0137
     1                  ' GEOINP, ITYPE = 2: SLANT PATH INTERSECTS',    GIN 0138
     2                  ' THE EARTH OR GNDALT, AND CANNOT REACH H2.'    GIN 0139
                      IERROR=1                                          GIN 0140
                  ENDIF                                                 GIN 0141
              ENDIF                                                     GIN 0142
          ELSEIF((LOGLOS .AND. ISLCT.EQ.24) .OR.                        GIN 0143
     1      (.NOT.LOGLOS .AND. BETA.GT.0.))THEN                         GIN 0144
C                                                                       GIN 0145
C             CASE 2D:  H1, H2, BETA                                    GIN 0146
              CALL FDBETA(H1,H2,BETA,ANGLE,PHI,LENN,HMIN,IERROR)        GIN 0147
          ELSEIF((LOGLOS .AND. ISLCT.EQ.22) .OR.                        GIN 0148
     1      (.NOT.LOGLOS .AND. ANGLE.GT.0.))THEN                        GIN 0149
C                                                                       GIN 0150
C             CASE 2B:  H1, ANGLE, RANGE                                GIN 0151
              CALL NEWH2(H1,H2,ANGLE,RANGE,BETA,LENN,HMIN,PHI)          GIN 0152
              IF(ANGLE.GT.90. .AND.PHI.GT.90.)LENN=1                    GIN 0153
              CALL DPFNMN(H1,ANGLE,H2,LENN,ITER,HMIN,PHI,IERROR)        GIN 0154
          ELSE                                                          GIN 0155
C                                                                       GIN 0156
C             CASE 2C:  H1, H2, RANGE                                   GIN 0157
              CALL FTRANG(H1,H2,RANGE,ANGLE,PHI,LENN,HMIN,IERROR)       GIN 0158
          ENDIF                                                         GIN 0159
      ELSE                                                              GIN 0160
          WRITE(IPR,'(/2A,I10)')' GEOINP:  ERROR IN INPUT DATA,',       GIN 0161
     1      ' ITYPE NOT EQUAL TO 2 OR 3.   ITYPE =',ITYPE               GIN 0162
          IERROR=1                                                      GIN 0163
      ENDIF                                                             GIN 0164
C                                                                       GIN 0165
C     TEST IERROR AND RECHECK LENN                                      GIN 0166
      IF(IERROR.EQ.0)THEN                                               GIN 0167
          LENN=0                                                        GIN 0168
          IF(ABS(HMIN-MIN(H1,H2)).GT..00005)LENN=1                      GIN 0169
C                                                                       GIN 0170
C         REDUCE PATH END POINTS ABOVE ZMAX TO ZMAX                     GIN 0171
          IF(HMIN.GE.ZMAX)THEN                                          GIN 0172
              WRITE(IPR,'(/2A,//4(4X,A,F11.5))')                        GIN 0173
     1          ' GEOINP:  THE ENTIRE PATH LIES ABOVE ZMAX,',           GIN 0174
     2          ' THE TOP OF THE ATMOSPHERIC PROFILE',                  GIN 0175
     3          'ZMAX =',ZMAX,'H1 =',H1,'H2 =',H2,'HMIN =',HMIN         GIN 0176
              IERROR=1                                                  GIN 0177
          ELSE                                                          GIN 0178
              IF(H1.GT.ZMAX .OR. H2.GT.ZMAX)CALL REDUCE(H1,H2,ANGLE,PHI)GIN 0179
              IF(NPR.LT.1)WRITE(IPR,'(//A,/5(/10X,A,F11.5,A),           GIN 0180
     1          /10X,A,I11)')' SLANT PATH PARAMETERS IN STANDARD FORM', GIN 0181
     2          'H1      =',H1   ,' KM' ,'H2      =',H2 ,' KM' ,        GIN 0182
     3          'ANGLE   =',ANGLE,' DEG','PHI     =',PHI,' DEG',        GIN 0183
     4          'HMIN    =',HMIN ,' KM' ,'LENN    =',LENN               GIN 0184
          ENDIF                                                         GIN 0185
      ENDIF                                                             GIN 0186
C                                                                       GIN 0187
C     DEFINE SINGLE PRECISION VALUES AND RETURN                         GIN 0188
      SPH1=REAL(H1)                                                     GIN 0189
      SPH2=REAL(H2)                                                     GIN 0190
      SPANGL=REAL(ANGLE)                                                GIN 0191
      SPRANG=REAL(RANGE)                                                GIN 0192
      SPBETA=REAL(BETA)                                                 GIN 0193
      SPHMIN=REAL(HMIN)                                                 GIN 0194
      SPPHI=REAL(PHI)                                                   GIN 0195
      RETURN                                                            GIN 0196
      END                                                               GIN 0197
