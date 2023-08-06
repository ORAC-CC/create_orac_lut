      SUBROUTINE FTRANG(H1,H2,RANGE,ANGLE,PHI,LENN,HMIN,IERROR)         FTR 0001
C                                                                       FTR 0002
C     GIVEN H1, H2 AND RANGE, THIS SUBROUTINE CALCULATES THE ZENITH     FTR 0003
C     ANGLE AT H1 (ANGLE) AND AT H2 (PHI).                              FTR 0004
C                                                                       FTR 0005
C     DECLARE ARGUMENTS:                                                FTR 0006
C       H1      OBSERVER ALTITUDE [KM].                                 FTR 0007
C       H2      FINAL ALTITUDE [KM].                                    FTR 0008
C       RANGE   SLANT PATH RANGE [KM].                                  FTR 0009
C       ANGLE   ZENITH ANGLE AT H1 TOWARDS H2 [DEG].                    FTR 0010
C       PHI     ZENITH ANGLE AT H2 TOWARDS H1 [DEG].                    FTR 0011
C       LENN    PATH LENGTH SWITCH (=0 SHORT PATH,                      FTR 0012
C                 =1 FOR PATH THROUGH TANGENT POINT).                   FTR 0013
C       HMIN    PATH MINIMUM ALTITUDE [KM].                             FTR 0014
C       IERROR  ERROR FLAG (=0 FOR SUCCESS, =1 FOR FAILURE).            FTR 0015
      DOUBLE PRECISION H1,H2,RANGE,ANGLE,PHI,HMIN                       FTR 0016
      INTEGER LENN,IERROR                                               FTR 0017
C                                                                       FTR 0018
C     LIST COMMONS:                                                     FTR 0019
C                                                                       FTR 0020
C     FILE UNIT NUMBERS                                                 FTR 0021
C       IRD      MODTRAN INPUT FILE, tape5, UNIT NUMBER (1).            FTR 0022
C       IPR      MODTRAN STANDARD OUTPUT FILE, tape6, UNIT NUMBER (2).  FTR 0023
C       IPU      MODTRAN SPECTRAL DATA FILE, tape7, UNIT NUMBER (7).    FTR 0024
C       NPR      PRINTOUT LEVEL SWITCH (1=small,0=normal,-1=large).     FTR 0025
C       IPR1     MODTRAN FLUX OUTPUT FILE, tape8, UNIT NUMBER (8).      FTR 0026
C       ISCRCH   MULTIPLE SCATTERING SCRATCH FILE UNIT NUMBER (10).     FTR 0027
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               FTR 0028
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           FTR 0029
      REAL RE,ZMAX                                                      FTR 0030
      INTEGER IMAX,IMOD,IPATH                                           FTR 0031
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             FTR 0032
      REAL GNDALT                                                       FTR 0033
      COMMON/GRAUND/GNDALT                                              FTR 0034
C                                                                       FTR 0035
C     DECLARE LOCAL VARIABLES:                                          FTR 0036
C       LOWER   FLAG, .TRUE. WHEN ANGLE INCREMENT NEEDS TO BE DECREASED.FTR 0037
C       DANGLE  ANGLE INCREMENT [DEG].                                  FTR 0038
C       HA      THE LOWER OF H1 AND H2 [KM].                            FTR 0039
C       HB      THE HIGHER OF H1 AND H2 [KM].                           FTR 0040
C       HBSAV   SAVED VALUE OF HB [KM].                                 FTR 0041
C       ANGLE1  CURRENT GUESS FOR ANGLE [DEG].                          FTR 0042
C       ITER    ITERATION COUNTER.                                      FTR 0043
      LOGICAL LOWER                                                     FTR 0044
      DOUBLE PRECISION DANGLE,HA,HB,HBSAV,RA,COEF,STORE,ANGLE1,         FTR 0045
     1  BETA,RANGE1,BENDNG,TERM1,TERM2,ADDRAN,ADDANG,ARG,DPRE           FTR 0046
      INTEGER ITER                                                      FTR 0047
C                                                                       FTR 0048
C     DECLARE DATA:                                                     FTR 0049
C       TOLRNC  MAXIMUM ALLOWED RANGE ERROR [KM].                       FTR 0050
C       ITERMX  MAXIMUM NUMBER OF ITERATIONS.                           FTR 0051
C       IAMT    COLUMN AMOUNT CALCULATION SWITCH.                       FTR 0052
      REAL TOLRNC                                                       FTR 0053
      INTEGER ITERMX,IAMT                                               FTR 0054
      DOUBLE PRECISION DPDEG                                            FTR 0055
      DATA TOLRNC/1.E-3/,ITERMX/30/,IAMT/2/,DPDEG/57.2957795131D0/      FTR 0056
C                                                                       FTR 0057
C     INITIALIZE PARAMETERS                                             FTR 0058
      LOWER=.FALSE.                                                     FTR 0059
      DANGLE=.1D0                                                       FTR 0060
      HA=H1                                                             FTR 0061
      HB=H2                                                             FTR 0062
      IF(H1.GT.H2)THEN                                                  FTR 0063
          HA=H2                                                         FTR 0064
          HB=H1                                                         FTR 0065
      ENDIF                                                             FTR 0066
C                                                                       FTR 0067
C     HBSAV IS DEFINED TO PROTECT AGAINST RPATH REDEFINING HB           FTR 0068
      HBSAV=HB                                                          FTR 0069
C                                                                       FTR 0070
C     GUESS AT ANGLE, INTEGRATE OVER PATH TO FIND RANGE, TEST           FTR 0071
C     FOR CONVERGENCE, AND ITERATE ANGLE IF NECESSARY.                  FTR 0072
      WRITE(IPR,'(///A,3F10.5,A,//A,//(2A))')                           FTR 0073
     1  ' CASE 2C: GIVEN H1, H2, RANGE:  (',H1,H2,RANGE,' )',           FTR 0074
     2  ' ITERATE AROUND ANGLE UNTIL RANGE CONVERGES',                  FTR 0075
     3  ' ITER    ANGLE      RANGE     DRANGE    ',                     FTR 0076
     4  '   BETA       HMIN        PHI    BENDING',                     FTR 0077
     5  '         (DEG)       (KM)      (KM)     ',                     FTR 0078
     6  '  (DEG)       (KM)      (DEG)     (DEG)'                       FTR 0079
C                                                                       FTR 0080
C     CALCULATE THE NO REFRACTION ZENITH ANGLE AS AN INITIAL GUESS      FTR 0081
      DPRE=DBLE(RE)                                                     FTR 0082
      RA=DPRE+HA                                                        FTR 0083
      COEF=.5D0/RA                                                      FTR 0084
      STORE=(HB-HA)*(DPRE+HB+RA)                                        FTR 0085
      ANGLE=DPDEG*ACOS(COEF*(STORE/RANGE-RANGE))                        FTR 0086
      ANGLE1=ANGLE                                                      FTR 0087
C                                                                       FTR 0088
C     BEGIN ITERATIVE PROCEDURE                                         FTR 0089
      ITER=0                                                            FTR 0090
   10 ITER=ITER+1                                                       FTR 0091
          IF(ITER.GT.ITERMX)THEN                                        FTR 0092
              WRITE(IPR,'(41h0FITRNG, CASE 2C (H1,H2,RANGE): SOLUTION , FTR 0093
     1          16hDID NOT CONVERGE,//10X,4hH1 =,F13.6,4x,4hH2 =,F13.6, FTR 0094
     2          4x,7hRANGE =,F13.6,4x,12hITERATIONS =,I5,//10X,5hLAST , FTR 0095
     3          9hITERATION,//10X,7hANGLE =,F16.9,4x,7hRANGE =,F16.9)') FTR 0096
     4          H1,H2,BETA,ITER,ANGLE1,RANGE1                           FTR 0097
              IERROR=1                                                  FTR 0098
              RETURN                                                    FTR 0099
          ENDIF                                                         FTR 0100
C                                                                       FTR 0101
C         DETERMINE RANGE1, THE RANGE CORRESPONDING TO ANGLE1           FTR 0102
          HB=HBSAV                                                      FTR 0103
          IF(ANGLE1.LE.90.)THEN                                         FTR 0104
C                                                                       FTR 0105
C             SHORT UPWARD PATH                                         FTR 0106
              LENN=0                                                    FTR 0107
              CALL DPFNMN(HA,ANGLE1,HB,LENN,ITER,HMIN,PHI,IERROR)       FTR 0108
              CALL DPRFPA(HA,HB,ANGLE1,PHI,LENN,                        FTR 0109
     1          HMIN,IAMT,BETA,RANGE1,BENDNG)                           FTR 0110
          ELSE                                                          FTR 0111
C                                                                       FTR 0112
C             LONG PATH                                                 FTR 0113
              LENN=1                                                    FTR 0114
              CALL DPFNMN(HA,ANGLE1,HB,LENN,ITER,HMIN,PHI,IERROR)       FTR 0115
              CALL DPRFPA(HA,HB,ANGLE1,PHI,LENN,                        FTR 0116
     1          HMIN,IAMT,BETA,RANGE1,BENDNG)                           FTR 0117
              IF(LENN.EQ.0)THEN                                         FTR 0118
C                                                                       FTR 0119
C                 THE PATH INTERSECTED THE EARTH.  DECREASE ANGLE BY .1 FTR 0120
                  WRITE(IPR,'(I3,7F11.5)')ITER,ANGLE1,                  FTR 0121
     1              RANGE1,RANGE-RANGE1,BETA,HMIN,PHI,BENDNG            FTR 0122
                  LOWER=.TRUE.                                          FTR 0123
                  ANGLE1=ANGLE1-DANGLE                                  FTR 0124
                  GOTO10                                                FTR 0125
              ENDIF                                                     FTR 0126
          ENDIF                                                         FTR 0127
C                                                                       FTR 0128
C         IF THE FINAL ALTITUDE HB HAS BEEN LOWERED (BECAUSE            FTR 0129
C         HB WAS ABOVE THE TOP OF THE ATMOSPHERE), ADD ON THE           FTR 0130
C         REMAINDER OF PATH LENGTH ASSUMING NO REFRACTION.              FTR 0131
          IF(HBSAV.GT.HB)THEN                                           FTR 0132
              TERM1=(DPRE+HB)*COS(PHI/DPDEG)                            FTR 0133
              TERM2=(HBSAV-HB)*(2*DPRE+HBSAV+HB)                        FTR 0134
              ADDRAN=TERM2/(SQRT(TERM1*TERM1+TERM2)-TERM1)              FTR 0135
              RANGE1=RANGE1+ADDRAN                                      FTR 0136
              ADDANG=DPDEG*ASIN(ADDRAN*SIN(PHI/DPDEG)/(DPRE+HBSAV))     FTR 0137
              BETA=BETA+ADDANG                                          FTR 0138
              PHI=PHI+ADDANG                                            FTR 0139
          ENDIF                                                         FTR 0140
          WRITE(IPR,'(I3,7F11.5)')                                      FTR 0141
     1      ITER,ANGLE1,RANGE1,RANGE-RANGE1,BETA,HMIN,PHI,BENDNG        FTR 0142
C                                                                       FTR 0143
C         CHECK FOR CONVERGENCE                                         FTR 0144
          IF(ABS(RANGE-RANGE1).LT.TOLRNC .OR.                           FTR 0145
     1      ABS(1.D0-RANGE1/RANGE).LT.2.D-6)THEN                        FTR 0146
              RANGE=RANGE1                                              FTR 0147
              IF(H1.LE.H2)THEN                                          FTR 0148
                  ANGLE=ANGLE1                                          FTR 0149
              ELSE                                                      FTR 0150
                  ANGLE=PHI                                             FTR 0151
                  PHI=ANGLE1                                            FTR 0152
              ENDIF                                                     FTR 0153
              IF(HMIN.GE.GNDALT)RETURN                                  FTR 0154
              WRITE(IPR,'(/A,//8X,2A)')                                 FTR 0155
     1          'FTRANG, CASE 2B (H1,H2,RANGE):  RANGE IS TOO LARGE.',  FTR 0156
     2          ' REFRACTED PATH TANGENT HEIGHT IS LESS THAN',          FTR 0157
     3          ' GROUND ALTITUDE; THE PATH INTERSECTS THE EARTH.'      FTR 0158
              IERROR=1                                                  FTR 0159
              RETURN                                                    FTR 0160
          ENDIF                                                         FTR 0161
C                                                                       FTR 0162
C         DETERMINE NEW VALUE FOR ANGLE (.6 IS A FUDGE FACTOR INTRODUCEDFTR 0163
C         TO AVOID OVER CORRECTING AND SPEED UP CONVERGENCE).           FTR 0164
          IF(RANGE1.LE.0.)STOP 'FTRANG error:  Range is zero or less.'  FTR 0165
          ARG=COEF*(STORE/RANGE1-RANGE1)                                FTR 0166
          IF(ARG.GE.1)THEN                                              FTR 0167
              ANGLE1=ANGLE1+.6D0*ANGLE                                  FTR 0168
          ELSEIF(ARG.GT.-1)THEN                                         FTR 0169
              ANGLE1=ANGLE1+.6D0*(ANGLE-DPDEG*ACOS(ARG))                FTR 0170
          ELSE                                                          FTR 0171
              ANGLE1=ANGLE1+.6D0*(ANGLE-180.D0)                         FTR 0172
          ENDIF                                                         FTR 0173
C                                                                       FTR 0174
C         CHECK IF ANGLE INCREMENT MUST BE LOWERED                      FTR 0175
          IF(LOWER)THEN                                                 FTR 0176
              DANGLE=.2D0*DANGLE                                        FTR 0177
              LOWER=.FALSE.                                             FTR 0178
          ENDIF                                                         FTR 0179
      GOTO10                                                            FTR 0180
      END                                                               FTR 0181

