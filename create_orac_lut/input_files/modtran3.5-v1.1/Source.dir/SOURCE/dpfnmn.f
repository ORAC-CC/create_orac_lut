      SUBROUTINE DPFNMN(H1,ANGLE,H2,LENN,ITER,HMIN,PHI,IERROR)          FMN 0001
      INCLUDE 'PARAM.LST'                                               FMN 0002
C***********************************************************************FMN 0003
C     DOUBLE PRECISION VERSION OF THE FORMER ROUTINE FNDHMN.            FMN 0004
C                                                                       FMN 0005
C     THIS SUBROUTINE CALCULATES THE MINIMUM ALTITUDE HMIN ALONG        FMN 0006
C     THE REFRACTED PATH AND THE FINAL ZENITH ANGLE PHI.                FMN 0007
C     THE PARAMETER LENN INDICATES WHETHER THE PATH GOES THROUGH        FMN 0008
C     A TANGENT HEIGHT (LENN=1) OR NOT (LENN=0).  IF ANGLE > 90 AND     FMN 0009
C     H1 > H2, THEN LENN CAN EITHER BE 1 OR 0, AND THE CHOICE IS        FMN 0010
C     LEFT TO THE USER.                                                 FMN 0011
C     THE (INDEX OF REFRACTION - 1.0) IS MODELED AS AN EXPONENTIAL      FMN 0012
C     BETWEEN THE LAYER BOUNDARIES, WITH A SCALE HEIGHT SH AND AN       FMN 0013
C     AMOUNT AT THE GROUND GAMMA.                                       FMN 0014
C     CPATH IS THE REFRACTIVE CONSTANT FOR THIS PATH AND                FMN 0015
C     EQUALS  INDEX(H1)*(RE+H1)*SIN(ANGLE).                             FMN 0016
C***********************************************************************FMN 0017
      INTEGER IRD,IPR,IPU,NPR,IPR1,LENN,ITER,IERROR,N                   FMN 0018
      DOUBLE PRECISION H1,ANGLE,H2,HMIN,PHI,CRFRCT,H,DPANDX,SH,TANHT,   FMN 0019
     1  GAMMA,CPATH,CH2,CMIN,HT1,CT1,HT2,CT2,HT3,CT3,DC,HT,HTAN,ZMIN,   FMN 0020
     2  DPRE,DPDH                                                       FMN 0021
      REAL GNDALT                                                       FMN 0022
      COMMON/GRAUND/GNDALT                                              FMN 0023
C                                                                       FMN 0024
C     HTAN IS USED CALCULATE TANGENT HEIGHT FOR PRINTING ERROR MESSAGES FMN 0025
C     HTAN IS NOT MEANT TO REPLACE HT WHICH HOLDS THE TANGENT HEIGHT IN FMN 0026
C     THIS ROUTINE AS FAR AS MODTRAN GOES.  FOR SOME INPUT ERRORS       FMN 0027
C     HT IS CALCULATED.                                                 FMN 0028
C                                                                       FMN 0029
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          FMN 0030
      REAL RE,ZMAX                                                      FMN 0031
      INTEGER IMAX,IMOD,IPATH                                           FMN 0032
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             FMN 0033
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     FMN 0034
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                FMN 0035
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     FMN 0036
      REAL TBOUND,SALB                                                  FMN 0037
      LOGICAL MODTRN                                                    FMN 0038
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   FMN 0039
     1  MODTRN                                                          FMN 0040
      REAL ETA                                                          FMN 0041
      DOUBLE PRECISION DPDEG                                            FMN 0042
      DATA DPDH/1./,ETA/5.0E-8/,DPDEG/57.2957795131D0/                  FMN 0043
C***  ETA MAY BE TOO SMALL FOR SOME COMPUTERS. TRY 1.0E-7 FOR 32 BIT    FMN 0044
C***  WORD MACHINES                                                     FMN 0045
C     CRFRCT IS REFRACTIVE CONSTANT FOR THE PATH                        FMN 0046
      CRFRCT(H) = (DPRE+H)*DPANDX(H,SH,GAMMA)                           FMN 0047
C                                                                       FMN 0048
      DPRE=DBLE(RE)                                                     FMN 0049
      N = 0                                                             FMN 0050
      CALL DPFISH(H1,SH,GAMMA)                                          FMN 0051
      CPATH = CRFRCT(H1)*SIN(ANGLE/DPDEG)                               FMN 0052
      CALL DPFISH(H2,SH,GAMMA)                                          FMN 0053
      CH2 = CRFRCT(H2)                                                  FMN 0054
      IF(ABS(CPATH/CH2).LE.1.0) THEN                                    FMN 0055
         IF(ANGLE.LE.90.0)  THEN                                        FMN 0056
            LENN = 0                                                    FMN 0057
            HMIN = H1                                                   FMN 0058
C***        CALCULATE THE ZENITH ANGLE PHI AT H2                        FMN 0059
            PHI = ASIN(CPATH/CH2)*DPDEG                                 FMN 0060
            IF(ANGLE.LE.90.0 .OR. LENN.EQ.1)PHI=DBLE(180.)-PHI          FMN 0061
            RETURN                                                      FMN 0062
         ENDIF                                                          FMN 0063
         IF(H1.LE.H2)  LENN = 1                                         FMN 0064
         IF(LENN.NE.1)  THEN                                            FMN 0065
            LENN = 0                                                    FMN 0066
            HMIN = H2                                                   FMN 0067
C***        CALCULATE THE ZENITH ANGLE PHI AT H2                        FMN 0068
            PHI = ASIN(CPATH/CH2)*DPDEG                                 FMN 0069
            IF(ANGLE.LE.90. .OR. LENN.EQ.1)PHI=DBLE(180.)-PHI           FMN 0070
            RETURN                                                      FMN 0071
         ENDIF                                                          FMN 0072
C***     LONG PATH THROUGH A TANGENT HEIGHT.                            FMN 0073
C***     SOLVE ITERATIVELY FOR THE TANGENT HEIGHT HT.                   FMN 0074
C***     HT IS THE HEIGHT FOR WHICH  INDEX(HT)*(DPRE+HT) = CPATH.       FMN 0075
         ZMIN=DBLE(GNDALT)                                              FMN 0076
         CALL DPFISH(ZMIN,SH,GAMMA)                                     FMN 0077
         CMIN = CRFRCT(ZMIN)                                            FMN 0078
C***     FOR BETA CASES (ITER>0), ALLOW FOR HT < 0.0                    FMN 0079
         IF(ITER.NE.0 .OR. CPATH.GE.CMIN)  THEN                         FMN 0080
            HT1 = (DPRE+H1)*SIN(ANGLE/DPDEG)-DPRE                       FMN 0081
            CALL DPFISH(HT1,SH,GAMMA)                                   FMN 0082
            CT1 = CRFRCT(HT1)                                           FMN 0083
            HT2 = HT1-DPDH                                              FMN 0084
            CALL DPFISH(HT2,SH,GAMMA)                                   FMN 0085
            CT2 = CRFRCT(HT2)                                           FMN 0086
C***        ITERATE TO FIND HT                                          FMN 0087
            N = 2                                                       FMN 0088
 120        CONTINUE                                                    FMN 0089
            IF(CT2 .NE. CT1) THEN                                       FMN 0090
               N = N+1                                                  FMN 0091
               HT3 = HT2+(HT2-HT1)*(CPATH-CT2)/(CT2-CT1)                FMN 0092
               CALL DPFISH(HT3,SH,GAMMA)                                FMN 0093
               CT3 = CRFRCT(HT3)                                        FMN 0094
               DC = CPATH-CT3                                           FMN 0095
               IF(ABS((CPATH-CT3)/CPATH).LT.ETA) THEN                   FMN 0096
                  HT = HT3                                              FMN 0097
                  HMIN = HT                                             FMN 0098
C***              CALCULATE THE ZENITH ANGLE PHI AT H2                  FMN 0099
                  PHI = ASIN(CPATH/CH2)*DPDEG                           FMN 0100
                  IF(ANGLE.LE.90. .OR. LENN.EQ.1)PHI=DBLE(180.)-PHI     FMN 0101
                  RETURN                                                FMN 0102
               ENDIF                                                    FMN 0103
               IF(N.GT.15) THEN                                         FMN 0104
                  DC = CPATH-CT3                                        FMN 0105
                  WRITE(IPR,24)  N,CPATH,CT3,DC,HT3                     FMN 0106
                  GO TO 25                                              FMN 0107
               ENDIF                                                    FMN 0108
               HT1 = HT2                                                FMN 0109
               CT1 = CT2                                                FMN 0110
               HT2 = HT3                                                FMN 0111
               CT2 = CT3                                                FMN 0112
               GO TO 120                                                FMN 0113
            ENDIF                                                       FMN 0114
            HT3 = HT2                                                   FMN 0115
            HT = HT3                                                    FMN 0116
            HMIN = HT                                                   FMN 0117
C***        CALCULATE THE ZENITH ANGLE PHI AT H2                        FMN 0118
            PHI = ASIN(CPATH/CH2)*DPDEG                                 FMN 0119
            IF(ANGLE.LE.90. .OR. LENN.EQ.1)PHI=DBLE(180.)-PHI           FMN 0120
            RETURN                                                      FMN 0121
         ENDIF                                                          FMN 0122
C                                                                       FMN 0123
C        RETURN WITH AN ERROR FOR PATHS TO THE SUN (OR MOON).           FMN 0124
         IF(ISSGEO.NE.0 .OR. IEMSCT.EQ.3)THEN                           FMN 0125
            IERROR=-5                                                   FMN 0126
            RETURN                                                      FMN 0127
         ENDIF                                                          FMN 0128
C***     TANGENT PATH INTERSECTS EARTH                                  FMN 0129
         WRITE(IPR,'(///2(A,F11.3),A,//9X,2(A,F11.3),A,I2,A)')          FMN 0130
     1     ' TANGENT PATH WITH H1 =',H1,' KM AND ANGLE =',ANGLE,        FMN 0131
     2     ' DEG INTERSECTS THE EARTH.',' H2 HAS BEEN RESET FROM',H2,   FMN 0132
     3     ' KM TO',ZMIN,' KM, AND LEN HAS BEEN RESET FROM',LENN,' TO 0'FMN 0133
         H2 = ZMIN                                                      FMN 0134
         HMIN = ZMIN                                                    FMN 0135
         LENN = 0                                                       FMN 0136
         CH2 = CMIN                                                     FMN 0137
C***     CALCULATE THE ZENITH ANGLE PHI AT H2                           FMN 0138
         PHI = ASIN(CPATH/CH2)*DPDEG                                    FMN 0139
         IF(ANGLE.LE.90. .OR. LENN.EQ.1)PHI=DBLE(180.)-PHI              FMN 0140
         RETURN                                                         FMN 0141
      ENDIF                                                             FMN 0142
C                                                                       FMN 0143
C***  H2 LT TANGENT HEIGHT FOR THIS H1 AND ANGLE                        FMN 0144
      IERROR = 2                                                        FMN 0145
 25   HTAN=TANHT(CPATH,H1)                                              FMN 0146
      WRITE(IPR,20)HTAN                                                 FMN 0147
C                                                                       FMN 0148
 20   FORMAT('0H2 IS LESS THAN THE TANGENT HEIGHT FOR THIS PATH AND ',  FMN 0149
     $'CAN''T BE REACHED;',/,' TANGENT HEIGHT = ',1X,F7.3,/)            FMN 0150
 24   FORMAT(///,'0FROM SUBROUTINE FNDHMN :',//,                        FMN 0151
     $     10X,'THE PROCEDURE TO FIND THE TANGENT HEIGHT DID NOT ',     FMN 0152
     $     'CONVERGE AFTER ',I3,'  ITERATIONS',//,                      FMN 0153
     $     10X,'CPATH   = ',F12.5,' KM',//,10X,'CT3     = ',F12.5,' KM',FMN 0154
     $     //,10X,'DC      = ',E12.3,' KM',//,                          FMN 0155
     $     10X,'HT3     = ',F12.5,' KM')                                FMN 0156
      END                                                               FMN 0157
