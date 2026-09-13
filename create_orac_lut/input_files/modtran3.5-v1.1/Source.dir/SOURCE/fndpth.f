      SUBROUTINE FNDPTH(CPATH,H1,HTAN,H2,RANGEI,BETA,LENN,ANGLE,PHI)    FPT 0001
C                                                                       FPT 0002
C     THIS ROUTINE DETERMINES H2, BETA AND LENN GIVEN H1, ANGLE, RANGE, FPT 0003
C     HTAN, AND CPATH (THE LOWTRAN 6 MANUAL, AFGL-TR-83-0187, PP 14-17).FPT 0004
C                                                                       FPT 0005
C     DECLARE ARGUMENTS:                                                FPT 0006
C       CPATH   REFRACTIVE PATH CONSTANT [KM].                          FPT 0007
C       H1      OBSERVER ALTITUDE [KM].                                 FPT 0008
C       HTAN    SLANT PATH TANGENT ALTITUDE [KM].                       FPT 0009
C       H2      FINAL ALTITUDE [KM].                                    FPT 0010
C       RANGEI  SLANT PATH RANGE [KM].                                  FPT 0011
C       BETA    EACH CENTER ANGLE [DEG].                                FPT 0012
C       LENN    LENGTH SWITCH (=0 FOR SHORT PATHS,                      FPT 0013
C                 =1 FOR PATHS THROUGH TANGENT POINT WITH H2<H1).       FPT 0014
C       ANGLE   PATH ZENITH ANGLE AT H1 TOWARDS H2 [DEG].               FPT 0015
C       PHI     PATH ZENITH ANGLE AT H2 TOWARDS H1 [DEG].               FPT 0016
      INTEGER LENN                                                      FPT 0017
      DOUBLE PRECISION CPATH,H1,HTAN,H2,RANGEI,BETA,ANGLE,PHI           FPT 0018
C                                                                       FPT 0019
C     LIST COMMONS:                                                     FPT 0020
      REAL GNDALT                                                       FPT 0021
      COMMON/GRAUND/GNDALT                                              FPT 0022
      REAL RE,ZMAX                                                      FPT 0023
      INTEGER IMAX,IMOD,IPATH                                           FPT 0024
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             FPT 0025
C                                                                       FPT 0026
C     DECLARE LOCAL VARIABLES                                           FPT 0027
C       CAPRJ IS FOR CAPITAL R WITH SUBSCRIPT J                         FPT 0028
C       PNTGRN IS THE INTEGRAND OF EQUATION 21.                         FPT 0029
      INTEGER I                                                         FPT 0030
      DOUBLE PRECISION R1,R2,R,DR,RPLDR,RANGEO,PNTGRN,SAVE,             FPT 0031
     1  RX,RATIO,CAPRJ,XJ,XJPL1,DX,CTHETA,STHETA,DBETA,                 FPT 0032
     2  Z,DZ,DPDEG,BASE,PERP,DRNG,DIFF,DPRE                             FPT 0033
      DATA DR,DPDEG/0.005,57.2957795131D0/                              FPT 0034
C                                                                       FPT 0035
C     (RANGEI .LT. DR) SHOULD NOT HAPPEN; SO THIS CHECK IS REDUNDANT.   FPT 0036
      IF (RANGEI .LT. DR) STOP'STOPPED IN FNDPTH'                       FPT 0037
      DPRE=DBLE(RE)                                                     FPT 0038
      RANGEO=0                                                          FPT 0039
      BETA=0                                                            FPT 0040
      DO 200 I=1,2                                                      FPT 0041
         IF (ANGLE.LE.90.D0 .AND. I.EQ.1)GOTO200                        FPT 0042
         IF (I .EQ. 1) THEN                                             FPT 0043
            R1=DPRE+H1                                                  FPT 0044
            R2=DPRE+HTAN                                                FPT 0045
         ELSEIF (I .EQ. 2) THEN                                         FPT 0046
            IF(HTAN.LT.DBLE(GNDALT+.001) .AND. ANGLE.GT.90.D0)GOTO200   FPT 0047
C           IF (HTAN APPROXIMATELY 0) THEN YOU ARE ABOUT TO HIT THE EARTFPT 0048
            R2=DPRE+DBLE(ZMAX)                                          FPT 0049
            IF (ANGLE.LE.90.D0) THEN                                    FPT 0050
               R1=DPRE + H1                                             FPT 0051
            ELSE                                                        FPT 0052
               R1 =DPRE + HTAN                                          FPT 0053
            ENDIF                                                       FPT 0054
         ENDIF                                                          FPT 0055
         IF (R2 .LT. R1) THEN                                           FPT 0056
            DZ=-DR                                                      FPT 0057
         ELSE                                                           FPT 0058
            DZ=DR                                                       FPT 0059
         ENDIF                                                          FPT 0060
         DO 100 R=R1, R2-DZ, DZ                                         FPT 0061
            Z=R-DPRE                                                    FPT 0062
            CALL IRFXN(Z,RX,RATIO)                                      FPT 0063
            STHETA=CPATH/(RX*R)                                         FPT 0064
            IF(STHETA.GT. DBLE(1.))STHETA= DBLE(1.)                     FPT 0065
            IF(STHETA.LT.-DBLE(1.))STHETA=-DBLE(1.)                     FPT 0066
            SAVE=STHETA                                                 FPT 0067
            CTHETA=SQRT(DBLE(1.)-STHETA**2)                             FPT 0068
C                                                                       FPT 0069
C           IF (R1 .GT. R2) THEN CTHETA IS NEGATIVE BECAUSE THETA .GT. 9FPT 0070
            IF (R1 .GT. R2) CTHETA=-CTHETA                              FPT 0071
            XJ=R*CTHETA                                                 FPT 0072
            CAPRJ=-R/RATIO                                              FPT 0073
            PNTGRN=1/(1-CAPRJ*STHETA*STHETA)                            FPT 0074
            RPLDR=R+DZ                                                  FPT 0075
            Z=RPLDR-DPRE                                                FPT 0076
            CALL IRFXN(Z,RX,RATIO)                                      FPT 0077
            STHETA=CPATH/(RX*RPLDR)                                     FPT 0078
            CTHETA=SQRT(1-STHETA**2)                                    FPT 0079
            IF (R1 .GT. R2) CTHETA=-CTHETA                              FPT 0080
            XJPL1=RPLDR*CTHETA                                          FPT 0081
            DX=XJPL1-XJ                                                 FPT 0082
            DRNG=PNTGRN*DX                                              FPT 0083
            RANGEO=RANGEO + DRNG                                        FPT 0084
            DBETA=(((SAVE+STHETA)/2)*(PNTGRN*DX))/(R-DZ/2)              FPT 0085
            BETA=BETA+DBETA                                             FPT 0086
            IF (RANGEO .GE. RANGEI) THEN                                FPT 0087
               DIFF=(RANGEI-(RANGEO-DRNG))                              FPT 0088
               H2=R -DPRE + (DZ/DRNG)*DIFF                              FPT 0089
               BETA=BETA*DPDEG                                          FPT 0090
               IF (I .EQ. 2) THEN                                       FPT 0091
                  LENN=1                                                FPT 0092
                  IF (ANGLE.LE.90.D0)LENN=0                             FPT 0093
                  IF(H2.LT.HTAN)THEN                                    FPT 0094
C                                                                       FPT 0095
C                    THIS WILL BE THE CASE IF I=2, AND YOU HAVE         FPT 0096
C                    GONE THROUGH THE R-LOOP BARELY (ONLY) ONCE.        FPT 0097
                     H2=HTAN                                            FPT 0098
                     LENN=0                                             FPT 0099
                  ENDIF                                                 FPT 0100
               ELSE                                                     FPT 0101
                  LENN=0                                                FPT 0102
               ENDIF                                                    FPT 0103
C                                                                       FPT 0104
C              CORRECTION FOR VERY SHORT PATHS; HERE IT IS ABOUT 5 KM ORFPT 0105
               IF (RANGEI .LT. 5.0 .AND. RANGEO/RANGEI .GT. 1.05) THEN  FPT 0106
C                 CALCULATE BETA BY STARIGHT LINE GEOMETRY.             FPT 0107
                  PERP =SIN(ANGLE/DPDEG)*RANGEI                         FPT 0108
                  BASE=COS(ANGLE/DPDEG)*RANGEI + DPRE+H1                FPT 0109
                  BETA=ATAN(PERP/BASE)*DPDEG                            FPT 0110
                  RANGEO=RANGEI                                         FPT 0111
                  H2=BASE-DPRE                                          FPT 0112
               ENDIF                                                    FPT 0113
               PHI=DBLE(180.)-ACOS(CTHETA)*DPDEG                        FPT 0114
               RETURN                                                   FPT 0115
            ENDIF                                                       FPT 0116
 100     CONTINUE                                                       FPT 0117
 200  CONTINUE                                                          FPT 0118
C                                                                       FPT 0119
C     COMES HERE IF YOU HAVE REACHED ZMAX, BUT YOUR RANGEI IS STILL     FPT 0120
C     NOT EQUAL TO OUTPUT VALUE.                                        FPT 0121
C     IN THIS CASE DO THE FOLLOWING.                                    FPT 0122
C                                                                       FPT 0123
      RANGEI=RANGEO                                                     FPT 0124
      H2=DBLE(ZMAX)                                                     FPT 0125
      IF (ANGLE.LE.90.D0)THEN                                           FPT 0126
         LENN=0                                                         FPT 0127
      ELSE                                                              FPT 0128
         LENN=1                                                         FPT 0129
      ENDIF                                                             FPT 0130
      IF(HTAN.LT.DBLE(GNDALT+.001) .AND. ANGLE.GT.90.D0)THEN            FPT 0131
C        YOU HAVE HIT THE EARTH IF YOU ARE AT THIS POINT OF THE CODE    FPT 0132
         LENN=0                                                         FPT 0133
         H2=DBLE(GNDALT)                                                FPT 0134
      ENDIF                                                             FPT 0135
      BETA=BETA*DPDEG                                                   FPT 0136
      PHI=DBLE(180.)-ACOS(CTHETA)*DPDEG                                 FPT 0137
      RETURN                                                            FPT 0138
      END                                                               FPT 0139
