      FUNCTION PSI(PSIO,DELO,BETA,IARB,IARBO)                           PSI 0001
C                                                                       PSI 0002
C     FUNCTION PSI RETURNS THE VALUE OF SOLAR AZIMUTH RELATIVE TO       PSI 0003
C     THE LINE OF SIGHT, AT THE CURRENT SCATTERING LOCATION             PSI 0004
C                                                                       PSI 0005
C       PI       THE CONSTANT PI                                        PSI 0006
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       PSI 0007
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       PSI 0008
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         PSI 0009
      REAL PI,DEG,BIGNUM,BIGEXP                                         PSI 0010
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                PSI 0011
      DATA  EPSILN/1.0E-5/                                              PSI 0012
      PSI=PSIO                                                          PSI 0013
      DELOR=DELO/DEG                                                    PSI 0014
      BETAR=BETA/DEG                                                    PSI 0015
      IF(IARBO.EQ.0) GO TO 5                                            PSI 0016
C     SPECIAL CASES WHEN PSIO IS ARBITRARY                              PSI 0017
      IARB=IARBO                                                        PSI 0018
      PSI = 0.                                                          PSI 0019
      IF(IARBO.EQ.1.OR.IARBO.EQ.3) RETURN                               PSI 0020
      IF(BETA.LE.EPSILN) RETURN                                         PSI 0021
C     PSI=180.0 (MOVED OUT FROM UNDER THE SUN)                          PSI 0022
      IARB=0                                                            PSI 0023
      PSI=180.0                                                         PSI 0024
      RETURN                                                            PSI 0025
5     CONTINUE                                                          PSI 0026
C     GENERAL CASE                                                      PSI 0027
      PSIOR=PSIO/DEG                                                    PSI 0028
      IARB=0                                                            PSI 0029
      ANUMER=SIN(DELOR)*SIN(PSIOR)                                      PSI 0030
      DENOM=COS(BETAR)*SIN(DELOR)*COS(PSIOR)-SIN(BETAR)*COS(DELOR)      PSI 0031
C     SPECIAL CASES                                                     PSI 0032
C     NUMERATOR GOES TO ZERO IN THE FOLLOWING 3 CASES                   PSI 0033
C 1)  DELO=0.0                                                          PSI 0034
      IF(DELO.GT.EPSILN) GO TO 20                                       PSI 0035
      IF(BETA.GT.EPSILN) GO TO 10                                       PSI 0036
      IARB=2                                                            PSI 0037
      RETURN                                                            PSI 0038
10    PSI=180.0                                                         PSI 0039
      RETURN                                                            PSI 0040
C 2)  PSIO=0.0                                                          PSI 0041
20    IF(ABS(PSIO).GT.EPSILN) GO TO 40                                  PSI 0042
      IF(ABS(BETA-DELO).GE.EPSILN) GO TO 30                             PSI 0043
C     SCATTERING POINT IS DIRECTLY UNDER THE SUN                        PSI 0044
      IARB=2                                                            PSI 0045
      RETURN                                                            PSI 0046
30    IF(BETA.LT.DELO) PSI=0.0                                          PSI 0047
      IF(BETA.GT.DELO) PSI=180.0                                        PSI 0048
      RETURN                                                            PSI 0049
C 3)  PSIO=180.0                                                        PSI 0050
40    IF(ABS(PSIO).LT.(180.0-EPSILN)) GO TO 60                          PSI 0051
      PSI=180.0                                                         PSI 0052
      RETURN                                                            PSI 0053
60    CONTINUE                                                          PSI 0054
C     DENOMINATOR CAN GO TO ZERO FOR THE FOLLOWING 2 CASES              PSI 0055
C 1)  BETA=DELO AND PSIO=0.0                                            PSI 0056
C     THIS CASE WAS HANDLED EARLIER                                     PSI 0057
C 2)  GENERAL CASE                                                      PSI 0058
      IF(ABS(DENOM).GT.EPSILN) GO TO 80                                 PSI 0059
      IF(PSIO.LT.0.0) PSI=-90.0                                         PSI 0060
      IF(PSIO.GT.0.0) PSI=90.0                                          PSI 0061
      RETURN                                                            PSI 0062
80    CONTINUE                                                          PSI 0063
      PSI=DEG*ATAN(ANUMER/DENOM)                                        PSI 0064
C     NOTE ATAN RETURNS ARGUMENTS BETWEEN -90 AND 90, PSI               PSI 0065
C     AND PSIO SHOULD BE OF THE SAME SIGN.                              PSI 0066
      IF(PSIO.GT.0.0.AND.PSI.LT.0.0) PSI=PSI+180.                       PSI 0067
      IF(PSIO.LT.0.0.AND.PSI.GT.0.0) PSI=PSI-180.                       PSI 0068
      RETURN                                                            PSI 0069
      END                                                               PSI 0070
