      FUNCTION   DEL(PSIO,DELO,BETA,IARBO)                              DEL 0001
C                                                                       DEL 0002
C     FUNCTION DEL RETURNS THE VALUE OF THE SUN'S ZENITH ANGLE          DEL 0003
C     AT ANY POINT ALONG THE OPTICAL PATH BASED UPON STRAIGHT           DEL 0004
C     LINE GEOMETRY (NO REFRACTION). THIS ANGLE IS USED TO SPECIFY      DEL 0005
C     THE SCATTERING POINT TO SUN PATHS. THE BENDING DUE TO REFRACTION  DEL 0006
C     ALONG THIS PATH IS DETERMINED BY THE GEO ROUTINES. IF THE BENDING DEL 0007
C     IS GREATER THAN ONE DEGREE THE ZENITH ANGLE IS CORRECTED ACCORDINGDEL 0008
C     AND THE PATH CALCULATION IS REPEATED.                             DEL 0009
C                                                                       DEL 0010
C     LIST COMMONS:                                                     DEL 0011
C       PI       THE CONSTANT PI                                        DEL 0012
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       DEL 0013
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       DEL 0014
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         DEL 0015
      REAL PI,DEG,BIGNUM,BIGEXP                                         DEL 0016
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                DEL 0017
      IF(IARBO.EQ.0) GO TO 10                                           DEL 0018
C     SPECIAL CASES IF PSIO IS ARBITRARY                                DEL 0019
      IF(IARBO.EQ.1) DEL=DELO                                           DEL 0020
      IF(IARBO.EQ.2) DEL=BETA                                           DEL 0021
      IF(IARBO.EQ.3) DEL=0.0                                            DEL 0022
      RETURN                                                            DEL 0023
10    CONTINUE                                                          DEL 0024
      PSIOR=PSIO/DEG                                                    DEL 0025
      DELOR=DELO/DEG                                                    DEL 0026
      BETAR=BETA/DEG                                                    DEL 0027
C     GENERAL CASE                                                      DEL 0028
      X=COS(DELOR)*COS(BETAR)+SIN(DELOR)*SIN(BETAR)*COS(PSIOR)          DEL 0029
      DEL=DEG*ACOS(X)                                                   DEL 0030
      RETURN                                                            DEL 0031
      END                                                               DEL 0032
