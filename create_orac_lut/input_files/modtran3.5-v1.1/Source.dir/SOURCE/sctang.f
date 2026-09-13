      FUNCTION   SCTANG(ANGLST,THTST,PSIST,IARB)                        SCT 0001
C                                                                       SCT 0002
C     FUNCTION SCTANG RETURNS THE SCATTERING ANGLE (THAT IS, THE        SCT 0003
C     ANGLE BETWEEN THE SUN'S RAYS AND THE LINE OF SIGHT) AT ANY        SCT 0004
C     POINT ALONG THE OPTICAL PATH.                                     SCT 0005
C                                                                       SCT 0006
C       PI       THE CONSTANT PI                                        SCT 0007
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       SCT 0008
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       SCT 0009
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         SCT 0010
      REAL PI,DEG,BIGNUM,BIGEXP                                         SCT 0011
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                SCT 0012
C                                                                       SCT 0013
      SUNZEN=ANGLST/DEG                                                 SCT 0014
      PTHZEN=THTST/DEG                                                  SCT 0015
      IF(IARB.EQ.0) GO TO 10                                            SCT 0016
C     SPECIAL CASES IF PSI IS ARBITRARY                                 SCT 0017
      COSANG= COS(SUNZEN)*COS(PTHZEN)                                   SCT 0018
      IF(COSANG.GT. 1.)COSANG= 1.                                       SCT 0019
      IF(COSANG.LT.-1.)COSANG=-1.                                       SCT 0020
      SCTANG=DEG*ACOS(COSANG)                                           SCT 0021
      RETURN                                                            SCT 0022
10    CONTINUE                                                          SCT 0023
      PSI=PSIST/DEG                                                     SCT 0024
C     GENERAL CASE                                                      SCT 0025
      X=SIN(SUNZEN)*SIN(PTHZEN)*COS(PSI)+COS(SUNZEN)*COS(PTHZEN)        SCT 0026
      IF(X.GT. 1.)X= 1.-1.E-6                                           SCT 0027
      IF(X.LT.-1.)X=-1.+1.E-6                                           SCT 0028
      SCTANG=DEG*ACOS(X)                                                SCT 0029
      RETURN                                                            SCT 0030
      END                                                               SCT 0031
