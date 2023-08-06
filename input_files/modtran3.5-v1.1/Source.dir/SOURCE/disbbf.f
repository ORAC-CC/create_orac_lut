      DOUBLE PRECISION FUNCTION DISBBF(T,V)                             DBB 0001
C                                                                       DBB 0002
C     PLANCK BLACK BODY FUNCTION [WATTS CM-2 STER-1 / CM-1].            DBB 0003
C                                                                       DBB 0004
C     DECLARE ARGUMENTS:                                                DBB 0005
C       T        TEMPERATURE [K].                                       DBB 0006
C       V        FREQUENCY [CM-1].                                      DBB 0007
      DOUBLE PRECISION T,V                                              DBB 0008
C                                                                       DBB 0009
C     LIST COMMONS:                                                     DBB 0010
C       PI       THE CONSTANT PI                                        DBB 0011
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       DBB 0012
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       DBB 0013
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         DBB 0014
      REAL PI,DEG,BIGNUM,BIGEXP                                         DBB 0015
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                DBB 0016
C                                                                       DBB 0017
C     DECLARE LOCAL VARIABLES:                                          DBB 0018
C       X        EXPONENT USED IN PLANCK FUNCTION.                      DBB 0019
      REAL X                                                            DBB 0020
      DISBBF=DBLE(0.)                                                   DBB 0021
      IF(V.LE.0.)RETURN                                                 DBB 0022
      X=1.43879*REAL(V/T)                                               DBB 0023
C                                                                       DBB 0024
C     PROTECT AGAINST EXPONENTIAL OVERFLOW                              DBB 0025
      IF(X.GT.BIGEXP)RETURN                                             DBB 0026
      DISBBF=DBLE(1.190956E-12/(EXP(X)-1.))*V**3                        DBB 0027
      RETURN                                                            DBB 0028
      END                                                               DBB 0029
