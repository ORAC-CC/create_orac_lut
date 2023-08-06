      REAL FUNCTION BBFN(T,V)                                           BBF 0001
C                                                                       BBF 0002
C     PLANCK BLACK BODY FUNCTION [WATTS CM-2 STER-1 / CM-1].            BBF 0003
C                                                                       BBF 0004
C     DECLARE ARGUMENTS:                                                BBF 0005
C       T        TEMPERATURE [K].                                       BBF 0006
C       V        FREQUENCY [CM-1].                                      BBF 0007
      REAL T,V                                                          BBF 0008
C                                                                       BBF 0009
C     LIST COMMONS:                                                     BBF 0010
C       PI       THE CONSTANT PI                                        BBF 0011
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       BBF 0012
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       BBF 0013
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         BBF 0014
      REAL PI,DEG,BIGNUM,BIGEXP                                         BBF 0015
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                BBF 0016
C                                                                       BBF 0017
C     DECLARE LOCAL VARIABLES:                                          BBF 0018
C       X        EXPONENT USED IN PLANCK FUNCTION.                      BBF 0019
      REAL X                                                            BBF 0020
      BBFN=0.                                                           BBF 0021
      IF(V.LE.0.)RETURN                                                 BBF 0022
      X=1.43879*V/T                                                     BBF 0023
C                                                                       BBF 0024
C     PROTECT AGAINST EXPONENTIAL OVERFLOW                              BBF 0025
      IF(X.GT.BIGEXP)RETURN                                             BBF 0026
      BBFN=1.190956E-12*V**3/(EXP(X)-1.)                                BBF 0027
      RETURN                                                            BBF 0028
      END                                                               BBF 0029
