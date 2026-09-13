      DOUBLE PRECISION FUNCTION DPRARF(H,SH,GAMMA)                      RAD 0001
C                                                                       RAD 0002
C     DOUBLE PRECISION (DP) VERSION OF THE ROUTINE PREVIOUSLY CALLED    RAD 0003
C     RADREF.  COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FORRAD 0004
C     HORIZONTAL PATH:  DPRARF = DPANDX/ D(DPANDX)/D(RADIUS)            RAD 0005
C     DPANDX = DOUBLE PRECISION (DP) VERSION OF ANDEX                   RAD 0006
      DOUBLE PRECISION H,SH,GAMMA,HSH                                   RAD 0007
C                                                                       RAD 0008
C       PI       THE CONSTANT PI                                        RAD 0009
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       RAD 0010
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       RAD 0011
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         RAD 0012
      REAL PI,DEG,BIGNUM,BIGEXP                                         RAD 0013
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                RAD 0014
      IF(SH.EQ.0.)GOTO20                                                RAD 0015
      HSH=H/SH                                                          RAD 0016
      IF(HSH .GT. BIGEXP) GO TO 20                                      RAD 0017
      DPRARF=SH*(DBLE(1.)+EXP(HSH )/GAMMA)                              RAD 0018
      RETURN                                                            RAD 0019
   20 DPRARF=DBLE(BIGNUM)                                               RAD 0020
      RETURN                                                            RAD 0021
      END                                                               RAD 0022
