      DOUBLE PRECISION FUNCTION DPANDX(H,SH,GAMMA)                      NDX 0001
C                                                                       NDX 0002
C     DOUBLE PRECISION (DP) VERSION OF THE ROUTINE PREVIOUSLY           NDX 0003
C     CALLED ANDEX.  COMPUTES THE INDEX OF REFRACTION AT                NDX 0004
C     HEIGHT H, SH IS THE SCALE HEIGHT, GAMMA IS THE VALUE              NDX 0005
C     AT H=0 OF THE REFRACTIVITY = INDEX OF REFRACTION - 1.             NDX 0006
      DOUBLE PRECISION H,SH,GAMMA,HSH                                   NDX 0007
C                                                                       NDX 0008
C     LIST COMMONS:                                                     NDX 0009
C       PI       THE CONSTANT PI                                        NDX 0010
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       NDX 0011
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       NDX 0012
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         NDX 0013
      REAL PI,DEG,BIGNUM,BIGEXP                                         NDX 0014
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                NDX 0015
      DPANDX=DBLE(1.)+GAMMA                                             NDX 0016
      IF(SH.EQ.0. .OR. H.LE.0.)RETURN                                   NDX 0017
      HSH=H/SH                                                          NDX 0018
      IF(HSH.GT.BIGEXP)THEN                                             NDX 0019
          DPANDX=DBLE(1.)                                               NDX 0020
      ELSE                                                              NDX 0021
          DPANDX=DBLE(1.)+GAMMA*EXP(-HSH)                               NDX 0022
      ENDIF                                                             NDX 0023
      RETURN                                                            NDX 0024
      END                                                               NDX 0025
