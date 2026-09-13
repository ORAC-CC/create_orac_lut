      FUNCTION TNRAIN(RR,V,TM,RADFLD)                                   TNR 0001
C                                                                       TNR 0002
C       PI       THE CONSTANT PI                                        TNR 0003
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       TNR 0004
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       TNR 0005
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         TNR 0006
      REAL PI,DEG,BIGNUM,BIGEXP                                         TNR 0007
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                TNR 0008
CCC   CALCULATES TRANSMISSION DUE TO RAIN AS A FUNCTION OF              TNR 0009
CCC   RR=RAIN RATE IN MM/HR                                             TNR 0010
CCC   OR WITHIN 350CM-1 USES THE MICROWAVE TABLE ROUTINE TO             TNR 0011
CCC   OBTAIN THE EXTINCTION DUE TO RAIN                                 TNR 0012
CCC   RANGE=SLANT RANGE KM                                              TNR 0013
CCC                                                                     TNR 0014
CCC   ASSUMES A MARSHALL-PALMER RAIN DROP SIZE DISTRIBUTION             TNR 0015
CCC   N(D)=NZERO*EXP(-A*D)                                              TNR 0016
CCC   NZERO=8.E3 (MM-1)  (M-3)                                          TNR 0017
CCC   A=41.*RR**(-0.21)                                                 TNR 0018
CCC   D=DROP DIAMETER (CM)                                              TNR 0019
CCC                                                                     TNR 0020
      REAL NZERO                                                        TNR 0021
      DATA NZERO /8000./                                                TNR 0022
CCC                                                                     TNR 0023
      A=41./RR**0.21                                                    TNR 0024
CCC                                                                     TNR 0025
      IF(RR.LE.0)TNRAIN=1.                                              TNR 0026
      IF(RR.LE.0)RETURN                                                 TNR 0027
CCC                                                                     TNR 0028
      IF(V.GE.350.0) THEN                                               TNR 0029
       TNRAIN=PI*NZERO/A**3                                             TNR 0030
      ELSE                                                              TNR 0031
       TNRAIN=GMRAIN(V,TM,RR)                                           TNR 0032
       TNRAIN = TNRAIN * RADFLD                                         TNR 0033
      END IF                                                            TNR 0034
      RETURN                                                            TNR 0035
      END                                                               TNR 0036
