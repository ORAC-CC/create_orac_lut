      REAL FUNCTION BMERFU(Y)                                           ERF 0001
C                                                                       ERF 0002
C     APPROXIMATION FOR EXP(-Y*Y) + SQRT(PI) * Y * ERF(Y) - 1           ERF 0003
      DATA P,A1,A2,A3,A4,A5,RTPI/.3275911,.451673692,-.504257335,       ERF 0004
     1  2.519390259,-2.575644906,1.881292140,1.772453851/               ERF 0005
C                                                                       ERF 0006
C     TO SINGLE PRECISION ACCURACY, THE FUNCTION SIMPLY                 ERF 0007
C     EQUALS THE EXPRESSION BELOW FOR Y > 3.2                           ERF 0008
      BMERFU=RTPI*Y-1.                                                  ERF 0009
      IF(Y.GE.3.2)RETURN                                                ERF 0010
C                                                                       ERF 0011
C     FOR Y > 2.33, THE ERROR FUNCTION IS CALCULATED FROM               ERF 0012
C     ITS CONTINUED FRACTION (ABRAMOWITZ & STEGUN, 7.1.14).             ERF 0013
C     FOR Y < 2.33, THE RATIONAL APPROXIMATION IS USED                  ERF 0014
C     (ABRAMOWITZ & STEGUN, 7.1.7.1.25).  THE TWO EXPRESSIONS           ERF 0015
C     ARE EQUAL (TO SINGLE PRECISION ACCURACY) AT Y = 2.33              ERF 0016
      IF(Y.GE.2.33)THEN                                                 ERF 0017
          T=Y*Y                                                         ERF 0018
          BMERFU=BMERFU+(1.75+T*.5)/(3.75+T*(5.+T))*EXP(-T)             ERF 0019
          RETURN                                                        ERF 0020
      ENDIF                                                             ERF 0021
      T=1./(1.+P*Y)                                                     ERF 0022
      BMERFU=BMERFU+(1.-Y*T*(A1+T*(A2+T*(A3+T*(A4+T*A5)))))*EXP(-Y*Y)   ERF 0023
      RETURN                                                            ERF 0024
      END                                                               ERF 0025
