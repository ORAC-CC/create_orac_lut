      REAL*8 FUNCTION  RATIO( A, B )                                    RAT 0001
C             INSERT FOR DOUBLE PRECISION MODIFICATION - NORTH          RAT 0002
                IMPLICIT DOUBLE PRECISION ( A-H, O-Z)                   RAT 0003
C        CALCULATE RATIO  A/B  WITH OVER- AND UNDER-FLOW PROTECTION     RAT 0004
                                                                        RAT 0005
         IF ( DABS(A).LT.1.0E-8 .AND. DABS(B).LT.1.0E-8 )  THEN         RAT 0006
            RATIO = 1.0                                                 RAT 0007
         ELSE IF ( B.EQ.0.0 )  THEN                                     RAT 0008
            RATIO = 1.E+20                                              RAT 0009
         ELSE                                                           RAT 0010
            RATIO = A / B                                               RAT 0011
         END IF                                                         RAT 0012
                                                                        RAT 0013
      RETURN                                                            RAT 0014
      END                                                               RAT 0015
