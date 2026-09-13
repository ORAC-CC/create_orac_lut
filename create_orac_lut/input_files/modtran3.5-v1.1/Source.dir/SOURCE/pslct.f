      INTEGER FUNCTION PSLCT(ANGLE,RANGE,BETA)                          SLT 0001
C                                                                       SLT 0002
C     THIS ROUTINE RETURNS THE AN INTEGER VALUE INDICATING THE TYPE     SLT 0003
C     OF SLANT (ITYPE=2) PATH.  THE FOLLOWING VALUES ARE RETURNED:      SLT 0004
C          PSLCT = 21 FOR CASE 2A (H1,H2,ANGLE)                         SLT 0005
C          PSLCT = 22 FOR CASE 2B (H1,ANGLE,RANGE)                      SLT 0006
C          PSLCT = 23 FOR CASE 2C (H1,H2,RANGE)                         SLT 0007
C          PSLCT = 24 FOR CASE 2D (H1,H2,BETA)                          SLT 0008
C                                                                       SLT 0009
C     H1     H2     ANGLE     RANGE      BETA        CASE  PSLCT        SLT 0010
C--------------------------------------------     ----------------      SLT 0011
C     X      X      X                                 2A     21         SLT 0012
C                                                                       SLT 0013
C     X             X         X                       2B     22         SLT 0014
C                                                                       SLT 0015
C     X      X                X                       2C     23         SLT 0016
C                                                                       SLT 0017
C     X      X                           X            2D     24         SLT 0018
C                                                                       SLT 0019
C     DECLARE INPUTS                                                    SLT 0020
      REAL ANGLE,RANGE,BETA                                             SLT 0021
      IF(BETA.GT.0.)THEN                                                SLT 0022
C                                                                       SLT 0023
C         BETA > 0 IMPLIES CASE 2D (H1,H2,BETA)                         SLT 0024
          PSLCT=24                                                      SLT 0025
      ELSEIF(RANGE.GT.0.)THEN                                           SLT 0026
C                                                                       SLT 0027
C         RANGE > 0 AND BETA = 0 IMPLIES CASE 2B OR 2C.                 SLT 0028
          IF(ANGLE.EQ.0.)THEN                                           SLT 0029
C                                                                       SLT 0030
C             RANGE > 0, BETA = 0 AND ANGLE = 0 IMPLIES CASE 2C         SLT 0031
              PSLCT=23                                                  SLT 0032
          ELSE                                                          SLT 0033
C                                                                       SLT 0034
C             RANGE > 0, BETA = 0 AND ANGLE NON-ZERO IMPLIES CASE 2B    SLT 0035
              PSLCT=22                                                  SLT 0036
          ENDIF                                                         SLT 0037
      ELSE                                                              SLT 0038
C                                                                       SLT 0039
C         RANGE =0 AND BETA = 0 IMPLIES CASE 2A                         SLT 0040
          PSLCT=21                                                      SLT 0041
      ENDIF                                                             SLT 0042
      RETURN                                                            SLT 0043
      END                                                               SLT 0044
