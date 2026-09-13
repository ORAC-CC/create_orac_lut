      SUBROUTINE INTERP(INTYPE,X,X1,X2,F,F1,F2)                         INT 0001
C     SUBROUTINE INTERP INTERPOLATES TO DETERMINE THE VALUE OF F        INT 0002
C     AT X, GIVEN F1 AT X1 AND F2 AT X2.                                INT 0003
C     INTYPE=1 FOR LINEAR INTERPOLATION                                 INT 0004
C     INTYPE=2 FOR LOGARITHMIC INTERPOLATION                            INT 0005
      ITYPE=INTYPE                                                      INT 0006
      IF(F1.LE.0.0.OR.F2.LE.0.0) ITYPE=1                                INT 0007
      IF(ITYPE.EQ.2) GO TO 100                                          INT 0008
C     LINEAR INTERPOLATION                                              INT 0009
      F=F1+(X-X1)*(F2-F1)/(X2-X1)                                       INT 0010
      RETURN                                                            INT 0011
100   CONTINUE                                                          INT 0012
      A1=ALOG(F1)                                                       INT 0013
      A2=ALOG(F2)                                                       INT 0014
      A=A1+(X-X1)*(A2-A1)/(X2-X1)                                       INT 0015
      F=EXP(A)                                                          INT 0016
      RETURN                                                            INT 0017
      END                                                               INT 0018
