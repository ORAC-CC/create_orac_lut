      SUBROUTINE DPEXNT(X,X1,X2,A)                                      XNT 0001
C                                                                       XNT 0002
C     DOUBLE PRECISION VERSION OF THE ROUTINE EXPINT                    XNT 0003
C                                                                       XNT 0004
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)                               XNT 0005
C     EXPONENTIAL INTERPOLATION                                         XNT 0006
CJ DB                                                                   XNT 0007
         IF(X1.LT.0.)X1=0.                                              XNT 0008
         IF(X2.LT.0.)X2=0.                                              XNT 0009
CJ                                                                      XNT 0010
      IF(X1.EQ.0.0 .OR. X2.EQ.0.0)  GO TO 100                           XNT 0011
      X = X1*(X2/X1)**A                                                 XNT 0012
      RETURN                                                            XNT 0013
  100 X = X1+(X2-X1)*A                                                  XNT 0014
      RETURN                                                            XNT 0015
      END                                                               XNT 0016
