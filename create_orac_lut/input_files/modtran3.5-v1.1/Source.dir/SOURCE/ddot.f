      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)                 DOT 0001
C                                                                       DOT 0002
C     FORMS THE DOT PRODUCT OF TWO VECTORS.                             DOT 0003
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  DOT 0004
C     JACK DONGARRA, LINPACK, 3/11/78.                                  DOT 0005
C                                                                       DOT 0006
      DOUBLE PRECISION DX(*),DY(*),DTEMP                                DOT 0007
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 DOT 0008
C                                                                       DOT 0009
      DDOT = 0.0D0                                                      DOT 0010
      DTEMP = 0.0D0                                                     DOT 0011
      IF(N.LE.0)RETURN                                                  DOT 0012
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               DOT 0013
C                                                                       DOT 0014
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                DOT 0015
C          NOT EQUAL TO 1                                               DOT 0016
C                                                                       DOT 0017
      IX = 1                                                            DOT 0018
      IY = 1                                                            DOT 0019
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 DOT 0020
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 DOT 0021
      DO 10 I = 1,N                                                     DOT 0022
        DTEMP = DTEMP + DX(IX)*DY(IY)                                   DOT 0023
        IX = IX + INCX                                                  DOT 0024
        IY = IY + INCY                                                  DOT 0025
   10 CONTINUE                                                          DOT 0026
      DDOT = DTEMP                                                      DOT 0027
      RETURN                                                            DOT 0028
C                                                                       DOT 0029
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            DOT 0030
C                                                                       DOT 0031
C                                                                       DOT 0032
C        CLEAN-UP LOOP                                                  DOT 0033
C                                                                       DOT 0034
   20 M = MOD(N,5)                                                      DOT 0035
      IF( M .EQ. 0 ) GO TO 40                                           DOT 0036
      DO 30 I = 1,M                                                     DOT 0037
        DTEMP = DTEMP + DX(I)*DY(I)                                     DOT 0038
   30 CONTINUE                                                          DOT 0039
      IF( N .LT. 5 ) GO TO 60                                           DOT 0040
   40 MP1 = M + 1                                                       DOT 0041
      DO 50 I = MP1,N,5                                                 DOT 0042
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +             DOT 0043
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)DOT 0044
   50 CONTINUE                                                          DOT 0045
   60 DDOT = DTEMP                                                      DOT 0046
      RETURN                                                            DOT 0047
      END                                                               DOT 0048
