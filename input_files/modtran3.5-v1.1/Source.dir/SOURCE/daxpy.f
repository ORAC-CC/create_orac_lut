      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)                            DAX 0001
C                                                                       DAX 0002
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            DAX 0003
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  DAX 0004
C     JACK DONGARRA, LINPACK, 3/11/78.                                  DAX 0005
C                                                                       DAX 0006
      DOUBLE PRECISION DX(*),DY(*),DA                                   DAX 0007
      INTEGER I,INCX,INCY,M,MP1,N                                       DAX 0008
C                                                                       DAX 0009
      IF(N.LE.0)RETURN                                                  DAX 0010
      IF (DA .EQ. 0.0D0) RETURN                                         DAX 0011
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               DAX 0012
C                                                                       DAX 0013
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                DAX 0014
C          NOT EQUAL TO 1                                               DAX 0015
C                                                                       DAX 0016
      IX = 1                                                            DAX 0017
      IY = 1                                                            DAX 0018
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 DAX 0019
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 DAX 0020
      DO 10 I = 1,N                                                     DAX 0021
        DY(IY) = DY(IY) + DA*DX(IX)                                     DAX 0022
        IX = IX + INCX                                                  DAX 0023
        IY = IY + INCY                                                  DAX 0024
   10 CONTINUE                                                          DAX 0025
      RETURN                                                            DAX 0026
C                                                                       DAX 0027
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            DAX 0028
C                                                                       DAX 0029
C                                                                       DAX 0030
C        CLEAN-UP LOOP                                                  DAX 0031
C                                                                       DAX 0032
   20 M = MOD(N,4)                                                      DAX 0033
      IF( M .EQ. 0 ) GO TO 40                                           DAX 0034
      DO 30 I = 1,M                                                     DAX 0035
        DY(I) = DY(I) + DA*DX(I)                                        DAX 0036
   30 CONTINUE                                                          DAX 0037
      IF( N .LT. 4 ) RETURN                                             DAX 0038
   40 MP1 = M + 1                                                       DAX 0039
      DO 50 I = MP1,N,4                                                 DAX 0040
        DY(I) = DY(I) + DA*DX(I)                                        DAX 0041
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)                            DAX 0042
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)                            DAX 0043
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)                            DAX 0044
   50 CONTINUE                                                          DAX 0045
      RETURN                                                            DAX 0046
      END                                                               DAX 0047
