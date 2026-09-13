      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)                             SWP 0001
C                                                                       SWP 0002
C     INTERCHANGES TWO VECTORS.                                         SWP 0003
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.                     SWP 0004
C     JACK DONGARRA, LINPACK, 3/11/78.                                  SWP 0005
C                                                                       SWP 0006
      DOUBLE PRECISION DX(*),DY(*),DTEMP                                SWP 0007
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 SWP 0008
C                                                                       SWP 0009
      IF(N.LE.0)RETURN                                                  SWP 0010
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               SWP 0011
C                                                                       SWP 0012
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL       SWP 0013
C         TO 1                                                          SWP 0014
C                                                                       SWP 0015
      IX = 1                                                            SWP 0016
      IY = 1                                                            SWP 0017
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 SWP 0018
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 SWP 0019
      DO 10 I = 1,N                                                     SWP 0020
        DTEMP = DX(IX)                                                  SWP 0021
        DX(IX) = DY(IY)                                                 SWP 0022
        DY(IY) = DTEMP                                                  SWP 0023
        IX = IX + INCX                                                  SWP 0024
        IY = IY + INCY                                                  SWP 0025
   10 CONTINUE                                                          SWP 0026
      RETURN                                                            SWP 0027
C                                                                       SWP 0028
C       CODE FOR BOTH INCREMENTS EQUAL TO 1                             SWP 0029
C                                                                       SWP 0030
C                                                                       SWP 0031
C       CLEAN-UP LOOP                                                   SWP 0032
C                                                                       SWP 0033
   20 M = MOD(N,3)                                                      SWP 0034
      IF( M .EQ. 0 ) GO TO 40                                           SWP 0035
      DO 30 I = 1,M                                                     SWP 0036
        DTEMP = DX(I)                                                   SWP 0037
        DX(I) = DY(I)                                                   SWP 0038
        DY(I) = DTEMP                                                   SWP 0039
   30 CONTINUE                                                          SWP 0040
      IF( N .LT. 3 ) RETURN                                             SWP 0041
   40 MP1 = M + 1                                                       SWP 0042
      DO 50 I = MP1,N,3                                                 SWP 0043
        DTEMP = DX(I)                                                   SWP 0044
        DX(I) = DY(I)                                                   SWP 0045
        DY(I) = DTEMP                                                   SWP 0046
        DTEMP = DX(I + 1)                                               SWP 0047
        DX(I + 1) = DY(I + 1)                                           SWP 0048
        DY(I + 1) = DTEMP                                               SWP 0049
        DTEMP = DX(I + 2)                                               SWP 0050
        DX(I + 2) = DY(I + 2)                                           SWP 0051
        DY(I + 2) = DTEMP                                               SWP 0052
   50 CONTINUE                                                          SWP 0053
      RETURN                                                            SWP 0054
      END                                                               SWP 0055
