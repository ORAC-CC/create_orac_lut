      SUBROUTINE  DSCAL(N,DA,DX,INCX)                                   SCL 0001
C                                                                       SCL 0002
C     SCALES A VECTOR BY A CONSTANT.                                    SCL 0003
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.                   SCL 0004
C     JACK DONGARRA, LINPACK, 3/11/78.                                  SCL 0005
C                                                                       SCL 0006
      DOUBLE PRECISION DA,DX(*)                                         SCL 0007
      INTEGER I,INCX,M,MP1,N,NINCX                                      SCL 0008
C                                                                       SCL 0009
      IF(N.LE.0)RETURN                                                  SCL 0010
      IF(INCX.EQ.1)GO TO 20                                             SCL 0011
C                                                                       SCL 0012
C        CODE FOR INCREMENT NOT EQUAL TO 1                              SCL 0013
C                                                                       SCL 0014
      NINCX = N*INCX                                                    SCL 0015
      DO 10 I = 1,NINCX,INCX                                            SCL 0016
        DX(I) = DA*DX(I)                                                SCL 0017
   10 CONTINUE                                                          SCL 0018
      RETURN                                                            SCL 0019
C                                                                       SCL 0020
C        CODE FOR INCREMENT EQUAL TO 1                                  SCL 0021
C                                                                       SCL 0022
C                                                                       SCL 0023
C        CLEAN-UP LOOP                                                  SCL 0024
C                                                                       SCL 0025
   20 M = MOD(N,5)                                                      SCL 0026
      IF( M .EQ. 0 ) GO TO 40                                           SCL 0027
      DO 30 I = 1,M                                                     SCL 0028
        DX(I) = DA*DX(I)                                                SCL 0029
   30 CONTINUE                                                          SCL 0030
      IF( N .LT. 5 ) RETURN                                             SCL 0031
   40 MP1 = M + 1                                                       SCL 0032
      DO 50 I = MP1,N,5                                                 SCL 0033
        DX(I) = DA*DX(I)                                                SCL 0034
        DX(I + 1) = DA*DX(I + 1)                                        SCL 0035
        DX(I + 2) = DA*DX(I + 2)                                        SCL 0036
        DX(I + 3) = DA*DX(I + 3)                                        SCL 0037
        DX(I + 4) = DA*DX(I + 4)                                        SCL 0038
   50 CONTINUE                                                          SCL 0039
      RETURN                                                            SCL 0040
      END                                                               SCL 0041
