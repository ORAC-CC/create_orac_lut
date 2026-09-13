      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)                        SUM 0001
C                                                                       SUM 0002
C     TAKES THE SUM OF THE ABSOLUTE VALUES.                             SUM 0003
C     JACK DONGARRA, LINPACK, 3/11/78.                                  SUM 0004
C                                                                       SUM 0005
      DOUBLE PRECISION DX(*),DTEMP                                      SUM 0006
      INTEGER I,INCX,M,MP1,N,NINCX                                      SUM 0007
                                                                        SUM 0008
C                                                                       SUM 0009
      DASUM = 0.0D0                                                     SUM 0010
      DTEMP = 0.0D0                                                     SUM 0011
      IF(N.LE.0)RETURN                                                  SUM 0012
      IF(INCX.EQ.1)GO TO 20                                             SUM 0013
C                                                                       SUM 0014
C        CODE FOR INCREMENT NOT EQUAL TO 1                              SUM 0015
C                                                                       SUM 0016
      NINCX = N*INCX                                                    SUM 0017
      DO 10 I = 1,NINCX,INCX                                            SUM 0018
        DTEMP = DTEMP + DABS(DX(I))                                     SUM 0019
   10 CONTINUE                                                          SUM 0020
      DASUM = DTEMP                                                     SUM 0021
      RETURN                                                            SUM 0022
C                                                                       SUM 0023
C        CODE FOR INCREMENT EQUAL TO 1                                  SUM 0024
C                                                                       SUM 0025
C                                                                       SUM 0026
C        CLEAN-UP LOOP                                                  SUM 0027
C                                                                       SUM 0028
   20 M = MOD(N,6)                                                      SUM 0029
                                                                        SUM 0030
                                                                        SUM 0031
                                                                        SUM 0032
                                                                        SUM 0033
      IF( M .EQ. 0 ) GO TO 40                                           SUM 0034
      DO 30 I = 1,M                                                     SUM 0035
                                                                        SUM 0036
                                                                        SUM 0037
        DTEMP = DTEMP + DABS(DX(I))                                     SUM 0038
   30 CONTINUE                                                          SUM 0039
      IF( N .LT. 6 ) GO TO 60                                           SUM 0040
   40 MP1 = M + 1                                                       SUM 0041
      DO 50 I = MP1,N,6                                                 SUM 0042
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2)) SUM 0043
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))           SUM 0044
   50 CONTINUE                                                          SUM 0045
   60 DASUM = DTEMP                                                     SUM 0046
      RETURN                                                            SUM 0047
      END                                                               SUM 0048
