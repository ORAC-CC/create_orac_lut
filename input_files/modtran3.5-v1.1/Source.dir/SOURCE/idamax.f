      INTEGER FUNCTION IDAMAX(N,DX,INCX)                                IMX 0001
C                                                                       IMX 0002
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            IMX 0003
C     JACK DONGARRA, LINPACK, 3/11/78.                                  IMX 0004
C                                                                       IMX 0005
      DOUBLE PRECISION DX(*),DMAX                                       IMX 0006
      INTEGER I,INCX,IX,N                                               IMX 0007
C                                                                       IMX 0008
      IDAMAX = 0                                                        IMX 0009
      IF( N .LT. 1 ) RETURN                                             IMX 0010
      IDAMAX = 1                                                        IMX 0011
      IF(N.EQ.1)RETURN                                                  IMX 0012
      IF(INCX.EQ.1)GO TO 20                                             IMX 0013
C                                                                       IMX 0014
C        CODE FOR INCREMENT NOT EQUAL TO 1                              IMX 0015
C                                                                       IMX 0016
      IX = 1                                                            IMX 0017
      DMAX = DABS(DX(1))                                                IMX 0018
      IX = IX + INCX                                                    IMX 0019
      DO 10 I = 2,N                                                     IMX 0020
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5                               IMX 0021
         IDAMAX = I                                                     IMX 0022
         DMAX = DABS(DX(IX))                                            IMX 0023
    5    IX = IX + INCX                                                 IMX 0024
   10 CONTINUE                                                          IMX 0025
      RETURN                                                            IMX 0026
C                                                                       IMX 0027
C        CODE FOR INCREMENT EQUAL TO 1                                  IMX 0028
C                                                                       IMX 0029
   20 DMAX = DABS(DX(1))                                                IMX 0030
      DO 30 I = 2,N                                                     IMX 0031
         IF(DABS(DX(I)).LE.DMAX) GO TO 30                               IMX 0032
         IDAMAX = I                                                     IMX 0033
         DMAX = DABS(DX(I))                                             IMX 0034
   30 CONTINUE                                                          IMX 0035
      RETURN                                                            IMX 0036
      END                                                               IMX 0037
