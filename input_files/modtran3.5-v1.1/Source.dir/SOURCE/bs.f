      SUBROUTINE BS(I,A,B,N,S)                                          BS  0001
C********************************************************************** BS  0002
      DIMENSION B(9)                                                    BS  0003
C                                                                       BS  0004
C             THIS SUBROUTINE DOES THE BINARY SEARCH FOR THE INDEX I    BS  0005
C             SUCH THAT A IS IN BETWEEN B(I) AND B(I+1)                 BS  0006
C             AND CALCULATES THE INTERPOLATION PARAMETER S              BS  0007
C             SUCH THAT A=S*B(I+1)+(1.-S)*B(I)                          BS  0008
C                                                                       BS  0009
      I=1                                                               BS  0010
      J=N                                                               BS  0011
   10 M=(I+J)/2                                                         BS  0012
      IF(A.LE.B(M)) THEN                                                BS  0013
      J=M                                                               BS  0014
      ELSE                                                              BS  0015
      I=M                                                               BS  0016
      END IF                                                            BS  0017
      IF(J.GT.I+1) GO TO 10                                             BS  0018
      S=(A-B(I))/(B(I+1)-B(I))                                          BS  0019
      RETURN                                                            BS  0020
      END                                                               BS  0021
