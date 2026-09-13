      FUNCTION   AITK(ARG,VAL,X,NDIM)                                   AIT 0001
C                                                                       AIT 0002
C      IBM SCIENTIFIC SUBROUTINE                                        AIT 0003
C     AITKEN INTERPOLATION ROUTINE                                      AIT 0004
C                                                                       AIT 0005
      DIMENSION ARG(NDIM),VAL(NDIM)                                     AIT 0006
      IF(NDIM-1)9,7,1                                                   AIT 0007
C                                                                       AIT 0008
C     START OF AITKEN-LOOP                                              AIT 0009
    1 DO 6 J=2,NDIM                                                     AIT 0010
      IEND=J-1                                                          AIT 0011
      DO 2 I=1,IEND                                                     AIT 0012
      H=ARG(I)-ARG(J)                                                   AIT 0013
      IF(H)2,13,2                                                       AIT 0014
    2 VAL(J)=(VAL(I)*(X-ARG(J))-VAL(J)*(X-ARG(I)))/H                    AIT 0015
    6 CONTINUE                                                          AIT 0016
C     END OF AITKEN-LOOP                                                AIT 0017
C                                                                       AIT 0018
    7 J=NDIM                                                            AIT 0019
    8 AITK=VAL(J)                                                       AIT 0020
    9 RETURN                                                            AIT 0021
C                                                                       AIT 0022
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             AIT 0023
   13 CONTINUE                                                          AIT 0024
CCC   IER=3                                                             AIT 0025
      J=IEND                                                            AIT 0026
      GO TO 8                                                           AIT 0027
      END                                                               AIT 0028
