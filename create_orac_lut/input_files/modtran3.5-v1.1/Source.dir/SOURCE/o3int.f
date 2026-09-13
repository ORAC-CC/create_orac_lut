      SUBROUTINE O3INT(V1C,V1,DV,NPT,CONTI,CONTO)                       O3I 0001
C                                                                       O3I 0002
C     INTERPOLATION  FOR  O3 CONTINUUM WITH LOWTRAN                     O3I 0003
C                                                                       O3I 0004
      DIMENSION CONTI(2687)                                             O3I 0005
      CONTO=0.                                                          O3I 0006
      I=(V1C-V1)/DV+1.00001                                             O3I 0007
      IF(I.LT.1  )GO TO 10                                              O3I 0008
      IF(I.GT.NPT)GO TO 10                                              O3I 0009
      CONTO=CONTI(I)                                                    O3I 0010
10    CONTINUE                                                          O3I 0011
      RETURN                                                            O3I 0012
      END                                                               O3I 0013
