      SUBROUTINE O2INT(V1C,V1,DV,NPT,C,CARRAY,A,AARRAY,B,BARRAY)        O2I 0001
C                                                                       O2I 0002
C     INTERPOLATION FOR O2 PRESSURE INDUCED CONTINUUM, NECESSARY FOR    O2I 0003
C          LOWTRAN7 FORMULATION  (MODELED AFTER THE LOWTRAN UV-O3 BANDS)O2I 0004
C                                                                       O2I 0005
      DIMENSION CARRAY(74),AARRAY(74),BARRAY(74)                        O2I 0006
      C=0.                                                              O2I 0007
      A=0.                                                              O2I 0008
      B=0.                                                              O2I 0009
      I=(V1C-V1)/DV+1.00001                                             O2I 0010
      IF(I.LT.1  )GO TO 10                                              O2I 0011
      IF(I.GT.NPT)GO TO 10                                              O2I 0012
      C=CARRAY(I)                                                       O2I 0013
      A=AARRAY(I)                                                       O2I 0014
      B=BARRAY(I)                                                       O2I 0015
10    CONTINUE                                                          O2I 0016
      RETURN                                                            O2I 0017
      END                                                               O2I 0018
