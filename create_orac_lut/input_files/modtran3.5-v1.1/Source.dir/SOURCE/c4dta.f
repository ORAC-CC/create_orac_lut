      SUBROUTINE C4DTA (C4L,V)                                          C4D 0001
C **  N2 CONTINUUM                                                      C4D 0002
      COMMON /C4C8/ C4(133),C8(102)                                     C4D 0003
      C4L=0.                                                            C4D 0004
      IF(V.LT.2080.) RETURN                                             C4D 0005
      IF(V.GT.2740.) RETURN                                             C4D 0006
      IV=V                                                              C4D 0007
      L=(IV-2080)/5+1                                                   C4D 0008
      C4L=C4(L)                                                         C4D 0009
      RETURN                                                            C4D 0010
      END                                                               C4D 0011
