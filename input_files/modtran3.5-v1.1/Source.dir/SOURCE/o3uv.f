      SUBROUTINE O3UV(V,C)                                              O3U 0001
      COMMON /O3UVF/ V1 ,V2 ,DV ,NPT ,S(133)                            O3U 0002
C                                                                       O3U 0003
C     INTERPOLATION  FOR  O3 CONTINUUM WITH LOWTRAN                     O3U 0004
C                                                                       O3U 0005
      C    =0.                                                          O3U 0006
      I=(V  -V1)/DV+1.00001                                             O3U 0007
      IF(I.LT.1   )GO TO 10                                             O3U 0008
      IF(I.GT.NPT )GO TO 10                                             O3U 0009
      VR = (I-1)*DV + V1                                                O3U 0010
      IF(VR. LE. (V+.1) .AND .VR.GE. (V-.1)) GO TO 5                    O3U 0011
      IF(I .EQ. NPT ) I=NPT-1                                           O3U 0012
      AM = (S(I+1) -S(I))/DV                                            O3U 0013
      C0 = S(I) - AM * VR                                               O3U 0014
      C  = AM * V + C0                                                  O3U 0015
      GO TO 10                                                          O3U 0016
5     C    =    S(I)                                                    O3U 0017
10    CONTINUE                                                          O3U 0018
C                                                                       O3U 0019
      RETURN                                                            O3U 0020
      END                                                               O3U 0021
