      SUBROUTINE O3HHT1(V,C)                                            1O3 0001
C     SUBROUTINE O3HHT1(V1C,V2C,DVC,NPTC,C)                             1O3 0002
      COMMON /O3HH1/ V1S,V2S,DVS,NPTS,S(2687)                           1O3 0003
C                                                                       1O3 0004
      CALL O3INT(V ,V1S,DVS,NPTS,S,C)                                   1O3 0005
C                                                                       1O3 0006
      RETURN                                                            1O3 0007
      END                                                               1O3 0008
