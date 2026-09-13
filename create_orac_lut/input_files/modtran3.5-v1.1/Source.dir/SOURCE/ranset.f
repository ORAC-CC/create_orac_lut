      SUBROUTINE RANSET(I)                                              RNS 0001
      COMMON/RST/II                                                     RNS 0002
CSSI  RANET=RAN(I)                                                      RNS 0003
      RANET=RANDOM(I)                                                   RNS 0004
      II=RANET                                                          RNS 0005
      RETURN                                                            RNS 0006
      END                                                               RNS 0007
CC    FUNCTION   RAYSCT(V)                                              RNS 0008
C     RADIATION FLD OUT                                                 RNS 0009
C **  MOLECULAR SCATTERING                                              RNS 0010
CC    RAYSCT=0.                                                         RNS 0011
CC    IF(V.LE.3000.) RETURN                                             RNS 0012
CC    RAYSCT=V**3/(9.26799E+18-1.07123E+09*V**2)                        RNS 0013
C     V**4 FOR RADIATION FLD IN                                         RNS 0014
CC    RETURN                                                            RNS 0015
CC    END                                                               RNS 0016
