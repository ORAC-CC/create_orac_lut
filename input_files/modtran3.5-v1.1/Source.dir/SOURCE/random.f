      FUNCTION RANDOM(I)                                                RAN 0001
      DOUBLE PRECISION DRND,Z                                           RAN 0002
      DATA DRND/0.574914365D0/                                          RAN 0003
      DRND=97301.D0*DRND                                                RAN 0004
      K=DRND                                                            RAN 0005
      Z=K                                                               RAN 0006
      DRND=DRND-Z                                                       RAN 0007
      RANDOM=DRND                                                       RAN 0008
      RETURN                                                            RAN 0009
      END                                                               RAN 0010
