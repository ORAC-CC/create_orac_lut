      FUNCTION   PF(NN,I,J)                                             PF  0001
C     CALL THE APPROPRIATE PHASE FUNCTION                               PF  0002
      COMMON/MNMPHS/ MNUM(27,26),PHSFNC(34,70)                          PF  0003
      M=MNUM(I,NN)                                                      PF  0004
      PF=PHSFNC(J,M)                                                    PF  0005
      RETURN                                                            PF  0006
      END                                                               PF  0007
