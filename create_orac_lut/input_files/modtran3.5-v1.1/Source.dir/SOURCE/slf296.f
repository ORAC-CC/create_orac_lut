      SUBROUTINE SLF296(V1C,SH2OT0)                                     S96 0001
C     LOADS SELF CONTINUUM  296K                                        S96 0002
      COMMON /SH2O/ V1,V2,DV,NPT,S296(2003)                             S96 0003
      CALL SINT(V1,V1C,DV,NPT,S296,SH2OT0)                              S96 0004
      RETURN                                                            S96 0005
      END                                                               S96 0006
