      SUBROUTINE SLF260(V1C,SH2OT1)                                     S60 0001
C     LOADS SELF CONTINUUM  260K                                        S60 0002
      COMMON /S260/ V1,V2,DV,NPT,S260(2003)                             S60 0003
      CALL SINT(V1,V1C,DV,NPT,S260,SH2OT1)                              S60 0004
      RETURN                                                            S60 0005
      END                                                               S60 0006
