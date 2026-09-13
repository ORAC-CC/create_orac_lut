      SUBROUTINE FRN296(V1C,FH2O)                                       FRN 0001
C     LOADS FOREIGN CONTINUUM  296K                                     FRN 0002
      COMMON /FH2O/ V1,V2,DV,NPT,F296(2003)                             FRN 0003
      CALL SINT(V1,V1C,DV,NPT,F296,FH2O)                                FRN 0004
      RETURN                                                            FRN 0005
      END                                                               FRN 0006
