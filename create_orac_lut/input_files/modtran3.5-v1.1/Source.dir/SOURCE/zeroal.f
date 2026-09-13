      SUBROUTINE  ZEROAL( ND1, XR0, XR1, TAUC,                          ZAL 0001
     $                    ND2, CMU, CWT, PSIX, EVAL, WK, Z0, Z1, ZJ,    ZAL 0002
     $                    ND3, HLPR, YLM0,                              ZAL 0003
     $                    ND4, ARRAY, CC, EVECC, YLMC,                  ZAL 0004
     $                    ND5, YLMU,                                    ZAL 0005
     $                    ND6, AMB, APB,                                ZAL 0006
     $                    ND7, KK, LL, ZZ, ZPLK0, ZPLK1,                ZAL 0007
     $                    ND8, Z0U, Z1U, ZBEAM,                         ZAL 0008
     $                    ND9, GC,                                      ZAL 0009
     $                    ND10, GU,                                     ZAL 0010
     $                    ND11, Z )                                     ZAL 0011
                                                                        ZAL 0012
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            ZAL 0013
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                ZAL 0014
C            ZERO ARRAYS                                                ZAL 0015
                                                                        ZAL 0016
      REAL*8  AMB(*), APB(*), ARRAY(*), CC(*), CMU(*), CWT(*),          ZAL 0017
     $      EVAL(*), EVECC(*), GC(*), GU(*), HLPR(0:*), KK(*),          ZAL 0018
     $      LL(*), PSIX(*), TAUC(*), WK(*), XR0(*), XR1(*),             ZAL 0019
     $      YLM0(0:*), YLMC(*), YLMU(*), Z(*), Z0(*), Z1(*),            ZAL 0020
     $      Z0U(*), Z1U(*), ZJ(*), ZZ(*), ZPLK0(*), ZPLK1(*),           ZAL 0021
     $      ZBEAM(*)                                                    ZAL 0022
                                                                        ZAL 0023
                                                                        ZAL 0024
      DO 2  N = 1, ND1                                                  ZAL 0025
         XR0(N) = 0.0                                                   ZAL 0026
         XR1(N) = 0.0                                                   ZAL 0027
         TAUC(N) = 0.0                                                  ZAL 0028
 2    CONTINUE                                                          ZAL 0029
                                                                        ZAL 0030
      DO 4  N = 1, ND2                                                  ZAL 0031
         CMU(N) = 0.0                                                   ZAL 0032
         CWT(N) = 0.0                                                   ZAL 0033
         PSIX(N) = 0.0                                                  ZAL 0034
         EVAL(N) = 0.0                                                  ZAL 0035
         WK(N) = 0.0                                                    ZAL 0036
         Z0(N) = 0.0                                                    ZAL 0037
         Z1(N) = 0.0                                                    ZAL 0038
         ZJ(N) = 0.0                                                    ZAL 0039
 4    CONTINUE                                                          ZAL 0040
                                                                        ZAL 0041
      DO 6  N = 1, ND3                                                  ZAL 0042
         HLPR(N) = 0.0                                                  ZAL 0043
         YLM0(N) = 0.0                                                  ZAL 0044
 6    CONTINUE                                                          ZAL 0045
                                                                        ZAL 0046
      DO 8  N = 1, ND4                                                  ZAL 0047
         ARRAY(N) = 0.0                                                 ZAL 0048
         CC(N) = 0.0                                                    ZAL 0049
         EVECC(N) = 0.0                                                 ZAL 0050
 8    CONTINUE                                                          ZAL 0051
                                                                        ZAL 0052
      DO 10  N = 1, ND5                                                 ZAL 0053
         YLMC(N) = 0.0                                                  ZAL 0054
         YLMU(N) = 0.0                                                  ZAL 0055
 10   CONTINUE                                                          ZAL 0056
                                                                        ZAL 0057
      DO 12  N = 1, ND6                                                 ZAL 0058
         AMB(N) = 0.0                                                   ZAL 0059
         APB(N) = 0.0                                                   ZAL 0060
 12   CONTINUE                                                          ZAL 0061
                                                                        ZAL 0062
      DO 14  N = 1, ND7                                                 ZAL 0063
         KK(N) = 0.0                                                    ZAL 0064
         LL(N) = 0.0                                                    ZAL 0065
         ZZ(N) = 0.0                                                    ZAL 0066
         ZPLK0(N) = 0.0                                                 ZAL 0067
         ZPLK1(N) = 0.0                                                 ZAL 0068
 14   CONTINUE                                                          ZAL 0069
                                                                        ZAL 0070
      DO 16  N = 1, ND8                                                 ZAL 0071
         Z0U(N) = 0.0                                                   ZAL 0072
         Z1U(N) = 0.0                                                   ZAL 0073
         ZBEAM(N) = 0.0                                                 ZAL 0074
 16   CONTINUE                                                          ZAL 0075
                                                                        ZAL 0076
      DO 18  N = 1, ND9                                                 ZAL 0077
         GC(N) = 0.0                                                    ZAL 0078
 18   CONTINUE                                                          ZAL 0079
                                                                        ZAL 0080
      DO 20  N = 1, ND10                                                ZAL 0081
         GU(N) = 0.0                                                    ZAL 0082
 20   CONTINUE                                                          ZAL 0083
                                                                        ZAL 0084
      DO 22  N = 1, ND11                                                ZAL 0085
         Z(N) = 0.0                                                     ZAL 0086
 22   CONTINUE                                                          ZAL 0087
                                                                        ZAL 0088
      RETURN                                                            ZAL 0089
      END                                                               ZAL 0090
