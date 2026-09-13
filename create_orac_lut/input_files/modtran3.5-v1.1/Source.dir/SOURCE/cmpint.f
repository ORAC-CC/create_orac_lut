      SUBROUTINE  CMPINT( FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM,      CMP 0001
     $                    MXCMU, MXULV, MXUMU, NCUT, NN, NSTR,          CMP 0002
     $                    PLANK, NTAU, TAUCPR, UMU0, UTAUPR,            CMP 0003
     $                    ZZ, ZPLK0, ZPLK1, UUM )                       CMP 0004
                                                                        CMP 0005
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            CMP 0006
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                CMP 0007
C       CALCULATES THE FOURIER INTENSITY COMPONENTS AT THE QUADRATURE   CMP 0008
C       ANGLES FOR AZIMUTHAL EXPANSION TERMS (MAZIM) IN EQ. SD(2)       CMP 0009
                                                                        CMP 0010
C    I N P U T    V A R I A B L E S:                                    CMP 0011
                                                                        CMP 0012
C       KK      :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)            CMP 0013
C       GC      :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)       CMP 0014
C       LL      :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED      CMP 0015
C                  BY SOLVING SCALED VERSION OF EQ. SC(5);              CMP 0016
C                  EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED          CMP 0017
C       LYRCUT  :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER         CMP 0018
C       MAZIM   :  ORDER OF AZIMUTHAL COMPONENT                         CMP 0019
C       NCUT    :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION       CMP 0020
C                  OPTICAL DEPTH EXCEEDS -ABSCUT-                       CMP 0021
C       NN      :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)            CMP 0022
C       TAUCPR  :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)            CMP 0023
C       UTAUPR  :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M      CMP 0024
C                  COORDINATES;  EQUAL TO -UTAU- IF NO DELTA-M          CMP 0025
C       ZZ      :  BEAM SOURCE VECTORS IN EQ. SS(19)                    CMP 0026
C       ZPLK0   :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)   CMP 0027
C       ZPLK1   :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)   CMP 0028
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        CMP 0029
                                                                        CMP 0030
C    O U T P U T   V A R I A B L E S:                                   CMP 0031
                                                                        CMP 0032
C       UUM     :  FOURIER COMPONENTS OF THE INTENSITY IN EQ.  SD(12)   CMP 0033
C                   ( AT POLAR QUADRATURE ANGLES )                      CMP 0034
                                                                        CMP 0035
C    I N T E R N A L   V A R I A B L E S:                               CMP 0036
                                                                        CMP 0037
C       FACT    :  EXP( - UTAUPR / UMU0 )                               CMP 0038
C       ZINT    :  INTENSITY OF M=0 CASE, IN EQ. SC(1)                  CMP 0039
C+----------------------------------------------------------------------CMP 0040
                                                                        CMP 0041
       LOGICAL  LYRCUT, PLANK                                           CMP 0042
       INTEGER  LAYRU(*)                                                CMP 0043
       REAL*8     UUM( MXUMU, MXULV, 0:* )                              CMP 0044
       REAL*8     GC( MXCMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),    CMP 0045
     $          TAUCPR( 0:* ), UTAUPR(*), ZZ( MXCMU, *),                CMP 0046
     $          ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* )                      CMP 0047
                                                                        CMP 0048
                                                                        CMP 0049
C                                                  ** ZERO OUTPUT ARRAY CMP 0050
       CALL ZEROIT( UUM, MXUMU*MXULV*(MXCMU + 1) )                      CMP 0051
                                                                        CMP 0052
C                                       ** LOOP OVER USER LEVELS        CMP 0053
       DO 100  LU = 1, NTAU                                             CMP 0054
                                                                        CMP 0055
          LYU = LAYRU(LU)                                               CMP 0056
          IF ( LYRCUT .AND. LYU.GT.NCUT )  GO TO 100                    CMP 0057
                                                                        CMP 0058
          DO 20  IQ = 1, NSTR                                           CMP 0059
             ZINT = 0.0                                                 CMP 0060
             DO 10  JQ = 1, NN                                          CMP 0061
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               CMP 0062
     $                 DEXP( - KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)) )  CMP 0063
10           CONTINUE                                                   CMP 0064
             DO 11  JQ = NN+1, NSTR                                     CMP 0065
                ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *              CMP 0066
     $              DEXP( - KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU-1)) )   CMP 0067
11           CONTINUE                                                   CMP 0068
                                                                        CMP 0069
             UUM(IQ,LU,MAZIM) = ZINT                                    CMP 0070
             IF ( FBEAM.GT.0.0 )                                        CMP 0071
     $            UUM(IQ,LU,MAZIM) = ZINT + ZZ(IQ,LYU)                  CMP 0072
     $                                    * DEXP( - UTAUPR(LU) / UMU0 ) CMP 0073
             IF ( PLANK .AND. MAZIM.EQ.0 )                              CMP 0074
     $            UUM(IQ,LU,MAZIM) = UUM(IQ,LU,MAZIM) + ZPLK0(IQ,LYU) + CMP 0075
     $                             ZPLK1(IQ,LYU) * UTAUPR(LU)           CMP 0076
20        CONTINUE                                                      CMP 0077
                                                                        CMP 0078
100   CONTINUE                                                          CMP 0079
                                                                        CMP 0080
      RETURN                                                            CMP 0081
      END                                                               CMP 0082
