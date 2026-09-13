      SUBROUTINE  SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,            SAT 0001
     $                    NN, NSTR, TAUCPR, SFLUP, SFLDN )              SAT 0002
                                                                        SAT 0003
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            SAT 0004
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                SAT 0005
C       CALCULATES SPHERICAL ALBEDO AND TRANSMISSIVITY FOR THE ENTIRE   SAT 0006
C       MEDIUM FROM THE M=0 INTENSITY COMPONENTS                        SAT 0007
C       (THIS IS A VERY SPECIALIZED VERSION OF 'FLUXES')                SAT 0008
                                                                        SAT 0009
C    I N P U T    V A R I A B L E S:                                    SAT 0010
                                                                        SAT 0011
C       CMU     :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE     SAT 0012
C       CWT     :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE       SAT 0013
C       KK      :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)            SAT 0014
C       GC      :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)       SAT 0015
C       LL      :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED      SAT 0016
C                  BY SOLVING SCALED VERSION OF EQ. SC(5);              SAT 0017
C                  EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED          SAT 0018
C       NN      :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)            SAT 0019
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        SAT 0020
                                                                        SAT 0021
C    O U T P U T   V A R I A B L E S:                                   SAT 0022
                                                                        SAT 0023
C       SFLUP   :  UP-FLUX AT TOP (EQUIVALENT TO SPHERICAL ALBEDO DUE TOSAT 0024
C                  RECIPROCITY).  FOR ILLUMINATION FROM BELOW IT GIVES  SAT 0025
C                  SPHERICAL TRANSMISSIVITY                             SAT 0026
C       SFLDN   :  DOWN-FLUX AT BOTTOM (FOR SINGLE LAYER                SAT 0027
C                  EQUIVALENT TO SPHERICAL TRANSMISSIVITY               SAT 0028
C                  DUE TO RECIPROCITY)                                  SAT 0029
                                                                        SAT 0030
C    I N T E R N A L   V A R I A B L E S:                               SAT 0031
                                                                        SAT 0032
C       ZINT    :  INTENSITY OF M=0 CASE, IN EQ. SC(1)                  SAT 0033
C+----------------------------------------------------------------------SAT 0034
                                                                        SAT 0035
      REAL*8  CMU(*), CWT(*), GC( MXCMU,MXCMU,* ), KK( MXCMU,* ),       SAT 0036
     $      LL( MXCMU,* ), TAUCPR( 0:* )                                SAT 0037
                                                                        SAT 0038
                                                                        SAT 0039
      SFLUP = 0.0                                                       SAT 0040
      DO 20  IQ = NN+1, NSTR                                            SAT 0041
         ZINT  = 0.0                                                    SAT 0042
         DO 10   JQ = 1, NN                                             SAT 0043
            ZINT = ZINT + GC(IQ,JQ,1) * LL(JQ,1) *                      SAT 0044
     $                    DEXP( KK(JQ,1) * TAUCPR(1) )                  SAT 0045
10       CONTINUE                                                       SAT 0046
         DO 11  JQ = NN+1, NSTR                                         SAT 0047
            ZINT = ZINT + GC(IQ,JQ,1) * LL(JQ,1)                        SAT 0048
11       CONTINUE                                                       SAT 0049
                                                                        SAT 0050
         SFLUP = SFLUP + CWT(IQ-NN) * CMU(IQ-NN) * ZINT                 SAT 0051
20    CONTINUE                                                          SAT 0052
                                                                        SAT 0053
      SFLDN  = 0.0                                                      SAT 0054
      DO 40  IQ = 1, NN                                                 SAT 0055
         ZINT   = 0.0                                                   SAT 0056
         DO 30  JQ = 1, NN                                              SAT 0057
             ZINT = ZINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)                 SAT 0058
30       CONTINUE                                                       SAT 0059
         DO 31  JQ = NN+1, NSTR                                         SAT 0060
             ZINT = ZINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR) *               SAT 0061
     $         DEXP( - KK(JQ,NLYR)*(TAUCPR(NLYR) - TAUCPR(NLYR-1)) )    SAT 0062
31       CONTINUE                                                       SAT 0063
                                                                        SAT 0064
         SFLDN = SFLDN + CWT(NN+1-IQ) * CMU(NN+1-IQ) * ZINT             SAT 0065
40    CONTINUE                                                          SAT 0066
                                                                        SAT 0067
      SFLUP = 2.0 * SFLUP                                               SAT 0068
      SFLDN = 2.0 * SFLDN                                               SAT 0069
                                                                        SAT 0070
      RETURN                                                            SAT 0071
      END                                                               SAT 0072
