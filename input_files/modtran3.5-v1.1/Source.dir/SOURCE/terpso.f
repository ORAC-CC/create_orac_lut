      SUBROUTINE  TERPSO( CWT, DELM0, FBEAM, GL, MAZIM, MXCMU,          TSO 0001
     $                      PLANK, NUMU, NSTR, OPRIM, PI, YLM0, YLMC,   TSO 0002
     $                      YLMU, PSIX, XR0, XR1, Z0, ZJ, ZBEAM, Z0U,   TSO 0003
     $                      Z1U, Z0UMS, Z1UMS, BEAMMS )                 TSO 0004
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            TSO 0005
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                TSO 0006
C       LOWER CASE VARIABLES ADDED                                      TSO 0007
                                                                        TSO 0008
C         INTERPOLATES SOURCE FUNCTIONS TO USER ANGLES                  TSO 0009
                                                                        TSO 0010
C    I N P U T      V A R I A B L E S:                                  TSO 0011
                                                                        TSO 0012
C       CWT    :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE        TSO 0013
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0                         TSO 0014
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTIONTSO 0015
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO) TSO 0016
C       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT                          TSO 0017
C       OPRIM  :  SINGLE SCATTERING ALBEDO                              TSO 0018
C       XR0    :  EXPANSION OF THERMAL SOURCE FUNCTION                  TSO 0019
C       XR1    :  EXPANSION OF THERMAL SOURCE FUNCTION EQS.SS(14-16)    TSO 0020
C       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             TSO 0021
C                 AT THE BEAM ANGLE                                     TSO 0022
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             TSO 0023
C                 AT THE QUADRATURE ANGLES                              TSO 0024
C       YLMU   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             TSO 0025
C                 AT THE USER ANGLES                                    TSO 0026
C       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)             TSO 0027
C       ZJ     :  SOLUTION VECTOR CAPITAL -Z-SUB-ZERO AFTER SOLVING     TSO 0028
C                 EQ. SS(19)                                            TSO 0029
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        TSO 0030
                                                                        TSO 0031
C    O U T P U T     V A R I A B L E S:                                 TSO 0032
                                                                        TSO 0033
C       ZBEAM  :  INCIDENT-BEAM SOURCE FUNCTION AT USER ANGLES          TSO 0034
C       Z0U,Z1U:  COMPONENTS OF A LINEAR-IN-OPTICAL-DEPTH-DEPENDENT     TSO 0035
C                    SOURCE (APPROXIMATING THE PLANCK EMISSION SOURCE)  TSO 0036
                                                                        TSO 0037
C   I N T E R N A L       V A R I A B L E S:                            TSO 0038
                                                                        TSO 0039
C       PSIX    :  SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)          TSO 0040
C+---------------------------------------------------------------------+TSO 0041
      LOGICAL  PLANK                                                    TSO 0042
      REAL*8   CWT(*), GL(0:*), PSIX(*),  YLM0(0:*), YLMC( 0:MXCMU,* ), TSO 0043
     $         YLMU( 0:MXCMU,*), Z0(*), ZJ(*), ZBEAM(*), Z0U(*),        TSO 0044
     $         Z1U(*), BEAMMS(*), Z1UMS(*), Z0UMS(*), SUM1              TSO 0045
C              LOWER CASE VARIABLES ADDED                               TSO 0046
                                                                        TSO 0047
                                                                        TSO 0048
      IF ( FBEAM.GT.0.0 )  THEN                                         TSO 0049
C                                  ** BEAM SOURCE TERMS; EQ. SD(9)      TSO 0050
         DO 20  IQ = MAZIM, NSTR-1                                      TSO 0051
            PSUM = 0.                                                   TSO 0052
            DO 10  JQ = 1, NSTR                                         TSO 0053
               PSUM = PSUM + CWT(JQ) * YLMC(IQ,JQ) * ZJ(JQ)             TSO 0054
10          CONTINUE                                                    TSO 0055
            PSIX(IQ+1) = 0.5 * GL(IQ) * PSUM                            TSO 0056
20       CONTINUE                                                       TSO 0057
                                                                        TSO 0058
         FACT = ( 2. - DELM0 ) * FBEAM / (4.0*PI)                       TSO 0059
         DO 40  IU = 1, NUMU                                            TSO 0060
            SUM = 0.                                                    TSO 0061
C             LOWER CASE VARIABLES ADDED                                TSO 0062
              SUM1= 0.                                                  TSO 0063
            DO 30 IQ = MAZIM, NSTR-1                                    TSO 0064
               SUM = SUM + YLMU(IQ,IU) *                                TSO 0065
     $                    ( PSIX(IQ+1) + FACT * GL(IQ) * YLM0(IQ) )     TSO 0066
C                LOWER CASE VARIABLES ADDED                             TSO 0067
                 SUM1 = SUM1 + YLMU(IQ,IU) * ( PSIX(IQ+1) )             TSO 0068
30            CONTINUE                                                  TSO 0069
            ZBEAM(IU) = SUM                                             TSO 0070
C             LOWER CASE VARIABLES ADDED                                TSO 0071
              BEAMMS(IU)=SUM1                                           TSO 0072
40         CONTINUE                                                     TSO 0073
      END IF                                                            TSO 0074
                                                                        TSO 0075
      IF ( PLANK .AND. MAZIM.EQ.0 )  THEN                               TSO 0076
                                                                        TSO 0077
C                                ** THERMAL SOURCE TERMS, STWJ(27C)     TSO 0078
         DO 80  IQ = MAZIM, NSTR-1                                      TSO 0079
            PSUM = 0.0                                                  TSO 0080
            DO 70  JQ = 1, NSTR                                         TSO 0081
               PSUM = PSUM + CWT(JQ) * YLMC(IQ,JQ) * Z0(JQ)             TSO 0082
 70           CONTINUE                                                  TSO 0083
            PSIX(IQ+1) = 0.5 * GL(IQ) * PSUM                            TSO 0084
 80        CONTINUE                                                     TSO 0085
                                                                        TSO 0086
         DO 100  IU = 1, NUMU                                           TSO 0087
            SUM = 0.0                                                   TSO 0088
            DO 90   IQ = MAZIM, NSTR-1                                  TSO 0089
               SUM = SUM + YLMU(IQ,IU) * PSIX(IQ+1)                     TSO 0090
90            CONTINUE                                                  TSO 0091
C             LOWER CASE VARIABLES ADDED                                TSO 0092
              Z0UMS(IU)=SUM                                             TSO 0093
              Z1UMS(IU)=OPRIM*XR1                                       TSO 0094
            Z0U(IU) = SUM + (1.-OPRIM) * XR0                            TSO 0095
            Z1U(IU) = XR1                                               TSO 0096
100        CONTINUE                                                     TSO 0097
                                                                        TSO 0098
      END IF                                                            TSO 0099
                                                                        TSO 0100
      RETURN                                                            TSO 0101
      END                                                               TSO 0102
