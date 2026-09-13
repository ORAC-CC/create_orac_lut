      SUBROUTINE  UPISOT( ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR, OPRIM, UPI 0001
     $                    WK, XR0, XR1, Z0, Z1, ZPLK0, ZPLK1 )          UPI 0002
                                                                        UPI 0003
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            UPI 0004
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                UPI 0005
C       FINDS THE PARTICULAR SOLUTION OF THERMAL RADIATION OF SS(15)    UPI 0006
                                                                        UPI 0007
C     ROUTINES CALLED:  DGECO, DGESL                                    UPI 0008
                                                                        UPI 0009
C   I N P U T     V A R I A B L E S:                                    UPI 0010
                                                                        UPI 0011
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)                         UPI 0012
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE      UPI 0013
C       OPRIM  :  DELTA-M SCALED SINGLE SCATTERING ALBEDO               UPI 0014
C       XR0    :  EXPANSION OF THERMAL SOURCE FUNCTION                  UPI 0015
C       XR1    :  EXPANSION OF THERMAL SOURCE FUNCTION EQS. SS(14-16)   UPI 0016
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        UPI 0017
                                                                        UPI 0018
C    O U T P U T    V A R I A B L E S:                                  UPI 0019
                                                                        UPI 0020
C       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)             UPI 0021
C       Z1     :  SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)             UPI 0022
C       ZPLK0, :  PERMANENT STORAGE FOR -Z0,Z1-, BUT RE-ORDERED         UPI 0023
C        ZPLK1                                                          UPI 0024
                                                                        UPI 0025
C   I N T E R N A L    V A R I A B L E S:                               UPI 0026
                                                                        UPI 0027
C       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(16)    UPI 0028
C       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK* UPI 0029
C       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*                   UPI 0030
C+---------------------------------------------------------------------+UPI 0031
                                                                        UPI 0032
      INTEGER IPVT(*)                                                   UPI 0033
      REAL*8    ARRAY( MXCMU,* ), CC( MXCMU,* ), CMU(*), WK(*),         UPI 0034
     $        Z0(*), Z1(*), ZPLK0(*), ZPLK1(*)                          UPI 0035
                                                                        UPI 0036
                                                                        UPI 0037
      DO 20 IQ = 1, NSTR                                                UPI 0038
                                                                        UPI 0039
         DO 10 JQ = 1, NSTR                                             UPI 0040
            ARRAY(IQ,JQ) = - CC(IQ,JQ)                                  UPI 0041
10       CONTINUE                                                       UPI 0042
         ARRAY(IQ,IQ) = 1.0 + ARRAY(IQ,IQ)                              UPI 0043
                                                                        UPI 0044
         Z1(IQ) = XR1                                                   UPI 0045
         Z0(IQ) = (1.-OPRIM) * XR0 + CMU(IQ) * Z1(IQ)                   UPI 0046
20    CONTINUE                                                          UPI 0047
C                       ** SOLVE LINEAR EQUATIONS: SAME AS IN *UPBEAM*, UPI 0048
C                       ** EXCEPT -ZJ- REPLACED BY -Z0-                 UPI 0049
      RCOND = 0.0                                                       UPI 0050
      CALL  DGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )                UPI 0051
      IF ( 1.0+RCOND .EQ. 1.0 )  CALL  ERRMSG                           UPI 0052
     $   ( 'UPISOT--dGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)           UPI 0053
                                                                        UPI 0054
      CALL  DGESL( ARRAY, MXCMU, NSTR, IPVT, Z0, 0 )                    UPI 0055
                                                                        UPI 0056
      DO 30  IQ = 1, NN                                                 UPI 0057
         ZPLK0( IQ+NN )   = Z0( IQ )                                    UPI 0058
         ZPLK1( IQ+NN )   = Z1( IQ )                                    UPI 0059
         ZPLK0( NN+1-IQ ) = Z0( IQ+NN )                                 UPI 0060
         ZPLK1( NN+1-IQ ) = Z1( IQ+NN )                                 UPI 0061
30    CONTINUE                                                          UPI 0062
                                                                        UPI 0063
      RETURN                                                            UPI 0064
      END                                                               UPI 0065
