      SUBROUTINE  UPBEAM( ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM,UPB 0001
     $                    MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ,UPB 0002
     $                    ZZ )                                          UPB 0003
                                                                        UPB 0004
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            UPB 0005
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                UPB 0006
C         FINDS THE INCIDENT-BEAM PARTICULAR SOLUTION  OF SS(18)        UPB 0007
                                                                        UPB 0008
C     ROUTINES CALLED:  DGECO, DGESL                                    UPB 0009
                                                                        UPB 0010
C   I N P U T    V A R I A B L E S:                                     UPB 0011
                                                                        UPB 0012
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)                         UPB 0013
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE      UPB 0014
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0                         UPB 0015
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTIONUPB 0016
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO) UPB 0017
C       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT                          UPB 0018
C       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             UPB 0019
C                 AT THE BEAM ANGLE                                     UPB 0020
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             UPB 0021
C                 AT THE QUADRATURE ANGLES                              UPB 0022
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        UPB 0023
                                                                        UPB 0024
C   O U T P U T    V A R I A B L E S:                                   UPB 0025
                                                                        UPB 0026
C       ZJ     :  RIGHT-HAND SIDE VECTOR CAPITAL-X-SUB-ZERO IN SS(19);  UPB 0027
C                 ALSO THE SOLUTION VECTOR CAPITAL-Z-SUB-ZERO           UPB 0028
C                 AFTER SOLVING THAT SYSTEM                             UPB 0029
C       ZZ     :  PERMANENT STORAGE FOR -ZJ-, BUT RE-ORDERED            UPB 0030
                                                                        UPB 0031
C   I N T E R N A L    V A R I A B L E S:                               UPB 0032
                                                                        UPB 0033
C       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(19)    UPB 0034
C       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK* UPB 0035
C       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*                   UPB 0036
C+---------------------------------------------------------------------+UPB 0037
                                                                        UPB 0038
      INTEGER  IPVT(*)                                                  UPB 0039
      REAL*8     ARRAY( MXCMU,* ), CC( MXCMU,* ), CMU(*), GL(0:*),      UPB 0040
     $         WK(*), YLM0(0:*), YLMC( 0:MXCMU,* ), ZJ(*), ZZ(*)        UPB 0041
                                                                        UPB 0042
      DO 40  IQ = 1, NSTR                                               UPB 0043
                                                                        UPB 0044
         DO 10  JQ = 1, NSTR                                            UPB 0045
            ARRAY(IQ,JQ) = - CC(IQ,JQ)                                  UPB 0046
10       CONTINUE                                                       UPB 0047
         ARRAY(IQ,IQ) = 1. + CMU(IQ) / UMU0 + ARRAY(IQ,IQ)              UPB 0048
                                                                        UPB 0049
         SUM = 0.                                                       UPB 0050
         DO 20  K = MAZIM, NSTR-1                                       UPB 0051
            SUM = SUM + GL(K) * YLMC(K,IQ) * YLM0(K)                    UPB 0052
20       CONTINUE                                                       UPB 0053
         ZJ(IQ) = ( 2. - DELM0 ) * FBEAM * SUM / (4.0*PI)               UPB 0054
40    CONTINUE                                                          UPB 0055
                                                                        UPB 0056
C                  ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION   UPB 0057
C                  ** OF -ARRAY- AND SEE IF IT IS NEARLY SINGULAR       UPB 0058
C                  ** (NOTE:  -ARRAY- IS DESTROYED)                     UPB 0059
      RCOND = 0.0                                                       UPB 0060
      CALL  DGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )                UPB 0061
      IF ( 1.0+RCOND .EQ. 1.0 )  CALL  ERRMSG                           UPB 0062
     $   ( 'UPBEAM--dGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)           UPB 0063
                                                                        UPB 0064
C                ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -ARRAY-       UPB 0065
C                ** (ASSUMED ALREADY L-U DECOMPOSED) AND R.H. SIDE(S)   UPB 0066
C                ** -ZJ-;  RETURN SOLUTION(S) IN -ZJ-                   UPB 0067
      JOB =   0                                                         UPB 0068
      CALL  DGESL( ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB )                  UPB 0069
      DO 50  IQ = 1, NN                                                 UPB 0070
         ZZ( IQ+NN )   = ZJ( IQ )                                       UPB 0071
         ZZ( NN+1-IQ ) = ZJ( IQ+NN )                                    UPB 0072
50    CONTINUE                                                          UPB 0073
      RETURN                                                            UPB 0074
      END                                                               UPB 0075
