      SUBROUTINE  SOLVE1( B, CBAND, FISOT, IHOM, IPVT, LL, MI9M2, MXCMU,SL1 0001
     $                    NCOL, NCUT, NN, NNLYRI, NSTR, Z )             SL1 0002
                                                                        SL1 0003
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            SL1 0004
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                SL1 0005
C        CONSTRUCT RIGHT-HAND SIDE VECTOR -B- FOR ISOTROPIC INCIDENCE   SL1 0006
C        (ONLY) ON EITHER TOP OR BOTTOM BOUNDARY AND SOLVE SYSTEM       SL1 0007
C        OF EQUATIONS OBTAINED FROM THE BOUNDARY CONDITIONS AND THE     SL1 0008
C        CONTINUITY-OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS           SL1 0009
                                                                        SL1 0010
C     ROUTINES CALLED:  DGBCO, DGBSL, ZEROIT                            SL1 0011
                                                                        SL1 0012
C     I N P U T      V A R I A B L E S:                                 SL1 0013
                                                                        SL1 0014
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),   SL1 0015
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED       SL1 0016
C                   BY LINPACK SOLUTION ROUTINES                        SL1 0017
C       IHOM     :  DIRECTION OF ILLUMINATION FLAG                      SL1 0018
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-                        SL1 0019
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)           SL1 0020
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        SL1 0021
                                                                        SL1 0022
C    O U T P U T     V A R I A B L E S:                                 SL1 0023
                                                                        SL1 0024
C       B        :  RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO      SL1 0025
C                   *DGBSL*; RETURNS AS SOLUTION VECTOR OF EQ.          SL1 0026
C                   SC(12), CONSTANTS OF INTEGRATION WITHOUT            SL1 0027
C                   EXPONENTIAL TERM                                    SL1 0028
C       LL      :   PERMANENT STORAGE FOR -B-, BUT RE-ORDERED           SL1 0029
                                                                        SL1 0030
C   I N T E R N A L    V A R I A B L E S:                               SL1 0031
                                                                        SL1 0032
C       IPVT     :  INTEGER VECTOR OF PIVOT INDICES                     SL1 0033
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL    SL1 0034
C       RCOND    :  INDICATOR OF SINGULARITY FOR -CBAND-                SL1 0035
C       Z        :  SCRATCH ARRAY REQUIRED BY *DGBCO*                   SL1 0036
C----------------------------------------------------------------------+SL1 0037
      INTEGER  IPVT(*)                                                  SL1 0038
      REAL*8     B( NNLYRI ), CBAND( MI9M2,NNLYRI ), LL( MXCMU,* ), Z(*)SL1 0039
                                                                        SL1 0040
                                                                        SL1 0041
      CALL  ZEROIT( B, NNLYRI )                                         SL1 0042
      NCD = 3*NN - 1                                                    SL1 0043
                                                                        SL1 0044
      IF ( IHOM.EQ.1 )  THEN                                            SL1 0045
C                             ** BECAUSE THERE ARE NO BEAM OR EMISSION  SL1 0046
C                             ** SOURCES, REMAINDER OF -B- ARRAY IS ZEROSL1 0047
         DO 10  I = 1, NN                                               SL1 0048
            B(I) = FISOT                                                SL1 0049
            B( NCOL-NN+I ) = 0.0                                        SL1 0050
10       CONTINUE                                                       SL1 0051
                                                                        SL1 0052
         RCOND = 0.0                                                    SL1 0053
         CALL  DGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )    SL1 0054
         IF ( 1.0+RCOND .EQ. 1.0 )  CALL  ERRMSG                        SL1 0055
     $         ( 'SOLVE1--dGBCO SAYS MATRIX NEAR SINGULAR', .FALSE. )   SL1 0056
                                                                        SL1 0057
      ELSE IF ( IHOM.EQ.2 )  THEN                                       SL1 0058
                                                                        SL1 0059
         DO 20 I = 1, NN                                                SL1 0060
            B(I) = 0.0                                                  SL1 0061
            B( NCOL-NN+I ) = FISOT                                      SL1 0062
20       CONTINUE                                                       SL1 0063
                                                                        SL1 0064
      END IF                                                            SL1 0065
                                                                        SL1 0066
      CALL  DGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )           SL1 0067
                                                                        SL1 0068
C                          ** ZERO -CBAND- TO GET RID OF 'FOREIGN'      SL1 0069
C                          ** ELEMENTS PUT IN BY LINPACK                SL1 0070
      DO 30  LC = 1, NCUT                                               SL1 0071
         IPNT = LC*NSTR - NN                                            SL1 0072
         DO 30  IQ = 1, NN                                              SL1 0073
            LL( NN+1-IQ, LC) = B( IPNT+1-IQ )                           SL1 0074
            LL( IQ+NN,   LC) = B( IQ+IPNT )                             SL1 0075
30    CONTINUE                                                          SL1 0076
                                                                        SL1 0077
      RETURN                                                            SL1 0078
      END                                                               SL1 0079
