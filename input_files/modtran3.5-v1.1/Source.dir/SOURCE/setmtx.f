      SUBROUTINE  SETMTX( BDR, CBAND, CMU, CWT, DELM0, GC, KK, LAMBER,  MTX 0001
     $                    LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI, MTX 0002
     $                    NN, NSTR, TAUCPR, WK )                        MTX 0003
                                                                        MTX 0004
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            MTX 0005
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                MTX 0006
C        CALCULATE COEFFICIENT MATRIX FOR THE SET OF EQUATIONS          MTX 0007
C        OBTAINED FROM THE BOUNDARY CONDITIONS AND THE CONTINUITY-      MTX 0008
C        OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS;  STORE IN THE       MTX 0009
C        SPECIAL BANDED-MATRIX FORMAT REQUIRED BY LINPACK ROUTINES      MTX 0010
                                                                        MTX 0011
C     ROUTINES CALLED:  ZEROIT                                          MTX 0012
                                                                        MTX 0013
C     I N P U T      V A R I A B L E S:                                 MTX 0014
                                                                        MTX 0015
C       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY                  MTX 0016
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE    MTX 0017
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE      MTX 0018
C       DELM0    :  KRONECKER DELTA, DELTA-SUB-M0                       MTX 0019
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)      MTX 0020
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)           MTX 0021
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER        MTX 0022
C       NN       :  NUMBER OF STREAMS IN A HEMISPHERE (NSTR/2)          MTX 0023
C       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED     MTX 0024
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)           MTX 0025
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        MTX 0026
                                                                        MTX 0027
C   O U T P U T     V A R I A B L E S:                                  MTX 0028
                                                                        MTX 0029
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),   MTX 0030
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED       MTX 0031
C                   BY LINPACK SOLUTION ROUTINES                        MTX 0032
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-                        MTX 0033
                                                                        MTX 0034
C   I N T E R N A L    V A R I A B L E S:                               MTX 0035
                                                                        MTX 0036
C       IROW     :  POINTS TO ROW IN  -CBAND-                           MTX 0037
C       JCOL     :  POINTS TO POSITION IN LAYER BLOCK                   MTX 0038
C       LDA      :  ROW DIMENSION OF -CBAND-                            MTX 0039
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL    MTX 0040
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-                        MTX 0041
C       NSHIFT   :  FOR POSITIONING NUMBER OF ROWS IN BAND STORAGE      MTX 0042
C       WK       :  TEMPORARY STORAGE FOR 'EXP' EVALUATIONS             MTX 0043
C ---------------------------------------------------------------------+MTX 0044
      LOGICAL LAMBER, LYRCUT                                            MTX 0045
      REAL*8    BDR( MI,0:* ), CBAND( MI9M2,NNLYRI ), CMU(*), CWT(*),   MTX 0046
     $        GC( MXCMU,MXCMU,* ), KK( MXCMU,* ), TAUCPR( 0:* ), WK(*)  MTX 0047
                                                                        MTX 0048
                                                                        MTX 0049
      CALL  ZEROIT( CBAND, MI9M2*NNLYRI )                               MTX 0050
      NCD    = 3*NN - 1                                                 MTX 0051
      LDA    = 3*NCD + 1                                                MTX 0052
      NSHIFT = LDA - 2*NSTR + 1                                         MTX 0053
      NCOL   = 0                                                        MTX 0054
C                         ** USE CONTINUITY CONDITIONS OF EQ. STWJ(17)  MTX 0055
C                         ** TO FORM COEFFICIENT MATRIX IN STWJ(20);    MTX 0056
C                         ** EMPLOY SCALING TRANSFORMATION STWJ(22)     MTX 0057
      DO 30  LC = 1, NCUT                                               MTX 0058
                                                                        MTX 0059
         DO 4  IQ = 1, NN                                               MTX 0060
            WK(IQ) = DEXP( KK(IQ,LC) * (TAUCPR(LC) - TAUCPR(LC-1)) )    MTX 0061
 4       CONTINUE                                                       MTX 0062
                                                                        MTX 0063
         JCOL = 0                                                       MTX 0064
         DO 10  IQ = 1, NN                                              MTX 0065
            NCOL = NCOL + 1                                             MTX 0066
            IROW = NSHIFT - JCOL                                        MTX 0067
            DO 5  JQ = 1, NSTR                                          MTX 0068
               CBAND(IROW+NSTR,NCOL) =   GC(JQ,IQ,LC)                   MTX 0069
               CBAND(IROW,     NCOL) = - GC(JQ,IQ,LC) * WK(IQ)          MTX 0070
               IROW = IROW + 1                                          MTX 0071
 5          CONTINUE                                                    MTX 0072
            JCOL = JCOL + 1                                             MTX 0073
10       CONTINUE                                                       MTX 0074
                                                                        MTX 0075
         DO 20  IQ = NN+1, NSTR                                         MTX 0076
            NCOL = NCOL + 1                                             MTX 0077
            IROW = NSHIFT - JCOL                                        MTX 0078
            DO 15  JQ = 1, NSTR                                         MTX 0079
               CBAND(IROW+NSTR,NCOL) =   GC(JQ,IQ,LC) * WK(NSTR+1-IQ)   MTX 0080
               CBAND(IROW,     NCOL) = - GC(JQ,IQ,LC)                   MTX 0081
               IROW = IROW + 1                                          MTX 0082
15          CONTINUE                                                    MTX 0083
            JCOL = JCOL + 1                                             MTX 0084
20       CONTINUE                                                       MTX 0085
                                                                        MTX 0086
30    CONTINUE                                                          MTX 0087
C                  ** USE TOP BOUNDARY CONDITION OF STWJ(20A) FOR       MTX 0088
C                  ** FIRST LAYER                                       MTX 0089
      JCOL = 0                                                          MTX 0090
      DO 40  IQ = 1, NN                                                 MTX 0091
         EXPA = DEXP( KK(IQ,1) * TAUCPR(1) )                            MTX 0092
         IROW = NSHIFT - JCOL + NN                                      MTX 0093
         DO 35  JQ = NN, 1, -1                                          MTX 0094
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1) * EXPA                     MTX 0095
            IROW = IROW+1                                               MTX 0096
35       CONTINUE                                                       MTX 0097
         JCOL = JCOL+1                                                  MTX 0098
40    CONTINUE                                                          MTX 0099
                                                                        MTX 0100
      DO 50  IQ = NN+1, NSTR                                            MTX 0101
         IROW = NSHIFT - JCOL + NN                                      MTX 0102
         DO 45  JQ = NN, 1, -1                                          MTX 0103
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1)                            MTX 0104
            IROW = IROW+1                                               MTX 0105
45       CONTINUE                                                       MTX 0106
         JCOL = JCOL+1                                                  MTX 0107
50    CONTINUE                                                          MTX 0108
C                           ** USE BOTTOM BOUNDARY CONDITION OF         MTX 0109
C                           ** STWJ(20C) FOR LAST LAYER                 MTX 0110
      NNCOL = NCOL - NSTR                                               MTX 0111
      JCOL  = 0                                                         MTX 0112
      DO 70  IQ = 1, NN                                                 MTX 0113
         NNCOL = NNCOL + 1                                              MTX 0114
         IROW  = NSHIFT - JCOL + NSTR                                   MTX 0115
                                                                        MTX 0116
         DO 60  JQ = NN+1, NSTR                                         MTX 0117
            IF ( LYRCUT .OR. (LAMBER .AND. DELM0.EQ.0) ) THEN           MTX 0118
                                                                        MTX 0119
C                          ** NO AZIMUTHAL-DEPENDENT INTENSITY IF LAM-  MTX 0120
C                          ** BERT SURFACE; NO INTENSITY COMPONENT IF   MTX 0121
C                          ** TRUNCATED BOTTOM LAYER                    MTX 0122
                                                                        MTX 0123
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT)                       MTX 0124
            ELSE                                                        MTX 0125
               SUM = 0.0                                                MTX 0126
               DO 55  K = 1, NN                                         MTX 0127
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)            MTX 0128
     $                        * GC(NN+1-K,IQ,NCUT)                      MTX 0129
55             CONTINUE                                                 MTX 0130
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) - (1.+DELM0) * SUM    MTX 0131
            END IF                                                      MTX 0132
                                                                        MTX 0133
            IROW = IROW + 1                                             MTX 0134
60       CONTINUE                                                       MTX 0135
         JCOL = JCOL + 1                                                MTX 0136
70    CONTINUE                                                          MTX 0137
                                                                        MTX 0138
      DO 90  IQ = NN+1, NSTR                                            MTX 0139
         NNCOL = NNCOL + 1                                              MTX 0140
         IROW  = NSHIFT - JCOL + NSTR                                   MTX 0141
         EXPA = WK(NSTR+1-IQ)                                           MTX 0142
                                                                        MTX 0143
         DO 80  JQ = NN+1, NSTR                                         MTX 0144
                                                                        MTX 0145
            IF ( LYRCUT .OR. (LAMBER .AND. DELM0.EQ.0) ) THEN           MTX 0146
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) * EXPA                MTX 0147
            ELSE                                                        MTX 0148
               SUM = 0.0                                                MTX 0149
               DO 75  K = 1, NN                                         MTX 0150
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)            MTX 0151
     $                        * GC(NN+1-K,IQ,NCUT)                      MTX 0152
75             CONTINUE                                                 MTX 0153
               CBAND(IROW,NNCOL) = ( GC(JQ,IQ,NCUT)                     MTX 0154
     $                               - (1.+DELM0) * SUM ) * EXPA        MTX 0155
            END IF                                                      MTX 0156
                                                                        MTX 0157
            IROW = IROW + 1                                             MTX 0158
80       CONTINUE                                                       MTX 0159
         JCOL = JCOL + 1                                                MTX 0160
90    CONTINUE                                                          MTX 0161
                                                                        MTX 0162
      RETURN                                                            MTX 0163
      END                                                               MTX 0164
