      SUBROUTINE  SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA, SL0 0001
     $                    FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,       SL0 0002
     $                    MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR,SL0 0003
     $                    NNLYRI, PI, TPLANK, TAUCPR, UMU0, Z, ZZ,      SL0 0004
     $                    ZPLK0, ZPLK1 )                                SL0 0005
                                                                        SL0 0006
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            SL0 0007
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                SL0 0008
C        CONSTRUCT RIGHT-HAND SIDE VECTOR -B- FOR GENERAL BOUNDARY      SL0 0009
C        CONDITIONS STWJ(17) AND SOLVE SYSTEM OF EQUATIONS OBTAINED     SL0 0010
C        FROM THE BOUNDARY CONDITIONS AND THE                           SL0 0011
C        CONTINUITY-OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS.          SL0 0012
C        THERMAL EMISSION CONTRIBUTES ONLY IN AZIMUTHAL INDEPENDENCE.   SL0 0013
                                                                        SL0 0014
C     ROUTINES CALLED:  DGBCO, DGBSL, ZEROIT                            SL0 0015
                                                                        SL0 0016
C     I N P U T      V A R I A B L E S:                                 SL0 0017
                                                                        SL0 0018
C       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY                  SL0 0019
C       BEM      :  SURFACE BIDIRECTIONAL EMISSIVITY                    SL0 0020
C       BPLANK   :  BOTTOM BOUNDARY THERMAL EMISSION                    SL0 0021
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),   SL0 0022
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED       SL0 0023
C                   BY LINPACK SOLUTION ROUTINES                        SL0 0024
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE    SL0 0025
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE      SL0 0026
C       EXPBEA   :  TRANSMISSION OF INCIDENT BEAM, EXP(-TAUCPR/UMU0)    SL0 0027
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER        SL0 0028
C       MAZIM    :  ORDER OF AZIMUTHAL COMPONENT                        SL0 0029
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-                        SL0 0030
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)           SL0 0031
C       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED     SL0 0032
C       TPLANK   :  TOP BOUNDARY THERMAL EMISSION                       SL0 0033
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)           SL0 0034
C       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)                   SL0 0035
C       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)  SL0 0036
C       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)  SL0 0037
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        SL0 0038
                                                                        SL0 0039
C   O U T P U T     V A R I A B L E S:                                  SL0 0040
                                                                        SL0 0041
C       B        :  RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO      SL0 0042
C                   *DGBSL*; RETURNS AS SOLUTION VECTOR OF EQ.          SL0 0043
C                   SC(12), CONSTANTS OF INTEGRATION WITHOUT            SL0 0044
C                   EXPONENTIAL TERM                                    SL0 0045
C      LL        :  PERMANENT STORAGE FOR -B-, BUT RE-ORDERED           SL0 0046
                                                                        SL0 0047
C   I N T E R N A L    V A R I A B L E S:                               SL0 0048
                                                                        SL0 0049
C       IPVT     :  INTEGER VECTOR OF PIVOT INDICES                     SL0 0050
C       IT       :  POINTER FOR POSITION IN  -B-                        SL0 0051
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL    SL0 0052
C       RCOND    :  INDICATOR OF SINGULARITY FOR -CBAND-                SL0 0053
C       Z        :  SCRATCH ARRAY REQUIRED BY *DGBCO*                   SL0 0054
C+---------------------------------------------------------------------+SL0 0055
                                                                        SL0 0056
      LOGICAL  LAMBER, LYRCUT                                           SL0 0057
      INTEGER  IPVT(*)                                                  SL0 0058
      REAL*8     B(*), BDR( MI,0:* ), BEM(*), CBAND( MI9M2,NNLYRI ),    SL0 0059
     $         CMU(*), CWT(*), EXPBEA(0:*), LL( MXCMU,* ),              SL0 0060
     $         TAUCPR( 0:* ), Z(*), ZZ( MXCMU,* ), ZPLK0( MXCMU,* ),    SL0 0061
     $         ZPLK1( MXCMU,* )                                         SL0 0062
                                                                        SL0 0063
      CALL  ZEROIT( B, NNLYRI )                                         SL0 0064
C                             ** CONSTRUCT -B-,  STWJ(20A,C) FOR        SL0 0065
C                             ** PARALLEL BEAM + BOTTOM REFLECTION +    SL0 0066
C                             ** THERMAL EMISSION AT TOP AND/OR BOTTOM  SL0 0067
      IF ( MAZIM.GT.0 .AND. FBEAM.GT.0.0 )  THEN                        SL0 0068
C                                         ** AZIMUTH-DEPENDENT CASE     SL0 0069
C                                         ** (NEVER CALLED IF FBEAM = 0)SL0 0070
         IF ( LYRCUT .OR. LAMBER ) THEN                                 SL0 0071
C               ** NO AZIMUTHAL-DEPENDENT INTENSITY FOR LAMBERT SURFACE;SL0 0072
C               ** NO INTENSITY COMPONENT FOR TRUNCATED BOTTOM LAYER    SL0 0073
                                                                        SL0 0074
            DO 10  IQ = 1, NN                                           SL0 0075
C                                                     ** TOP BOUNDARY   SL0 0076
               B(IQ) = - ZZ(NN+1-IQ,1)                                  SL0 0077
C                                                  ** BOTTOM BOUNDARY   SL0 0078
               B(NCOL-NN+IQ) = - ZZ(IQ+NN,NCUT) * EXPBEA(NCUT)          SL0 0079
10          CONTINUE                                                    SL0 0080
                                                                        SL0 0081
         ELSE                                                           SL0 0082
                                                                        SL0 0083
            DO 20  IQ = 1, NN                                           SL0 0084
               B(IQ) = - ZZ(NN+1-IQ,1)                                  SL0 0085
                                                                        SL0 0086
               SUM = 0.                                                 SL0 0087
               DO 15  JQ = 1, NN                                        SL0 0088
                  SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)            SL0 0089
     $                        * ZZ(NN+1-JQ,NCUT) * EXPBEA(NCUT)         SL0 0090
15             CONTINUE                                                 SL0 0091
               B(NCOL-NN+IQ) = SUM                                      SL0 0092
               IF ( FBEAM.GT.0.0 )                                      SL0 0093
     $              B(NCOL-NN+IQ) = SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI   SL0 0094
     $                                 - ZZ(IQ+NN,NCUT) ) * EXPBEA(NCUT)SL0 0095
20          CONTINUE                                                    SL0 0096
                                                                        SL0 0097
                                                                        SL0 0098
         END IF                                                         SL0 0099
                                                                        SL0 0100
                                                                        SL0 0101
C                             ** CONTINUITY CONDITION FOR LAYER         SL0 0102
C                             ** INTERFACES OF EQ. STWJ(20B)            SL0 0103
         IT = NN                                                        SL0 0104
         DO 40  LC = 1, NCUT-1                                          SL0 0105
            DO 30  IQ = 1, NSTR                                         SL0 0106
                                                                        SL0 0107
               IT    = IT + 1                                           SL0 0108
               B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)         SL0 0109
30          CONTINUE                                                    SL0 0110
40       CONTINUE                                                       SL0 0111
      ELSE                                                              SL0 0112
C                                   ** AZIMUTH-INDEPENDENT CASE         SL0 0113
         IF ( FBEAM.EQ.0.0 )  THEN                                      SL0 0114
                                                                        SL0 0115
            DO 50 IQ = 1, NN                                            SL0 0116
C                                      ** TOP BOUNDARY                  SL0 0117
               B(IQ) = - ZPLK0(NN+1-IQ,1) + FISOT + TPLANK              SL0 0118
50          CONTINUE                                                    SL0 0119
                                                                        SL0 0120
            IF ( LYRCUT ) THEN                                          SL0 0121
C                               ** NO INTENSITY COMPONENT FOR TRUNCATED SL0 0122
C                               ** BOTTOM LAYER                         SL0 0123
               DO 60 IQ = 1, NN                                         SL0 0124
C                                      ** BOTTOM BOUNDARY               SL0 0125
                  B(NCOL-NN+IQ) = - ZPLK0(IQ+NN,NCUT)                   SL0 0126
     $                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)    SL0 0127
60             CONTINUE                                                 SL0 0128
                                                                        SL0 0129
            ELSE                                                        SL0 0130
                                                                        SL0 0131
               DO 80 IQ = 1, NN                                         SL0 0132
                  SUM = 0.                                              SL0 0133
                 DO 70 JQ = 1, NN                                       SL0 0134
                                                                        SL0 0135
                                                                        SL0 0136
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)         SL0 0137
     $                          * ( ZPLK0(NN+1-JQ,NCUT)                 SL0 0138
     $                            + ZPLK1(NN+1-JQ,NCUT) * TAUCPR(NCUT) )SL0 0139
70                CONTINUE                                              SL0 0140
                  B(NCOL-NN+IQ) = 2.*SUM + BEM(IQ) * BPLANK             SL0 0141
     $                            - ZPLK0(IQ+NN,NCUT)                   SL0 0142
     $                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)    SL0 0143
80             CONTINUE                                                 SL0 0144
            END IF                                                      SL0 0145
C                             ** CONTINUITY CONDITION FOR LAYER         SL0 0146
C                             ** INTERFACES, STWJ(20B)                  SL0 0147
                                                                        SL0 0148
            IT = NN                                                     SL0 0149
            DO 100  LC = 1, NCUT-1                                      SL0 0150
               DO 90  IQ = 1, NSTR                                      SL0 0151
                  IT    = IT + 1                                        SL0 0152
                  B(IT) = ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +               SL0 0153
     $                  ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)  SL0 0154
90             CONTINUE                                                 SL0 0155
100         CONTINUE                                                    SL0 0156
                                                                        SL0 0157
         ELSE                                                           SL0 0158
            DO 150 IQ = 1, NN                                           SL0 0159
                                                                        SL0 0160
                                                                        SL0 0161
                                                                        SL0 0162
               B(IQ) = - ZZ(NN+1-IQ,1) - ZPLK0(NN+1-IQ,1) +FISOT +TPLANKSL0 0163
                                                                        SL0 0164
150         CONTINUE                                                    SL0 0165
                                                                        SL0 0166
            IF ( LYRCUT ) THEN                                          SL0 0167
               DO 160 IQ = 1, NN                                        SL0 0168
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN,NCUT) * EXPBEA(NCUT)       SL0 0169
     $                            - ZPLK0(IQ+NN,NCUT)                   SL0 0170
     $                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)    SL0 0171
160            CONTINUE                                                 SL0 0172
                                                                        SL0 0173
            ELSE                                                        SL0 0174
                                                                        SL0 0175
               DO 180 IQ = 1, NN                                        SL0 0176
                  SUM = 0.                                              SL0 0177
                  DO 170 JQ = 1, NN                                     SL0 0178
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)         SL0 0179
     $                          * ( ZZ(NN+1-JQ,NCUT) * EXPBEA(NCUT)     SL0 0180
     $                            + ZPLK0(NN+1-JQ,NCUT)                 SL0 0181
     $                            + ZPLK1(NN+1-JQ,NCUT) * TAUCPR(NCUT) )SL0 0182
170               CONTINUE                                              SL0 0183
                  B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI  SL0 0184
     $                                 - ZZ(IQ+NN,NCUT) ) * EXPBEA(NCUT)SL0 0185
     $                            + BEM(IQ) * BPLANK                    SL0 0186
     $                            - ZPLK0(IQ+NN,NCUT)                   SL0 0187
     $                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)    SL0 0188
180            CONTINUE                                                 SL0 0189
            END IF                                                      SL0 0190
            IT = NN                                                     SL0 0191
            DO 200  LC = 1, NCUT-1                                      SL0 0192
               DO 190  IQ = 1, NSTR                                     SL0 0193
                  IT    = IT + 1                                        SL0 0194
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)      SL0 0195
     $                    + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +             SL0 0196
     $                    ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)SL0 0197
190            CONTINUE                                                 SL0 0198
200         CONTINUE                                                    SL0 0199
         END IF                                                         SL0 0200
                                                                        SL0 0201
      END IF                                                            SL0 0202
C                     ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITIONSL0 0203
C                     ** OF BAND MATRIX -CBAND- AND TEST IF IT IS NEARLYSL0 0204
C                     ** SINGULAR (NOTE: -CBAND- IS DESTROYED)          SL0 0205
C                     ** (-CBAND- IS IN LINPACK PACKED FORMAT)          SL0 0206
      RCOND = 0.0                                                       SL0 0207
      NCD = 3*NN - 1                                                    SL0 0208
      CALL  DGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )       SL0 0209
      IF ( 1.0+RCOND .EQ. 1.0 )  CALL  ERRMSG                           SL0 0210
     $   ( 'SOLVE0--dGBCO SAYS MATRIX NEAR SINGULAR',.FALSE.)           SL0 0211
C                   ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -CBAND-    SL0 0212
C                   ** AND R.H. SIDE(S) -B- AFTER -CBAND- HAS BEEN L-U  SL0 0213
C                   ** DECOMPOSED.  SOLUTION IS RETURNED IN -B-.        SL0 0214
      CALL  DGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )           SL0 0215
C                   ** ZERO -CBAND- (IT MAY CONTAIN 'FOREIGN'           SL0 0216
C                   ** ELEMENTS UPON RETURNING FROM LINPACK);           SL0 0217
C                   ** NECESSARY TO PREVENT ERRORS                      SL0 0218
                                                                        SL0 0219
      CALL  ZEROIT( CBAND, MI9M2*NNLYRI )                               SL0 0220
      DO 220  LC = 1, NCUT                                              SL0 0221
         IPNT = LC*NSTR - NN                                            SL0 0222
         DO 220  IQ = 1, NN                                             SL0 0223
            LL(NN+1-IQ,LC) = B(IPNT+1-IQ)                               SL0 0224
            LL(IQ+NN,  LC) = B(IQ+IPNT)                                 SL0 0225
220   CONTINUE                                                          SL0 0226
      RETURN                                                            SL0 0227
      END                                                               SL0 0228
