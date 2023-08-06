      SUBROUTINE DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    BCO 0001
      INTEGER LDA,N,ML,MU,IPVT(*)                                       BCO 0002
      DOUBLE PRECISION ABD(LDA,*),Z(*)                                  BCO 0003
      DOUBLE PRECISION RCOND                                            BCO 0004
C                                                                       BCO 0005
C     DGBCO FACTORS A DOUBLE PRECISION BAND MATRIX BY GAUSSIAN          BCO 0006
C     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.            BCO 0007
C                                                                       BCO 0008
C     IF  RCOND  IS NOT NEEDED, DGBFA IS SLIGHTLY FASTER.               BCO 0009
C     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL.                        BCO 0010
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL.                 BCO 0011
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI.               BCO 0012
C                                                                       BCO 0013
C     ON ENTRY                                                          BCO 0014
C                                                                       BCO 0015
C        ABD     DOUBLE PRECISION(LDA, N)                               BCO 0016
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      BCO 0017
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   BCO 0018
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         BCO 0019
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       BCO 0020
C                SEE THE COMMENTS BELOW FOR DETAILS.                    BCO 0021
C                                                                       BCO 0022
C        LDA     INTEGER                                                BCO 0023
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              BCO 0024
C                LDA MUST BE .GE. 2*ML + MU + 1 .                       BCO 0025
C                                                                       BCO 0026
C        N       INTEGER                                                BCO 0027
C                THE ORDER OF THE ORIGINAL MATRIX.                      BCO 0028
C                                                                       BCO 0029
C        ML      INTEGER                                                BCO 0030
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           BCO 0031
C                0 .LE. ML .LT. N .                                     BCO 0032
C                                                                       BCO 0033
C        MU      INTEGER                                                BCO 0034
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           BCO 0035
C                0 .LE. MU .LT. N .                                     BCO 0036
C                MORE EFFICIENT IF  ML .LE. MU .                        BCO 0037
C                                                                       BCO 0038
C     ON RETURN                                                         BCO 0039
C                                                                       BCO 0040
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         BCO 0041
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          BCO 0042
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       BCO 0043
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          BCO 0044
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       BCO 0045
C                                                                       BCO 0046
C        IPVT    INTEGER(N)                                             BCO 0047
C                AN INTEGER VECTOR OF PIVOT INDICES.                    BCO 0048
C                                                                       BCO 0049
C        RCOND   DOUBLE PRECISION                                       BCO 0050
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        BCO 0051
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       BCO 0052
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             BCO 0053
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . BCO 0054
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     BCO 0055
C                           1.0 + RCOND .EQ. 1.0                        BCO 0056
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           BCO 0057
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         BCO 0058
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          BCO 0059
C                UNDERFLOWS.                                            BCO 0060
C                                                                       BCO 0061
C        Z       DOUBLE PRECISION(N)                                    BCO 0062
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  BCO 0063
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      BCO 0064
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           BCO 0065
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    BCO 0066
C                                                                       BCO 0067
C     BAND STORAGE                                                      BCO 0068
C                                                                       BCO 0069
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      BCO 0070
C           WILL SET UP THE INPUT.                                      BCO 0071
C                                                                       BCO 0072
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)                BCO 0073
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                BCO 0074
C                   M = ML + MU + 1                                     BCO 0075
C                   DO 20 J = 1, N                                      BCO 0076
C                      I1 = MAX0(1, J-MU)                               BCO 0077
C                      I2 = MIN0(N, J+ML)                               BCO 0078
C                      DO 10 I = I1, I2                                 BCO 0079
C                         K = I - J + M                                 BCO 0080
C                         ABD(K,J) = A(I,J)                             BCO 0081
C                10    CONTINUE                                         BCO 0082
C                20 CONTINUE                                            BCO 0083
C                                                                       BCO 0084
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         BCO 0085
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      BCO 0086
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            BCO 0087
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    BCO 0088
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            BCO 0089
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          BCO 0090
C                                                                       BCO 0091
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS                              BCO 0092
C                                                                       BCO 0093
C           11 12 13  0  0  0                                           BCO 0094
C           21 22 23 24  0  0                                           BCO 0095
C            0 32 33 34 35  0                                           BCO 0096
C            0  0 43 44 45 46                                           BCO 0097
C            0  0  0 54 55 56                                           BCO 0098
C            0  0  0  0 65 66                                           BCO 0099
C                                                                       BCO 0100
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN  BCO 0101
C                                                                       BCO 0102
C            *  *  *  +  +  +  , * = NOT USED                           BCO 0103
C            *  * 13 24 35 46  , + = USED FOR PIVOTING                  BCO 0104
C            * 12 23 34 45 56                                           BCO 0105
C           11 22 33 44 55 66                                           BCO 0106
C           21 32 43 54 65  *                                           BCO 0107
C                                                                       BCO 0108
C     LINPACK. THIS VERSION DATED 08/14/78 .                            BCO 0109
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      BCO 0110
C                                                                       BCO 0111
C     SUBROUTINES AND FUNCTIONS                                         BCO 0112
C                                                                       BCO 0113
C     LINPACK DGBFA                                                     BCO 0114
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       BCO 0115
C     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN                                BCO 0116
C                                                                       BCO 0117
C     INTERNAL VARIABLES                                                BCO 0118
C                                                                       BCO 0119
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 BCO 0120
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           BCO 0121
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM                     BCO 0122
C                                                                       BCO 0123
C                                                                       BCO 0124
C     COMPUTE 1-NORM OF A                                               BCO 0125
C                                                                       BCO 0126
      ANORM = 0.0D0                                                     BCO 0127
      L = ML + 1                                                        BCO 0128
      IS = L + MU                                                       BCO 0129
      DO 10 J = 1, N                                                    BCO 0130
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1))                      BCO 0131
         IF (IS .GT. ML + 1) IS = IS - 1                                BCO 0132
         IF (J .LE. MU) L = L + 1                                       BCO 0133
         IF (J .GE. N - ML) L = L - 1                                   BCO 0134
   10 CONTINUE                                                          BCO 0135
C                                                                       BCO 0136
C     FACTOR                                                            BCO 0137
C                                                                       BCO 0138
      CALL DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                             BCO 0139
C                                                                       BCO 0140
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              BCO 0141
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .  BCO 0142
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE      BCO 0143
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE  BCO 0144
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID    BCO 0145
C     OVERFLOW.                                                         BCO 0146
C                                                                       BCO 0147
C     SOLVE TRANS(U)*W = E                                              BCO 0148
C                                                                       BCO 0149
      EK = 1.0D0                                                        BCO 0150
      DO 20 J = 1, N                                                    BCO 0151
         Z(J) = 0.0D0                                                   BCO 0152
   20 CONTINUE                                                          BCO 0153
      M = ML + MU + 1                                                   BCO 0154
      JU = 0                                                            BCO 0155
      DO 100 K = 1, N                                                   BCO 0156
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      BCO 0157
         IF (DABS(EK-Z(K)) .LE. DABS(ABD(M,K))) GO TO 30                BCO 0158
            S = DABS(ABD(M,K))/DABS(EK-Z(K))                            BCO 0159
            CALL DSCAL(N,S,Z,1)                                         BCO 0160
            EK = S*EK                                                   BCO 0161
   30    CONTINUE                                                       BCO 0162
         WK = EK - Z(K)                                                 BCO 0163
         WKM = -EK - Z(K)                                               BCO 0164
         S = DABS(WK)                                                   BCO 0165
         SM = DABS(WKM)                                                 BCO 0166
         IF (ABD(M,K) .EQ. 0.0D0) GO TO 40                              BCO 0167
            WK = WK/ABD(M,K)                                            BCO 0168
            WKM = WKM/ABD(M,K)                                          BCO 0169
         GO TO 50                                                       BCO 0170
   40    CONTINUE                                                       BCO 0171
            WK = 1.0D0                                                  BCO 0172
            WKM = 1.0D0                                                 BCO 0173
   50    CONTINUE                                                       BCO 0174
         KP1 = K + 1                                                    BCO 0175
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                               BCO 0176
         MM = M                                                         BCO 0177
         IF (KP1 .GT. JU) GO TO 90                                      BCO 0178
            DO 60 J = KP1, JU                                           BCO 0179
               MM = MM - 1                                              BCO 0180
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J))                       BCO 0181
               Z(J) = Z(J) + WK*ABD(MM,J)                               BCO 0182
               S = S + DABS(Z(J))                                       BCO 0183
   60       CONTINUE                                                    BCO 0184
            IF (S .GE. SM) GO TO 80                                     BCO 0185
               T = WKM - WK                                             BCO 0186
               WK = WKM                                                 BCO 0187
               MM = M                                                   BCO 0188
               DO 70 J = KP1, JU                                        BCO 0189
                  MM = MM - 1                                           BCO 0190
                  Z(J) = Z(J) + T*ABD(MM,J)                             BCO 0191
   70          CONTINUE                                                 BCO 0192
   80       CONTINUE                                                    BCO 0193
   90    CONTINUE                                                       BCO 0194
         Z(K) = WK                                                      BCO 0195
  100 CONTINUE                                                          BCO 0196
      S = 1.0D0/DASUM(N,Z,1)                                            BCO 0197
      CALL DSCAL(N,S,Z,1)                                               BCO 0198
C                                                                       BCO 0199
C     SOLVE TRANS(L)*Y = W                                              BCO 0200
C                                                                       BCO 0201
      DO 120 KB = 1, N                                                  BCO 0202
         K = N + 1 - KB                                                 BCO 0203
         LM = MIN0(ML,N-K)                                              BCO 0204
         IF (K .LT. N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1)     BCO 0205
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110                           BCO 0206
            S = 1.0D0/DABS(Z(K))                                        BCO 0207
            CALL DSCAL(N,S,Z,1)                                         BCO 0208
  110    CONTINUE                                                       BCO 0209
         L = IPVT(K)                                                    BCO 0210
         T = Z(L)                                                       BCO 0211
         Z(L) = Z(K)                                                    BCO 0212
         Z(K) = T                                                       BCO 0213
  120 CONTINUE                                                          BCO 0214
      S = 1.0D0/DASUM(N,Z,1)                                            BCO 0215
      CALL DSCAL(N,S,Z,1)                                               BCO 0216
C                                                                       BCO 0217
      YNORM = 1.0D0                                                     BCO 0218
C                                                                       BCO 0219
C     SOLVE L*V = Y                                                     BCO 0220
C                                                                       BCO 0221
      DO 140 K = 1, N                                                   BCO 0222
         L = IPVT(K)                                                    BCO 0223
         T = Z(L)                                                       BCO 0224
         Z(L) = Z(K)                                                    BCO 0225
         Z(K) = T                                                       BCO 0226
         LM = MIN0(ML,N-K)                                              BCO 0227
         IF (K .LT. N) CALL DAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)           BCO 0228
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130                           BCO 0229
            S = 1.0D0/DABS(Z(K))                                        BCO 0230
            CALL DSCAL(N,S,Z,1)                                         BCO 0231
            YNORM = S*YNORM                                             BCO 0232
  130    CONTINUE                                                       BCO 0233
  140 CONTINUE                                                          BCO 0234
      S = 1.0D0/DASUM(N,Z,1)                                            BCO 0235
      CALL DSCAL(N,S,Z,1)                                               BCO 0236
      YNORM = S*YNORM                                                   BCO 0237
C                                                                       BCO 0238
C     SOLVE  U*Z = W                                                    BCO 0239
C                                                                       BCO 0240
      DO 160 KB = 1, N                                                  BCO 0241
         K = N + 1 - KB                                                 BCO 0242
         IF (DABS(Z(K)) .LE. DABS(ABD(M,K))) GO TO 150                  BCO 0243
            S = DABS(ABD(M,K))/DABS(Z(K))                               BCO 0244
            CALL DSCAL(N,S,Z,1)                                         BCO 0245
            YNORM = S*YNORM                                             BCO 0246
  150    CONTINUE                                                       BCO 0247
         IF (ABD(M,K) .NE. 0.0D0) Z(K) = Z(K)/ABD(M,K)                  BCO 0248
         IF (ABD(M,K) .EQ. 0.0D0) Z(K) = 1.0D0                          BCO 0249
         LM = MIN0(K,M) - 1                                             BCO 0250
         LA = M - LM                                                    BCO 0251
         LZ = K - LM                                                    BCO 0252
         T = -Z(K)                                                      BCO 0253
         CALL DAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)                           BCO 0254
  160 CONTINUE                                                          BCO 0255
C     MAKE ZNORM = 1.0                                                  BCO 0256
      S = 1.0D0/DASUM(N,Z,1)                                            BCO 0257
      CALL DSCAL(N,S,Z,1)                                               BCO 0258
      YNORM = S*YNORM                                                   BCO 0259
C                                                                       BCO 0260
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         BCO 0261
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               BCO 0262
      RETURN                                                            BCO 0263
      END                                                               BCO 0264
