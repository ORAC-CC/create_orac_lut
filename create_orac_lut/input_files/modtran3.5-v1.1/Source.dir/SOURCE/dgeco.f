      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)                            ECO 0001
      INTEGER LDA,N,IPVT(*)                                             ECO 0002
      DOUBLE PRECISION A(LDA,*),Z(*)                                    ECO 0003
      DOUBLE PRECISION RCOND                                            ECO 0004
C                                                                       ECO 0005
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION   ECO 0006
C     AND ESTIMATES THE CONDITION OF THE MATRIX.                        ECO 0007
C                                                                       ECO 0008
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.               ECO 0009
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.                        ECO 0010
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.                 ECO 0011
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.               ECO 0012
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.                   ECO 0013
C                                                                       ECO 0014
C     ON ENTRY                                                          ECO 0015
C                                                                       ECO 0016
C        A       DOUBLE PRECISION(LDA, N)                               ECO 0017
C                THE MATRIX TO BE FACTORED.                             ECO 0018
C                                                                       ECO 0019
C        LDA     INTEGER                                                ECO 0020
C                THE LEADING DIMENSION OF THE ARRAY  A .                ECO 0021
C                                                                       ECO 0022
C        N       INTEGER                                                ECO 0023
C                THE ORDER OF THE MATRIX  A .                           ECO 0024
C                                                                       ECO 0025
C     ON RETURN                                                         ECO 0026
C                                                                       ECO 0027
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         ECO 0028
C                WHICH WERE USED TO OBTAIN IT.                          ECO 0029
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       ECO 0030
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          ECO 0031
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       ECO 0032
C                                                                       ECO 0033
C        IPVT    INTEGER(N)                                             ECO 0034
C                AN INTEGER VECTOR OF PIVOT INDICES.                    ECO 0035
C                                                                       ECO 0036
C        RCOND   DOUBLE PRECISION                                       ECO 0037
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        ECO 0038
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       ECO 0039
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             ECO 0040
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . ECO 0041
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     ECO 0042
C                           1.0 + RCOND .EQ. 1.0                        ECO 0043
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           ECO 0044
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         ECO 0045
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          ECO 0046
C                UNDERFLOWS.                                            ECO 0047
C                                                                       ECO 0048
C        Z       DOUBLE PRECISION(N)                                    ECO 0049
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  ECO 0050
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      ECO 0051
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           ECO 0052
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    ECO 0053
C                                                                       ECO 0054
C     LINPACK. THIS VERSION DATED 08/14/78 .                            ECO 0055
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      ECO 0056
C                                                                       ECO 0057
C     SUBROUTINES AND FUNCTIONS                                         ECO 0058
C                                                                       ECO 0059
C     LINPACK DGEFA                                                     ECO 0060
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       ECO 0061
C     FORTRAN DABS,DMAX1,DSIGN                                          ECO 0062
C                                                                       ECO 0063
C     INTERNAL VARIABLES                                                ECO 0064
C                                                                       ECO 0065
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 ECO 0066
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           ECO 0067
      INTEGER INFO,J,K,KB,KP1,L                                         ECO 0068
C                                                                       ECO 0069
C                                                                       ECO 0070
C     COMPUTE 1-NORM OF A                                               ECO 0071
                                                                        ECO 0072
C                                                                       ECO 0073
      ANORM = 0.0D0                                                     ECO 0074
      DO 10 J = 1, N                                                    ECO 0075
                                                                        ECO 0076
           ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))                       ECO 0077
                                                                        ECO 0078
   10 CONTINUE                                                          ECO 0079
                                                                        ECO 0080
C                                                                       ECO 0081
C     FACTOR                                                            ECO 0082
C                                                                       ECO 0083
      CALL DGEFA(A,LDA,N,IPVT,INFO)                                     ECO 0084
                                                                        ECO 0085
                                                                        ECO 0086
                                                                        ECO 0087
                                                                        ECO 0088
                                                                        ECO 0089
C                                                                       ECO 0090
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              ECO 0091
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .  ECO 0092
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE      ECO 0093
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE  ECO 0094
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID    ECO 0095
C     OVERFLOW.                                                         ECO 0096
C                                                                       ECO 0097
C     SOLVE TRANS(U)*W = E                                              ECO 0098
C                                                                       ECO 0099
      EK = 1.0D0                                                        ECO 0100
      DO 20 J = 1, N                                                    ECO 0101
         Z(J) = 0.0D0                                                   ECO 0102
   20 CONTINUE                                                          ECO 0103
      DO 100 K = 1, N                                                   ECO 0104
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      ECO 0105
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30                  ECO 0106
            S = DABS(A(K,K))/DABS(EK-Z(K))                              ECO 0107
            CALL DSCAL(N,S,Z,1)                                         ECO 0108
            EK = S*EK                                                   ECO 0109
   30    CONTINUE                                                       ECO 0110
                                                                        ECO 0111
         WK = EK - Z(K)                                                 ECO 0112
         WKM = -EK - Z(K)                                               ECO 0113
         S = DABS(WK)                                                   ECO 0114
         SM = DABS(WKM)                                                 ECO 0115
         IF (A(K,K) .EQ. 0.0D0) GO TO 40                                ECO 0116
            WK = WK/A(K,K)                                              ECO 0117
            WKM = WKM/A(K,K)                                            ECO 0118
         GO TO 50                                                       ECO 0119
   40    CONTINUE                                                       ECO 0120
            WK = 1.0D0                                                  ECO 0121
            WKM = 1.0D0                                                 ECO 0122
   50    CONTINUE                                                       ECO 0123
         KP1 = K + 1                                                    ECO 0124
         IF (KP1 .GT. N) GO TO 90                                       ECO 0125
            DO 60 J = KP1, N                                            ECO 0126
               SM = SM + DABS(Z(J)+WKM*A(K,J))                          ECO 0127
               Z(J) = Z(J) + WK*A(K,J)                                  ECO 0128
               S = S + DABS(Z(J))                                       ECO 0129
   60       CONTINUE                                                    ECO 0130
            IF (S .GE. SM) GO TO 80                                     ECO 0131
               T = WKM - WK                                             ECO 0132
               WK = WKM                                                 ECO 0133
               DO 70 J = KP1, N                                         ECO 0134
                  Z(J) = Z(J) + T*A(K,J)                                ECO 0135
   70          CONTINUE                                                 ECO 0136
   80       CONTINUE                                                    ECO 0137
   90    CONTINUE                                                       ECO 0138
         Z(K) = WK                                                      ECO 0139
  100 CONTINUE                                                          ECO 0140
      S = 1.0D0/DASUM(N,Z,1)                                            ECO 0141
      CALL DSCAL(N,S,Z,1)                                               ECO 0142
C                                                                       ECO 0143
C     SOLVE TRANS(L)*Y = W                                              ECO 0144
C                                                                       ECO 0145
      DO 120 KB = 1, N                                                  ECO 0146
         K = N + 1 - KB                                                 ECO 0147
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)      ECO 0148
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110                           ECO 0149
            S = 1.0D0/DABS(Z(K))                                        ECO 0150
            CALL DSCAL(N,S,Z,1)                                         ECO 0151
  110    CONTINUE                                                       ECO 0152
         L = IPVT(K)                                                    ECO 0153
         T = Z(L)                                                       ECO 0154
         Z(L) = Z(K)                                                    ECO 0155
         Z(K) = T                                                       ECO 0156
  120 CONTINUE                                                          ECO 0157
      S = 1.0D0/DASUM(N,Z,1)                                            ECO 0158
      CALL DSCAL(N,S,Z,1)                                               ECO 0159
C                                                                       ECO 0160
      YNORM = 1.0D0                                                     ECO 0161
C                                                                       ECO 0162
C     SOLVE L*V = Y                                                     ECO 0163
C                                                                       ECO 0164
      DO 140 K = 1, N                                                   ECO 0165
         L = IPVT(K)                                                    ECO 0166
         T = Z(L)                                                       ECO 0167
         Z(L) = Z(K)                                                    ECO 0168
         Z(K) = T                                                       ECO 0169
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)            ECO 0170
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130                           ECO 0171
            S = 1.0D0/DABS(Z(K))                                        ECO 0172
            CALL DSCAL(N,S,Z,1)                                         ECO 0173
            YNORM = S*YNORM                                             ECO 0174
  130    CONTINUE                                                       ECO 0175
  140 CONTINUE                                                          ECO 0176
      S = 1.0D0/DASUM(N,Z,1)                                            ECO 0177
      CALL DSCAL(N,S,Z,1)                                               ECO 0178
      YNORM = S*YNORM                                                   ECO 0179
C                                                                       ECO 0180
C     SOLVE  U*Z = V                                                    ECO 0181
C                                                                       ECO 0182
      DO 160 KB = 1, N                                                  ECO 0183
         K = N + 1 - KB                                                 ECO 0184
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150                    ECO 0185
            S = DABS(A(K,K))/DABS(Z(K))                                 ECO 0186
            CALL DSCAL(N,S,Z,1)                                         ECO 0187
            YNORM = S*YNORM                                             ECO 0188
  150    CONTINUE                                                       ECO 0189
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                      ECO 0190
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                            ECO 0191
         T = -Z(K)                                                      ECO 0192
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                              ECO 0193
  160 CONTINUE                                                          ECO 0194
C     MAKE ZNORM = 1.0                                                  ECO 0195
      S = 1.0D0/DASUM(N,Z,1)                                            ECO 0196
      CALL DSCAL(N,S,Z,1)                                               ECO 0197
      YNORM = S*YNORM                                                   ECO 0198
C                                                                       ECO 0199
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         ECO 0200
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               ECO 0201
      RETURN                                                            ECO 0202
      END                                                               ECO 0203
