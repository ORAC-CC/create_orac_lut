      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)                      BSL 0001
      INTEGER LDA,N,ML,MU,IPVT(*),JOB                                   BSL 0002
      DOUBLE PRECISION ABD(LDA,*),B(*)                                  BSL 0003
C                                                                       BSL 0004
C     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM                     BSL 0005
C     A * X = B  OR  TRANS(A) * X = B                                   BSL 0006
C     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.                     BSL 0007
C                                                                       BSL 0008
C     ON ENTRY                                                          BSL 0009
C                                                                       BSL 0010
C        ABD     DOUBLE PRECISION(LDA, N)                               BSL 0011
C                THE OUTPUT FROM DGBCO OR DGBFA.                        BSL 0012
C                                                                       BSL 0013
C        LDA     INTEGER                                                BSL 0014
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              BSL 0015
C                                                                       BSL 0016
C        N       INTEGER                                                BSL 0017
C                THE ORDER OF THE ORIGINAL MATRIX.                      BSL 0018
C                                                                       BSL 0019
C        ML      INTEGER                                                BSL 0020
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           BSL 0021
C                                                                       BSL 0022
C        MU      INTEGER                                                BSL 0023
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           BSL 0024
C                                                                       BSL 0025
C        IPVT    INTEGER(N)                                             BSL 0026
C                THE PIVOT VECTOR FROM DGBCO OR DGBFA.                  BSL 0027
C                                                                       BSL 0028
C        B       DOUBLE PRECISION(N)                                    BSL 0029
C                THE RIGHT HAND SIDE VECTOR.                            BSL 0030
C                                                                       BSL 0031
C        JOB     INTEGER                                                BSL 0032
C                = 0         TO SOLVE  A*X = B ,                        BSL 0033
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE           BSL 0034
C                            TRANS(A)  IS THE TRANSPOSE.                BSL 0035
C                                                                       BSL 0036
C     ON RETURN                                                         BSL 0037
C                                                                       BSL 0038
C        B       THE SOLUTION VECTOR  X .                               BSL 0039
C                                                                       BSL 0040
C     ERROR CONDITION                                                   BSL 0041
C                                                                       BSL 0042
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   BSL 0043
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  BSL 0044
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       BSL 0045
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     BSL 0046
C        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0           BSL 0047
C        OR DGBFA HAS SET INFO .EQ. 0 .                                 BSL 0048
C                                                                       BSL 0049
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 BSL 0050
C     WITH  P  COLUMNS                                                  BSL 0051
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    BSL 0052
C           IF (RCOND IS TOO SMALL) GO TO ...                           BSL 0053
C           DO 10 J = 1, P                                              BSL 0054
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)                BSL 0055
C        10 CONTINUE                                                    BSL 0056
C                                                                       BSL 0057
C     LINPACK. THIS VERSION DATED 08/14/78 .                            BSL 0058
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      BSL 0059
C                                                                       BSL 0060
C     SUBROUTINES AND FUNCTIONS                                         BSL 0061
C                                                                       BSL 0062
C     BLAS DAXPY,DDOT                                                   BSL 0063
C     FORTRAN MIN0                                                      BSL 0064
C                                                                       BSL 0065
C     INTERNAL VARIABLES                                                BSL 0066
C                                                                       BSL 0067
      DOUBLE PRECISION DDOT,T                                           BSL 0068
      INTEGER K,KB,L,LA,LB,LM,M,NM1                                     BSL 0069
C                                                                       BSL 0070
      M = MU + ML + 1                                                   BSL 0071
      NM1 = N - 1                                                       BSL 0072
      IF (JOB .NE. 0) GO TO 50                                          BSL 0073
C                                                                       BSL 0074
C        JOB = 0 , SOLVE  A * X = B                                     BSL 0075
C        FIRST SOLVE L*Y = B                                            BSL 0076
C                                                                       BSL 0077
         IF (ML .EQ. 0) GO TO 30                                        BSL 0078
         IF (NM1 .LT. 1) GO TO 30                                       BSL 0079
            DO 20 K = 1, NM1                                            BSL 0080
               LM = MIN0(ML,N-K)                                        BSL 0081
               L = IPVT(K)                                              BSL 0082
               T = B(L)                                                 BSL 0083
               IF (L .EQ. K) GO TO 10                                   BSL 0084
                  B(L) = B(K)                                           BSL 0085
                  B(K) = T                                              BSL 0086
   10          CONTINUE                                                 BSL 0087
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)                   BSL 0088
   20       CONTINUE                                                    BSL 0089
   30    CONTINUE                                                       BSL 0090
C                                                                       BSL 0091
C        NOW SOLVE  U*X = Y                                             BSL 0092
C                                                                       BSL 0093
         DO 40 KB = 1, N                                                BSL 0094
            K = N + 1 - KB                                              BSL 0095
            B(K) = B(K)/ABD(M,K)                                        BSL 0096
            LM = MIN0(K,M) - 1                                          BSL 0097
            LA = M - LM                                                 BSL 0098
            LB = K - LM                                                 BSL 0099
            T = -B(K)                                                   BSL 0100
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)                        BSL 0101
   40    CONTINUE                                                       BSL 0102
      GO TO 100                                                         BSL 0103
   50 CONTINUE                                                          BSL 0104
C                                                                       BSL 0105
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         BSL 0106
C        FIRST SOLVE  TRANS(U)*Y = B                                    BSL 0107
C                                                                       BSL 0108
         DO 60 K = 1, N                                                 BSL 0109
            LM = MIN0(K,M) - 1                                          BSL 0110
            LA = M - LM                                                 BSL 0111
            LB = K - LM                                                 BSL 0112
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)                            BSL 0113
            B(K) = (B(K) - T)/ABD(M,K)                                  BSL 0114
   60    CONTINUE                                                       BSL 0115
C                                                                       BSL 0116
C        NOW SOLVE TRANS(L)*X = Y                                       BSL 0117
C                                                                       BSL 0118
         IF (ML .EQ. 0) GO TO 90                                        BSL 0119
         IF (NM1 .LT. 1) GO TO 90                                       BSL 0120
            DO 80 KB = 1, NM1                                           BSL 0121
               K = N - KB                                               BSL 0122
               LM = MIN0(ML,N-K)                                        BSL 0123
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)             BSL 0124
               L = IPVT(K)                                              BSL 0125
               IF (L .EQ. K) GO TO 70                                   BSL 0126
                  T = B(L)                                              BSL 0127
                  B(L) = B(K)                                           BSL 0128
                  B(K) = T                                              BSL 0129
   70          CONTINUE                                                 BSL 0130
   80       CONTINUE                                                    BSL 0131
   90    CONTINUE                                                       BSL 0132
  100 CONTINUE                                                          BSL 0133
      RETURN                                                            BSL 0134
      END                                                               BSL 0135
