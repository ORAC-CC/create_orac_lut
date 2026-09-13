      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)                              ESL 0001
      INTEGER LDA,N,IPVT(*),JOB                                         ESL 0002
      DOUBLE PRECISION A(LDA,*),B(*)                                    ESL 0003
C                                                                       ESL 0004
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM                          ESL 0005
C     A * X = B  OR  TRANS(A) * X = B                                   ESL 0006
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     ESL 0007
C                                                                       ESL 0008
C     ON ENTRY                                                          ESL 0009
C                                                                       ESL 0010
C        A       DOUBLE PRECISION(LDA, N)                               ESL 0011
C                THE OUTPUT FROM DGECO OR DGEFA.                        ESL 0012
C                                                                       ESL 0013
C        LDA     INTEGER                                                ESL 0014
C                THE LEADING DIMENSION OF THE ARRAY  A .                ESL 0015
C                                                                       ESL 0016
C        N       INTEGER                                                ESL 0017
C                THE ORDER OF THE MATRIX  A .                           ESL 0018
C                                                                       ESL 0019
C        IPVT    INTEGER(N)                                             ESL 0020
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  ESL 0021
C                                                                       ESL 0022
C        B       DOUBLE PRECISION(N)                                    ESL 0023
C                THE RIGHT HAND SIDE VECTOR.                            ESL 0024
C                                                                       ESL 0025
C        JOB     INTEGER                                                ESL 0026
C                = 0         TO SOLVE  A*X = B ,                        ESL 0027
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE            ESL 0028
C                            TRANS(A)  IS THE TRANSPOSE.                ESL 0029
C                                                                       ESL 0030
C     ON RETURN                                                         ESL 0031
C                                                                       ESL 0032
C        B       THE SOLUTION VECTOR  X .                               ESL 0033
C                                                                       ESL 0034
C     ERROR CONDITION                                                   ESL 0035
C                                                                       ESL 0036
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   ESL 0037
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  ESL 0038
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       ESL 0039
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     ESL 0040
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0           ESL 0041
C        OR DGEFA HAS SET INFO .EQ. 0 .                                 ESL 0042
C                                                                       ESL 0043
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 ESL 0044
C     WITH  P  COLUMNS                                                  ESL 0045
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            ESL 0046
C           IF (RCOND IS TOO SMALL) GO TO ...                           ESL 0047
C           DO 10 J = 1, P                                              ESL 0048
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        ESL 0049
C        10 CONTINUE                                                    ESL 0050
C                                                                       ESL 0051
C     LINPACK. THIS VERSION DATED 08/14/78 .                            ESL 0052
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      ESL 0053
C                                                                       ESL 0054
C     SUBROUTINES AND FUNCTIONS                                         ESL 0055
C                                                                       ESL 0056
C     BLAS DAXPY,DDOT                                                   ESL 0057
C                                                                       ESL 0058
C     INTERNAL VARIABLES                                                ESL 0059
C                                                                       ESL 0060
      DOUBLE PRECISION DDOT,T                                           ESL 0061
      INTEGER K,KB,L,NM1                                                ESL 0062
C                                                                       ESL 0063
      NM1 = N - 1                                                       ESL 0064
      IF (JOB .NE. 0) GO TO 50                                          ESL 0065
C                                                                       ESL 0066
C        JOB = 0 , SOLVE  A * X = B                                     ESL 0067
C        FIRST SOLVE  L*Y = B                                           ESL 0068
C                                                                       ESL 0069
         IF (NM1 .LT. 1) GO TO 30                                       ESL 0070
         DO 20 K = 1, NM1                                               ESL 0071
            L = IPVT(K)                                                 ESL 0072
            T = B(L)                                                    ESL 0073
            IF (L .EQ. K) GO TO 10                                      ESL 0074
               B(L) = B(K)                                              ESL 0075
               B(K) = T                                                 ESL 0076
   10       CONTINUE                                                    ESL 0077
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)                       ESL 0078
   20    CONTINUE                                                       ESL 0079
   30    CONTINUE                                                       ESL 0080
C                                                                       ESL 0081
C        NOW SOLVE  U*X = Y                                             ESL 0082
C                                                                       ESL 0083
         DO 40 KB = 1, N                                                ESL 0084
            K = N + 1 - KB                                              ESL 0085
            B(K) = B(K)/A(K,K)                                          ESL 0086
            T = -B(K)                                                   ESL 0087
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                           ESL 0088
   40    CONTINUE                                                       ESL 0089
      GO TO 100                                                         ESL 0090
   50 CONTINUE                                                          ESL 0091
C                                                                       ESL 0092
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         ESL 0093
C        FIRST SOLVE  TRANS(U)*Y = B                                    ESL 0094
C                                                                       ESL 0095
         DO 60 K = 1, N                                                 ESL 0096
            T = DDOT(K-1,A(1,K),1,B(1),1)                               ESL 0097
            B(K) = (B(K) - T)/A(K,K)                                    ESL 0098
   60    CONTINUE                                                       ESL 0099
C                                                                       ESL 0100
C        NOW SOLVE TRANS(L)*X = Y                                       ESL 0101
C                                                                       ESL 0102
         IF (NM1 .LT. 1) GO TO 90                                       ESL 0103
         DO 80 KB = 1, NM1                                              ESL 0104
            K = N - KB                                                  ESL 0105
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)                 ESL 0106
            L = IPVT(K)                                                 ESL 0107
            IF (L .EQ. K) GO TO 70                                      ESL 0108
               T = B(L)                                                 ESL 0109
               B(L) = B(K)                                              ESL 0110
               B(K) = T                                                 ESL 0111
   70       CONTINUE                                                    ESL 0112
   80    CONTINUE                                                       ESL 0113
   90    CONTINUE                                                       ESL 0114
  100 CONTINUE                                                          ESL 0115
      RETURN                                                            ESL 0116
      END                                                               ESL 0117
