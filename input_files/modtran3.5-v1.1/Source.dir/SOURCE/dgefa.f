      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)                               EFA 0001
      INTEGER LDA,N,IPVT(*),INFO                                        EFA 0002
      DOUBLE PRECISION A(LDA,*)                                         EFA 0003
C                                                                       EFA 0004
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.  EFA 0005
C                                                                       EFA 0006
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED            EFA 0007
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          EFA 0008
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .                   EFA 0009
C                                                                       EFA 0010
C     ON ENTRY                                                          EFA 0011
C                                                                       EFA 0012
C        A       DOUBLE PRECISION(LDA, N)                               EFA 0013
C                THE MATRIX TO BE FACTORED.                             EFA 0014
C                                                                       EFA 0015
C        LDA     INTEGER                                                EFA 0016
C                THE LEADING DIMENSION OF THE ARRAY  A .                EFA 0017
C                                                                       EFA 0018
C        N       INTEGER                                                EFA 0019
C                THE ORDER OF THE MATRIX  A .                           EFA 0020
C                                                                       EFA 0021
C     ON RETURN                                                         EFA 0022
C                                                                       EFA 0023
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         EFA 0024
C                WHICH WERE USED TO OBTAIN IT.                          EFA 0025
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       EFA 0026
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          EFA 0027
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       EFA 0028
C                                                                       EFA 0029
C        IPVT    INTEGER(N)                                             EFA 0030
C                AN INTEGER VECTOR OF PIVOT INDICES.                    EFA 0031
C                                                                       EFA 0032
C        INFO    INTEGER                                                EFA 0033
C                = 0  NORMAL VALUE.                                     EFA 0034
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       EFA 0035
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        EFA 0036
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO  EFA 0037
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE   EFA 0038
C                     INDICATION OF SINGULARITY.                        EFA 0039
C                                                                       EFA 0040
C     LINPACK. THIS VERSION DATED 08/14/78 .                            EFA 0041
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      EFA 0042
C                                                                       EFA 0043
C     SUBROUTINES AND FUNCTIONS                                         EFA 0044
C                                                                       EFA 0045
C     BLAS DAXPY,DSCAL,IDAMAX                                           EFA 0046
C                                                                       EFA 0047
C     INTERNAL VARIABLES                                                EFA 0048
C                                                                       EFA 0049
      DOUBLE PRECISION T                                                EFA 0050
      INTEGER IDAMAX,J,K,KP1,L,NM1                                      EFA 0051
C                                                                       EFA 0052
C                                                                       EFA 0053
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        EFA 0054
C                                                                       EFA 0055
      INFO = 0                                                          EFA 0056
      NM1 = N - 1                                                       EFA 0057
      IF (NM1 .LT. 1) GO TO 70                                          EFA 0058
      DO 60 K = 1, NM1                                                  EFA 0059
         KP1 = K + 1                                                    EFA 0060
C                                                                       EFA 0061
C        FIND L = PIVOT INDEX                                           EFA 0062
C                                                                       EFA 0063
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1                             EFA 0064
         IPVT(K) = L                                                    EFA 0065
C                                                                       EFA 0066
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          EFA 0067
C                                                                       EFA 0068
         IF (A(L,K) .EQ. 0.0D0) GO TO 40                                EFA 0069
C                                                                       EFA 0070
C           INTERCHANGE IF NECESSARY                                    EFA 0071
C                                                                       EFA 0072
            IF (L .EQ. K) GO TO 10                                      EFA 0073
               T = A(L,K)                                               EFA 0074
               A(L,K) = A(K,K)                                          EFA 0075
               A(K,K) = T                                               EFA 0076
   10       CONTINUE                                                    EFA 0077
C                                                                       EFA 0078
C           COMPUTE MULTIPLIERS                                         EFA 0079
C                                                                       EFA 0080
            T = -1.0D0/A(K,K)                                           EFA 0081
            CALL DSCAL(N-K,T,A(K+1,K),1)                                EFA 0082
C                                                                       EFA 0083
C           ROW ELIMINATION WITH COLUMN INDEXING                        EFA 0084
C                                                                       EFA 0085
            DO 30 J = KP1, N                                            EFA 0086
               T = A(L,J)                                               EFA 0087
               IF (L .EQ. K) GO TO 20                                   EFA 0088
                  A(L,J) = A(K,J)                                       EFA 0089
                  A(K,J) = T                                            EFA 0090
   20          CONTINUE                                                 EFA 0091
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)                  EFA 0092
   30       CONTINUE                                                    EFA 0093
         GO TO 50                                                       EFA 0094
   40    CONTINUE                                                       EFA 0095
            INFO = K                                                    EFA 0096
   50    CONTINUE                                                       EFA 0097
   60 CONTINUE                                                          EFA 0098
   70 CONTINUE                                                          EFA 0099
      IPVT(N) = N                                                       EFA 0100
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                   EFA 0101
      RETURN                                                            EFA 0102
      END                                                               EFA 0103
