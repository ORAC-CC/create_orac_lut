      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                       BFA 0001
      INTEGER LDA,N,ML,MU,IPVT(*),INFO                                  BFA 0002
      DOUBLE PRECISION ABD(LDA,*)                                       BFA 0003
C                                                                       BFA 0004
C     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.      BFA 0005
C                                                                       BFA 0006
C     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED            BFA 0007
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          BFA 0008
C                                                                       BFA 0009
C     ON ENTRY                                                          BFA 0010
C                                                                       BFA 0011
C        ABD     DOUBLE PRECISION(LDA, N)                               BFA 0012
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      BFA 0013
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   BFA 0014
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         BFA 0015
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       BFA 0016
C                SEE THE COMMENTS BELOW FOR DETAILS.                    BFA 0017
C                                                                       BFA 0018
C        LDA     INTEGER                                                BFA 0019
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              BFA 0020
C                LDA MUST BE .GE. 2*ML + MU + 1 .                       BFA 0021
C                                                                       BFA 0022
C        N       INTEGER                                                BFA 0023
C                THE ORDER OF THE ORIGINAL MATRIX.                      BFA 0024
C                                                                       BFA 0025
C        ML      INTEGER                                                BFA 0026
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           BFA 0027
C                0 .LE. ML .LT. N .                                     BFA 0028
C                                                                       BFA 0029
C        MU      INTEGER                                                BFA 0030
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           BFA 0031
C                0 .LE. MU .LT. N .                                     BFA 0032
C                MORE EFFICIENT IF  ML .LE. MU .                        BFA 0033
C     ON RETURN                                                         BFA 0034
C                                                                       BFA 0035
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         BFA 0036
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          BFA 0037
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       BFA 0038
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          BFA 0039
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       BFA 0040
C                                                                       BFA 0041
C        IPVT    INTEGER(N)                                             BFA 0042
C                AN INTEGER VECTOR OF PIVOT INDICES.                    BFA 0043
C                                                                       BFA 0044
C        INFO    INTEGER                                                BFA 0045
C                = 0  NORMAL VALUE.                                     BFA 0046
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       BFA 0047
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        BFA 0048
C                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF        BFA 0049
C                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE      BFA 0050
C                     INDICATION OF SINGULARITY.                        BFA 0051
C                                                                       BFA 0052
C     BAND STORAGE                                                      BFA 0053
C                                                                       BFA 0054
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      BFA 0055
C           WILL SET UP THE INPUT.                                      BFA 0056
C                                                                       BFA 0057
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)                BFA 0058
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                BFA 0059
C                   M = ML + MU + 1                                     BFA 0060
C                   DO 20 J = 1, N                                      BFA 0061
C                      I1 = MAX0(1, J-MU)                               BFA 0062
C                      I2 = MIN0(N, J+ML)                               BFA 0063
C                      DO 10 I = I1, I2                                 BFA 0064
C                         K = I - J + M                                 BFA 0065
C                         ABD(K,J) = A(I,J)                             BFA 0066
C                10    CONTINUE                                         BFA 0067
C                20 CONTINUE                                            BFA 0068
C                                                                       BFA 0069
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         BFA 0070
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      BFA 0071
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            BFA 0072
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    BFA 0073
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            BFA 0074
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          BFA 0075
C                                                                       BFA 0076
C     LINPACK. THIS VERSION DATED 08/14/78 .                            BFA 0077
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      BFA 0078
C                                                                       BFA 0079
C     SUBROUTINES AND FUNCTIONS                                         BFA 0080
C                                                                       BFA 0081
C     BLAS DAXPY,DSCAL,IDAMAX                                           BFA 0082
C     FORTRAN MAX0,MIN0                                                 BFA 0083
C                                                                       BFA 0084
C     INTERNAL VARIABLES                                                BFA 0085
C                                                                       BFA 0086
      DOUBLE PRECISION T                                                BFA 0087
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1             BFA 0088
C                                                                       BFA 0089
C                                                                       BFA 0090
      M = ML + MU + 1                                                   BFA 0091
      INFO = 0                                                          BFA 0092
C                                                                       BFA 0093
C     ZERO INITIAL FILL-IN COLUMNS                                      BFA 0094
C                                                                       BFA 0095
      J0 = MU + 2                                                       BFA 0096
      J1 = MIN0(N,M) - 1                                                BFA 0097
      IF (J1 .LT. J0) GO TO 30                                          BFA 0098
      DO 20 JZ = J0, J1                                                 BFA 0099
         I0 = M + 1 - JZ                                                BFA 0100
         DO 10 I = I0, ML                                               BFA 0101
            ABD(I,JZ) = 0.0D0                                           BFA 0102
   10    CONTINUE                                                       BFA 0103
   20 CONTINUE                                                          BFA 0104
   30 CONTINUE                                                          BFA 0105
      JZ = J1                                                           BFA 0106
      JU = 0                                                            BFA 0107
C                                                                       BFA 0108
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        BFA 0109
C                                                                       BFA 0110
      NM1 = N - 1                                                       BFA 0111
      IF (NM1 .LT. 1) GO TO 130                                         BFA 0112
      DO 120 K = 1, NM1                                                 BFA 0113
         KP1 = K + 1                                                    BFA 0114
C                                                                       BFA 0115
C        ZERO NEXT FILL-IN COLUMN                                       BFA 0116
C                                                                       BFA 0117
         JZ = JZ + 1                                                    BFA 0118
         IF (JZ .GT. N) GO TO 50                                        BFA 0119
         IF (ML .LT. 1) GO TO 50                                        BFA 0120
            DO 40 I = 1, ML                                             BFA 0121
               ABD(I,JZ) = 0.0D0                                        BFA 0122
   40       CONTINUE                                                    BFA 0123
   50    CONTINUE                                                       BFA 0124
C                                                                       BFA 0125
C        FIND L = PIVOT INDEX                                           BFA 0126
C                                                                       BFA 0127
         LM = MIN0(ML,N-K)                                              BFA 0128
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1                            BFA 0129
         IPVT(K) = L + K - M                                            BFA 0130
C                                                                       BFA 0131
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          BFA 0132
C                                                                       BFA 0133
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100                             BFA 0134
C                                                                       BFA 0135
C           INTERCHANGE IF NECESSARY                                    BFA 0136
C                                                                       BFA 0137
            IF (L .EQ. M) GO TO 60                                      BFA 0138
               T = ABD(L,K)                                             BFA 0139
               ABD(L,K) = ABD(M,K)                                      BFA 0140
               ABD(M,K) = T                                             BFA 0141
   60       CONTINUE                                                    BFA 0142
C                                                                       BFA 0143
C           COMPUTE MULTIPLIERS                                         BFA 0144
C                                                                       BFA 0145
            T = -1.0D0/ABD(M,K)                                         BFA 0146
            CALL DSCAL(LM,T,ABD(M+1,K),1)                               BFA 0147
C                                                                       BFA 0148
C           ROW ELIMINATION WITH COLUMN INDEXING                        BFA 0149
C                                                                       BFA 0150
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                            BFA 0151
            MM = M                                                      BFA 0152
            IF (JU .LT. KP1) GO TO 90                                   BFA 0153
            DO 80 J = KP1, JU                                           BFA 0154
               L = L - 1                                                BFA 0155
               MM = MM - 1                                              BFA 0156
               T = ABD(L,J)                                             BFA 0157
               IF (L .EQ. MM) GO TO 70                                  BFA 0158
                  ABD(L,J) = ABD(MM,J)                                  BFA 0159
                  ABD(MM,J) = T                                         BFA 0160
   70          CONTINUE                                                 BFA 0161
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)              BFA 0162
   80       CONTINUE                                                    BFA 0163
   90       CONTINUE                                                    BFA 0164
         GO TO 110                                                      BFA 0165
  100    CONTINUE                                                       BFA 0166
            INFO = K                                                    BFA 0167
  110    CONTINUE                                                       BFA 0168
  120 CONTINUE                                                          BFA 0169
  130 CONTINUE                                                          BFA 0170
      IPVT(N) = N                                                       BFA 0171
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N                                 BFA 0172
      RETURN                                                            BFA 0173
      END                                                               BFA 0174
