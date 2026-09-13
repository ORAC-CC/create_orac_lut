      SUBROUTINE DTRDI(T,LDT,N,DET,JOB,INFO)                            DET 0001
      INTEGER LDT,N,JOB,INFO                                            DET 0002
      DOUBLE PRECISION T(LDT,*),DET(2)                                  DET 0003
C                                                                       DET 0004
C     DTRDI COMPUTES THE DETERMINANT AND INVERSE OF A DOUBLE PRECISION  DET 0005
C     TRIANGULAR MATRIX.                                                DET 0006
C                                                                       DET 0007
C     ON ENTRY                                                          DET 0008
C                                                                       DET 0009
C        T       DOUBLE PRECISION(LDT,N)                                DET 0010
C                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO             DET 0011
C                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND         DET 0012
C                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE         DET 0013
C                USED TO STORE OTHER INFORMATION.                       DET 0014
C                                                                       DET 0015
C        LDT     INTEGER                                                DET 0016
C                LDT IS THE LEADING DIMENSION OF THE ARRAY T.           DET 0017
C                                                                       DET 0018
C        N       INTEGER                                                DET 0019
C                N IS THE ORDER OF THE SYSTEM.                          DET 0020
C                                                                       DET 0021
C        JOB     INTEGER                                                DET 0022
C                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.       DET 0023
C                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.       DET 0024
C                = 100       DET, NO INVERSE.                           DET 0025
C                = 110       DET, INVERSE OF LOWER TRIANGULAR.          DET 0026
C                = 111       DET, INVERSE OF UPPER TRIANGULAR.          DET 0027
C                                                                       DET 0028
C     ON RETURN                                                         DET 0029
C                                                                       DET 0030
C        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.               DET 0031
C                OTHERWISE UNCHANGED.                                   DET 0032
C                                                                       DET 0033
C        DET     DOUBLE PRECISION(2)                                    DET 0034
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.           DET 0035
C                OTHERWISE NOT REFERENCED.                              DET 0036
C                DETERMINANT = DET(1) * 10.0**DET(2)                    DET 0037
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0                  DET 0038
C                OR  DET(1) .EQ. 0.0 .                                  DET 0039
C                                                                       DET 0040
C        INFO    INTEGER                                                DET 0041
C                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR        DET 0042
C                AND THE INVERSE IS REQUESTED.                          DET 0043
C                OTHERWISE INFO CONTAINS THE INDEX OF                   DET 0044
C                A ZERO DIAGONAL ELEMENT OF T.                          DET 0045
C                                                                       DET 0046
C                                                                       DET 0047
C     LINPACK. THIS VERSION DATED 08/14/78 .                            DET 0048
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      DET 0049
C                                                                       DET 0050
C     SUBROUTINES AND FUNCTIONS                                         DET 0051
C                                                                       DET 0052
C     BLAS DAXPY,DSCAL                                                  DET 0053
C     FORTRAN DABS,MOD                                                  DET 0054
C                                                                       DET 0055
C     INTERNAL VARIABLES                                                DET 0056
C                                                                       DET 0057
      DOUBLE PRECISION TEMP                                             DET 0058
      DOUBLE PRECISION TEN                                              DET 0059
      INTEGER I,J,K,KB,KM1,KP1                                          DET 0060
C                                                                       DET 0061
C     BEGIN BLOCK PERMITTING ...EXITS TO 180                            DET 0062
C                                                                       DET 0063
C        COMPUTE DETERMINANT                                            DET 0064
C                                                                       DET 0065
         IF (JOB/100 .EQ. 0) GO TO 70                                   DET 0066
            DET(1) = 1.0D0                                              DET 0067
            DET(2) = 0.0D0                                              DET 0068
            TEN = 10.0D0                                                DET 0069
            DO 50 I = 1, N                                              DET 0070
               DET(1) = T(I,I)*DET(1)                                   DET 0071
C           ...EXIT                                                     DET 0072
               IF (DET(1) .EQ. 0.0D0) GO TO 60                          DET 0073
   10          IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20                    DET 0074
                  DET(1) = TEN*DET(1)                                   DET 0075
                  DET(2) = DET(2) - 1.0D0                               DET 0076
               GO TO 10                                                 DET 0077
   20          CONTINUE                                                 DET 0078
   30          IF (DABS(DET(1)) .LT. TEN) GO TO 40                      DET 0079
                  DET(1) = DET(1)/TEN                                   DET 0080
                  DET(2) = DET(2) + 1.0D0                               DET 0081
               GO TO 30                                                 DET 0082
   40          CONTINUE                                                 DET 0083
   50       CONTINUE                                                    DET 0084
   60       CONTINUE                                                    DET 0085
   70    CONTINUE                                                       DET 0086
C                                                                       DET 0087
C        COMPUTE INVERSE OF UPPER TRIANGULAR                            DET 0088
C                                                                       DET 0089
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170                           DET 0090
            IF (MOD(JOB,10) .EQ. 0) GO TO 120                           DET 0091
C              BEGIN BLOCK PERMITTING ...EXITS TO 110                   DET 0092
                  DO 100 K = 1, N                                       DET 0093
                     INFO = K                                           DET 0094
C              ......EXIT                                               DET 0095
                     IF (T(K,K) .EQ. 0.0D0) GO TO 110                   DET 0096
                     T(K,K) = 1.0D0/T(K,K)                              DET 0097
                     TEMP = -T(K,K)                                     DET 0098
                     CALL DSCAL(K-1,TEMP,T(1,K),1)                      DET 0099
                     KP1 = K + 1                                        DET 0100
                     IF (N .LT. KP1) GO TO 90                           DET 0101
                     DO 80 J = KP1, N                                   DET 0102
                        TEMP = T(K,J)                                   DET 0103
                        T(K,J) = 0.0D0                                  DET 0104
                        CALL DAXPY(K,TEMP,T(1,K),1,T(1,J),1)            DET 0105
   80                CONTINUE                                           DET 0106
   90                CONTINUE                                           DET 0107
  100             CONTINUE                                              DET 0108
                  INFO = 0                                              DET 0109
  110          CONTINUE                                                 DET 0110
            GO TO 160                                                   DET 0111
  120       CONTINUE                                                    DET 0112
C                                                                       DET 0113
C              COMPUTE INVERSE OF LOWER TRIANGULAR                      DET 0114
C                                                                       DET 0115
               DO 150 KB = 1, N                                         DET 0116
                  K = N + 1 - KB                                        DET 0117
                  INFO = K                                              DET 0118
C     ............EXIT                                                  DET 0119
                  IF (T(K,K) .EQ. 0.0D0) GO TO 180                      DET 0120
                  T(K,K) = 1.0D0/T(K,K)                                 DET 0121
                  TEMP = -T(K,K)                                        DET 0122
                  IF (K .NE. N) CALL DSCAL(N-K,TEMP,T(K+1,K),1)         DET 0123
                  KM1 = K - 1                                           DET 0124
                  IF (KM1 .LT. 1) GO TO 140                             DET 0125
                  DO 130 J = 1, KM1                                     DET 0126
                     TEMP = T(K,J)                                      DET 0127
                     T(K,J) = 0.0D0                                     DET 0128
                     CALL DAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)           DET 0129
  130             CONTINUE                                              DET 0130
  140             CONTINUE                                              DET 0131
  150          CONTINUE                                                 DET 0132
               INFO = 0                                                 DET 0133
  160       CONTINUE                                                    DET 0134
  170    CONTINUE                                                       DET 0135
  180 CONTINUE                                                          DET 0136
      RETURN                                                            DET 0137
      END                                                               DET 0138
