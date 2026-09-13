      SUBROUTINE DCHDC(A,LDA,P,WORK,JPVT,JOB,INFO)                      DCH 0001
      INTEGER LDA,P,JPVT(*),JOB,INFO                                    DCH 0002
      DOUBLE PRECISION A(LDA,*),WORK(*)                                 DCH 0003
C                                                                       DCH 0004
C     DCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE  DCH 0005
C     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE        DCH 0006
C     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK     DCH 0007
C     OF A POSITIVE SEMIDEFINITE MATRIX.                                DCH 0008
C                                                                       DCH 0009
C     ON ENTRY                                                          DCH 0010
C                                                                       DCH 0011
C         A      DOUBLE PRECISION(LDA,P).                               DCH 0012
C                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO        DCH 0013
C                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED. DCH 0014
C                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED.       DCH 0015
C                                                                       DCH 0016
C         LDA    INTEGER.                                               DCH 0017
C                LDA IS THE LEADING DIMENSION OF THE ARRAY A.           DCH 0018
C                                                                       DCH 0019
C         P      INTEGER.                                               DCH 0020
C                P IS THE ORDER OF THE MATRIX.                          DCH 0021
C                                                                       DCH 0022
C         WORK   DOUBLE PRECISION.                                      DCH 0023
C                WORK IS A WORK ARRAY.                                  DCH 0024
C                                                                       DCH 0025
C         JPVT   INTEGER(P).                                            DCH 0026
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION      DCH 0027
C                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED. DCH 0028
C                EACH DIAGONAL ELEMENT A(K,K)                           DCH 0029
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE     DCH 0030
C                VALUE OF JPVT(K).                                      DCH 0031
C                                                                       DCH 0032
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL          DCH 0033
C                                      ELEMENT.                         DCH 0034
C                                                                       DCH 0035
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT.     DCH 0036
C                                                                       DCH 0037
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT.    DCH 0038
C                                                                       DCH 0039
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS DCH 0040
C                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO  DCH 0041
C                THE BEGINNING OF THE ARRAY A AND FINAL                 DCH 0042
C                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS  DCH 0043
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY    DCH 0044
C                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE     DCH 0045
C                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT     DCH 0046
C                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT       DCH 0047
C                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF       DCH 0048
C                JOB .EQ. 0.                                            DCH 0049
C                                                                       DCH 0050
C        JOB     INTEGER.                                               DCH 0051
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.      DCH 0052
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.                    DCH 0053
C                IF JOB .NE. 0, PIVOTING IS DONE.                       DCH 0054
C                                                                       DCH 0055
C     ON RETURN                                                         DCH 0056
C                                                                       DCH 0057
C         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR       DCH 0058
C                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING.   DCH 0059
C                                                                       DCH 0060
C         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT     DCH 0061
C                OF A THAT WAS MOVED INTO THE J-TH POSITION,            DCH 0062
C                PROVIDED PIVOTING WAS REQUESTED.                       DCH 0063
C                                                                       DCH 0064
C         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL       DCH 0065
C                ELEMENT OF THE CHOLESKY FACTOR.                        DCH 0066
C                                                                       DCH 0067
C     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN.     DCH 0068
C     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL        DCH 0069
C     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN     DCH 0070
C     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO   DCH 0071
C     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE      DCH 0072
C     INFO TO BE LESS THAN P.                                           DCH 0073
C                                                                       DCH 0074
C     LINPACK. THIS VERSION DATED 03/19/79 .                            DCH 0075
C     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND   DCH 0076
C     UNIVERSITY OF MARYLAND.                                           DCH 0077
C                                                                       DCH 0078
C                                                                       DCH 0079
C     BLAS DAXPY,DSWAP                                                  DCH 0080
C     FORTRAN DSQRT                                                     DCH 0081
C                                                                       DCH 0082
C     INTERNAL VARIABLES                                                DCH 0083
C                                                                       DCH 0084
      INTEGER PU,PL,PLP1,J,JP,JT,K,KB,KM1,KP1,L,MAXL                    DCH 0085
      DOUBLE PRECISION TEMP                                             DCH 0086
      DOUBLE PRECISION MAXDIA                                           DCH 0087
      LOGICAL SWAPK,NEGK                                                DCH 0088
C                                                                       DCH 0089
      PL = 1                                                            DCH 0090
      PU = 0                                                            DCH 0091
      INFO = P                                                          DCH 0092
      IF (JOB .EQ. 0) GO TO 160                                         DCH 0093
C                                                                       DCH 0094
C        PIVOTING HAS BEEN REQUESTED. REARRANGE THE                     DCH 0095
C        THE ELEMENTS ACCORDING TO JPVT.                                DCH 0096
C                                                                       DCH 0097
         DO 70 K = 1, P                                                 DCH 0098
            SWAPK = JPVT(K) .GT. 0                                      DCH 0099
            NEGK = JPVT(K) .LT. 0                                       DCH 0100
            JPVT(K) = K                                                 DCH 0101
            IF (NEGK) JPVT(K) = -JPVT(K)                                DCH 0102
            IF (.NOT.SWAPK) GO TO 60                                    DCH 0103
               IF (K .EQ. PL) GO TO 50                                  DCH 0104
                  CALL DSWAP(PL-1,A(1,K),1,A(1,PL),1)                   DCH 0105
                  TEMP = A(K,K)                                         DCH 0106
                  A(K,K) = A(PL,PL)                                     DCH 0107
                  A(PL,PL) = TEMP                                       DCH 0108
                  PLP1 = PL + 1                                         DCH 0109
                  IF (P .LT. PLP1) GO TO 40                             DCH 0110
                  DO 30 J = PLP1, P                                     DCH 0111
                     IF (J .GE. K) GO TO 10                             DCH 0112
                        TEMP = A(PL,J)                                  DCH 0113
                        A(PL,J) = A(J,K)                                DCH 0114
                        A(J,K) = TEMP                                   DCH 0115
                     GO TO 20                                           DCH 0116
   10                CONTINUE                                           DCH 0117
                     IF (J .EQ. K) GO TO 20                             DCH 0118
                        TEMP = A(K,J)                                   DCH 0119
                        A(K,J) = A(PL,J)                                DCH 0120
                        A(PL,J) = TEMP                                  DCH 0121
   20                CONTINUE                                           DCH 0122
   30             CONTINUE                                              DCH 0123
   40             CONTINUE                                              DCH 0124
                  JPVT(K) = JPVT(PL)                                    DCH 0125
                  JPVT(PL) = K                                          DCH 0126
   50          CONTINUE                                                 DCH 0127
               PL = PL + 1                                              DCH 0128
   60       CONTINUE                                                    DCH 0129
   70    CONTINUE                                                       DCH 0130
         PU = P                                                         DCH 0131
         IF (P .LT. PL) GO TO 150                                       DCH 0132
         DO 140 KB = PL, P                                              DCH 0133
            K = P - KB + PL                                             DCH 0134
            IF (JPVT(K) .GE. 0) GO TO 130                               DCH 0135
               JPVT(K) = -JPVT(K)                                       DCH 0136
               IF (PU .EQ. K) GO TO 120                                 DCH 0137
                  CALL DSWAP(K-1,A(1,K),1,A(1,PU),1)                    DCH 0138
                  TEMP = A(K,K)                                         DCH 0139
                  A(K,K) = A(PU,PU)                                     DCH 0140
                  A(PU,PU) = TEMP                                       DCH 0141
                  KP1 = K + 1                                           DCH 0142
                  IF (P .LT. KP1) GO TO 110                             DCH 0143
                  DO 100 J = KP1, P                                     DCH 0144
                     IF (J .GE. PU) GO TO 80                            DCH 0145
                        TEMP = A(K,J)                                   DCH 0146
                        A(K,J) = A(J,PU)                                DCH 0147
                        A(J,PU) = TEMP                                  DCH 0148
                     GO TO 90                                           DCH 0149
   80                CONTINUE                                           DCH 0150
                     IF (J .EQ. PU) GO TO 90                            DCH 0151
                        TEMP = A(K,J)                                   DCH 0152
                        A(K,J) = A(PU,J)                                DCH 0153
                        A(PU,J) = TEMP                                  DCH 0154
   90                CONTINUE                                           DCH 0155
  100             CONTINUE                                              DCH 0156
  110             CONTINUE                                              DCH 0157
                  JT = JPVT(K)                                          DCH 0158
                  JPVT(K) = JPVT(PU)                                    DCH 0159
                  JPVT(PU) = JT                                         DCH 0160
  120          CONTINUE                                                 DCH 0161
               PU = PU - 1                                              DCH 0162
  130       CONTINUE                                                    DCH 0163
  140    CONTINUE                                                       DCH 0164
  150    CONTINUE                                                       DCH 0165
  160 CONTINUE                                                          DCH 0166
      DO 270 K = 1, P                                                   DCH 0167
C                                                                       DCH 0168
C        REDUCTION LOOP.                                                DCH 0169
C                                                                       DCH 0170
         MAXDIA = A(K,K)                                                DCH 0171
         KP1 = K + 1                                                    DCH 0172
         MAXL = K                                                       DCH 0173
C                                                                       DCH 0174
C        DETERMINE THE PIVOT ELEMENT.                                   DCH 0175
C                                                                       DCH 0176
         IF (K .LT. PL .OR. K .GE. PU) GO TO 190                        DCH 0177
            DO 180 L = KP1, PU                                          DCH 0178
               IF (A(L,L) .LE. MAXDIA) GO TO 170                        DCH 0179
                  MAXDIA = A(L,L)                                       DCH 0180
                  MAXL = L                                              DCH 0181
  170          CONTINUE                                                 DCH 0182
  180       CONTINUE                                                    DCH 0183
  190    CONTINUE                                                       DCH 0184
C                                                                       DCH 0185
C        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.                     DCH 0186
C                                                                       DCH 0187
         IF (MAXDIA .GT. 0.0D0) GO TO 200                               DCH 0188
            INFO = K - 1                                                DCH 0189
C     ......EXIT                                                        DCH 0190
            GO TO 280                                                   DCH 0191
  200    CONTINUE                                                       DCH 0192
         IF (K .EQ. MAXL) GO TO 210                                     DCH 0193
C                                                                       DCH 0194
C           START THE PIVOTING AND UPDATE JPVT.                         DCH 0195
C                                                                       DCH 0196
            KM1 = K - 1                                                 DCH 0197
            CALL DSWAP(KM1,A(1,K),1,A(1,MAXL),1)                        DCH 0198
            A(MAXL,MAXL) = A(K,K)                                       DCH 0199
            A(K,K) = MAXDIA                                             DCH 0200
            JP = JPVT(MAXL)                                             DCH 0201
            JPVT(MAXL) = JPVT(K)                                        DCH 0202
            JPVT(K) = JP                                                DCH 0203
  210    CONTINUE                                                       DCH 0204
C                                                                       DCH 0205
C        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.         DCH 0206
C                                                                       DCH 0207
         WORK(K) = DSQRT(A(K,K))                                        DCH 0208
         A(K,K) = WORK(K)                                               DCH 0209
         IF (P .LT. KP1) GO TO 260                                      DCH 0210
         DO 250 J = KP1, P                                              DCH 0211
            IF (K .EQ. MAXL) GO TO 240                                  DCH 0212
               IF (J .GE. MAXL) GO TO 220                               DCH 0213
                  TEMP = A(K,J)                                         DCH 0214
                  A(K,J) = A(J,MAXL)                                    DCH 0215
                  A(J,MAXL) = TEMP                                      DCH 0216
               GO TO 230                                                DCH 0217
  220          CONTINUE                                                 DCH 0218
               IF (J .EQ. MAXL) GO TO 230                               DCH 0219
                  TEMP = A(K,J)                                         DCH 0220
                  A(K,J) = A(MAXL,J)                                    DCH 0221
                  A(MAXL,J) = TEMP                                      DCH 0222
  230          CONTINUE                                                 DCH 0223
  240       CONTINUE                                                    DCH 0224
            A(K,J) = A(K,J)/WORK(K)                                     DCH 0225
            WORK(J) = A(K,J)                                            DCH 0226
            TEMP = -A(K,J)                                              DCH 0227
            CALL DAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)                 DCH 0228
  250    CONTINUE                                                       DCH 0229
  260    CONTINUE                                                       DCH 0230
  270 CONTINUE                                                          DCH 0231
  280 CONTINUE                                                          DCH 0232
      RETURN                                                            DCH 0233
      END                                                               DCH 0234
