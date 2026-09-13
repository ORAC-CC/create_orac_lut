      SUBROUTINE  ASYMTX( A, EVEC, EVAL, M, IA, IEVEC, IER, WK,         ASM 0001
     $                    AAD, EVECD, EVALD, WKD )                      ASM 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            ASM 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                ASM 0004
                                                                        ASM 0005
C    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ====== ASM 0006
                                                                        ASM 0007
C       SOLVES EIGENFUNCTION PROBLEM FOR REAL ASYMMETRIC MATRIX         ASM 0008
C       FOR WHICH IT IS KNOWN A PRIORI THAT THE EIGENVALUES ARE REAL.   ASM 0009
                                                                        ASM 0010
C       THIS IS AN ADAPTATION OF A SUBROUTINE EIGRF IN THE IMSL         ASM 0011
C       LIBRARY TO USE REAL INSTEAD OF COMPLEX ARITHMETIC, ACCOUNTING   ASM 0012
C       FOR THE KNOWN FACT THAT THE EIGENVALUES AND EIGENVECTORS IN     ASM 0013
C       THE DISCRETE ORDINATE SOLUTION ARE REAL.  OTHER CHANGES INCLUDE ASM 0014
C       PUTTING ALL THE CALLED SUBROUTINES IN-LINE, DELETING THE        ASM 0015
C       PERFORMANCE INDEX CALCULATION, UPDATING MANY DO-LOOPS           ASM 0016
C       TO FORTRAN77, AND IN CALCULATING THE MACHINE PRECISION          ASM 0017
C       TOL INSTEAD OF SPECIFYING IT IN A DATA STATEMENT.               ASM 0018
                                                                        ASM 0019
C       EIGRF IS BASED PRIMARILY ON EISPACK ROUTINES.  THE MATRIX IS    ASM 0020
C       FIRST BALANCED USING THE PARLETT-REINSCH ALGORITHM.  THEN       ASM 0021
C       THE MARTIN-WILKINSON ALGORITHM IS APPLIED.                      ASM 0022
                                                                        ASM 0023
C       REFERENCES:                                                     ASM 0024
C          DONGARRA, J. AND C. MOLER, EISPACK -- A PACKAGE FOR SOLVING  ASM 0025
C             MATRIX EIGENVALUE PROBLEMS, IN COWELL, ED., 1984:         ASM 0026
C             SOURCES AND DEVELOPMENT OF MATHEMATICAL SOFTWARE,         ASM 0027
C             PRENTICE-HALL, ENGLEWOOD CLIFFS, NJ                       ASM 0028
C         PARLETT AND REINSCH, 1969: BALANCING A MATRIX FOR CALCULATION ASM 0029
C             OF EIGENVALUES AND EIGENVECTORS, NUM. MATH. 13, 293-304   ASM 0030
C         WILKINSON, J., 1965: THE ALGEBRAIC EIGENVALUE PROBLEM,        ASM 0031
C             CLARENDON PRESS, OXFORD                                   ASM 0032
                                                                        ASM 0033
C   I N P U T    V A R I A B L E S:                                     ASM 0034
                                                                        ASM 0035
C        A    :  INPUT ASYMMETRIC MATRIX, DESTROYED AFTER SOLVED        ASM 0036
C        M    :  ORDER OF  A                                            ASM 0037
C       IA    :  FIRST DIMENSION OF  A                                  ASM 0038
C    IEVEC    :  FIRST DIMENSION OF  EVEC                               ASM 0039
                                                                        ASM 0040
C   O U T P U T    V A R I A B L E S:                                   ASM 0041
                                                                        ASM 0042
C       EVEC  :  (UNNORMALIZED) EIGENVECTORS OF  A                      ASM 0043
C                   ( COLUMN J CORRESPONDS TO EVAL(J) )                 ASM 0044
                                                                        ASM 0045
C       EVAL  :  (UNORDERED) EIGENVALUES OF  A ( DIMENSION AT LEAST M ) ASM 0046
                                                                        ASM 0047
C       IER   :  IF .NE. 0, SIGNALS THAT EVAL(IER) FAILED TO CONVERGE;  ASM 0048
C                   IN THAT CASE EIGENVALUES IER+1,IER+2,...,M  ARE     ASM 0049
C                   CORRECT BUT EIGENVALUES 1,...,IER ARE SET TO ZERO.  ASM 0050
                                                                        ASM 0051
C   S C R A T C H   V A R I A B L E S:                                  ASM 0052
                                                                        ASM 0053
C       WK    :  WORK AREA ( DIMENSION AT LEAST 2*M )                   ASM 0054
C       AAD   :  DOUBLE PRECISION STAND-IN FOR -A-                      ASM 0055
C       EVECD :  DOUBLE PRECISION STAND-IN FOR -EVEC-                   ASM 0056
C       EVALD :  DOUBLE PRECISION STAND-IN FOR -EVAL-                   ASM 0057
C       WKD   :  DOUBLE PRECISION STAND-IN FOR -WK-                     ASM 0058
C+---------------------------------------------------------------------+ASM 0059
                                                                        ASM 0060
C      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                           ASM 0061
      REAL*8              A( IA,* ),   WK(*),  EVAL(*),  EVEC( IEVEC,* )ASM 0062
      DOUBLE PRECISION  AAD( IA,* ), WKD(*), EVALD(*), EVECD( IA,* )    ASM 0063
      DOUBLE PRECISION  D1MACH                                          ASM 0064
      LOGICAL           NOCONV, NOTLAS                                  ASM 0065
      DATA     C1 / 0.4375D0 /, C2/ 0.5D0 /, C3/ 0.75D0 /, C4/ 0.95D0 /,ASM 0066
     $         C5/ 16.D0 /, C6/ 256.D0 /, ZERO / 0.D0 /, ONE / 1.D0 /   ASM 0067
                                                                        ASM 0068
                                                                        ASM 0069
      IER = 0                                                           ASM 0070
      TOL = D1MACH(3)                                                   ASM 0071
      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )                        ASM 0072
     $     CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )       ASM 0073
                                                                        ASM 0074
C                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES         ASM 0075
      IF ( M.EQ.1 )  THEN                                               ASM 0076
         EVAL(1) = A(1,1)                                               ASM 0077
         EVEC(1,1) = 1.0                                                ASM 0078
         RETURN                                                         ASM 0079
                                                                        ASM 0080
      ELSE IF ( M.EQ.2 )  THEN                                          ASM 0081
         DISCRI = ( A(1,1) - A(2,2) )**2 + 4. * A(1,2) * A(2,1)         ASM 0082
         IF ( DISCRI.LT.0.0 )                                           ASM 0083
     $        CALL ERRMSG( 'ASYMTX--COMPLEX EVALS IN 2X2 CASE', .TRUE. )ASM 0084
         SGN = 1.0                                                      ASM 0085
         IF ( A(1,1).LT.A(2,2) )  SGN = - 1.0                           ASM 0086
         EVAL(1) = 0.5 * ( A(1,1) + A(2,2) + SGN*DSQRT(DISCRI) )        ASM 0087
         EVAL(2) = 0.5 * ( A(1,1) + A(2,2) - SGN*DSQRT(DISCRI) )        ASM 0088
         EVEC(1,1) = 1.0                                                ASM 0089
         EVEC(2,2) = 1.0                                                ASM 0090
         IF ( A(1,1).EQ.A(2,2) .AND. (A(2,1).EQ.0.0.OR.A(1,2).EQ.0.0) ) ASM 0091
     $        THEN                                                      ASM 0092
       RNORM=DABS(A(1,1))+DABS(A(1,2))+DABS(A(2,1))+DABS(A(2,2))        ASM 0093
            W = TOL * RNORM                                             ASM 0094
            EVEC(2,1) = A(2,1) / W                                      ASM 0095
            EVEC(1,2) = - A(1,2) / W                                    ASM 0096
         ELSE                                                           ASM 0097
            EVEC(2,1) = A(2,1) / ( EVAL(1) - A(2,2) )                   ASM 0098
            EVEC(1,2) = A(1,2) / ( EVAL(2) - A(1,1) )                   ASM 0099
         ENDIF                                                          ASM 0100
         RETURN                                                         ASM 0101
      END IF                                                            ASM 0102
C                               ** PUT S.P. MATRIX INTO D.P. MATRIX     ASM 0103
      DO 1  J = 1, M                                                    ASM 0104
         DO 1  K = 1, M                                                 ASM 0105
            AAD( J,K ) = DBLE( A(J,K) )                                 ASM 0106
    1 CONTINUE                                                          ASM 0107
C                                        ** INITIALIZE OUTPUT VARIABLES ASM 0108
      DO 20 I = 1, M                                                    ASM 0109
         EVALD(I) = ZERO                                                ASM 0110
         DO 10 J = 1, M                                                 ASM 0111
            EVECD(I,J) = ZERO                                           ASM 0112
10       CONTINUE                                                       ASM 0113
         EVECD(I,I) = ONE                                               ASM 0114
20    CONTINUE                                                          ASM 0115
C                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY   ASM 0116
C                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;  ASM 0117
C                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE      ASM 0118
C                  ** AND PUSH THEM DOWN                                ASM 0119
      RNORM = ZERO                                                      ASM 0120
      L  = 1                                                            ASM 0121
      K  = M                                                            ASM 0122
                                                                        ASM 0123
30    KKK = K                                                           ASM 0124
         DO 70  J = KKK, 1, -1                                          ASM 0125
            ROW = ZERO                                                  ASM 0126
            DO 40 I = 1, K                                              ASM 0127
               IF ( I.NE.J ) ROW = ROW + DABS( AAD(J,I) )               ASM 0128
40          CONTINUE                                                    ASM 0129
            IF ( ROW.EQ.ZERO ) THEN                                     ASM 0130
               WKD(K) = J                                               ASM 0131
               IF ( J.NE.K ) THEN                                       ASM 0132
                  DO 50 I = 1, K                                        ASM 0133
                     REPL   = AAD(I,J)                                  ASM 0134
                     AAD(I,J) = AAD(I,K)                                ASM 0135
                     AAD(I,K) = REPL                                    ASM 0136
50                CONTINUE                                              ASM 0137
                  DO 60 I = L, M                                        ASM 0138
                     REPL   = AAD(J,I)                                  ASM 0139
                     AAD(J,I) = AAD(K,I)                                ASM 0140
                     AAD(K,I) = REPL                                    ASM 0141
60                CONTINUE                                              ASM 0142
               END IF                                                   ASM 0143
               K = K - 1                                                ASM 0144
               GO TO 30                                                 ASM 0145
            END IF                                                      ASM 0146
70       CONTINUE                                                       ASM 0147
C                                     ** SEARCH FOR COLUMNS ISOLATING ANASM 0148
C                                       ** EIGENVALUE AND PUSH THEM LEFTASM 0149
80    LLL = L                                                           ASM 0150
         DO 120 J = LLL, K                                              ASM 0151
            COL = ZERO                                                  ASM 0152
            DO 90 I = L, K                                              ASM 0153
               IF ( I.NE.J ) COL = COL + DABS( AAD(I,J) )               ASM 0154
90          CONTINUE                                                    ASM 0155
            IF ( COL.EQ.ZERO ) THEN                                     ASM 0156
               WKD(L) = J                                               ASM 0157
               IF ( J.NE.L ) THEN                                       ASM 0158
                  DO 100 I = 1, K                                       ASM 0159
                     REPL   = AAD(I,J)                                  ASM 0160
                     AAD(I,J) = AAD(I,L)                                ASM 0161
                     AAD(I,L) = REPL                                    ASM 0162
100               CONTINUE                                              ASM 0163
                  DO 110 I = L, M                                       ASM 0164
                     REPL   = AAD(J,I)                                  ASM 0165
                     AAD(J,I) = AAD(L,I)                                ASM 0166
                     AAD(L,I) = REPL                                    ASM 0167
110               CONTINUE                                              ASM 0168
               END IF                                                   ASM 0169
               L = L + 1                                                ASM 0170
               GO TO 80                                                 ASM 0171
            END IF                                                      ASM 0172
120      CONTINUE                                                       ASM 0173
C                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH KASM 0174
      DO 130 I = L, K                                                   ASM 0175
         WKD(I) = ONE                                                   ASM 0176
130   CONTINUE                                                          ASM 0177
                                                                        ASM 0178
140   NOCONV = .FALSE.                                                  ASM 0179
         DO 200 I = L, K                                                ASM 0180
            COL = ZERO                                                  ASM 0181
            ROW = ZERO                                                  ASM 0182
            DO 150 J = L, K                                             ASM 0183
               IF ( J.NE.I ) THEN                                       ASM 0184
                  COL = COL + DABS( AAD(J,I) )                          ASM 0185
                  ROW = ROW + DABS( AAD(I,J) )                          ASM 0186
               END IF                                                   ASM 0187
150         CONTINUE                                                    ASM 0188
            F = ONE                                                     ASM 0189
            G = ROW / C5                                                ASM 0190
            H = COL + ROW                                               ASM 0191
160         IF ( COL.LT.G ) THEN                                        ASM 0192
               F   = F * C5                                             ASM 0193
               COL = COL * C6                                           ASM 0194
               GO TO 160                                                ASM 0195
            END IF                                                      ASM 0196
            G = ROW * C5                                                ASM 0197
170         IF ( COL.GE.G ) THEN                                        ASM 0198
               F   = F / C5                                             ASM 0199
               COL = COL / C6                                           ASM 0200
               GO TO 170                                                ASM 0201
            END IF                                                      ASM 0202
C                                                         ** NOW BALANCEASM 0203
            IF ( (COL+ROW)/F .LT. C4*H ) THEN                           ASM 0204
               WKD(I)  = WKD(I) * F                                     ASM 0205
               NOCONV = .TRUE.                                          ASM 0206
               DO 180 J = L, M                                          ASM 0207
                  AAD(I,J) = AAD(I,J) / F                               ASM 0208
180            CONTINUE                                                 ASM 0209
               DO 190 J = 1, K                                          ASM 0210
                  AAD(J,I) = AAD(J,I) * F                               ASM 0211
190            CONTINUE                                                 ASM 0212
            END IF                                                      ASM 0213
200      CONTINUE                                                       ASM 0214
                                                                        ASM 0215
      IF ( NOCONV ) GO TO 140                                           ASM 0216
C                                  ** IS -A- ALREADY IN HESSENBERG FORM?ASM 0217
      IF ( K-1 .LT. L+1 ) GO TO 350                                     ASM 0218
C                                   ** TRANSFER -A- TO A HESSENBERG FORMASM 0219
      DO 290 N = L+1, K-1                                               ASM 0220
         H        = ZERO                                                ASM 0221
         WKD(N+M) = ZERO                                                ASM 0222
         SCALE    = ZERO                                                ASM 0223
C                                                        ** SCALE COLUMNASM 0224
         DO 210 I = N, K                                                ASM 0225
            SCALE = SCALE + DABS(AAD(I,N-1))                            ASM 0226
210      CONTINUE                                                       ASM 0227
         IF ( SCALE.NE.ZERO ) THEN                                      ASM 0228
            DO 220 I = K, N, -1                                         ASM 0229
               WKD(I+M) = AAD(I,N-1) / SCALE                            ASM 0230
               H = H + WKD(I+M)**2                                      ASM 0231
220         CONTINUE                                                    ASM 0232
            G = - SIGN( DSQRT(H), WKD(N+M) )                            ASM 0233
            H = H - WKD(N+M) * G                                        ASM 0234
            WKD(N+M) = WKD(N+M) - G                                     ASM 0235
C                                                 ** FORM (I-(U*UT)/H)*AASM 0236
            DO 250 J = N, M                                             ASM 0237
               F = ZERO                                                 ASM 0238
               DO 230  I = K, N, -1                                     ASM 0239
                  F = F + WKD(I+M) * AAD(I,J)                           ASM 0240
230            CONTINUE                                                 ASM 0241
               DO 240 I = N, K                                          ASM 0242
                  AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H                ASM 0243
240            CONTINUE                                                 ASM 0244
250         CONTINUE                                                    ASM 0245
C                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)ASM 0246
            DO 280 I = 1, K                                             ASM 0247
               F = ZERO                                                 ASM 0248
               DO 260  J = K, N, -1                                     ASM 0249
                  F = F + WKD(J+M) * AAD(I,J)                           ASM 0250
260            CONTINUE                                                 ASM 0251
               DO 270 J = N, K                                          ASM 0252
                  AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H                ASM 0253
270            CONTINUE                                                 ASM 0254
280         CONTINUE                                                    ASM 0255
            WKD(N+M)  = SCALE * WKD(N+M)                                ASM 0256
            AAD(N,N-1) = SCALE * G                                      ASM 0257
         END IF                                                         ASM 0258
290   CONTINUE                                                          ASM 0259
                                                                        ASM 0260
      DO 340  N = K-2, L, -1                                            ASM 0261
         N1 = N + 1                                                     ASM 0262
         N2 = N + 2                                                     ASM 0263
         F  = AAD(N1,N)                                                 ASM 0264
         IF ( F.NE.ZERO ) THEN                                          ASM 0265
            F  = F * WKD(N1+M)                                          ASM 0266
            DO 300 I = N2, K                                            ASM 0267
               WKD(I+M) = AAD(I,N)                                      ASM 0268
300         CONTINUE                                                    ASM 0269
            IF ( N1.LE.K ) THEN                                         ASM 0270
               DO 330 J = 1, M                                          ASM 0271
                  G = ZERO                                              ASM 0272
                  DO 310 I = N1, K                                      ASM 0273
                     G = G + WKD(I+M) * EVECD(I,J)                      ASM 0274
310               CONTINUE                                              ASM 0275
                  G = G / F                                             ASM 0276
                  DO 320 I = N1, K                                      ASM 0277
                     EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)             ASM 0278
320               CONTINUE                                              ASM 0279
330            CONTINUE                                                 ASM 0280
            END IF                                                      ASM 0281
         END IF                                                         ASM 0282
340   CONTINUE                                                          ASM 0283
                                                                        ASM 0284
350   CONTINUE                                                          ASM 0285
      N = 1                                                             ASM 0286
      DO 370 I = 1, M                                                   ASM 0287
         DO 360 J = N, M                                                ASM 0288
            RNORM = RNORM + DABS(AAD(I,J))                              ASM 0289
360      CONTINUE                                                       ASM 0290
         N = I                                                          ASM 0291
         IF ( I.LT.L .OR. I.GT.K ) EVALD(I) = AAD(I,I)                  ASM 0292
370   CONTINUE                                                          ASM 0293
      N = K                                                             ASM 0294
      T = ZERO                                                          ASM 0295
C                                         ** SEARCH FOR NEXT EIGENVALUESASM 0296
380   IF ( N.LT.L ) GO TO 530                                           ASM 0297
      IN = 0                                                            ASM 0298
      N1 = N - 1                                                        ASM 0299
      N2 = N - 2                                                        ASM 0300
C                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENTASM 0301
390   CONTINUE                                                          ASM 0302
      DO 400 I = L, N                                                   ASM 0303
         LB = N+L - I                                                   ASM 0304
         IF ( LB.EQ.L ) GO TO 410                                       ASM 0305
         S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )                ASM 0306
         IF ( S.EQ.ZERO ) S = RNORM                                     ASM 0307
         IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) GO TO 410                 ASM 0308
400   CONTINUE                                                          ASM 0309
                                                                        ASM 0310
410   X = AAD(N,N)                                                      ASM 0311
      IF ( LB.EQ.N ) THEN                                               ASM 0312
C                                        ** ONE EIGENVALUE FOUND        ASM 0313
         AAD(N,N)  = X + T                                              ASM 0314
         EVALD(N) = AAD(N,N)                                            ASM 0315
         N = N1                                                         ASM 0316
         GO TO 380                                                      ASM 0317
      END IF                                                            ASM 0318
                                                                        ASM 0319
      Y = AAD(N1,N1)                                                    ASM 0320
      W = AAD(N,N1) * AAD(N1,N)                                         ASM 0321
      IF ( LB.EQ.N1 ) THEN                                              ASM 0322
C                                        ** TWO EIGENVALUES FOUND       ASM 0323
         P = (Y-X) * C2                                                 ASM 0324
         Q = P**2 + W                                                   ASM 0325
         Z = DSQRT( DABS(Q) )                                           ASM 0326
         AAD(N,N) = X + T                                               ASM 0327
         X = AAD(N,N)                                                   ASM 0328
         AAD(N1,N1) = Y + T                                             ASM 0329
C                                        ** REAL PAIR                   ASM 0330
         Z = P + DSIGN(Z,P)                                             ASM 0331
         EVALD(N1) = X + Z                                              ASM 0332
         EVALD(N)  = EVALD(N1)                                          ASM 0333
         IF ( Z.NE.ZERO ) EVALD(N) = X - W / Z                          ASM 0334
         X = AAD(N,N1)                                                  ASM 0335
C                                  ** EMPLOY SCALE FACTOR IN CASE       ASM 0336
C                                  ** X AND Z ARE VERY SMALL            ASM 0337
         R = DSQRT( X*X + Z*Z )                                         ASM 0338
         P = X / R                                                      ASM 0339
         Q = Z / R                                                      ASM 0340
C                                             ** ROW MODIFICATION       ASM 0341
         DO 420 J = N1, M                                               ASM 0342
            Z = AAD(N1,J)                                               ASM 0343
            AAD(N1,J) = Q * Z + P * AAD(N,J)                            ASM 0344
            AAD(N,J)  = Q * AAD(N,J) - P * Z                            ASM 0345
420      CONTINUE                                                       ASM 0346
C                                             ** COLUMN MODIFICATION    ASM 0347
         DO 430 I = 1, N                                                ASM 0348
            Z = AAD(I,N1)                                               ASM 0349
            AAD(I,N1) = Q * Z + P * AAD(I,N)                            ASM 0350
            AAD(I,N)  = Q * AAD(I,N) - P * Z                            ASM 0351
430      CONTINUE                                                       ASM 0352
C                                          ** ACCUMULATE TRANSFORMATIONSASM 0353
         DO 440 I = L, K                                                ASM 0354
            Z = EVECD(I,N1)                                             ASM 0355
            EVECD(I,N1) = Q * Z + P * EVECD(I,N)                        ASM 0356
            EVECD(I,N)  = Q * EVECD(I,N) - P * Z                        ASM 0357
440      CONTINUE                                                       ASM 0358
                                                                        ASM 0359
         N = N2                                                         ASM 0360
         GO TO 380                                                      ASM 0361
      END IF                                                            ASM 0362
                                                                        ASM 0363
      IF ( IN.EQ.30 ) THEN                                              ASM 0364
C                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR   ASM 0365
C                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUEASM 0366
         IER = N                                                        ASM 0367
         GO TO 670                                                      ASM 0368
      END IF                                                            ASM 0369
C                                                          ** FORM SHIFTASM 0370
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN                                ASM 0371
         T = T + X                                                      ASM 0372
         DO 450 I = L, N                                                ASM 0373
            AAD(I,I) = AAD(I,I) - X                                     ASM 0374
450      CONTINUE                                                       ASM 0375
         S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))                         ASM 0376
         X = C3 * S                                                     ASM 0377
         Y = X                                                          ASM 0378
         W = - C1 * S**2                                                ASM 0379
      END IF                                                            ASM 0380
                                                                        ASM 0381
      IN = IN + 1                                                       ASM 0382
C                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTSASM 0383
                                                                        ASM 0384
      DO 460 J = LB, N2                                                 ASM 0385
         I = N2+LB - J                                                  ASM 0386
         Z = AAD(I,I)                                                   ASM 0387
         R = X - Z                                                      ASM 0388
         S = Y - Z                                                      ASM 0389
         P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)                    ASM 0390
         Q = AAD(I+1,I+1) - Z - R - S                                   ASM 0391
         R = AAD(I+2,I+1)                                               ASM 0392
         S = DABS(P) + DABS(Q) + DABS(R)                                ASM 0393
         P = P / S                                                      ASM 0394
         Q = Q / S                                                      ASM 0395
         R = R / S                                                      ASM 0396
         IF ( I.EQ.LB ) GO TO 470                                       ASM 0397
         UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )                ASM 0398
         VV=DABS(P)*(DABS(AAD(I-1,I-1))+DABS(Z)+DABS(AAD(I+1,I+1)))     ASM 0399
         IF ( UU .LE. TOL*VV ) GO TO 470                                ASM 0400
460   CONTINUE                                                          ASM 0401
                                                                        ASM 0402
470   CONTINUE                                                          ASM 0403
      AAD(I+2,I) = ZERO                                                 ASM 0404
      DO 480 J = I+3, N                                                 ASM 0405
         AAD(J,J-2) = ZERO                                              ASM 0406
         AAD(J,J-3) = ZERO                                              ASM 0407
480   CONTINUE                                                          ASM 0408
                                                                        ASM 0409
C             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO NASM 0410
                                                                        ASM 0411
      DO 520 KA = I, N1                                                 ASM 0412
         NOTLAS = KA.NE.N1                                              ASM 0413
         IF ( KA.EQ.I ) THEN                                            ASM 0414
            S = DSIGN( DSQRT( P*P + Q*Q + R*R ), P )                    ASM 0415
            IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)                ASM 0416
         ELSE                                                           ASM 0417
            P = AAD(KA,KA-1)                                            ASM 0418
            Q = AAD(KA+1,KA-1)                                          ASM 0419
            R = ZERO                                                    ASM 0420
            IF ( NOTLAS ) R = AAD(KA+2,KA-1)                            ASM 0421
            X = DABS(P) + DABS(Q) + DABS(R)                             ASM 0422
            IF ( X.EQ.ZERO ) GO TO 520                                  ASM 0423
            P = P / X                                                   ASM 0424
            Q = Q / X                                                   ASM 0425
            R = R / X                                                   ASM 0426
            S = DSIGN( DSQRT( P*P + Q*Q + R*R ), P )                    ASM 0427
            AAD(KA,KA-1) = - S * X                                      ASM 0428
         END IF                                                         ASM 0429
         P = P + S                                                      ASM 0430
         X = P / S                                                      ASM 0431
         Y = Q / S                                                      ASM 0432
         Z = R / S                                                      ASM 0433
         Q = Q / P                                                      ASM 0434
         R = R / P                                                      ASM 0435
C                                                    ** ROW MODIFICATIONASM 0436
         DO 490 J = KA, M                                               ASM 0437
            P = AAD(KA,J) + Q * AAD(KA+1,J)                             ASM 0438
            IF ( NOTLAS ) THEN                                          ASM 0439
               P = P + R * AAD(KA+2,J)                                  ASM 0440
               AAD(KA+2,J) = AAD(KA+2,J) - P * Z                        ASM 0441
            END IF                                                      ASM 0442
            AAD(KA+1,J) = AAD(KA+1,J) - P * Y                           ASM 0443
            AAD(KA,J)   = AAD(KA,J)   - P * X                           ASM 0444
490      CONTINUE                                                       ASM 0445
C                                                 ** COLUMN MODIFICATIONASM 0446
         DO 500 II = 1, MIN0(N,KA+3)                                    ASM 0447
            P = X * AAD(II,KA) + Y * AAD(II,KA+1)                       ASM 0448
            IF ( NOTLAS ) THEN                                          ASM 0449
               P = P + Z * AAD(II,KA+2)                                 ASM 0450
               AAD(II,KA+2) = AAD(II,KA+2) - P * R                      ASM 0451
            END IF                                                      ASM 0452
            AAD(II,KA+1) = AAD(II,KA+1) - P * Q                         ASM 0453
            AAD(II,KA)   = AAD(II,KA) - P                               ASM 0454
500      CONTINUE                                                       ASM 0455
C                                          ** ACCUMULATE TRANSFORMATIONSASM 0456
         DO 510 II = L, K                                               ASM 0457
            P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)                   ASM 0458
            IF ( NOTLAS ) THEN                                          ASM 0459
               P = P + Z * EVECD(II,KA+2)                               ASM 0460
               EVECD(II,KA+2) = EVECD(II,KA+2) - P * R                  ASM 0461
            END IF                                                      ASM 0462
            EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q                     ASM 0463
            EVECD(II,KA)   = EVECD(II,KA) - P                           ASM 0464
510      CONTINUE                                                       ASM 0465
                                                                        ASM 0466
520   CONTINUE                                                          ASM 0467
      GO TO 390                                                         ASM 0468
C                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTORASM 0469
530   CONTINUE                                                          ASM 0470
      IF ( RNORM.NE.ZERO ) THEN                                         ASM 0471
         DO 560  N = M, 1, -1                                           ASM 0472
            N2 = N                                                      ASM 0473
            AAD(N,N) = ONE                                              ASM 0474
            DO 550  I = N-1, 1, -1                                      ASM 0475
               W = AAD(I,I) - EVALD(N)                                  ASM 0476
               IF ( W.EQ.ZERO ) W = TOL * RNORM                         ASM 0477
               R = AAD(I,N)                                             ASM 0478
               DO 540 J = N2, N-1                                       ASM 0479
                  R = R + AAD(I,J) * AAD(J,N)                           ASM 0480
540            CONTINUE                                                 ASM 0481
               AAD(I,N) = - R / W                                       ASM 0482
               N2 = I                                                   ASM 0483
550         CONTINUE                                                    ASM 0484
560      CONTINUE                                                       ASM 0485
C                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALSASM 0486
                                                                        ASM 0487
         DO 580 I = 1, M                                                ASM 0488
            IF ( I.LT.L .OR. I.GT.K ) THEN                              ASM 0489
               DO 570 J = I, M                                          ASM 0490
                  EVECD(I,J) = AAD(I,J)                                 ASM 0491
570            CONTINUE                                                 ASM 0492
            END IF                                                      ASM 0493
580      CONTINUE                                                       ASM 0494
C                                   ** MULTIPLY BY TRANSFORMATION MATRIXASM 0495
         IF ( K.NE.0 ) THEN                                             ASM 0496
            DO 600  J = M, L, -1                                        ASM 0497
               DO 600 I = L, K                                          ASM 0498
                  Z = ZERO                                              ASM 0499
                  DO 590 N = L, MIN0(J,K)                               ASM 0500
                     Z = Z + EVECD(I,N) * AAD(N,J)                      ASM 0501
590               CONTINUE                                              ASM 0502
                  EVECD(I,J) = Z                                        ASM 0503
600         CONTINUE                                                    ASM 0504
         END IF                                                         ASM 0505
                                                                        ASM 0506
      END IF                                                            ASM 0507
                                                                        ASM 0508
      DO 620 I = L, K                                                   ASM 0509
         DO 620 J = 1, M                                                ASM 0510
            EVECD(I,J) = EVECD(I,J) * WKD(I)                            ASM 0511
620   CONTINUE                                                          ASM 0512
C                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURREDASM 0513
      DO 640  I = L-1, 1, -1                                            ASM 0514
         J = WKD(I)                                                     ASM 0515
         IF ( I.NE.J ) THEN                                             ASM 0516
            DO 630 N = 1, M                                             ASM 0517
               REPL       = EVECD(I,N)                                  ASM 0518
               EVECD(I,N) = EVECD(J,N)                                  ASM 0519
               EVECD(J,N) = REPL                                        ASM 0520
630         CONTINUE                                                    ASM 0521
         END IF                                                         ASM 0522
640   CONTINUE                                                          ASM 0523
                                                                        ASM 0524
      DO 660 I = K+1, M                                                 ASM 0525
         J = WKD(I)                                                     ASM 0526
         IF ( I.NE.J ) THEN                                             ASM 0527
            DO 650 N = 1, M                                             ASM 0528
               REPL       = EVECD(I,N)                                  ASM 0529
               EVECD(I,N) = EVECD(J,N)                                  ASM 0530
               EVECD(J,N) = REPL                                        ASM 0531
650         CONTINUE                                                    ASM 0532
         END IF                                                         ASM 0533
660   CONTINUE                                                          ASM 0534
C                         ** PUT RESULTS INTO OUTPUT ARRAYS             ASM 0535
  670 CONTINUE                                                          ASM 0536
      DO 680 J = 1, M                                                   ASM 0537
         EVAL( J ) = EVALD(J)                                           ASM 0538
         DO 680 K = 1, M                                                ASM 0539
            EVEC( J,K ) = EVECD(J,K)                                    ASM 0540
680   CONTINUE                                                          ASM 0541
                                                                        ASM 0542
      RETURN                                                            ASM 0543
      END                                                               ASM 0544
