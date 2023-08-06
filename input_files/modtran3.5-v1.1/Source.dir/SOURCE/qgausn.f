      SUBROUTINE  QGAUSN( M, GMU, GWT )                                 QGS 0001
                                                                        QGS 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            QGS 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                QGS 0004
C       COMPUTE WEIGHTS AND ABSCISSAE FOR ORDINARY GAUSSIAN QUADRATURE  QGS 0005
C       (NO WEIGHT FUNCTION INSIDE INTEGRAL) ON THE INTERVAL (0,1)      QGS 0006
                                                                        QGS 0007
C   INPUT :    M                     ORDER OF QUADRATURE RULE           QGS 0008
                                                                        QGS 0009
C   OUTPUT :  GMU(I)  I = 1 TO M,    ARRAY OF ABSCISSAE                 QGS 0010
C             GWT(I)  I = 1 TO M,    ARRAY OF WEIGHTS                   QGS 0011
                                                                        QGS 0012
C   REFERENCE:  DAVIS, P.J. AND P. RABINOWITZ, METHODS OF NUMERICAL     QGS 0013
C                   INTEGRATION, ACADEMIC PRESS, NEW YORK, PP. 87, 1975.QGS 0014
                                                                        QGS 0015
C   METHOD:  COMPUTE THE ABSCISSAE AS ROOTS OF THE LEGENDRE             QGS 0016
C            POLYNOMIAL P-SUB-M USING A CUBICALLY CONVERGENT            QGS 0017
C            REFINEMENT OF NEWTON'S METHOD.  COMPUTE THE                QGS 0018
C            WEIGHTS FROM EQ. 2.7.3.8 OF DAVIS/RABINOWITZ.  NOTE        QGS 0019
C            THAT NEWTON'S METHOD CAN VERY EASILY DIVERGE; ONLY A       QGS 0020
C            VERY GOOD INITIAL GUESS CAN GUARANTEE CONVERGENCE.         QGS 0021
C            THE INITIAL GUESS USED HERE HAS NEVER LED TO DIVERGENCE    QGS 0022
C            EVEN FOR M UP TO 1000.                                     QGS 0023
                                                                        QGS 0024
C   ACCURACY:  AT LEAST 13 SIGNIFICANT DIGITS                           QGS 0025
                                                                        QGS 0026
C   INTERNAL VARIABLES:                                                 QGS 0027
                                                                        QGS 0028
C    ITER      : NUMBER OF NEWTON METHOD ITERATIONS                     QGS 0029
C    MAXIT     : MAXIMUM ALLOWED ITERATIONS OF NEWTON METHOD            QGS 0030
C    PM2,PM1,P : 3 SUCCESSIVE LEGENDRE POLYNOMIALS                      QGS 0031
C    PPR       : DERIVATIVE OF LEGENDRE POLYNOMIAL                      QGS 0032
C    P2PRI     : 2ND DERIVATIVE OF LEGENDRE POLYNOMIAL                  QGS 0033
C    TOL       : CONVERGENCE CRITERION FOR LEGENDRE POLY ROOT ITERATION QGS 0034
C    X,XI      : SUCCESSIVE ITERATES IN CUBICALLY-CONVERGENT VERSION    QGS 0035
C                OF NEWTONS METHOD (SEEKING ROOTS OF LEGENDRE POLY.)    QGS 0036
C+---------------------------------------------------------------------+QGS 0037
      REAL*8     CONA, GMU(*), GWT(*), PI, T                            QGS 0038
      INTEGER  ITER, LIM, M, MAXIT, NP1                                 QGS 0039
      DOUBLE   PRECISION  D1MACH                                        QGS 0040
      DOUBLE   PRECISION  EN, NNP1, ONE, P, PM1, PM2, PPR, P2PRI, PROD, QGS 0041
     $                    TMP, TOL, TWO, X, XI                          QGS 0042
      SAVE     PI, TOL                                                  QGS 0043
      DATA     PI / 0.0 /,  MAXIT / 1000 /,  ONE / 1.D0 /,  TWO / 2.D0 /QGS 0044
                                                                        QGS 0045
                                                                        QGS 0046
      IF ( PI.EQ.0.0 )  THEN                                            QGS 0047
         PI = 2. * DASIN(1.0D0)                                         QGS 0048
         TOL = 10. * D1MACH(3)                                          QGS 0049
      END IF                                                            QGS 0050
                                                                        QGS 0051
      IF ( M.LT.1 )  CALL ERRMSG( 'QGAUSN--Bad value of M', .TRUE. )    QGS 0052
      IF ( M.EQ.1 )  THEN                                               QGS 0053
         GMU( 1 ) = 0.5                                                 QGS 0054
         GWT( 1 ) = 1.0                                                 QGS 0055
         RETURN                                                         QGS 0056
      END IF                                                            QGS 0057
                                                                        QGS 0058
      EN   = M                                                          QGS 0059
      NP1  = M + 1                                                      QGS 0060
      NNP1 = M * NP1                                                    QGS 0061
      CONA = DBLE( M-1 ) / ( 8 * M**3 )                                 QGS 0062
                                                                        QGS 0063
      LIM  = M / 2                                                      QGS 0064
      DO 30  K = 1, LIM                                                 QGS 0065
C                                        ** INITIAL GUESS FOR K-TH ROOT QGS 0066
C                                        ** OF LEGENDRE POLYNOMIAL, FROMQGS 0067
C                                        ** DAVIS/RABINOWITZ (2.7.3.3A) QGS 0068
         T = ( 4*K - 1 ) * PI / ( 4*M + 2 )                             QGS 0069
         X = DCOS ( T + CONA / DTAN( T ) )                              QGS 0070
         ITER = 0                                                       QGS 0071
C                                        ** UPWARD RECURRENCE FOR       QGS 0072
C                                        ** LEGENDRE POLYNOMIALS        QGS 0073
   10    ITER = ITER + 1                                                QGS 0074
         PM2 = ONE                                                      QGS 0075
         PM1 = X                                                        QGS 0076
         DO 20 NN = 2, M                                                QGS 0077
            P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN      QGS 0078
            PM2 = PM1                                                   QGS 0079
            PM1 = P                                                     QGS 0080
   20    CONTINUE                                                       QGS 0081
C                                              ** NEWTON METHOD         QGS 0082
         TMP   = ONE / ( ONE - X**2 )                                   QGS 0083
         PPR   = EN * ( PM2 - X * P ) * TMP                             QGS 0084
         P2PRI = ( TWO * X * PPR - NNP1 * P ) * TMP                     QGS 0085
         XI    = X - ( P / PPR ) * ( ONE +                              QGS 0086
     $               ( P / PPR ) * P2PRI / ( TWO * PPR ) )              QGS 0087
                                                                        QGS 0088
C                                              ** CHECK FOR CONVERGENCE QGS 0089
         IF ( DABS(XI-X) .GT. TOL ) THEN                                QGS 0090
            IF( ITER.GT.MAXIT )                                         QGS 0091
     $          CALL ERRMSG( 'QGAUSN--MAX ITERATION COUNT', .TRUE. )    QGS 0092
            X = XI                                                      QGS 0093
            GO TO 10                                                    QGS 0094
         END IF                                                         QGS 0095
C                             ** ITERATION FINISHED--CALCULATE WEIGHTS, QGS 0096
C                             ** ABSCISSAE FOR (-1,1)                   QGS 0097
         GMU( K ) = - X                                                 QGS 0098
         GWT( K ) = TWO / ( TMP * ( EN * PM2 )**2 )                     QGS 0099
         GMU( NP1 - K ) = - GMU( K )                                    QGS 0100
         GWT( NP1 - K ) =   GWT( K )                                    QGS 0101
  30  CONTINUE                                                          QGS 0102
C                                    ** SET MIDDLE ABSCISSA AND WEIGHT  QGS 0103
C                                    ** FOR RULES OF ODD ORDER          QGS 0104
      IF ( MOD( M,2 ) .NE. 0 )  THEN                                    QGS 0105
         GMU( LIM + 1 ) = 0.0                                           QGS 0106
         PROD = ONE                                                     QGS 0107
         DO 40 K = 3, M, 2                                              QGS 0108
            PROD = PROD * K / ( K-1 )                                   QGS 0109
  40     CONTINUE                                                       QGS 0110
         GWT( LIM + 1 ) = TWO / PROD**2                                 QGS 0111
      END IF                                                            QGS 0112
C                                        ** CONVERT FROM (-1,1) TO (0,1)QGS 0113
      DO 50  K = 1, M                                                   QGS 0114
         GMU( K ) = 0.5 * GMU( K ) + 0.5                                QGS 0115
         GWT( K ) = 0.5 * GWT( K )                                      QGS 0116
  50  CONTINUE                                                          QGS 0117
                                                                        QGS 0118
      RETURN                                                            QGS 0119
      END                                                               QGS 0120
