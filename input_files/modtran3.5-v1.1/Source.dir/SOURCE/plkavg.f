      REAL*8 FUNCTION  PLKAVG ( WNUMLO, WNUMHI, T )                     PLK 0001
                                                                        PLK 0002
C        COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS    PLK 0003
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            PLK 0004
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                PLK 0005
                                                                        PLK 0006
C  NOTE ** CHANGE 'R1MACH' TO 'D1MACH' TO RUN IN DOUBLE PRECISION       PLK 0007
                                                                        PLK 0008
C  I N P U T :  WNUMLO : LOWER WAVENUMBER ( INV CM ) OF SPECTRAL        PLK 0009
C                           INTERVAL                                    PLK 0010
C               WNUMHI : UPPER WAVENUMBER                               PLK 0011
C               T      : TEMPERATURE (K)                                PLK 0012
                                                                        PLK 0013
C  O U T P U T :  PLKAVG : INTEGRATED PLANCK FUNCTION ( WATTS/SQ M )    PLK 0014
C                           = INTEGRAL (WNUMLO TO WNUMHI) OF            PLK 0015
C                              2H C**2  NU**3 / ( EXP(HC NU/KT) - 1)    PLK 0016
C                              (WHERE H=PLANCKS CONSTANT, C=SPEED OF    PLK 0017
C                              LIGHT, NU=WAVENUMBER, T=TEMPERATURE,     PLK 0018
C                              AND K = BOLTZMANN CONSTANT)              PLK 0019
                                                                        PLK 0020
C  REFERENCE : SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE          PLK 0021
C                 OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.,      PLK 0022
C                 JAN. 1974                                             PLK 0023
                                                                        PLK 0024
C  METHOD :  FOR  -WNUMLO-  CLOSE TO  -WNUMHI-, A SIMPSON-RULE          PLK 0025
C            QUADRATURE IS DONE TO AVOID ILL-CONDITIONING; OTHERWISE    PLK 0026
                                                                        PLK 0027
C            (1)  FOR 'WNUMLO' OR 'WNUMHI' SMALL,                       PLK 0028
C                 INTEGRAL(0 TO WNUMLO/HI) IS CALCULATED BY EXPANDING   PLK 0029
C                 THE INTEGRAND IN A POWER SERIES AND INTEGRATING       PLK 0030
C                 TERM BY TERM;                                         PLK 0031
                                                                        PLK 0032
C            (2)  OTHERWISE, INTEGRAL(WNUMLO/HI TO INFINITY) IS         PLK 0033
C                 CALCULATED BY EXPANDING THE DENOMINATOR OF THE        PLK 0034
C                 INTEGRAND IN POWERS OF THE EXPONENTIAL AND            PLK 0035
C                 INTEGRATING TERM BY TERM.                             PLK 0036
                                                                        PLK 0037
C  ACCURACY :  AT LEAST 6 SIGNIFICANT DIGITS, ASSUMING THE              PLK 0038
C              PHYSICAL CONSTANTS ARE INFINITELY ACCURATE               PLK 0039
                                                                        PLK 0040
C  ERRORS WHICH ARE NOT TRAPPED:                                        PLK 0041
                                                                        PLK 0042
C      * POWER OR EXPONENTIAL SERIES MAY UNDERFLOW, GIVING NO           PLK 0043
C        SIGNIFICANT DIGITS.  THIS MAY OR MAY NOT BE OF CONCERN,        PLK 0044
C        DEPENDING ON THE APPLICATION.                                  PLK 0045
                                                                        PLK 0046
C      * SIMPSON-RULE SPECIAL CASE IS SKIPPED WHEN DENOMINATOR OF       PLK 0047
C        INTEGRAND WILL CAUSE OVERFLOW.  IN THAT CASE THE NORMAL        PLK 0048
C        PROCEDURE IS USED, WHICH MAY BE INACCURATE IF THE              PLK 0049
C        WAVENUMBER LIMITS (WNUMLO, WNUMHI) ARE CLOSE TOGETHER.         PLK 0050
C ----------------------------------------------------------------------PLK 0051
C                                   *** ARGUMENTS                       PLK 0052
      REAL*8     T, WNUMLO, WNUMHI                                      PLK 0053
C                                   *** LOCAL VARIABLES                 PLK 0054
                                                                        PLK 0055
C        A1,2,... :  POWER SERIES COEFFICIENTS                          PLK 0056
C        C2       :  H * C / K, IN UNITS CM*K (H = PLANCKS CONSTANT,    PLK 0057
C                      C = SPEED OF LIGHT, K = BOLTZMANN CONSTANT)      PLK 0058
C        D(I)     :  EXPONENTIAL SERIES EXPANSION OF INTEGRAL OF        PLK 0059
C                       PLANCK FUNCTION FROM WNUMLO (I=1) OR WNUMHI     PLK 0060
C                       (I=2) TO INFINITY                               PLK 0061
C        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 ON        PLK 0062
C                       COMPUTER                                        PLK 0063
C        EX       :  EXP( - V(I) )                                      PLK 0064
C        EXM      :  EX**M                                              PLK 0065
C        MMAX     :  NO. OF TERMS TO TAKE IN EXPONENTIAL SERIES         PLK 0066
C        MV       :  MULTIPLES OF 'V(I)'                                PLK 0067
C        P(I)     :  POWER SERIES EXPANSION OF INTEGRAL OF              PLK 0068
C                       PLANCK FUNCTION FROM ZERO TO WNUMLO (I=1) OR    PLK 0069
C                       WNUMHI (I=2)                                    PLK 0070
C        PI       :  3.14159...                                         PLK 0071
C        SIGMA    :  STEFAN-BOLTZMANN CONSTANT (W/M**2/K**4)            PLK 0072
C        SIGDPI   :  SIGMA / PI                                         PLK 0073
C        SMALLV   :  NUMBER OF TIMES THE POWER SERIES IS USED (0,1,2)   PLK 0074
C        V(I)     :  C2 * (WNUMLO(I=1) OR WNUMHI(I=2)) / TEMPERATURE    PLK 0075
C        VCUT     :  POWER-SERIES CUTOFF POINT                          PLK 0076
C        VCP      :  EXPONENTIAL SERIES CUTOFF POINTS                   PLK 0077
C        VMAX     :  LARGEST ALLOWABLE ARGUMENT OF 'EXP' FUNCTION       PLK 0078
                                                                        PLK 0079
      PARAMETER  ( A1 = 1./3., A2 = -1./8., A3 = 1./60., A4 = -1./5040.,PLK 0080
     $             A5 = 1./272160., A6 = -1./13305600. )                PLK 0081
      INTEGER  SMALLV                                                   PLK 0082
      REAL*8     C2, CONC, D(2), EPSIL, EX, MV, P(2), SIGMA, SIGDPI,    PLK 0083
     $         V(2), VCUT, VCP(7), VSQ                                  PLK 0084
      DOUBLE PRECISION   D1MACH                                         PLK 0085
      SAVE     PI, CONC, VMAX, EPSIL, SIGDPI                            PLK 0086
      DATA     C2 / 1.438786 /,  SIGMA / 5.67032E-8 /,                  PLK 0087
     $         VCUT / 1.5 /, VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /PLK 0088
      DATA     PI / 0.0 /                                               PLK 0089
      F(X) = X**3 / ( DEXP(X) - 1 )                                     PLK 0090
                                                                        PLK 0091
                                                                        PLK 0092
      IF ( PI.EQ.0.0 )  THEN                                            PLK 0093
         PI = 2. * DASIN( 1.0D0 )                                       PLK 0094
         VMAX = DLOG( D1MACH(2) )                                       PLK 0095
         EPSIL = D1MACH(3)                                              PLK 0096
         SIGDPI = SIGMA / PI                                            PLK 0097
         CONC = 15. / PI**4                                             PLK 0098
      END IF                                                            PLK 0099
                                                                        PLK 0100
      IF( T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. )            PLK 0101
     $    CALL ERRMSG( 'PLKAVG--TEMPERATURE OR WAVENUMS. WRONG', .TRUE.)PLK 0102
                                                                        PLK 0103
      IF ( T.LT.1.E-4 )  THEN                                           PLK 0104
         PLKAVG = 0.0                                                   PLK 0105
         RETURN                                                         PLK 0106
      ENDIF                                                             PLK 0107
                                                                        PLK 0108
      V(1) = C2 * WNUMLO / T                                            PLK 0109
      V(2) = C2 * WNUMHI / T                                            PLK 0110
      IF ( V(1).GT.EPSIL .AND. V(2).LT.VMAX .AND.                       PLK 0111
     $     (WNUMHI-WNUMLO)/WNUMHI .LT. 1.E-2 )  THEN                    PLK 0112
                                                                        PLK 0113
C                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL PLK 0114
C                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE. PLK 0115
         HH = V(2) - V(1)                                               PLK 0116
         OLDVAL = 0.0                                                   PLK 0117
         VAL0 = F( V(1) ) + F( V(2) )                                   PLK 0118
                                                                        PLK 0119
         DO  2  N = 1, 10                                               PLK 0120
            DEL = HH / (2*N)                                            PLK 0121
            VAL = VAL0                                                  PLK 0122
            DO  1  K = 1, 2*N-1                                         PLK 0123
               VAL = VAL + 2*(1+MOD(K,2)) * F( V(1) + K*DEL )           PLK 0124
    1       CONTINUE                                                    PLK 0125
            VAL = DEL/3. * VAL                                          PLK 0126
            IF ( DABS( (VAL-OLDVAL)/VAL ) .LE. 1.E-6 )  GO TO 3         PLK 0127
            OLDVAL = VAL                                                PLK 0128
    2    CONTINUE                                                       PLK 0129
         CALL ERRMSG( 'PLKAVG--SIMPSON RULE DIDNT CONVERGE', .FALSE. )  PLK 0130
                                                                        PLK 0131
    3    PLKAVG = SIGDPI * T**4 * CONC * VAL                            PLK 0132
         RETURN                                                         PLK 0133
      END IF                                                            PLK 0134
                                                                        PLK 0135
      SMALLV = 0                                                        PLK 0136
      DO  50  I = 1, 2                                                  PLK 0137
                                                                        PLK 0138
         IF( V(I).LT.VCUT )  THEN                                       PLK 0139
C                                   ** USE POWER SERIES                 PLK 0140
            SMALLV = SMALLV + 1                                         PLK 0141
            VSQ = V(I)**2                                               PLK 0142
            P(I) =  CONC * VSQ * V(I) * ( A1 + V(I) * ( A2 + V(I) *     PLK 0143
     $                ( A3 + VSQ * ( A4 + VSQ * ( A5 + VSQ*A6 ) ) ) ) ) PLK 0144
         ELSE                                                           PLK 0145
C                    ** USE EXPONENTIAL SERIES                          PLK 0146
            MMAX = 0                                                    PLK 0147
C                                ** FIND UPPER LIMIT OF SERIES          PLK 0148
   20       MMAX = MMAX + 1                                             PLK 0149
               IF ( V(I).LT.VCP( MMAX ) )  GO TO 20                     PLK 0150
                                                                        PLK 0151
            EX = DEXP( - V(I) )                                         PLK 0152
            EXM = 1.0                                                   PLK 0153
            D(I) = 0.0                                                  PLK 0154
                                                                        PLK 0155
            DO  30  M = 1, MMAX                                         PLK 0156
               MV = M * V(I)                                            PLK 0157
               EXM = EX * EXM                                           PLK 0158
               D(I) = D(I) +                                            PLK 0159
     $                EXM * ( 6. + MV*( 6. + MV*( 3. + MV ) ) ) / M**4  PLK 0160
   30       CONTINUE                                                    PLK 0161
                                                                        PLK 0162
            D(I) = CONC * D(I)                                          PLK 0163
         END IF                                                         PLK 0164
                                                                        PLK 0165
   50 CONTINUE                                                          PLK 0166
                                                                        PLK 0167
      IF ( SMALLV .EQ. 2 ) THEN                                         PLK 0168
C                                    ** WNUMLO AND WNUMHI BOTH SMALL    PLK 0169
         PLKAVG = P(2) - P(1)                                           PLK 0170
                                                                        PLK 0171
      ELSE IF ( SMALLV .EQ. 1 ) THEN                                    PLK 0172
C                                    ** WNUMLO SMALL, WNUMHI LARGE      PLK 0173
         PLKAVG = 1. - P(1) - D(2)                                      PLK 0174
                                                                        PLK 0175
      ELSE                                                              PLK 0176
C                                    ** WNUMLO AND WNUMHI BOTH LARGE    PLK 0177
         PLKAVG = D(1) - D(2)                                           PLK 0178
                                                                        PLK 0179
      END IF                                                            PLK 0180
                                                                        PLK 0181
      PLKAVG = SIGDPI * T**4 * PLKAVG                                   PLK 0182
      IF( PLKAVG.EQ.0.0 )                                               PLK 0183
     $    CALL ERRMSG( 'PLKAVG--RETURNS ZERO; POSSIBLE UNDERFLOW',      PLK 0184
     $                 .FALSE. )                                        PLK 0185
                                                                        PLK 0186
      RETURN                                                            PLK 0187
      END                                                               PLK 0188
