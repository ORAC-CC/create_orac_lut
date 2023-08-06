      SUBROUTINE  ALTRIN( GU, KK, LL, MXCMU, MXUMU, MAXUMU, NLYR,       ALT 0001
     $                    NN, NSTR, NUMU, TAUCPR, UMU, U0U, WK )        ALT 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            ALT 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                ALT 0004
                                                                        ALT 0005
C       COMPUTES AZIMUTHALLY-AVERAGED INTENSITY AT TOP AND BOTTOM       ALT 0006
C       OF MEDIUM (RELATED TO ALBEDO AND TRANSMISSION OF MEDIUM BY      ALT 0007
C       RECIPROCITY PRINCIPLES;  SEE REF S2).  USER POLAR ANGLES ARE    ALT 0008
C       USED AS INCIDENT BEAM ANGLES. (THIS IS A VERY SPECIALIZED       ALT 0009
C       VERSION OF 'USRINT')                                            ALT 0010
                                                                        ALT 0011
C       ** NOTE **  USER INPUT VALUES OF -UMU- (ASSUMED POSITIVE) ARE   ALT 0012
C                   TEMPORARILY IN UPPER LOCATIONS OF  -UMU-  AND       ALT 0013
C                   CORRESPONDING NEGATIVES ARE IN LOWER LOCATIONS      ALT 0014
C                   (THIS MAKES -GU- COME OUT RIGHT).  I.E. THE CONTENTSALT 0015
C                   OF THE TEMPORARY -UMU- ARRAY ARE:                   ALT 0016
                                                                        ALT 0017
C                     -UMU(NUMU),..., -UMU(1), UMU(1),..., UMU(NUMU)    ALT 0018
                                                                        ALT 0019
C   I N P U T    V A R I A B L E S:                                     ALT 0020
                                                                        ALT 0021
C       GU     :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES        ALT 0022
C                 (I.E., G IN EQ. SC(1) )                               ALT 0023
C       KK     :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)             ALT 0024
C       LL     :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED       ALT 0025
C                 BY SOLVING SCALED VERSION OF EQ. SC(5);               ALT 0026
C                 EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED           ALT 0027
C       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)             ALT 0028
C       TAUCPR :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)             ALT 0029
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        ALT 0030
                                                                        ALT 0031
C   O U T P U T    V A R I A B L E:                                     ALT 0032
                                                                        ALT 0033
C       U0U  :    DIFFUSE AZIMUTHALLY-AVERAGED INTENSITY AT TOP AND     ALT 0034
C                 BOTTOM OF MEDIUM (DIRECTLY TRANSMITTED COMPONENT,     ALT 0035
C                 CORRESPONDING TO -BNDINT- IN 'USRINT', IS OMITTED).   ALT 0036
                                                                        ALT 0037
C   I N T E R N A L    V A R I A B L E S:                               ALT 0038
                                                                        ALT 0039
C       DTAU   :  OPTICAL DEPTH OF A COMPUTATIONAL LAYER                ALT 0040
C       PALINT :  NON-BOUNDARY-FORCED INTENSITY COMPONENT               ALT 0041
C       UTAUPR :  OPTICAL DEPTHS OF USER OUTPUT LEVELS (DELTA-M SCALED) ALT 0042
C       WK     :  SCRATCH VECTOR FOR SAVING 'EXP' EVALUATIONS           ALT 0043
C       ALL THE EXPONENTIAL FACTORS (I.E., EXP1, EXPN,... ETC.)         ALT 0044
C       COME FROM THE SUBSTITUTION OF CONSTANTS OF INTEGRATION IN       ALT 0045
C       EQ. SC(12) INTO EQS. S1(8-9).  ALL HAVE NEGATIVE ARGUMENTS.     ALT 0046
C+---------------------------------------------------------------------+ALT 0047
                                                                        ALT 0048
      REAL*8     UTAUPR( 2 )                                            ALT 0049
      REAL*8     GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ), MU, ALT 0050
     $         TAUCPR( 0:* ), UMU(*), U0U( MAXUMU,* ), WK(*)            ALT 0051
                                                                        ALT 0052
                                                                        ALT 0053
      UTAUPR(1) = 0.0                                                   ALT 0054
      UTAUPR(2) = TAUCPR( NLYR )                                        ALT 0055
      DO 100  LU = 1, 2                                                 ALT 0056
         IF ( LU.EQ.1 )  THEN                                           ALT 0057
            IUMIN = NUMU / 2 + 1                                        ALT 0058
            IUMAX = NUMU                                                ALT 0059
            SGN = 1.0                                                   ALT 0060
         ELSE                                                           ALT 0061
            IUMIN = 1                                                   ALT 0062
            IUMAX = NUMU / 2                                            ALT 0063
            SGN = - 1.0                                                 ALT 0064
         END IF                                                         ALT 0065
C                                   ** LOOP OVER POLAR ANGLES AT WHICH  ALT 0066
C                                   ** ALBEDOS/TRANSMISSIVITIES DESIRED ALT 0067
C                                   ** ( UPWARD ANGLES AT TOP BOUNDARY, ALT 0068
C                                   ** DOWNWARD ANGLES AT BOTTOM )      ALT 0069
         DO 50  IU = IUMIN, IUMAX                                       ALT 0070
            MU = UMU(IU)                                                ALT 0071
C                                     ** INTEGRATE FROM TOP TO BOTTOM   ALT 0072
C                                     ** COMPUTATIONAL LAYER            ALT 0073
            PALINT = 0.0                                                ALT 0074
            DO 30  LC = 1, NLYR                                         ALT 0075
                                                                        ALT 0076
               DTAU = TAUCPR(LC) - TAUCPR(LC-1)                         ALT 0077
               EXP1 =  DEXP( (UTAUPR(LU) - TAUCPR(LC-1)) / MU )         ALT 0078
               EXP2 =  DEXP( (UTAUPR(LU) - TAUCPR( LC )) / MU )         ALT 0079
                                                                        ALT 0080
C                                      ** -KK- IS NEGATIVE              ALT 0081
               DO 20  IQ = 1, NN                                        ALT 0082
                  WK(IQ) = DEXP( KK(IQ,LC) * DTAU )                     ALT 0083
                  DENOM = 1.0 + MU * KK(IQ,LC)                          ALT 0084
                  IF ( DABS(DENOM).LT.0.0001 ) THEN                     ALT 0085
C                                                   ** L'HOSPITAL LIMIT ALT 0086
                     EXPN = DTAU / MU * EXP2                            ALT 0087
                  ELSE                                                  ALT 0088
                     EXPN = ( EXP1 * WK(IQ) - EXP2 ) * SGN / DENOM      ALT 0089
                  END IF                                                ALT 0090
                  PALINT = PALINT + GU(IU,IQ,LC) * LL(IQ,LC) * EXPN     ALT 0091
20             CONTINUE                                                 ALT 0092
C                                      ** -KK- IS POSITIVE              ALT 0093
               DO 21  IQ = NN+1, NSTR                                   ALT 0094
                  DENOM = 1.0 + MU * KK(IQ,LC)                          ALT 0095
                  IF ( DABS(DENOM).LT.0.0001 ) THEN                     ALT 0096
                     EXPN = - DTAU / MU * EXP1                          ALT 0097
                  ELSE                                                  ALT 0098
                     EXPN = ( EXP1 - EXP2 * WK(NSTR+1-IQ) ) *SGN / DENOMALT 0099
                  END IF                                                ALT 0100
                  PALINT = PALINT + GU(IU,IQ,LC) * LL(IQ,LC) * EXPN     ALT 0101
21             CONTINUE                                                 ALT 0102
                                                                        ALT 0103
30          CONTINUE                                                    ALT 0104
                                                                        ALT 0105
            U0U( IU, LU ) = PALINT                                      ALT 0106
                                                                        ALT 0107
 50      CONTINUE                                                       ALT 0108
100   CONTINUE                                                          ALT 0109
                                                                        ALT 0110
      RETURN                                                            ALT 0111
      END                                                               ALT 0112
