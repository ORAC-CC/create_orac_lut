      SUBROUTINE  SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZIM,     EIG 0001
     $                    MXCMU, NN, NSTR, WK, YLMC, CC, EVECC, EVAL,   EIG 0002
     $                    KK, GC, AAD, WKD, EVECCD, EVALD )             EIG 0003
                                                                        EIG 0004
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            EIG 0005
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                EIG 0006
C         SOLVES EIGENVALUE/VECTOR PROBLEM NECESSARY TO CONSTRUCT       EIG 0007
C         HOMOGENEOUS PART OF DISCRETE ORDINATE SOLUTION; STWJ(8B)      EIG 0008
C         ** NOTE ** EIGENVALUE PROBLEM IS DEGENERATE WHEN SINGLE       EIG 0009
C                    SCATTERING ALBEDO = 1;  PRESENT WAY OF DOING IT    EIG 0010
C                    SEEMS NUMERICALLY MORE STABLE THAN ALTERNATIVE     EIG 0011
C                    METHODS THAT WE TRIED                              EIG 0012
                                                                        EIG 0013
C     ROUTINES CALLED:  ASYMTX                                          EIG 0014
                                                                        EIG 0015
C   I N P U T     V A R I A B L E S:                                    EIG 0016
                                                                        EIG 0017
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTIONEIG 0018
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO) EIG 0019
C       CMU    :  COMPUTATIONAL POLAR ANGLE COSINES                     EIG 0020
C       CWT    :  WEIGHTS FOR QUADRATURE OVER POLAR ANGLE COSINE        EIG 0021
C       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT                          EIG 0022
C       NN     :  HALF THE TOTAL NUMBER OF STREAMS                      EIG 0023
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             EIG 0024
C                 AT THE QUADRATURE ANGLES -CMU-                        EIG 0025
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        EIG 0026
                                                                        EIG 0027
C   O U T P U T    V A R I A B L E S:                                   EIG 0028
                                                                        EIG 0029
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5); NEEDED IN SS(15&18)    EIG 0030
C       EVAL   :  -NN- EIGENVALUES OF EQ. SS(12) ON RETURN FROM 'ASYMTX'EIG 0031
C                    BUT THEN SQUARE ROOTS TAKEN                        EIG 0032
C       EVECC  :  -NN- EIGENVECTORS  (G+) - (G-)  ON RETURN             EIG 0033
C                    FROM 'ASYMTX' ( COLUMN J CORRESPONDS TO -EVAL(J)- )EIG 0034
C                    BUT THEN  (G+) + (G-)  IS CALCULATED FROM SS(10),  EIG 0035
C                    G+  AND  G-  ARE SEPARATED, AND  G+  IS STACKED ON EIG 0036
C                    TOP OF  G-  TO FORM -NSTR- EIGENVECTORS OF SS(7)   EIG 0037
C       GC     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVECTORS, BUT    EIG 0038
C                    IN AN ORDER CORRESPONDING TO -KK-                  EIG 0039
C       KK     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVALUES OF SS(7),EIG 0040
C                    BUT RE-ORDERED WITH NEGATIVE VALUES FIRST ( SQUARE EIG 0041
C                    ROOTS OF -EVAL- TAKEN AND NEGATIVES ADDED )        EIG 0042
                                                                        EIG 0043
C   I N T E R N A L   V A R I A B L E S:                                EIG 0044
                                                                        EIG 0045
C       AMB,APB :  MATRICES (ALPHA-BETA), (ALPHA+BETA) IN REDUCED       EIG 0046
C                    EIGENVALUE PROBLEM                                 EIG 0047
C       ARRAY   :  COMPLETE COEFFICIENT MATRIX OF REDUCED EIGENVALUE    EIG 0048
C                    PROBLEM: (ALFA+BETA)*(ALFA-BETA)                   EIG 0049
C       GPPLGM  :  (G+) + (G-) (CF. EQS. SS(10-11))                     EIG 0050
C       GPMIGM  :  (G+) - (G-) (CF. EQS. SS(10-11))                     EIG 0051
C       WK      :  SCRATCH ARRAY REQUIRED BY 'ASYMTX'                   EIG 0052
C+---------------------------------------------------------------------+EIG 0053
      REAL*8    AMB( MI,* ), APB( MI,* ), ARRAY( MI,* ), CC( MXCMU,* ), EIG 0054
     $        CMU(*), CWT(*), EVAL(*), EVECC( MXCMU,* ), GC( MXCMU,* ), EIG 0055
     $        GL(0:*), KK(*), WK(*), YLMC( 0:MXCMU,* )                  EIG 0056
      DOUBLE PRECISION   EVECCD( MI,* ), EVALD(*), WKD(*), AAD( MI,* )  EIG 0057
                                                                        EIG 0058
                                                                        EIG 0059
C                             ** CALCULATE QUANTITIES IN EQS. SS(5-6)   EIG 0060
      DO 40 IQ  = 1, NN                                                 EIG 0061
                                                                        EIG 0062
         DO 20  JQ = 1, NSTR                                            EIG 0063
            SUM = 0.0                                                   EIG 0064
            DO 10  L = MAZIM, NSTR-1                                    EIG 0065
               SUM = SUM + GL(L) * YLMC(L,IQ) * YLMC(L,JQ)              EIG 0066
10          CONTINUE                                                    EIG 0067
            CC(IQ,JQ) = 0.5 * SUM * CWT(JQ)                             EIG 0068
20       CONTINUE                                                       EIG 0069
                                                                        EIG 0070
         DO 30  JQ = 1, NN                                              EIG 0071
C                             ** FILL REMAINDER OF ARRAY USING SYMMETRY EIG 0072
C                             ** RELATIONS  C(-MUI,MUJ) = C(MUI,-MUJ)   EIG 0073
C                             ** AND        C(-MUI,-MUJ) = C(MUI,MUJ)   EIG 0074
                                                                        EIG 0075
            CC(IQ+NN,JQ) = CC(IQ,JQ+NN)                                 EIG 0076
            CC(IQ+NN,JQ+NN) = CC(IQ,JQ)                                 EIG 0077
C                                      ** GET FACTORS OF COEFF. MATRIX  EIG 0078
C                                      ** OF REDUCED EIGENVALUE PROBLEM EIG 0079
            ALPHA =   CC(IQ,JQ) / CMU(IQ)                               EIG 0080
            BETA = CC(IQ,JQ+NN) / CMU(IQ)                               EIG 0081
            AMB(IQ,JQ) = ALPHA - BETA                                   EIG 0082
            APB(IQ,JQ) = ALPHA + BETA                                   EIG 0083
30       CONTINUE                                                       EIG 0084
         AMB(IQ,IQ) = AMB(IQ,IQ) - 1.0 / CMU(IQ)                        EIG 0085
         APB(IQ,IQ) = APB(IQ,IQ) - 1.0 / CMU(IQ)                        EIG 0086
                                                                        EIG 0087
C                                INSERT TO ELIMINATE SINGULARITY, IF ANYEIG 0088
         IF ( AMB(IQ,IQ).EQ.0. ) AMB(IQ,IQ) = 1.E-6                     EIG 0089
         IF ( APB(IQ,IQ).EQ.0. ) APB(IQ,IQ) = 1.E-6                     EIG 0090
                                                                        EIG 0091
40    CONTINUE                                                          EIG 0092
C                      ** FINISH CALCULATION OF COEFFICIENT MATRIX OF   EIG 0093
C                      ** REDUCED EIGENVALUE PROBLEM:  GET MATRIX       EIG 0094
C                      ** PRODUCT (ALFA+BETA)*(ALFA-BETA); SS(12)       EIG 0095
      DO 70  IQ = 1, NN                                                 EIG 0096
         DO 70  JQ = 1, NN                                              EIG 0097
            SUM = 0.                                                    EIG 0098
            DO 60  KQ = 1, NN                                           EIG 0099
               SUM = SUM + APB(IQ,KQ) * AMB(KQ,JQ)                      EIG 0100
60          CONTINUE                                                    EIG 0101
            ARRAY(IQ,JQ) = SUM                                          EIG 0102
70    CONTINUE                                                          EIG 0103
C                      ** FIND (REAL) EIGENVALUES AND EIGENVECTORS      EIG 0104
                                                                        EIG 0105
      CALL  ASYMTX( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WK,         EIG 0106
     $              AAD, EVECCD, EVALD, WKD )                           EIG 0107
                                                                        EIG 0108
      IF ( IER.GT.0 )  THEN                                             EIG 0109
         WRITE( *, '(//,A,I4,A)' )  ' ASYMTX--EIGENVALUE NO. ', IER,    EIG 0110
     $     '  DIDNT CONVERGE.  LOWER-NUMBERED EIGENVALUES WRONG.'       EIG 0111
         CALL  ERRMSG( 'ASYMTX--CONVERGENCE PROBLEMS', .TRUE. )         EIG 0112
      END IF                                                            EIG 0113
                                                                        EIG 0114
CDIR$ IVDEP                                                             EIG 0115
      DO 75  IQ = 1, NN                                                 EIG 0116
C                                INSERT TO ELIMINATE SINGULARITY, IF ANYEIG 0117
           IF ( EVAL(IQ).EQ.0. ) EVAL(IQ) = 1.E-8                       EIG 0118
         EVAL(IQ) = DSQRT( DABS( EVAL(IQ) ) )                           EIG 0119
         KK( IQ+NN ) = EVAL(IQ)                                         EIG 0120
C                                             ** ADD NEGATIVE EIGENVALUEEIG 0121
         KK( NN+1-IQ ) = - EVAL(IQ)                                     EIG 0122
75    CONTINUE                                                          EIG 0123
C                          ** FIND EIGENVECTORS (G+) + (G-) FROM SS(10) EIG 0124
C                          ** AND STORE TEMPORARILY IN -APB- ARRAY      EIG 0125
      DO 90  JQ = 1, NN                                                 EIG 0126
         DO 90  IQ = 1, NN                                              EIG 0127
            SUM = 0.                                                    EIG 0128
            DO 80  KQ = 1,NN                                            EIG 0129
               SUM = SUM + AMB(IQ,KQ) * EVECC(KQ,JQ)                    EIG 0130
80          CONTINUE                                                    EIG 0131
            APB(IQ,JQ) = SUM / EVAL(JQ)                                 EIG 0132
90    CONTINUE                                                          EIG 0133
                                                                        EIG 0134
      DO 100  JQ = 1, NN                                                EIG 0135
CDIR$ IVDEP                                                             EIG 0136
         DO 100  IQ = 1, NN                                             EIG 0137
            GPPLGM = APB(IQ,JQ)                                         EIG 0138
            GPMIGM = EVECC(IQ,JQ)                                       EIG 0139
C                                ** RECOVER EIGENVECTORS G+,G- FROM     EIG 0140
C                                   THEIR SUM AND DIFFERENCE; STACK THEMEIG 0141
C                                   TO GET EIGENVECTORS OF FULL SYSTEM  EIG 0142
C                                   SS(7) (JQ = EIGENVECTOR NUMBER)     EIG 0143
                                                                        EIG 0144
            EVECC(IQ,      JQ) = 0.5 * ( GPPLGM + GPMIGM )              EIG 0145
            EVECC(IQ+NN,   JQ) = 0.5 * ( GPPLGM - GPMIGM )              EIG 0146
                                                                        EIG 0147
C                                ** EIGENVECTORS CORRESPONDING TO       EIG 0148
C                                ** NEGATIVE EIGENVALUES (CORRESP. TO   EIG 0149
C                                ** REVERSING SIGN OF 'K' IN SS(10) )   EIG 0150
            GPPLGM = - GPPLGM                                           EIG 0151
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )              EIG 0152
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )              EIG 0153
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )                 EIG 0154
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )                 EIG 0155
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )              EIG 0156
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )              EIG 0157
100   CONTINUE                                                          EIG 0158
                                                                        EIG 0159
      RETURN                                                            EIG 0160
      END                                                               EIG 0161
