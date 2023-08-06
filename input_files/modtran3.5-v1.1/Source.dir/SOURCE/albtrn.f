      SUBROUTINE  ALBTRN( ALBEDO, AMB, APB, ARRAY, B, BDR, CBAND, CC,   ALB 0001
     $                    CMU, CWT, EVAL, EVECC, GL, GC, GU, IPVT,      ALB 0002
     $                    KK, LL, NLYR, NN, NSTR, NUMU, PRNT, TAUCPR,   ALB 0003
     $                    UMU, U0U, WK, YLMC, YLMU, Z, AAD, EVALD,      ALB 0004
     $                    EVECCD, WKD, MI, MI9M2, MAXULV, MAXUMU,       ALB 0005
     $                    MXCMU, MXUMU, NNLYRI, ALBMED, TRNMED )        ALB 0006
                                                                        ALB 0007
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            ALB 0008
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                ALB 0009
                                                                        ALB 0010
                                                                        ALB 0011
C        SPECIAL CASE TO GET ONLY ALBEDO AND TRANSMISSIVITY             ALB 0012
C        OF ENTIRE MEDIUM AS A FUNCTION OF INCIDENT BEAM ANGLE          ALB 0013
C        (MANY SIMPLIFICATIONS BECAUSE BOUNDARY CONDITION IS JUST       ALB 0014
C        ISOTROPIC ILLUMINATION, THERE ARE NO THERMAL SOURCES, AND      ALB 0015
C        PARTICULAR SOLUTIONS DO NOT NEED TO BE COMPUTED).  SEE         ALB 0016
C        REF. S2 AND REFERENCES THEREIN FOR THEORY.                     ALB 0017
                                                                        ALB 0018
C        ROUTINES CALLED:  ALTRIN, LEPOLY, PRALTR, SETMTX, SOLVE1,      ALB 0019
C                          SOLEIG, ZEROIT                               ALB 0020
                                                                        ALB 0021
      LOGICAL  PRNT(*)                                                  ALB 0022
      INTEGER  NLYR, NUMU, NSTR                                         ALB 0023
      REAL*8     UMU(*), U0U( MAXUMU,* )                                ALB 0024
                                                                        ALB 0025
      INTEGER IPVT(*)                                                   ALB 0026
      REAL*8    ALBMED(*), AMB( MI,* ), APB( MI,*), ARRAY( MXCMU,* ),   ALB 0027
     $        B(*), BDR( MI,0:* ), CBAND( MI9M2,* ), CC( MXCMU,* ),     ALB 0028
     $        CMU(*), CWT(*), EVAL(*), EVECC( MXCMU,* ),                ALB 0029
     $        GL( 0:MXCMU,* ), GC( MXCMU,MXCMU,* ), GU( MXUMU,MXCMU,* ),ALB 0030
     $        KK( MXCMU,* ), LL( MXCMU,* ), TAUCPR( 0:* ), TRNMED(*),   ALB 0031
     $        WK(*), YLMC( 0:MXCMU,* ), YLMU( 0:MXCMU,* ), Z(*)         ALB 0032
      DOUBLE PRECISION   AAD( MI,* ), EVALD(*) , EVECCD( MI,* ), WKD(*) ALB 0033
                                                                        ALB 0034
      LOGICAL  LAMBER, LYRCUT                                           ALB 0035
                                                                        ALB 0036
                                                                        ALB 0037
C                    ** SET DISORT VARIABLES THAT ARE IGNORED IN THIS   ALB 0038
C                    ** SPECIAL CASE BUT ARE NEEDED BELOW IN ARGUMENT   ALB 0039
C                    ** LISTS OF SUBROUTINES SHARED WITH GENERAL CASE   ALB 0040
      NCUT = NLYR                                                       ALB 0041
      LYRCUT = .FALSE.                                                  ALB 0042
      FISOT = 1.0                                                       ALB 0043
      LAMBER = .TRUE.                                                   ALB 0044
                                                                        ALB 0045
      MAZIM = 0                                                         ALB 0046
      DELM0 = 1.0                                                       ALB 0047
C                          ** GET LEGENDRE POLYNOMIALS FOR COMPUTATIONALALB 0048
C                          ** AND USER POLAR ANGLE COSINES              ALB 0049
                                                                        ALB 0050
      CALL  LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, YLMU )             ALB 0051
      CALL  LEPOLY( NN,   MAZIM, MXCMU, NSTR-1, CMU, YLMC )             ALB 0052
                                                                        ALB 0053
C                       ** EVALUATE LEGENDRE POLYNOMIALS WITH NEGATIVE  ALB 0054
C                       ** -CMU- FROM THOSE WITH POSITIVE -CMU-;        ALB 0055
C                       ** DAVE/ARMSTRONG EQ. (15)                      ALB 0056
      SGN  = -1.0                                                       ALB 0057
      DO  5  L = MAZIM, NSTR-1                                          ALB 0058
         SGN = - SGN                                                    ALB 0059
         DO  5  IQ = NN+1, NSTR                                         ALB 0060
            YLMC( L,IQ ) = SGN * YLMC( L,IQ-NN )                        ALB 0061
    5 CONTINUE                                                          ALB 0062
C                                  ** ZERO BOTTOM REFLECTIVITY          ALB 0063
C                                  ** (-ALBEDO- IS USED ONLY IN ANALYTICALB 0064
C                                  ** FORMULAE INVOLVING ALBEDO = 0     ALB 0065
C                                  ** SOLUTIONS; EQS 16-17 OF REF S2)   ALB 0066
      CALL  ZEROIT( BDR, MI*(MI+1) )                                    ALB 0067
                                                                        ALB 0068
                                                                        ALB 0069
C ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============ALB 0070
                                                                        ALB 0071
      DO 100  LC = 1, NLYR                                              ALB 0072
                                                                        ALB 0073
C                        ** SOLVE EIGENFUNCTION PROBLEM IN EQ. STWJ(8B) ALB 0074
                                                                        ALB 0075
         CALL  SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL(0,LC), MI, MAZIM,  ALB 0076
     $                 MXCMU, NN, NSTR, WK, YLMC, CC, EVECC, EVAL,      ALB 0077
     $                 KK(1,LC), GC(1,1,LC), AAD, WKD, EVECCD, EVALD)   ALB 0078
                                                                        ALB 0079
C                          ** INTERPOLATE EIGENVECTORS TO USER ANGLES   ALB 0080
                                                                        ALB 0081
         CALL  TERPEV( CWT, EVECC, GL(0,LC), GU(1,1,LC), MAZIM, MXCMU,  ALB 0082
     $                 MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU )          ALB 0083
100   CONTINUE                                                          ALB 0084
                                                                        ALB 0085
C ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============ALB 0086
                                                                        ALB 0087
C                      ** SET COEFFICIENT MATRIX OF EQUATIONS COMBINING ALB 0088
C                      ** BOUNDARY AND LAYER INTERFACE CONDITIONS       ALB 0089
                                                                        ALB 0090
      CALL  SETMTX( BDR, CBAND, CMU, CWT, DELM0, GC, KK, LAMBER,        ALB 0091
     $              LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,       ALB 0092
     $              NN, NSTR, TAUCPR, WK )                              ALB 0093
                                                                        ALB 0094
      CALL  ZEROIT( U0U, MAXUMU*MAXULV )                                ALB 0095
                                                                        ALB 0096
      NHOM = 2                                                          ALB 0097
      IF( NLYR.EQ.1 )  NHOM = 1                                         ALB 0098
      SPHALB = 0.0                                                      ALB 0099
      SPHTRN = 0.0                                                      ALB 0100
      DO 200  IHOM = 1, NHOM                                            ALB 0101
C                             ** SOLVE FOR CONSTANTS OF INTEGRATION IN  ALB 0102
C                             ** HOMOGENEOUS SOLUTION FOR ILLUMINATION  ALB 0103
C                             ** FROM TOP (IHOM=1), THEN BOTTOM (IHOM=2)ALB 0104
                                                                        ALB 0105
         CALL  SOLVE1( B, CBAND, FISOT, IHOM, IPVT, LL, MI9M2, MXCMU,   ALB 0106
     $                 NCOL, NLYR, NN, NNLYRI, NSTR, Z )                ALB 0107
                                                                        ALB 0108
C                             ** COMPUTE AZIMUTHALLY-AVERAGED INTENSITY ALB 0109
C                             ** AT USER ANGLES; GIVES ALBEDO IF MULTI- ALB 0110
C                             ** LAYER (EQ. 9 OF REF S2); GIVES BOTH    ALB 0111
C                             ** ALBEDO AND TRANSMISSIVITY IF SINGLE    ALB 0112
C                             ** LAYER (EQS. 3-4 OF REF S2)             ALB 0113
                                                                        ALB 0114
         CALL  ALTRIN( GU, KK, LL, MXCMU, MXUMU, MAXUMU, NLYR,          ALB 0115
     $                 NN, NSTR, NUMU, TAUCPR, UMU, U0U, WK )           ALB 0116
                                                                        ALB 0117
         IF ( IHOM.EQ.1 )  THEN                                         ALB 0118
C                                   ** SAVE ALBEDOS;  FLIP TRANSMISSIV. ALB 0119
C                                   ** END OVER END TO CORRESPOND TO    ALB 0120
C                                   ** POSITIVE -UMU- INST. OF NEGATIVE ALB 0121
            DO 120  IU = 1, NUMU/2                                      ALB 0122
               ALBMED(IU) = U0U( IU + NUMU/2, 1 )                       ALB 0123
               IF( NLYR.EQ.1 )  TRNMED(IU) = U0U( NUMU/2+1-IU, 2 )      ALB 0124
     $                          + EXP( - TAUCPR(NLYR) / UMU(IU+NUMU/2) )ALB 0125
120         CONTINUE                                                    ALB 0126
C                                    ** GET SPHERICAL ALBEDO AND, FOR 1 ALB 0127
C                                    ** LAYER, SPHERICAL TRANSMISSIVITY ALB 0128
            IF( ALBEDO.GT.0.0 )                                         ALB 0129
     $          CALL SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,         ALB 0130
     $                       NN, NSTR, TAUCPR, SPHALB, SPHTRN )         ALB 0131
                                                                        ALB 0132
         ELSE IF ( IHOM.EQ.2 )  THEN                                    ALB 0133
C                                      ** SAVE TRANSMISSIVITIES         ALB 0134
            DO 140  IU = 1, NUMU/2                                      ALB 0135
               TRNMED(IU) = U0U( IU + NUMU/2, 1 )                       ALB 0136
     $                      + EXP( - TAUCPR(NLYR) / UMU(IU+NUMU/2) )    ALB 0137
140         CONTINUE                                                    ALB 0138
C                             ** GET SPHERICAL ALBEDO AND TRANSMISSIVITYALB 0139
            IF( ALBEDO.GT.0.0 )                                         ALB 0140
     $          CALL SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,         ALB 0141
     $                       NN, NSTR, TAUCPR, SPHTRN, SPHALB )         ALB 0142
         END IF                                                         ALB 0143
200   CONTINUE                                                          ALB 0144
                                                                        ALB 0145
      IF ( ALBEDO.GT.0.0 )  THEN                                        ALB 0146
C                                ** REF. S2, EQS. 16-17 (THESE EQS. HAVEALB 0147
C                                ** A SIMPLE PHYSICAL INTERPRETATION    ALB 0148
C                                ** LIKE THAT OF THE DOUBLING EQS.)     ALB 0149
         DO 220  IU = 1, NUMU                                           ALB 0150
            ALBMED(IU) = ALBMED(IU) + ( ALBEDO / (1.-ALBEDO*SPHALB) )   ALB 0151
     $                                * SPHTRN * TRNMED(IU)             ALB 0152
            TRNMED(IU) = TRNMED(IU) + ( ALBEDO / (1.-ALBEDO*SPHALB) )   ALB 0153
     $                                * SPHALB * TRNMED(IU)             ALB 0154
220      CONTINUE                                                       ALB 0155
      END IF                                                            ALB 0156
C                          ** RETURN -UMU- TO ALL POSITIVE VALUES, TO   ALB 0157
C                          ** AGREE WITH ORDERING IN -ALBMED,TRNMED-    ALB 0158
      NUMU = NUMU / 2                                                   ALB 0159
      DO 230  IU = 1, NUMU                                              ALB 0160
        UMU(IU) = UMU(IU+NUMU)                                          ALB 0161
 230  CONTINUE                                                          ALB 0162
                                                                        ALB 0163
      IF ( PRNT(6) )  CALL  PRALTR( UMU, NUMU, ALBMED, TRNMED )         ALB 0164
                                                                        ALB 0165
      RETURN                                                            ALB 0166
      END                                                               ALB 0167
