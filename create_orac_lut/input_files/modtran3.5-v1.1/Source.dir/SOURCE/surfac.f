      SUBROUTINE  SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER,           SRF 0001
     $                    MI, MAZIM, MXCMU, MXUMU, NN, NUMU, NSTR,      SRF 0002
     $                    ONLYFL, UMU, USRANG, YLM0, YLMC, YLMU, BDR,   SRF 0003
     $                    EMU, BEM, RMU )                               SRF 0004
                                                                        SRF 0005
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            SRF 0006
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                SRF 0007
C       SPECIFIES USER'S SURFACE BIDIRECTIONAL PROPERTIES, STWJ(21)     SRF 0008
                                                                        SRF 0009
C   I N P U T     V A R I A B L E S:                                    SRF 0010
                                                                        SRF 0011
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0                         SRF 0012
C       HLPR   :  LEGENDRE MOMENTS OF SURFACE BIDIRECTIONAL REFLECTIVITYSRF 0013
C                    (WITH 2K+1 FACTOR INCLUDED)                        SRF 0014
C       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT                          SRF 0015
C       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)             SRF 0016
C       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL             SRF 0017
C                 AT THE BEAM ANGLE                                     SRF 0018
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS            SRF 0019
C                 AT THE QUADRATURE ANGLES                              SRF 0020
C       YLMU   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS            SRF 0021
C                 AT THE USER ANGLES                                    SRF 0022
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)                        SRF 0023
                                                                        SRF 0024
C    O U T P U T     V A R I A B L E S:                                 SRF 0025
                                                                        SRF 0026
C       BDR :  SURFACE BIDIRECTIONAL REFLECTIVITY (COMPUTATIONAL ANGLES)SRF 0027
C       RMU :  SURFACE BIDIRECTIONAL REFLECTIVITY (USER ANGLES)         SRF 0028
C       BEM :  SURFACE DIRECTIONAL EMISSIVITY (COMPUTATIONAL ANGLES)    SRF 0029
C       EMU :  SURFACE DIRECTIONAL EMISSIVITY (USER ANGLES)             SRF 0030
                                                                        SRF 0031
C    I N T E R N A L     V A R I A B L E S:                             SRF 0032
                                                                        SRF 0033
C       DREF      DIRECTIONAL REFLECTIVITY                              SRF 0034
C       NMUG   :  NUMBER OF ANGLE COSINE QUADRATURE POINTS              SRF 0035
C                 ON (0,1) FOR INTEGRATING BIDIRECTIONAL REFLECTIVITY   SRF 0036
C                 TO GET DIRECTIONAL EMISSIVITY (IT IS NECESSARY TO USE SRF 0037
C                 A QUADRATURE SET DISTINCT FROM THE COMPUTATIONAL      SRF 0038
C                 ANGLES, BECAUSE THE COMPUTATIONAL ANGLES MAY NOT BE   SRF 0039
C                 DENSE ENOUGH -- I.E. 'NSTR' MAY BE TOO SMALL-- TO GIVESRF 0040
C                 AN ACCURATE APPROXIMATION FOR THE INTEGRATION).       SRF 0041
C       GMU    :  THE 'NMUG' ANGLE COSINE QUADRATURE POINTS ON (0,1)    SRF 0042
C       GWT    :  THE 'NMUG' ANGLE COSINE QUADRATURE WEIGHTS ON (0,1)   SRF 0043
C       YLMG   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS            SRF 0044
C                 AT THE 'NMUG' QUADRATURE ANGLES                       SRF 0045
C+---------------------------------------------------------------------+SRF 0046
      LOGICAL  LAMBER, ONLYFL, USRANG                                   SRF 0047
      REAL*8     BDR( MI,0:* ), BEM(*), EMU(*),                         SRF 0048
     $         HLPR(0:*), RMU( MXUMU,0:* ), UMU(*),                     SRF 0049
     $         YLM0(0:*), YLMC( 0:MXCMU,* ), YLMU( 0:MXCMU,* )          SRF 0050
      PARAMETER  ( NMUG = 10, MAXSTR = 100 )                            SRF 0051
      LOGICAL  PASS1                                                    SRF 0052
      REAL*8     GMU( NMUG ), GWT( NMUG ), YLMG( 0:MAXSTR, NMUG )       SRF 0053
      SAVE  PASS1, GMU, GWT, YLMG                                       SRF 0054
      DATA  PASS1 / .TRUE. /                                            SRF 0055
                                                                        SRF 0056
                                                                        SRF 0057
      IF ( PASS1 )  THEN                                                SRF 0058
         PASS1 = .FALSE.                                                SRF 0059
         CALL QGAUSN( NMUG, GMU, GWT )                                  SRF 0060
                                                                        SRF 0061
         CALL LEPOLY( NMUG, 0, MAXSTR, MAXSTR, GMU, YLMG )              SRF 0062
C                       ** CONVERT LEGENDRE POLYS. TO NEGATIVE -GMU-    SRF 0063
         SGN  = - 1.0                                                   SRF 0064
         DO 1  K = 0, MAXSTR                                            SRF 0065
            SGN = - SGN                                                 SRF 0066
            DO 1  JG = 1, NMUG                                          SRF 0067
               YLMG( K,JG ) = SGN * YLMG( K,JG )                        SRF 0068
 1       CONTINUE                                                       SRF 0069
                                                                        SRF 0070
      END IF                                                            SRF 0071
                                                                        SRF 0072
      CALL  ZEROIT( BDR, MI*(MI+1) )                                    SRF 0073
      CALL  ZEROIT( BEM, MI )                                           SRF 0074
                                                                        SRF 0075
      IF ( LAMBER .AND. MAZIM.EQ.0 ) THEN                               SRF 0076
                                                                        SRF 0077
         DO 20 IQ = 1, NN                                               SRF 0078
            BEM(IQ) = 1.0 - ALBEDO                                      SRF 0079
            DO 20 JQ = 0, NN                                            SRF 0080
               BDR(IQ,JQ) = ALBEDO                                      SRF 0081
20       CONTINUE                                                       SRF 0082
                                                                        SRF 0083
      ELSE IF ( .NOT.LAMBER ) THEN                                      SRF 0084
C                                  ** COMPUTE SURFACE BIDIRECTIONAL     SRF 0085
C                                  ** PROPERTIES AT COMPUTATIONAL ANGLESSRF 0086
         DO 60 IQ = 1, NN                                               SRF 0087
                                                                        SRF 0088
            DO 40 JQ = 1, NN                                            SRF 0089
              SUM = 0.0                                                 SRF 0090
              DO 30 K = MAZIM, NSTR-1                                   SRF 0091
                 SUM = SUM + HLPR(K) * YLMC(K,IQ) * YLMC(K,JQ+NN)       SRF 0092
30            CONTINUE                                                  SRF 0093
              BDR(IQ,JQ) = (2.-DELM0) * SUM                             SRF 0094
40          CONTINUE                                                    SRF 0095
                                                                        SRF 0096
            IF ( FBEAM.GT.0.0 )  THEN                                   SRF 0097
               SUM = 0.0                                                SRF 0098
               DO 50 K = MAZIM, NSTR-1                                  SRF 0099
                  SUM = SUM + HLPR(K) * YLMC(K,IQ) * YLM0(K)            SRF 0100
50             CONTINUE                                                 SRF 0101
               BDR(IQ,0) = (2.-DELM0) * SUM                             SRF 0102
            ENDIF                                                       SRF 0103
                                                                        SRF 0104
60       CONTINUE                                                       SRF 0105
                                                                        SRF 0106
         IF ( MAZIM.EQ.0 ) THEN                                         SRF 0107
                                                                        SRF 0108
            IF ( NSTR.GT.MAXSTR )  CALL                                 SRF 0109
     $           ERRMSG( 'SURFAC--PARAMETER MAXSTR TOO SMALL', .TRUE. ) SRF 0110
                                                                        SRF 0111
C                              ** INTEGRATE BIDIRECTIONAL REFLECTIVITY  SRF 0112
C                              ** AT REFLECTION POLAR ANGLES -CMU- AND  SRF 0113
C                              ** INCIDENT ANGLES -GMU- TO GET          SRF 0114
C                              ** DIRECTIONAL EMISSIVITY AT             SRF 0115
C                              ** COMPUTATIONAL ANGLES -CMU-.           SRF 0116
            DO 100  IQ = 1, NN                                          SRF 0117
               DREF = 0.0                                               SRF 0118
               DO 90  JG = 1, NMUG                                      SRF 0119
                  SUM = 0.0                                             SRF 0120
                  DO 80  K = 0, NSTR-1                                  SRF 0121
                     SUM = SUM + HLPR(K) * YLMC(K,IQ) * YLMG(K,JG)      SRF 0122
80                CONTINUE                                              SRF 0123
                  DREF = DREF + 2.* GWT(JG) * GMU(JG) * SUM             SRF 0124
90             CONTINUE                                                 SRF 0125
               BEM(IQ) = 1.0 - DREF                                     SRF 0126
100         CONTINUE                                                    SRF 0127
                                                                        SRF 0128
         END IF                                                         SRF 0129
                                                                        SRF 0130
      END IF                                                            SRF 0131
C                                       ** COMPUTE SURFACE BIDIRECTIONALSRF 0132
C                                       ** PROPERTIES AT USER ANGLES    SRF 0133
                                                                        SRF 0134
      IF ( .NOT.ONLYFL .AND. USRANG )  THEN                             SRF 0135
                                                                        SRF 0136
         CALL  ZEROIT( EMU, MXUMU )                                     SRF 0137
         CALL  ZEROIT( RMU, MXUMU*(MI+1) )                              SRF 0138
                                                                        SRF 0139
         DO 170 IU = 1, NUMU                                            SRF 0140
            IF ( UMU(IU).GT.0.0 )  THEN                                 SRF 0141
                                                                        SRF 0142
               IF ( LAMBER .AND. MAZIM.EQ.0 )  THEN                     SRF 0143
                  DO 110 IQ = 0, NN                                     SRF 0144
                     RMU(IU,IQ) = ALBEDO                                SRF 0145
110               CONTINUE                                              SRF 0146
                  EMU(IU) = 1.0 - ALBEDO                                SRF 0147
                                                                        SRF 0148
               ELSE IF ( .NOT.LAMBER ) THEN                             SRF 0149
                  DO 130 IQ = 1, NN                                     SRF 0150
                     SUM = 0.0                                          SRF 0151
                     DO 120 K = MAZIM, NSTR-1                           SRF 0152
                        SUM = SUM + HLPR(K) * YLMU(K,IU) * YLMC(K,IQ+NN)SRF 0153
120                  CONTINUE                                           SRF 0154
                     RMU(IU,IQ) = (2.-DELM0) * SUM                      SRF 0155
130               CONTINUE                                              SRF 0156
                                                                        SRF 0157
                  IF ( FBEAM.GT.0.0 )  THEN                             SRF 0158
                     SUM = 0.0                                          SRF 0159
                     DO 140 K = MAZIM, NSTR-1                           SRF 0160
                        SUM = SUM + HLPR(K) * YLMU(K,IU) * YLM0(K)      SRF 0161
140                  CONTINUE                                           SRF 0162
                     RMU(IU,0) = (2.-DELM0) * SUM                       SRF 0163
                  END IF                                                SRF 0164
                                                                        SRF 0165
                  IF ( MAZIM.EQ.0 ) THEN                                SRF 0166
                                                                        SRF 0167
C                               ** INTEGRATE BIDIRECTIONAL REFLECTIVITY SRF 0168
C                               ** AT REFLECTION ANGLES -UMU- AND       SRF 0169
C                               ** INCIDENT ANGLES -GMU- TO GET         SRF 0170
C                               ** DIRECTIONAL EMISSIVITY AT            SRF 0171
C                               ** USER ANGLES -UMU-.                   SRF 0172
                     DREF = 0.0                                         SRF 0173
                     DO 160 JG = 1, NMUG                                SRF 0174
                        SUM = 0.0                                       SRF 0175
                        DO 150 K = 0, NSTR-1                            SRF 0176
                           SUM = SUM + HLPR(K) * YLMU(K,IU) * YLMG(K,JG)SRF 0177
150                     CONTINUE                                        SRF 0178
                        DREF = DREF + 2.* GWT(JG) * GMU(JG) * SUM       SRF 0179
160                  CONTINUE                                           SRF 0180
                                                                        SRF 0181
                     EMU(IU) = 1.0 - DREF                               SRF 0182
                  END IF                                                SRF 0183
                                                                        SRF 0184
               END IF                                                   SRF 0185
            END IF                                                      SRF 0186
170      CONTINUE                                                       SRF 0187
                                                                        SRF 0188
      END IF                                                            SRF 0189
                                                                        SRF 0190
      RETURN                                                            SRF 0191
      END                                                               SRF 0192
