      SUBROUTINE  SETDIS( CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR, SDS 0001
     $                    GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT,   SDS 0002
     $                    MAXUMU, MAXCMU, MXCMU, NCUT, NLYR, NTAU, NN,  SDS 0003
     $                    NSTR, PLANK, NUMU, ONLYFL, OPRIM, PMOM,SSALB, SDS 0004
     $                    TAUC, TAUCPR, UTAU, UTAUPR, UMU, UMU0, USRTAU,SDS 0005
     $                    USRANG )                                      SDS 0006
                                                                        SDS 0007
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            SDS 0008
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                SDS 0009
C          PERFORM MISCELLANEOUS SETTING-UP OPERATIONS                  SDS 0010
                                                                        SDS 0011
C       ROUTINES CALLED:  ERRMSG, QGAUSN, ZEROIT                        SDS 0012
                                                                        SDS 0013
C       INPUT :  ALL ARE DISORT INPUT VARIABLES (SEE DOC FILE)          SDS 0014
                                                                        SDS 0015
C       OUTPUT:  NTAU,UTAU   IF USRTAU = FALSE                          SDS 0016
C                NUMU,UMU    IF USRANG = FALSE                          SDS 0017
C                CMU,CWT     COMPUTATIONAL POLAR ANGLES AND             SDS 0018
C                               CORRESPONDING QUADRATURE WEIGHTS        SDS 0019
C                EXPBEA      TRANSMISSION OF DIRECT BEAM                SDS 0020
C                FLYR        TRUNCATED FRACTION IN DELTA-M METHOD       SDS 0021
C                GL          PHASE FUNCTION LEGENDRE COEFFICIENTS MULTI-SDS 0022
C                              PLIED BY (2L+1) AND SINGLE-SCATTER ALBEDOSDS 0023
C                HLPR        LEGENDRE MOMENTS OF SURFACE BIDIRECTIONAL  SDS 0024
C                              REFLECTIVITY, TIMES 2K+1                 SDS 0025
C                LAYRU       COMPUTATIONAL LAYER IN WHICH UTAU FALLS    SDS 0026
C                LYRCUT      FLAG AS TO WHETHER RADIATION WILL BE ZEROEDSDS 0027
C                              BELOW LAYER NCUT                         SDS 0028
C                NCUT        COMPUTATIONAL LAYER WHERE ABSORPTION       SDS 0029
C                              OPTICAL DEPTH FIRST EXCEEDS  ABSCUT      SDS 0030
C                NN          NSTR / 2                                   SDS 0031
C                OPRIM       DELTA-M-SCALED SINGLE-SCATTER ALBEDO       SDS 0032
C                TAUCPR      DELTA-M-SCALED OPTICAL DEPTH               SDS 0033
C                UTAUPR      DELTA-M-SCALED VERSION OF  UTAU            SDS 0034
                                                                        SDS 0035
      LOGICAL  DELTAM, LAMBER, LYRCUT, PLANK, ONLYFL, USRTAU, USRANG    SDS 0036
      INTEGER  LAYRU(*)                                                 SDS 0037
      REAL*8     CMU(*), CWT(*), DTAUC(*), EXPBEA(0:*), FLYR(*),        SDS 0038
     $         GL(0:MXCMU,*), HL(0:*), HLPR(0:*), OPRIM(*),             SDS 0039
     $         PMOM(0:MAXCMU,*), SSALB(*), TAUC(0:*), TAUCPR(0:*),      SDS 0040
     $         UTAU(*), UTAUPR(*), UMU(*)                               SDS 0041
      DATA  ABSCUT / 10. /                                              SDS 0042
      IF ( .NOT.USRTAU ) THEN                                           SDS 0043
C                              ** SET OUTPUT LEVELS AT COMPUTATIONAL    SDS 0044
C                              ** LAYER BOUNDARIES                      SDS 0045
         NTAU = NLYR + 1                                                SDS 0046
         DO 30  LC = 0, NTAU-1                                          SDS 0047
            UTAU(LC+1) = TAUC(LC)                                       SDS 0048
30       CONTINUE                                                       SDS 0049
      END IF                                                            SDS 0050
C                        ** APPLY DELTA-M SCALING AND MOVE DESCRIPTION  SDS 0051
C                        ** OF COMPUTATIONAL LAYERS TO LOCAL VARIABLES  SDS 0052
      EXPBEA( 0 ) = 1.0                                                 SDS 0053
      ABSTAU = 0.0                                                      SDS 0054
      DO  60  LC = 1, NLYR                                              SDS 0055
         PMOM(0,LC) = 1.0                                               SDS 0056
         IF ( ABSTAU.LT.ABSCUT )  NCUT = LC                             SDS 0057
         ABSTAU = ABSTAU + ( 1. - SSALB(LC) ) * DTAUC(LC)               SDS 0058
                                                                        SDS 0059
         IF ( .NOT.DELTAM )  THEN                                       SDS 0060
            OPRIM(LC) = SSALB(LC)                                       SDS 0061
            TAUCPR(LC) = TAUC(LC)                                       SDS 0062
            DO 40  K = 0, NSTR-1                                        SDS 0063
               GL(K,LC) = (2*K+1) * OPRIM(LC) * PMOM(K,LC)              SDS 0064
 40         CONTINUE                                                    SDS 0065
            F = 0.0                                                     SDS 0066
         ELSE                                                           SDS 0067
C                                    ** DO DELTA-M TRANSFORMATION       SDS 0068
            F = PMOM( NSTR,LC )                                         SDS 0069
            OPRIM(LC) = SSALB(LC) * ( 1. - F ) / ( 1. - F * SSALB(LC) ) SDS 0070
            TAUCPR(LC) = TAUCPR(LC-1) + ( 1. - F*SSALB(LC) ) * DTAUC(LC)SDS 0071
            DO 50  K = 0, NSTR-1                                        SDS 0072
               GL(K,LC) = (2*K+1) * OPRIM(LC) * (PMOM(K,LC)-F) / (1.-F) SDS 0073
 50         CONTINUE                                                    SDS 0074
         ENDIF                                                          SDS 0075
                                                                        SDS 0076
         FLYR(LC) = F                                                   SDS 0077
         EXPBEA(LC) = 0.0                                               SDS 0078
         IF ( FBEAM.GT.0.0 )  EXPBEA(LC) = DEXP( - TAUCPR(LC) / UMU0 )  SDS 0079
60    CONTINUE                                                          SDS 0080
C                      ** IF NO THERMAL EMISSION, CUT OFF MEDIUM BELOW  SDS 0081
C                      ** ABSORPTION OPTICAL DEPTH = ABSCUT ( NOTE THAT SDS 0082
C                      ** DELTA-M TRANSFORMATION LEAVES ABSORPTION      SDS 0083
C                      ** OPTICAL DEPTH INVARIANT ).  NOT WORTH THE     SDS 0084
C                      ** TROUBLE FOR ONE-LAYER PROBLEMS, THOUGH.       SDS 0085
      LYRCUT = .FALSE.                                                  SDS 0086
      IF ( ABSTAU.GE.ABSCUT .AND. .NOT.PLANK .AND. IBCND.NE.1           SDS 0087
     $     .AND. NLYR.GT.1 )  LYRCUT =.TRUE.                            SDS 0088
      IF ( .NOT.LYRCUT )  NCUT = NLYR                                   SDS 0089
                                                                        SDS 0090
C                             ** SET ARRAYS DEFINING LOCATION OF USER   SDS 0091
C                             ** OUTPUT LEVELS WITHIN DELTA-M-SCALED    SDS 0092
C                             ** COMPUTATIONAL MESH                     SDS 0093
      DO 90  LU = 1, NTAU                                               SDS 0094
         DO 70 LC = 1, NLYR                                             SDS 0095
            IF ( UTAU(LU).GE.TAUC(LC-1) .AND. UTAU(LU).LE.TAUC(LC) )    SDS 0096
     $           GO TO 80                                               SDS 0097
70       CONTINUE                                                       SDS 0098
         LC = NLYR                                                      SDS 0099
                                                                        SDS 0100
80       UTAUPR(LU) = UTAU(LU)                                          SDS 0101
         IF(DELTAM) UTAUPR(LU) = TAUCPR(LC-1) + (1.-SSALB(LC)*FLYR(LC)) SDS 0102
     $                                        * (UTAU(LU) - TAUC(LC-1)) SDS 0103
         LAYRU(LU) = LC                                                 SDS 0104
90    CONTINUE                                                          SDS 0105
C                      ** CALCULATE COMPUTATIONAL POLAR ANGLE COSINES   SDS 0106
C                      ** AND ASSOCIATED QUADRATURE WEIGHTS FOR GAUSSIANSDS 0107
C                      ** QUADRATURE ON THE INTERVAL (0,1) (UPWARD)     SDS 0108
      NN = NSTR / 2                                                     SDS 0109
      CALL  QGAUSN( NN, CMU, CWT )                                      SDS 0110
C                                  ** DOWNWARD (NEG) ANGLES AND WEIGHTS SDS 0111
      DO 100  IQ = 1, NN                                                SDS 0112
         CMU(IQ+NN) = - CMU(IQ)                                         SDS 0113
         CWT(IQ+NN) =   CWT(IQ)                                         SDS 0114
100   CONTINUE                                                          SDS 0115
                                                                        SDS 0116
      IF ( FBEAM.GT.0.0 )  THEN                                         SDS 0117
C                               ** COMPARE BEAM ANGLE TO COMPU'L ANGLES SDS 0118
         DO 110  IQ = 1, NN                                             SDS 0119
            IF ( DABS(UMU0-CMU(IQ))/UMU0 .LT. 1.E-4 )  CALL ERRMSG      SDS 0120
     $         ( 'SETDIS--BEAM ANGLE=COMPUTATIONAL ANGLE; CHANGE NSTR', SDS 0121
     $            .TRUE. )                                              SDS 0122
  110    CONTINUE                                                       SDS 0123
      END IF                                                            SDS 0124
      IF ( .NOT.USRANG .OR. (ONLYFL .AND. MAXUMU.GE.NSTR) )  THEN       SDS 0125
                                                                        SDS 0126
C                                   ** SET OUTPUT POLAR ANGLES TO       SDS 0127
C                                   ** COMPUTATIONAL POLAR ANGLES       SDS 0128
            NUMU = NSTR                                                 SDS 0129
            DO 120  IU = 1, NN                                          SDS 0130
               UMU(IU) = - CMU(NN+1-IU)                                 SDS 0131
120         CONTINUE                                                    SDS 0132
            DO 121  IU = NN+1, NSTR                                     SDS 0133
               UMU(IU) = CMU(IU-NN)                                     SDS 0134
121         CONTINUE                                                    SDS 0135
      END IF                                                            SDS 0136
             IF ( USRANG .AND. IBCND.EQ.1 )  THEN                       SDS 0137
                                                                        SDS 0138
C                               ** SHIFT POSITIVE USER ANGLE COSINES TO SDS 0139
C                               ** UPPER LOCATIONS AND PUT NEGATIVES    SDS 0140
C                               ** IN LOWER LOCATIONS                   SDS 0141
         DO 140  IU = 1, NUMU                                           SDS 0142
            UMU(IU+NUMU) = UMU(IU)                                      SDS 0143
140      CONTINUE                                                       SDS 0144
         DO 141  IU = 1, NUMU                                           SDS 0145
            UMU(IU) = - UMU( 2*NUMU+1-IU)                               SDS 0146
141      CONTINUE                                                       SDS 0147
         NUMU = 2*NUMU                                                  SDS 0148
      END IF                                                            SDS 0149
                                                                        SDS 0150
      IF ( .NOT.LYRCUT .AND. .NOT.LAMBER )  THEN                        SDS 0151
         DO 160  K = 0, NSTR                                            SDS 0152
            HLPR(K) = (2*K+1) * HL(K)                                   SDS 0153
160      CONTINUE                                                       SDS 0154
      END IF                                                            SDS 0155
                                                                        SDS 0156
      RETURN                                                            SDS 0157
      END                                                               SDS 0158
