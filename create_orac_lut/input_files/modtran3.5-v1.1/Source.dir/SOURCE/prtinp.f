      SUBROUTINE  PRTINP( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,     PIN 0001
     $                    WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,          PIN 0002
     $                    NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,          PIN 0003
     $                    FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,      PIN 0004
     $                    TEMIS, DELTAM, PLANK, ONLYFL, ACCUR,          PIN 0005
     $                    FLYR, LYRCUT, OPRIM, TAUC, TAUCPR,            PIN 0006
     $                    MAXCMU, PRTMOM )                              PIN 0007
                                                                        PIN 0008
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            PIN 0009
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                PIN 0010
C        PRINT VALUES OF INPUT VARIABLES                                PIN 0011
                                                                        PIN 0012
      LOGICAL  DELTAM, LAMBER, LYRCUT, PLANK, ONLYFL, PRTMOM            PIN 0013
      REAL*8     UMU(*), FLYR(*), DTAUC(*), OPRIM(*), PHI(*),           PIN 0014
     $         PMOM( 0:MAXCMU,* ), SSALB(*), UTAU(*), TAUC( 0:* ),      PIN 0015
     $         TAUCPR( 0:* ), TEMPER( 0:* ), HL( 0:MAXCMU )             PIN 0016
                                                                        PIN 0017
                                                                        PIN 0018
      WRITE( *,1010 )  NSTR, NLYR                                       PIN 0019
      IF ( IBCND.NE.1 ) WRITE( *,1030 )  NTAU, (UTAU(LU), LU = 1, NTAU) PIN 0020
      IF ( .NOT.ONLYFL )                                                PIN 0021
     $      WRITE( *,1040 )  NUMU, ( UMU(IU), IU = 1, NUMU )            PIN 0022
      IF ( .NOT.ONLYFL .AND. IBCND.NE.1 )                               PIN 0023
     $      WRITE( *,1050 )  NPHI, ( PHI(J), J = 1, NPHI )              PIN 0024
      IF ( .NOT.PLANK .OR. IBCND.EQ.1 )  WRITE( *,1100 )                PIN 0025
      WRITE( *,1055 )  IBCND                                            PIN 0026
      IF ( IBCND.EQ.0 )  THEN                                           PIN 0027
         WRITE( *,1060 ) FBEAM, UMU0, PHI0, FISOT                       PIN 0028
         IF ( LAMBER )   WRITE( *,1080 ) ALBEDO                         PIN 0029
         IF ( .NOT.LAMBER )  WRITE( *,1090 ) ( HL(K), K = 0, NSTR )     PIN 0030
         IF ( PLANK )  WRITE( *,1110 ) WVNMLO, WVNMHI, BTEMP,           PIN 0031
     $                                 TTEMP, TEMIS                     PIN 0032
      ELSE IF ( IBCND.EQ.1 )  THEN                                      PIN 0033
         WRITE( *,1070 )                                                PIN 0034
         WRITE( *,1080 ) ALBEDO                                         PIN 0035
      ENDIF                                                             PIN 0036
      IF ( DELTAM )      WRITE( *,1120 )                                PIN 0037
      IF ( .NOT.DELTAM ) WRITE( *,1130 )                                PIN 0038
      IF ( IBCND.EQ.1 )  THEN                                           PIN 0039
         WRITE( *,1135 )                                                PIN 0040
      ELSE IF ( ONLYFL )  THEN                                          PIN 0041
         WRITE( *,1140 )                                                PIN 0042
      ELSE                                                              PIN 0043
         WRITE( *,1150 )                                                PIN 0044
      ENDIF                                                             PIN 0045
      WRITE( *,1160 )  ACCUR                                            PIN 0046
      IF ( LYRCUT )  WRITE( *,1170 )                                    PIN 0047
      IF ( PLANK )  WRITE ( *,1190 )                                    PIN 0048
      IF ( .NOT.PLANK )  WRITE ( *,1191 )                               PIN 0049
      YESSCT = 0.0                                                      PIN 0050
      DO 10 LC = 1, NLYR                                                PIN 0051
         YESSCT = YESSCT + SSALB(LC)                                    PIN 0052
         IF( PLANK )                                                    PIN 0053
     $       WRITE( *,1200 )  LC, DTAUC(LC), TAUC(LC), SSALB(LC),       PIN 0054
     $                    FLYR(LC), TAUCPR(LC)-TAUCPR(LC-1), TAUCPR(LC),PIN 0055
     $                    OPRIM(LC), PMOM(1,LC), TEMPER(LC-1)           PIN 0056
         IF( .NOT.PLANK )                                               PIN 0057
     $       WRITE( *,1200 )  LC, DTAUC(LC), TAUC(LC), SSALB(LC),       PIN 0058
     $                    FLYR(LC), TAUCPR(LC)-TAUCPR(LC-1), TAUCPR(LC),PIN 0059
     $                    OPRIM(LC), PMOM(1,LC)                         PIN 0060
 10   CONTINUE                                                          PIN 0061
      IF( PLANK )  WRITE( *,1210 ) TEMPER(NLYR)                         PIN 0062
                                                                        PIN 0063
      IF( PRTMOM .AND. YESSCT.GT.0.0 )  THEN                            PIN 0064
         WRITE( *, '(/,A)' )  ' LAYER   PHASE FUNCTION MOMENTS'         PIN 0065
         DO 20 LC = 1, NLYR                                             PIN 0066
            IF( SSALB(LC).GT.0.0 )                                      PIN 0067
     $          WRITE( *,1300 )  LC, ( PMOM(K,LC), K = 0, NSTR )        PIN 0068
 20      CONTINUE                                                       PIN 0069
      ENDIF                                                             PIN 0070
                                                                        PIN 0071
      RETURN                                                            PIN 0072
                                                                        PIN 0073
1010  FORMAT ( /, ' NO. STREAMS =', I4,                                 PIN 0074
     $  '     NO. COMPUTATIONAL LAYERS =', I4 )                         PIN 0075
1030  FORMAT( I4,' USER OPTICAL DEPTHS :',10F10.4, /, (26X,10F10.4) )   PIN 0076
1040  FORMAT( I4,' USER POLAR ANGLE COSINES :',10F9.5,/,(31X,10F9.5) )  PIN 0077
1050  FORMAT( I4,' USER AZIMUTHAL ANGLES :', 10F9.2, /, (28X,10F9.2) )  PIN 0078
1055  FORMAT( ' BOUNDARY CONDITION FLAG: IBCND =', I2 )                 PIN 0079
1060  FORMAT( '    INCIDENT BEAM WITH INTENSITY =', 1P,E11.3, ' AND',   PIN 0080
     $ ' POLAR ANGLE COSINE = ', 0P,F8.5,'  AND AZIMUTH ANGLE =', F7.2, PIN 0081
     $ /,'    PLUS ISOTROPIC INCIDENT INTENSITY =', 1P,E11.3 )          PIN 0082
1070  FORMAT( '    ISOTROPIC ILLUMINATION FROM TOP AND BOTTOM' )        PIN 0083
1080  FORMAT( '    BOTTOM ALBEDO (LAMBERTIAN) =', 0P,F8.4 )             PIN 0084
1090  FORMAT( '    LEGENDRE COEFFS OF BOTTOM BIDIRECTIONAL',            PIN 0085
     $ ' REFLECTIVITY :', /, (10X,10F9.5) )                             PIN 0086
1100  FORMAT( ' NO THERMAL EMISSION' )                                  PIN 0087
1110  FORMAT( '    THERMAL EMISSION IN WAVENUMBER INTERVAL :', 2F14.4,/,PIN 0088
     $   '    BOTTOM TEMPERATURE =', F10.2, '     TOP TEMPERATURE =',   PIN 0089
     $   F10.2,'    TOP EMISSIVITY =', F8.4 )                           PIN 0090
1120  FORMAT( ' USES DELTA-M METHOD' )                                  PIN 0091
1130  FORMAT( ' DOES NOT USE DELTA-M METHOD' )                          PIN 0092
1135  FORMAT( ' CALCULATE ALBEDO AND TRANSMISSIVITY OF MEDIUM',         PIN 0093
     $   ' VS. INCIDENT BEAM ANGLE' )                                   PIN 0094
1140  FORMAT( ' CALCULATE FLUXES AND AZIM-AVERAGED INTENSITIES ONLY' )  PIN 0095
1150  FORMAT( ' CALCULATE FLUXES AND INTENSITIES' )                     PIN 0096
1160  FORMAT( ' RELATIVE CONVERGENCE CRITERION FOR AZIMUTH SERIES =',   PIN 0097
     $   1P,E11.2 )                                                     PIN 0098
1170  FORMAT( ' SETS RADIATION = 0 BELOW ABSORPTION OPTICAL DEPTH 10' ) PIN 0099
1190  FORMAT( /, 37X, '<------------- DELTA-M --------------->', /,     PIN 0100
     $'                   TOTAL    SINGLE                           ',  PIN 0101
     $               'TOTAL    SINGLE', /,                              PIN 0102
     $'       OPTICAL   OPTICAL   SCATTER   TRUNCATED   ',              PIN 0103
     $   'OPTICAL   OPTICAL   SCATTER    ASYMM', /,                     PIN 0104
     $'         DEPTH     DEPTH    ALBEDO    FRACTION     ',            PIN 0105
     $     'DEPTH     DEPTH    ALBEDO   FACTOR   TEMPERATURE' )         PIN 0106
1191  FORMAT( /, 37X, '<------------- DELTA-M --------------->', /,     PIN 0107
     $'                   TOTAL    SINGLE                           ',  PIN 0108
     $               'TOTAL    SINGLE', /,                              PIN 0109
     $'       OPTICAL   OPTICAL   SCATTER   TRUNCATED   ',              PIN 0110
     $   'OPTICAL   OPTICAL   SCATTER    ASYMM', /,                     PIN 0111
     $'         DEPTH     DEPTH    ALBEDO    FRACTION     ',            PIN 0112
     $     'DEPTH     DEPTH    ALBEDO   FACTOR' )                       PIN 0113
1200  FORMAT( I4, 2F10.4, F10.5, F12.5, 2F10.4, F10.5, F9.4,F14.3 )     PIN 0114
1210  FORMAT( 85X, F14.3 )                                              PIN 0115
1300  FORMAT( I6, 10F11.6, /, (6X,10F11.6) )                            PIN 0116
                                                                        PIN 0117
      END                                                               PIN 0118
