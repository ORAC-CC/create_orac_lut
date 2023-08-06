      SUBROUTINE  CHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,     CKN 0001
     $                    WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,     CKN 0002
     $                    NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,     CKN 0003
     $                    PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP,       CKN 0004
     $                    TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, MAXCLY,   CKN 0005
     $                    MAXULV, MAXUMU, MAXCMU, MAXPHI, MXCLY,        CKN 0006
     $                    MXULV,  MXUMU,  MXCMU,  MXPHI, TAUC )         CKN 0007
                                                                        CKN 0008
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            CKN 0009
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                CKN 0010
C           CHECKS THE INPUT DIMENSIONS AND VARIABLES                   CKN 0011
                                                                        CKN 0012
      LOGICAL  WRTBAD, WRTDIM                                           CKN 0013
      LOGICAL  LAMBER, PLANK, ONLYFL, USRANG, USRTAU, INPERR            CKN 0014
      INTEGER  IBCND, MAXCLY, MAXUMU, MAXULV, MAXCMU, MAXPHI, NLYR,     CKN 0015
     $         NUMU, NSTR, NPHI, NTAU, MXCMU, MXUMU, MXPHI, MXCLY,      CKN 0016
     $         MXULV                                                    CKN 0017
      REAL*8     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,   CKN 0018
     $         HL( 0:MAXCMU ), PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ), CKN 0019
     $         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP, CKN 0020
     $         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV ),     CKN 0021
     $         TAUC( 0:* )                                              CKN 0022
                                                                        CKN 0023
                                                                        CKN 0024
      INPERR = .FALSE.                                                  CKN 0025
      IF ( NLYR.LT.1 ) INPERR = WRTBAD( 'NLYR' )                        CKN 0026
      IF ( NLYR.GT.MAXCLY ) INPERR = WRTBAD( 'MAXCLY' )                 CKN 0027
                                                                        CKN 0028
      DO 10  LC = 1, NLYR                                               CKN 0029
         IF ( DTAUC(LC).LT.0.0 ) INPERR = WRTBAD( 'DTAUC' )             CKN 0030
         IF ( SSALB(LC).LT.0.0 .OR. SSALB(LC).GT.1.0 )                  CKN 0031
     $        INPERR = WRTBAD( 'SSALB' )                                CKN 0032
         IF ( PLANK .AND. IBCND.NE.1 )  THEN                            CKN 0033
            IF( LC.EQ.1 .AND. TEMPER(0).LT.0.0 )                        CKN 0034
     $          INPERR = WRTBAD( 'TEMPER' )                             CKN 0035
            IF( TEMPER(LC).LT.0.0 ) INPERR = WRTBAD( 'TEMPER' )         CKN 0036
         ENDIF                                                          CKN 0037
         DO 5  K = 0, NSTR                                              CKN 0038
            IF( PMOM(K,LC).LT.-1.0 .OR. PMOM(K,LC).GT.1.0 )             CKN 0039
     $          INPERR = WRTBAD( 'PMOM' )                               CKN 0040
 5       CONTINUE                                                       CKN 0041
10    CONTINUE                                                          CKN 0042
                                                                        CKN 0043
      IF ( IBCND.EQ.1 )  THEN                                           CKN 0044
         IF ( MAXULV.LT.2 ) INPERR = WRTBAD( 'MAXULV' )                 CKN 0045
      ELSE IF ( USRTAU )  THEN                                          CKN 0046
         IF ( NTAU.LT.1 ) INPERR = WRTBAD( 'NTAU' )                     CKN 0047
         IF ( MAXULV.LT.NTAU ) INPERR = WRTBAD( 'MAXULV' )              CKN 0048
         DO 20  LU = 1, NTAU                                            CKN 0049
         IF(DABS(UTAU(LU)-TAUC(NLYR)).LE.1.E-4)UTAU(LU)=TAUC(NLYR)      CKN 0050
            IF( UTAU(LU).LT.0.0 .OR. UTAU(LU).GT.TAUC(NLYR) )           CKN 0051
     $          INPERR = WRTBAD( 'UTAU' )                               CKN 0052
20       CONTINUE                                                       CKN 0053
      ELSE                                                              CKN 0054
         IF ( MAXULV.LT.NLYR+1 ) INPERR = WRTBAD( 'MAXULV' )            CKN 0055
      END IF                                                            CKN 0056
                                                                        CKN 0057
      IF ( NSTR.LT.2 .OR. MOD(NSTR,2).NE.0 ) INPERR = WRTBAD( 'NSTR' )  CKN 0058
      IF ( NSTR.GT.MAXCMU ) INPERR = WRTBAD( 'MAXCMU' )                 CKN 0059
                                                                        CKN 0060
      IF ( USRANG )  THEN                                               CKN 0061
         IF ( NUMU.LT.0 ) INPERR = WRTBAD( 'NUMU' )                     CKN 0062
         IF ( .NOT.ONLYFL .AND. NUMU.EQ.0 ) INPERR = WRTBAD( 'NUMU'  )  CKN 0063
         IF ( NUMU.GT.MAXUMU ) INPERR = WRTBAD( 'MAXUMU' )              CKN 0064
         IF ( IBCND.EQ.1 .AND. 2*NUMU.GT.MAXUMU )                       CKN 0065
     $        INPERR = WRTBAD( 'MAXUMU' )                               CKN 0066
         DO 30  IU = 1, NUMU                                            CKN 0067
            IF( UMU(IU).LT.-1.0 .OR. UMU(IU).GT.1.0 .OR. UMU(IU).EQ.0.0)CKN 0068
     $           INPERR = WRTBAD( 'UMU' )                               CKN 0069
            IF( IBCND.EQ.1 .AND. UMU(IU).LT.0.0 )                       CKN 0070
     $           INPERR = WRTBAD( 'UMU' )                               CKN 0071
            IF( IU.GT.1 .AND. UMU(IU).LT.UMU(IU-1) )                    CKN 0072
     $           INPERR = WRTBAD( 'UMU' )                               CKN 0073
30       CONTINUE                                                       CKN 0074
      ELSE                                                              CKN 0075
         IF( MAXUMU.LT.NSTR ) INPERR = WRTBAD( 'MAXUMU' )               CKN 0076
      END IF                                                            CKN 0077
                                                                        CKN 0078
      IF ( .NOT.ONLYFL .AND. IBCND.NE.1 )  THEN                         CKN 0079
         IF ( NPHI.LE.0 ) INPERR = WRTBAD( 'NPHI' )                     CKN 0080
         IF ( NPHI.GT.MAXPHI ) INPERR = WRTBAD( 'MAXPHI' )              CKN 0081
         DO 40  J = 1, NPHI                                             CKN 0082
            IF ( PHI(J).LT.0.0 .OR. PHI(J).GT.360.0 )                   CKN 0083
     $           INPERR = WRTBAD( 'PHI' )                               CKN 0084
40       CONTINUE                                                       CKN 0085
      END IF                                                            CKN 0086
                                                                        CKN 0087
      IF ( IBCND.LT.0 .OR. IBCND.GT.1 ) INPERR = WRTBAD( 'IBCND' )      CKN 0088
      IF ( IBCND.EQ.0 )  THEN                                           CKN 0089
         IF ( FBEAM.LT.0.0 ) INPERR = WRTBAD( 'FBEAM' )                 CKN 0090
         IF ( FBEAM.GT.0.0 .AND. ( UMU0.LE.0.0 .OR. UMU0.GT.1.0 ) )     CKN 0091
     $        INPERR = WRTBAD( 'UMU0' )                                 CKN 0092
         IF ( FBEAM.GT.0.0 .AND. ( PHI0.LT.0.0 .OR. PHI0.GT.360.0 ) )   CKN 0093
     $        INPERR = WRTBAD( 'PHI0' )                                 CKN 0094
         IF ( FISOT.LT.0.0 ) INPERR = WRTBAD( 'FISOT' )                 CKN 0095
         IF ( LAMBER )  THEN                                            CKN 0096
            IF ( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )                     CKN 0097
     $           INPERR = WRTBAD( 'ALBEDO' )                            CKN 0098
         ELSE                                                           CKN 0099
C                    ** MAKE SURE FLUX ALBEDO AT DENSE MESH OF INCIDENT CKN 0100
C                       ANGLES DOES NOT ASSUME UNPHYSICAL VALUES        CKN 0101
                                                                        CKN 0102
            DO 50  RMU = 0.0, 1.0, 0.01                                 CKN 0103
               FLXALB = DREF( RMU, HL, NSTR )                           CKN 0104
               IF ( FLXALB.LT.0.0 .OR. FLXALB.GT.1.0 )                  CKN 0105
     $              INPERR = WRTBAD( 'HL' )                             CKN 0106
50          CONTINUE                                                    CKN 0107
         ENDIF                                                          CKN 0108
                                                                        CKN 0109
      ELSE IF ( IBCND.EQ.1 )  THEN                                      CKN 0110
         IF ( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )                        CKN 0111
     $        INPERR = WRTBAD( 'ALBEDO' )                               CKN 0112
      END IF                                                            CKN 0113
                                                                        CKN 0114
      IF ( PLANK .AND. IBCND.NE.1 )  THEN                               CKN 0115
         IF ( WVNMLO.LT.0.0 .OR. WVNMHI.LE.WVNMLO )                     CKN 0116
     $        INPERR = WRTBAD( 'WVNMLO,HI' )                            CKN 0117
         IF ( TEMIS.LT.0.0 .OR. TEMIS.GT.1.0 )                          CKN 0118
     $        INPERR = WRTBAD( 'TEMIS' )                                CKN 0119
         IF ( BTEMP.LT.0.0 ) INPERR = WRTBAD( 'BTEMP' )                 CKN 0120
         IF ( TTEMP.LT.0.0 ) INPERR = WRTBAD( 'TTEMP' )                 CKN 0121
      END IF                                                            CKN 0122
                                                                        CKN 0123
      IF ( ACCUR.LT.0.0 .OR. ACCUR.GT.1.E-2 )                           CKN 0124
     $     INPERR = WRTBAD( 'ACCUR' )                                   CKN 0125
                                                                        CKN 0126
      IF ( MXCLY.LT.NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR )             CKN 0127
      IF ( IBCND.NE.1 )  THEN                                           CKN 0128
         IF ( USRTAU .AND. MXULV.LT.NTAU )                              CKN 0129
     $        INPERR = WRTDIM( 'MXULV', NTAU )                          CKN 0130
         IF ( .NOT.USRTAU .AND. MXULV.LT.NLYR+1 )                       CKN 0131
     $        INPERR = WRTDIM( 'MXULV', NLYR+1 )                        CKN 0132
      ELSE                                                              CKN 0133
         IF ( MXULV.LT.2 ) INPERR = WRTDIM( 'MXULV', 2 )                CKN 0134
      END IF                                                            CKN 0135
      IF ( MXCMU.LT.NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR )             CKN 0136
      IF ( USRANG .AND. MXUMU.LT.NUMU )                                 CKN 0137
     $     INPERR = WRTDIM( 'MXUMU', NUMU )                             CKN 0138
      IF ( USRANG .AND. IBCND.EQ.1 .AND. MXUMU.LT.2*NUMU )              CKN 0139
     $     INPERR = WRTDIM( 'MXUMU', NUMU )                             CKN 0140
      IF ( .NOT.USRANG .AND. MXUMU.LT.NSTR )                            CKN 0141
     $      INPERR = WRTDIM( 'MXUMU', NSTR )                            CKN 0142
      IF ( .NOT.ONLYFL .AND. IBCND.NE.1 .AND. MXPHI.LT.NPHI )           CKN 0143
     $      INPERR = WRTDIM( 'MXPHI', NPHI )                            CKN 0144
                                                                        CKN 0145
      IF ( INPERR )                                                     CKN 0146
     $   CALL ERRMSG( 'DISORT--input and/or dimension errors', .TRUE. ) CKN 0147
C                                                                       CKN 0148
C     COMMENT OUT THIS WARNING MESSAGE (J. VAIL, OCT95).                CKN 0149
C     IF ( PLANK )  THEN                                                CKN 0150
C        DO 100  LC = 1, NLYR                                           CKN 0151
C        IF ( DABS(TEMPER(LC)-TEMPER(LC-1)) .GT. 20.0 )                 CKN 0152
C    $        CALL ERRMSG( 'CHEKIN--vertical temperature step may'      CKN 0153
C    $                  // ' be too large for good accuracy', .FALSE. ) CKN 0154
C 100    CONTINUE                                                       CKN 0155
C     END IF                                                            CKN 0156
                                                                        CKN 0157
      RETURN                                                            CKN 0158
      END                                                               CKN 0159
