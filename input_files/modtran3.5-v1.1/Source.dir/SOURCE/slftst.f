      SUBROUTINE  SLFTST( ACCUR, ALBEDO, BTEMP, DELTAM, DTAUC, FBEAM,   SLF 0001
     $                    FISOT, IBCND, LAMBER, NLYR, PLANK, NPHI,      SLF 0002
     $                    NUMU, NSTR, NTAU, ONLYFL, PHI, PHI0, PMOM,    SLF 0003
     $                    PRNT, SSALB, TEMIS, TEMPER, TTEMP, UMU,       SLF 0004
     $                    USRANG, USRTAU, UTAU, UMU0, WVNMHI, WVNMLO,   SLF 0005
     $                    COMPAR, FLUP, RFLDIR, RFLDN, UU )             SLF 0006
                                                                        SLF 0007
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            SLF 0008
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                SLF 0009
C       IF  COMPAR = FALSE, SAVE USER INPUT VALUES THAT WOULD OTHERWISE SLF 0010
C       BE DESTROYED AND REPLACE THEM WITH INPUT VALUES FOR SELF-TEST.  SLF 0011
C       IF  COMPAR = TRUE, COMPARE SELF-TEST CASE RESULTS WITH CORRECT  SLF 0012
C       ANSWERS AND RESTORE USER INPUT VALUES IF TEST IS PASSED.        SLF 0013
                                                                        SLF 0014
C       (SEE FILE 'DISORT.DOC' FOR VARIABLE DEFINITIONS.)               SLF 0015
                                                                        SLF 0016
C     I N T E R N A L    V A R I A B L E S:                             SLF 0017
                                                                        SLF 0018
C         ACC     RELATIVE ACCURACY REQUIRED FOR PASSING SELF-TEST      SLF 0019
C         ERRORN  RELATIVE ERRORS IN 'DISORT' OUTPUT VARIABLES          SLF 0020
C         OK      LOGICAL VARIABLE FOR DETERMINING FAILURE OF SELF-TEST SLF 0021
C         ALL VARIABLES ENDING IN 'S':  TEMPORARY 'S'TORAGE FOR INPUT   SLF 0022
C+---------------------------------------------------------------------+SLF 0023
      REAL*8     PMOM( 0:* ), TEMPER( 0:* )                             SLF 0024
      LOGICAL  COMPAR, DELTAM, LAMBER, PLANK, OK, ONLYFL, PRNT(*),      SLF 0025
     $         USRANG, USRTAU                                           SLF 0026
      REAL*8     PMOMS( 0:4 ), TEMPES (0:1 )                            SLF 0027
      LOGICAL  DELTAS, LAMBES, NOPLNS, ONLYFS, PRNTS ( 7 ), USRANS,     SLF 0028
     $         USRTAS, TSTBAD                                           SLF 0029
      SAVE     NLYRS, DTAUCS, SSALBS, PMOMS, NSTRS, USRANS, NUMUS,      SLF 0030
     $         UMUS, USRTAS, NTAUS, UTAUS, NPHIS, PHIS, IBCNDS,         SLF 0031
     $         FBEAMS, UMU0S, PHI0S, FISOTS, LAMBES, ALBEDS, DELTAS,    SLF 0032
     $         ONLYFS, ACCURS, NOPLNS, WVNMLS, WVNMHS, BTEMPS, TTEMPS,  SLF 0033
     $         TEMISS, TEMPES, PRNTS                                    SLF 0034
      DATA     ACC / 1.E-4 /                                            SLF 0035
                                                                        SLF 0036
                                                                        SLF 0037
      IF  ( .NOT.COMPAR )  THEN                                         SLF 0038
C                                              ** SAVE USER INPUT VALUESSLF 0039
         NLYRS = NLYR                                                   SLF 0040
         DTAUCS = DTAUC                                                 SLF 0041
         SSALBS = SSALB                                                 SLF 0042
         DO 1  N = 0, 4                                                 SLF 0043
            PMOMS(N) = PMOM(N)                                          SLF 0044
 1       CONTINUE                                                       SLF 0045
         NSTRS = NSTR                                                   SLF 0046
         USRANS = USRANG                                                SLF 0047
         NUMUS  = NUMU                                                  SLF 0048
         UMUS  = UMU                                                    SLF 0049
         USRTAS = USRTAU                                                SLF 0050
         NTAUS = NTAU                                                   SLF 0051
         UTAUS  = UTAU                                                  SLF 0052
         NPHIS = NPHI                                                   SLF 0053
         PHIS  = PHI                                                    SLF 0054
         IBCNDS = IBCND                                                 SLF 0055
         FBEAMS = FBEAM                                                 SLF 0056
         UMU0S = UMU0                                                   SLF 0057
         PHI0S = PHI0                                                   SLF 0058
         FISOTS  = FISOT                                                SLF 0059
         LAMBES = LAMBER                                                SLF 0060
         ALBEDS = ALBEDO                                                SLF 0061
         DELTAS = DELTAM                                                SLF 0062
         ONLYFS = ONLYFL                                                SLF 0063
         ACCURS = ACCUR                                                 SLF 0064
         NOPLNS = PLANK                                                 SLF 0065
         WVNMLS = WVNMLO                                                SLF 0066
         WVNMHS = WVNMHI                                                SLF 0067
         BTEMPS = BTEMP                                                 SLF 0068
         TTEMPS = TTEMP                                                 SLF 0069
         TEMISS = TEMIS                                                 SLF 0070
         TEMPES( 0 ) = TEMPER( 0 )                                      SLF 0071
         TEMPES( 1 ) = TEMPER( 1 )                                      SLF 0072
         DO 3 I = 1, 7                                                  SLF 0073
            PRNTS( I ) = PRNT( I )                                      SLF 0074
    3    CONTINUE                                                       SLF 0075
C                                     ** SET INPUT VALUES FOR SELF-TEST SLF 0076
         NLYR = 1                                                       SLF 0077
         DTAUC = 1.0                                                    SLF 0078
         SSALB = 0.9                                                    SLF 0079
C                          ** HAZE L MOMENTS                            SLF 0080
         PMOM(0) = 1.0                                                  SLF 0081
         PMOM(1) = 0.8042                                               SLF 0082
         PMOM(2) = 0.646094                                             SLF 0083
         PMOM(3) = 0.481851                                             SLF 0084
         PMOM(4) = 0.359056                                             SLF 0085
         NSTR = 4                                                       SLF 0086
         USRANG = .TRUE.                                                SLF 0087
         NUMU  = 1                                                      SLF 0088
         UMU  = 0.5                                                     SLF 0089
         USRTAU = .TRUE.                                                SLF 0090
         NTAU = 1                                                       SLF 0091
         UTAU  = 0.5                                                    SLF 0092
         NPHI = 1                                                       SLF 0093
         PHI  = 90.0                                                    SLF 0094
         IBCND = 0                                                      SLF 0095
         FBEAM = 3.14159265                                             SLF 0096
         UMU0 = 0.866                                                   SLF 0097
         PHI0 = 0.0                                                     SLF 0098
         FISOT  = 1.0                                                   SLF 0099
         LAMBER = .TRUE.                                                SLF 0100
         ALBEDO = 0.7                                                   SLF 0101
         DELTAM = .TRUE.                                                SLF 0102
         ONLYFL = .FALSE.                                               SLF 0103
         ACCUR = 1.E-4                                                  SLF 0104
         PLANK = .TRUE.                                                 SLF 0105
         WVNMLO = 0.0                                                   SLF 0106
         WVNMHI = 50000.                                                SLF 0107
         BTEMP = 300.0                                                  SLF 0108
         TTEMP = 100.0                                                  SLF 0109
         TEMIS = 0.8                                                    SLF 0110
         TEMPER( 0 ) = 210.0                                            SLF 0111
         TEMPER( 1 ) = 200.0                                            SLF 0112
         DO 5 I = 1, 7                                                  SLF 0113
            PRNT( I ) = .FALSE.                                         SLF 0114
    5    CONTINUE                                                       SLF 0115
                                                                        SLF 0116
      ELSE                                                              SLF 0117
C                                    ** COMPARE TEST CASE RESULTS WITH  SLF 0118
C                                    ** CORRECT ANSWERS AND ABORT IF BADSLF 0119
         OK = .TRUE.                                                    SLF 0120
         ERROR1 = ( UU  - 47.86005 ) / 47.86005                         SLF 0121
         ERROR2 = ( RFLDIR - 1.527286 ) / 1.527286                      SLF 0122
         ERROR3 = ( RFLDN - 28.37223 ) / 28.37223                       SLF 0123
         ERROR4 = ( FLUP   - 152.5853 ) / 152.5853                      SLF 0124
         IF( DABS(ERROR1).GT.ACC ) OK = TSTBAD( 'UU',     ERROR1 )      SLF 0125
         IF( DABS(ERROR2).GT.ACC ) OK = TSTBAD( 'RFLDIR', ERROR2 )      SLF 0126
         IF( DABS(ERROR3).GT.ACC ) OK = TSTBAD( 'RFLDN',  ERROR3 )      SLF 0127
         IF( DABS(ERROR4).GT.ACC ) OK = TSTBAD( 'FLUP',   ERROR4 )      SLF 0128
                                                                        SLF 0129
         IF( .NOT. OK )                                                 SLF 0130
     $       CALL ERRMSG( 'DISORT--SELF-TEST FAILED', .TRUE. )          SLF 0131
                                                                        SLF 0132
C                                           ** RESTORE USER INPUT VALUESSLF 0133
         NLYR = NLYRS                                                   SLF 0134
         DTAUC = DTAUCS                                                 SLF 0135
         SSALB = SSALBS                                                 SLF 0136
         DO 11  N = 0, 4                                                SLF 0137
            PMOM(N) = PMOMS(N)                                          SLF 0138
 11      CONTINUE                                                       SLF 0139
         NSTR = NSTRS                                                   SLF 0140
         USRANG = USRANS                                                SLF 0141
         NUMU  = NUMUS                                                  SLF 0142
         UMU  = UMUS                                                    SLF 0143
         USRTAU = USRTAS                                                SLF 0144
         NTAU = NTAUS                                                   SLF 0145
         UTAU  = UTAUS                                                  SLF 0146
         NPHI = NPHIS                                                   SLF 0147
         PHI  = PHIS                                                    SLF 0148
         IBCND = IBCNDS                                                 SLF 0149
         FBEAM = FBEAMS                                                 SLF 0150
         UMU0 = UMU0S                                                   SLF 0151
         PHI0 = PHI0S                                                   SLF 0152
         FISOT  = FISOTS                                                SLF 0153
         LAMBER = LAMBES                                                SLF 0154
         ALBEDO = ALBEDS                                                SLF 0155
         DELTAM = DELTAS                                                SLF 0156
         ONLYFL = ONLYFS                                                SLF 0157
         ACCUR = ACCURS                                                 SLF 0158
         PLANK = NOPLNS                                                 SLF 0159
         WVNMLO = WVNMLS                                                SLF 0160
         WVNMHI = WVNMHS                                                SLF 0161
         BTEMP = BTEMPS                                                 SLF 0162
         TTEMP = TTEMPS                                                 SLF 0163
         TEMIS = TEMISS                                                 SLF 0164
         TEMPER( 0 ) = TEMPES( 0 )                                      SLF 0165
         TEMPER( 1 ) = TEMPES( 1 )                                      SLF 0166
         DO 13  I = 1, 7                                                SLF 0167
            PRNT( I ) = PRNTS( I )                                      SLF 0168
   13    CONTINUE                                                       SLF 0169
      END IF                                                            SLF 0170
                                                                        SLF 0171
      RETURN                                                            SLF 0172
      END                                                               SLF 0173
