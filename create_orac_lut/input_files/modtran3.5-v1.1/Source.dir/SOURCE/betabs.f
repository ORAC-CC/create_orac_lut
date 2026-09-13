      REAL FUNCTION BETABS(COSANG,G)                                    BET 0001
C                                                                       BET 0002
C     FUNCTION BETABS SUPPLIES THE BACK SCATTER FRACTION FOR A          BET 0003
C     GIVEN ASYMMETRY FACTOR AND COSINE OF ANGLE (A.E.R. 1986).         BET 0004
C                                                                       BET 0005
C                         1                                             BET 0006
C                      1  /                                             BET 0007
C     BETA(COSANG)  =  -  |  P (-M,COSANG) dM                           BET 0008
C                      2  /   G                                         BET 0009
C                         0                                             BET 0010
C                                                                       BET 0011
C     WHERE P (-M,COSANG) IS THE HENYEY-GREENSTEIN PHASE FUNCTION.      BET 0012
C            G                                                          BET 0013
C                                                                       BET 0014
C     DECLARE INPUTS                                                    BET 0015
      REAL COSANG,G                                                     BET 0016
C                                                                       BET 0017
C     LIST COMMONS                                                      BET 0018
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               BET 0019
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           BET 0020
C                                                                       BET 0021
C     DECLARE LOCAL VARIABLES                                           BET 0022
      INTEGER N                                                         BET 0023
      REAL ABSCOS,ABSG,COSN,BMAX,BMIN,COSINE(9),A(10,5)                 BET 0024
      DATA COSINE/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0/          BET 0025
      DATA A/                                                           BET 0026
     1  0.5,    0.5,               0.5,               0.5,              BET 0027
     2          0.5,               0.5,               0.5,              BET 0028
     3          0.5,               0.5,               0.5,              BET 0029
     4  0.0,    0.13979292889,    -0.12019000482,    -0.46017123188,    BET 0030
     5         -0.406828796532,   -0.3001541656,     -0.553474414782,   BET 0031
     6         -0.626794663911,   -0.84678101,       -0.406823676662,   BET 0032
     7  0.0,   -1.5989873995,     -0.2724219928,      1.18747390274,    BET 0033
     8          0.49409050834,    -0.35928947292,     0.37397957743,    BET 0034
     9          0.18057740986,     0.50718036245,     0.01406832224,    BET 0035
     &  0.0,    3.5184116963,      0.6385960123,     -2.081230849,      BET 0036
     1         -1.0144699491,      0.1589475781,     -0.74761782865,    BET 0037
     2         -0.37416958202,    -0.374040109,       0.1055607702,     BET 0038
     3  0.0,   -2.5592172257,     -0.7459840146,      0.85392041407,    BET 0039
     4          0.4272082373,      0.00049606046,     0.42711266606,    BET 0040
     5          0.32038683614,     0.2136746594,     -0.2128054157/     BET 0041
C                                                                       BET 0042
C     SPECIAL CASES                                                     BET 0043
      ABSCOS=ABS(COSANG)                                                BET 0044
      ABSG=ABS(G)                                                       BET 0045
      IF(ABSCOS.LT..000001 .OR. ABSG.LT..000001)THEN                    BET 0046
          BETABS=.5                                                     BET 0047
          RETURN                                                        BET 0048
      ELSEIF(ABSG.GT..999999)THEN                                       BET 0049
          BETABS=0.                                                     BET 0050
          IF(COSANG.LT.0.)BETABS=1.                                     BET 0051
          IF(G.LT.0.)BETABS=1.-BETABS                                   BET 0052
          RETURN                                                        BET 0053
      ENDIF                                                             BET 0054
C                                                                       BET 0055
C     BACKSCATTERING INTERPOLATION                                      BET 0056
      IF(ABSCOS.LT..7)THEN                                              BET 0057
          N=10*ABSCOS+1                                                 BET 0058
      ELSEIF(ABSCOS.LT..8)THEN                                          BET 0059
          N=7                                                           BET 0060
      ELSEIF(ABSCOS.LE.1.)THEN                                          BET 0061
          N=8                                                           BET 0062
      ELSE                                                              BET 0063
          WRITE(IPR,'(A,F15.7)')                                        BET 0064
     1      ' FATAL ERROR in BETABS:  COSANG =',COSANG                  BET 0065
          STOP                                                          BET 0066
      ENDIF                                                             BET 0067
      COSN=COSINE(N)                                                    BET 0068
      BMAX=A(N,1)+ABSG*(A(N,2)+ABSG*(A(N,3)+ABSG*(A(N,4)+ABSG*A(N,5)))) BET 0069
      N=N+1                                                             BET 0070
      BMIN=A(N,1)+ABSG*(A(N,2)+ABSG*(A(N,3)+ABSG*(A(N,4)+ABSG*A(N,5)))) BET 0071
      BETABS=BMAX+(BMIN-BMAX)*(ABSCOS-COSN)/(COSINE(N)-COSN)            BET 0072
C                                                                       BET 0073
C     CHECK FOR ROUND-OFF ERROR PROBLEMS                                BET 0074
      IF(BETABS.GE..5)THEN                                              BET 0075
          BETABS=.5                                                     BET 0076
          RETURN                                                        BET 0077
      ENDIF                                                             BET 0078
      IF(BETABS.LT.0.)BETABS=0.                                         BET 0079
C                                                                       BET 0080
C     IF G IS NEGATIVE, THE COMPLEMENT IS REQUIRED.                     BET 0081
      IF(G.LT.0.)BETABS=1.-BETABS                                       BET 0082
      IF(COSANG.LT.0.)BETABS=1.-BETABS                                  BET 0083
      RETURN                                                            BET 0084
      END                                                               BET 0085
