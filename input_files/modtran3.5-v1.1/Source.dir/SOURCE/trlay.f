      SUBROUTINE TRLAY(TAU,TSCAT,G,CSZEN,S0DEP,DEPRAT,EX,TDFS,REFS)     TRL 0001
C                                                                       TRL 0002
C     CALCULATE PARAMETERS FOR SOLAR HYBRID MODIFIED DELTA EDDINGTON    TRL 0003
C     2-STREAM APPROXIMATION.  THE EQUATIONS ARE DERIVED FROM           TRL 0004
C     W. E. MEADOR AND W. R. WEAVER, J. ATMOS. SCI. 37, 630-642 (1980). TRL 0005
C                                                                       TRL 0006
C     INPUTS                                                            TRL 0007
C       TAU     VERTICAL EXTINCTION OPTICAL DEPTH.                      TRL 0008
C       TSCAT   VERTICAL SCATTERING OPTICAL DEPTH.                      TRL 0009
C       G       HENYEY-GREENSTEIN ASYMMETRY FACTOR.                     TRL 0010
C       CSZEN   COSINE OF THE SOLAR ZENITH ANGLE.                       TRL 0011
C       S0DEP   EXTINCTION OPTICAL DEPTH FROM LAYER BOTTOM TO SUN/MOON. TRL 0012
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN   TRL 0013
C               ACROSS THE CURRENT LAYER.                               TRL 0014
C                                                                       TRL 0015
C     OUTPUTS                                                           TRL 0016
C       EX      SOLAR ATTENUATION ACROSS LAYER                          TRL 0017
C       TDFS    LAYER TRANSMITTANCE MINUS EX                            TRL 0018
C       REFS    LAYER REFLECTANCE                                       TRL 0019
      REAL TAU,TSCAT,G,S0DEP,DEPRAT,CSZEN,EX,TDFS,REFS                  TRL 0020
C                                                                       TRL 0021
C     LIST PARAMETERS                                                   TRL 0022
      REAL CUTOFF,PT3CUT                                                TRL 0023
      PARAMETER(CUTOFF=0.01,PT3CUT=.3*CUTOFF)                           TRL 0024
C                                                                       TRL 0025
C     LIST COMMONS                                                      TRL 0026
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               TRL 0027
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           TRL 0028
C                                                                       TRL 0029
C     DECLARE FUNCTIONS                                                 TRL 0030
      REAL BETABS                                                       TRL 0031
C                                                                       TRL 0032
C     DECLARE LOCAL VARIABLES                                           TRL 0033
      REAL OMEGA,MU0,GAMMA1,G1MG2,GAMMA3,GAMMA4,ONEMG2,DENOM,           TRL 0034
     1  A1,A2,K,KT,KM,X,ONEPKM,ONEMKM,G3M,G4M,XHALF,TERM1,              TRL 0035
     2  TERM2,TERM3,COEF,EKT,TWOEKT,RATP,RATM,XMKT,G3K,G4K              TRL 0036
      LOGICAL LWARN                                                     TRL 0037
C                                                                       TRL 0038
C     LIST DATA                                                         TRL 0039
      SAVE LWARN                                                        TRL 0040
      DATA LWARN/.TRUE./                                                TRL 0041
      IF(TAU.LE.0.)THEN                                                 TRL 0042
C                                                                       TRL 0043
C         CASE 1:  TAU = 0.                                             TRL 0044
          EX=1.                                                         TRL 0045
          REFS=0.                                                       TRL 0046
          TDFS=0.                                                       TRL 0047
      ELSE                                                              TRL 0048
          OMEGA=TSCAT/TAU                                               TRL 0049
          GAMMA3=BETABS(CSZEN,G)                                        TRL 0050
          GAMMA4=1.-GAMMA3                                              TRL 0051
          X=S0DEP*DEPRAT                                                TRL 0052
          IF(TAU.GT.X)THEN                                              TRL 0053
C                                                                       TRL 0054
C             ONLY WARN USER IF NOT PREVIOUSLY WARNED                   TRL 0055
C             AND TAU EXCEEDS X BY MORE THAN 5%.                        TRL 0056
              IF(LWARN .AND. TAU.GT.1.05*X)THEN                         TRL 0057
                  IF(NPR.EQ.-1)WRITE(IPR1,'(/2A,1PE12.4,A,              TRL 0058
     1              /10X,2A,E12.4,A,/10X,2A,/(10X,A,0PF7.4,A))')        TRL 0059
     2              ' WARNING:  THE ESTIMATED LAYER OPTICAL',           TRL 0060
     3              ' DEPTH TOWARDS THE SUN IS',X,' AND',               TRL 0061
     4              ' THE VERTICAL LAYER OPTICAL DEPTH',                TRL 0062
     5              ' IS',TAU,'.  THE LAYER OPTICAL',                   TRL 0063
     6              ' DEPTH TOWARDS THE SUN WAS RESET',                 TRL 0064
     7              ' TO THE VERTICAL LAYER EXTINCTION',                TRL 0065
     8              ' OPTICAL DEPTH OVER THE COSINE (=',                TRL 0066
     9              CSZEN,') OF THE SOLAR ZENITH.',                     TRL 0067
     &              ' ***  THIS WARNING WILL NOT BE REPEATED  ***'      TRL 0068
                  LWARN=.FALSE.                                         TRL 0069
              ENDIF                                                     TRL 0070
              MU0=ABS(CSZEN)+1.E-10                                     TRL 0071
              X=TAU/MU0                                                 TRL 0072
          ELSE                                                          TRL 0073
              MU0=TAU/X                                                 TRL 0074
          ENDIF                                                         TRL 0075
          ONEMG2=(1.-G)*(1.+G)                                          TRL 0076
          DENOM=ONEMG2*(1.-MU0)+MU0                                     TRL 0077
          GAMMA1=((1.-OMEGA)+.75*ONEMG2*(1.-G*OMEGA)                    TRL 0078
     1      +GAMMA3*OMEGA*(1.-ONEMG2))/DENOM                            TRL 0079
          G1MG2=(1.+ONEMG2)*(1.-OMEGA)/DENOM                            TRL 0080
          A1=GAMMA1-GAMMA3*G1MG2                                        TRL 0081
          A2=GAMMA1-GAMMA4*G1MG2                                        TRL 0082
          K=SQRT(G1MG2*(2*GAMMA1-G1MG2))                                TRL 0083
          KT=K*TAU                                                      TRL 0084
          KM=K*MU0                                                      TRL 0085
          ONEPKM=1.+KM                                                  TRL 0086
          IF(X+KT.LT.CUTOFF)THEN                                        TRL 0087
C                                                                       TRL 0088
C             CASE 2:  (K TAU) AND (TAU/MU0) SMALL                      TRL 0089
C                      [(K + 1/MU0) TAU < CUTOFF]                       TRL 0090
              G3M=GAMMA3/MU0                                            TRL 0091
              G4M=GAMMA4/MU0                                            TRL 0092
              XHALF=.5*X                                                TRL 0093
              TERM1=ONEPKM+KM                                           TRL 0094
              EX=1.-X*(1.-XHALF)                                        TRL 0095
              COEF=TSCAT/(1.-TAU*(K-GAMMA1))                            TRL 0096
              REFS=COEF*(G3M-XHALF*(G3M*TERM1-A2))                      TRL 0097
              TDFS=COEF*(G4M-XHALF*(G4M*TERM1-A1))                      TRL 0098
CNEXT                                                                   TRL 0099
CNEXT         FOR GREATER ACCURACY, ADD THE NEXT HIGHER                 TRL 0100
CNEXT         ORDER TERM BY REPLACING THE EXPRESSIONS ABOVE             TRL 0101
CNEXT         FOR EX, COEF, REFS AND TDFS WITH THOSE BELOW.             TRL 0102
CNEXT         XTHIRD=X/3.                                               TRL 0103
CNEXT         TERM2=1.+KM*(3.+4*KM)                                     TRL 0104
CNEXT         TERM3=TERM1+KM                                            TRL 0105
CNEXT         EX=1.-X*(1.-XHALF*(1.-XTHIRD))                            TRL 0106
CNEXT         COEF=TSCAT/(1.-TAU*(K-GAMMA1)*(1.-KT))                    TRL 0107
CNEXT         REFS=COEF*(G3M-XHALF*(G3M*TERM1-A2                        TRL 0108
CNEXT1          -XTHIRD*(G3M*TERM2-A2*TERM3)))                          TRL 0109
CNEXT         TDFS=COEF*(G4M-XHALF*(G4M*TERM1-A1                        TRL 0110
CNEXT1          -XTHIRD*(G4M*TERM2-A1*(1.+TERM3))))                     TRL 0111
          ELSEIF(KT.LT.PT3CUT)THEN                                      TRL 0112
C                                                                       TRL 0113
C             CASE 3:  (K TAU) SMALL.                                   TRL 0114
C                      [K TAU < .3 CUTOFF]                              TRL 0115
              COEF=.5*OMEGA/(1.-TAU*(K-GAMMA1)*(1.-KT))                 TRL 0116
              EX=EXP(-X)                                                TRL 0117
              TERM1=TAU*(1.-.5*KT*(1.-KT/3.))                           TRL 0118
              EKT=1.-K*TERM1                                            TRL 0119
              TERM2=(1.+EKT)*TERM1                                      TRL 0120
              TERM3=1.+EKT*EKT                                          TRL 0121
              TWOEKT=2*EKT                                              TRL 0122
              RATP=(1.-EKT*EX)/ONEPKM                                   TRL 0123
              ONEMKM=1.-KM                                              TRL 0124
              RATM=(EKT-EX)/ONEMKM                                      TRL 0125
              DENOM=ONEPKM*ONEMKM                                       TRL 0126
              REFS=COEF*(GAMMA3*(RATP+EKT*RATM)                         TRL 0127
     1          +A2*(MU0*(TWOEKT*EX-TERM3)+TERM2)/DENOM)                TRL 0128
              TDFS=COEF*(GAMMA4*(RATM+EKT*RATP)                         TRL 0129
     1          +A1*(MU0*(TWOEKT-EX*TERM3)-EX*TERM2)/DENOM)             TRL 0130
          ELSE                                                          TRL 0131
C                                                                       TRL 0132
C             CASE 4:  (K TAU) NOT SMALL.                               TRL 0133
C                      [K TAU > .3 CUTOFF]                              TRL 0134
              EKT=EXP(-KT)                                              TRL 0135
              EX=EXP(-X)                                                TRL 0136
              RATP=(1.-EKT*EX)/ONEPKM                                   TRL 0137
              XMKT=X-KT                                                 TRL 0138
              IF(ABS(XMKT).LT..004)THEN                                 TRL 0139
C                                                                       TRL 0140
C                 CASE 4A:  (K - 1/MU0) TAU NEAR ZERO.                  TRL 0141
                  RATM=X*EX*(1.+.5*XMKT*(1.+XMKT/3.))                   TRL 0142
              ELSE                                                      TRL 0143
C                                                                       TRL 0144
C                 CASE 4B:  (K - 1/MU0) TAU NOT NEAR ZERO.              TRL 0145
                  RATM=(EKT-EX)/(1.-KM)                                 TRL 0146
              ENDIF                                                     TRL 0147
              COEF=OMEGA/(K+GAMMA1+(K-GAMMA1)*EKT**2)                   TRL 0148
              G3K=GAMMA3*K                                              TRL 0149
              REFS=COEF*((G3K+A2)*RATP+(G3K-A2)*RATM*EKT)               TRL 0150
              G4K=GAMMA4*K                                              TRL 0151
              TDFS=COEF*((G4K+A1)*RATM+(G4K-A1)*RATP*EKT)               TRL 0152
          ENDIF                                                         TRL 0153
      ENDIF                                                             TRL 0154
      RETURN                                                            TRL 0155
      END                                                               TRL 0156
