      SUBROUTINE RNSCAT(V,R,TT,PHASE,DIST,CSSA,ASYMR)                   RNS 0001
CC    SUBROUTINE RNSCAT(V,R,TT,PHASE,DIST,IK,CSSA,ASYMR,IENT)           RNS 0002
C********************************************************************** RNS 0003
      INTEGER PHASE,DIST                                                RNS 0004
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          RNS 0005
      DIMENSION SC(3,4)                                                 RNS 0006
C                                                                       RNS 0007
C          ARGUMENTS:                                                   RNS 0008
C                                                                       RNS 0009
C          F = FREQUENCY (GHZ)                                          RNS 0010
C          R = RAINFALL RATE (MM/HR)                                    RNS 0011
C          TCEL = TEMPERATURE (DEGREES CELSIUS)                         RNS 0012
C          PHASE = PHASE PARAMETER (1=WATER, 2=ICE)                     RNS 0013
C          DIST = DROP SIZE DISTRIBUTION PARAMETER                      RNS 0014
C                     (1=MARSHALL-PALMER, 2=BEST)                       RNS 0015
C                                                                       RNS 0016
C          RESULTS:                                                     RNS 0017
C                                                                       RNS 0018
C          SC(1) = ABSORPTION COEFFICIENT (1/KM)                        RNS 0019
C          SC(2) = EXTINCTION COEFFICIENT (1/KM)                        RNS 0020
C          SC(I),I=3,NSC = LEGENDRE COEFFICIENTS #I-3  (NSC=10)         RNS 0021
C          ERR = ERROR RETURN CODE: 0=NO ERROR, 1=BAD FREQUENCY,        RNS 0022
C                2=BAD RAINFALL RATE, 3=BAD TEMPERATURE,                RNS 0023
C                4=BAD PHASE PARAMETER, 5=BAD DROP SIZE DISTRIBUTION    RNS 0024
C                                                                       RNS 0025
C          THE INTERNAL DATA:                                           RNS 0026
C                                                                       RNS 0027
      DIMENSION FR(9),TEMP(9)                                           RNS 0028
C                                                                       RNS 0029
C          FR(I),I=1,NF = TABULATED FREQUENCIES (GHZ)  (NF=9)           RNS 0030
C          TEMP(I),I=1,NT = TABULATED TEMPERATURES  (NT=3)              RNS 0031
C                                                                       RNS 0032
C          THE BLOCK-DATA SECTION                                       RNS 0033
C                                                                       RNS 0034
      DATA RMIN,RMAX/0.,50./,NF/9/,NT/3/,NSC/4/,MAXI/3/                 RNS 0035
      DATA TK/273.15/,CMT0/1.0/,C7500/0.5/,G0/0.0/,G7500/0.85/          RNS 0036
      DATA (TEMP(I),I=1,3)/-10.,0.,10./                                 RNS 0037
      DATA (FR(I),I=1,9)/19.35,37.,50.3,89.5,100.,118.,130.,183.,231./  RNS 0038
C                                                                       RNS 0039
C      THIS SUBROUTINE REQUIRES FREQUENCIES IN GHZ                      RNS 0040
      NOPR = 0                                                          RNS 0041
CC    IF(IK .EQ. 1) NOPR = 1                                            RNS 0042
CC    IF(IENT .GT.1) NOPR = 0                                           RNS 0043
      F= V *29.97925                                                    RNS 0044
      FSAV=F                                                            RNS 0045
      RSAV=R                                                            RNS 0046
      INT=0                                                             RNS 0047
C      CONVERT TEMP TO DEGREES CELSIUS                                  RNS 0048
      TCEL=TT-TK                                                        RNS 0049
      TSAV=TCEL                                                         RNS 0050
C      FREQ RANGE OF DATA 19.35-231 GHZ IF LESS THAN 19.35              RNS 0051
C      SET UP PARAMETERS FOR INTERPOLATION                              RNS 0052
      IF(F.LT.FR(1)) THEN                                               RNS 0053
        FL=0.0                                                          RNS 0054
        FM=FR(1)                                                        RNS 0055
        INT=1                                                           RNS 0056
        IF(NOPR .GT. 0) WRITE (IPR,801)                                 RNS 0057
      END IF                                                            RNS 0058
C      IF MORE THAN 231 GHZ SET UP PARAMETERS FOR EXTRAPOLATION         RNS 0059
       IF(F.GT.FR(NF)) THEN                                             RNS 0060
         FL=FR(NF)                                                      RNS 0061
         FM=7500.                                                       RNS 0062
         INT=2                                                          RNS 0063
         IF(NOPR .GT. 0) WRITE (IPR,801)                                RNS 0064
      END IF                                                            RNS 0065
C      TEMP RANGE OF DATA IS -10 TO +10 DEGREES CELCIUS                 RNS 0066
C      IF BELOW OR ABOVE EXTREME SET AND DO CALCULATIONS AT EXTREME     RNS 0067
      IF (TCEL.LT.TEMP(1)) THEN                                         RNS 0068
        TCEL=TEMP(1)                                                    RNS 0069
        IF(NOPR .GT. 0) WRITE (IPR,802)                                 RNS 0070
      END IF                                                            RNS 0071
C                                                                       RNS 0072
      IF (TCEL.GT.TEMP(3)) THEN                                         RNS 0073
        TCEL=TEMP(3)                                                    RNS 0074
        IF(NOPR .GT. 0) WRITE (IPR,802)                                 RNS 0075
      END IF                                                            RNS 0076
C                                                                       RNS 0077
C      RAIN RATE OF DATA IS FOR 0-50 MM/HR                              RNS 0078
C      IF GT 50 TREAT CALCULATIONS AS IF 50 MM/HR WAS INPUT             RNS 0079
      IF(R.GT.50) THEN                                                  RNS 0080
        R=50.                                                           RNS 0081
        IF(NOPR .GT. 0) WRITE (IPR,803)                                 RNS 0082
      END IF                                                            RNS 0083
C                                                                       RNS 0084
      KI=1                                                              RNS 0085
C             FIGURE OUT THE SECOND INDEX                               RNS 0086
   10 J=PHASE+2*DIST                                                    RNS 0087
C                                                                       RNS 0088
C                                                                       RNS 0089
C             GET THE TEMPERATURE INTERPOLATION PARAMETER ST            RNS 0090
C             IF NEEDED AND AMEND THE SECOND INDEX                      RNS 0091
      CALL BS(J,TCEL,TEMP,NT,ST)                                        RNS 0092
C                                                                       RNS 0093
C             FIGURE OUT THE THIRD INDEX AND THE FREQUENCY INTERPOLATIONRNS 0094
C             PARAMETER SF                                              RNS 0095
      CALL BS(K,F,FR,NF,SF)                                             RNS 0096
C                                                                       RNS 0097
C             INITIALIZE SC                                             RNS 0098
      DO 11 I=1,NSC                                                     RNS 0099
      SC(KI,I)=0.                                                       RNS 0100
   11 CONTINUE                                                          RNS 0101
      SC(KI,3)=1.                                                       RNS 0102
C                                                                       RNS 0103
C             NOW DO THE CALCULATIONS                                   RNS 0104
C                                                                       RNS 0105
C             THE WATER CONTENT IS                                      RNS 0106
      IF(DIST.EQ.1) THEN                                                RNS 0107
      WC=.0889*R**.84                                                   RNS 0108
      ELSE                                                              RNS 0109
      WC=.067*R**.846                                                   RNS 0110
      END IF                                                            RNS 0111
C                                                                       RNS 0112
C             FOR A TEMPERATURE DEPENDENT CASE, I.E.                    RNS 0113
      IF(J.LT.3) THEN                                                   RNS 0114
      S1=(1.-SF)*(1.-ST)                                                RNS 0115
      S2=(1.-SF)*ST                                                     RNS 0116
      S3=SF*(1.-ST)                                                     RNS 0117
      S4=SF*ST                                                          RNS 0118
      DO 14 I=1,MAXI                                                    RNS 0119
      IF(I.LE.2) THEN                                                   RNS 0120
      ISC=I                                                             RNS 0121
      ELSE                                                              RNS 0122
      ISC=I+1                                                           RNS 0123
      END IF                                                            RNS 0124
      SC(KI,ISC)=S1*TAB(I,J,K,WC)+S2*TAB(I,J+1,K,WC)+                   RNS 0125
     *             S3*TAB(I,J,K+1,WC)+S4*TAB(I,J+1,K+1,WC)              RNS 0126
   14 CONTINUE                                                          RNS 0127
C                                                                       RNS 0128
C             FOR A TEMPERATURE INDEPENDENT CASE                        RNS 0129
      ELSE                                                              RNS 0130
      S1=1.-SF                                                          RNS 0131
      S2=SF                                                             RNS 0132
      DO 17 I=1,MAXI                                                    RNS 0133
      IF(I.LE.2) THEN                                                   RNS 0134
      ISC=I                                                             RNS 0135
      ELSE                                                              RNS 0136
      ISC=I+1                                                           RNS 0137
      END IF                                                            RNS 0138
      SC(KI,ISC)=S1*TAB(I,J,K,WC)+S2*TAB(I,J,K+1,WC)                    RNS 0139
   17 CONTINUE                                                          RNS 0140
      END IF                                                            RNS 0141
      F=FSAV                                                            RNS 0142
      IF(INT.EQ.3) GO TO 20                                             RNS 0143
      IF(INT.EQ.4) GO TO 30                                             RNS 0144
      IF(INT.EQ.0) THEN                                                 RNS 0145
        CSSA=SC(KI,1)/SC(KI,2)                                          RNS 0146
        IF(CSSA.GT.1.0) CSSA=1.0                                        RNS 0147
        ASYMR=SC(KI,4)/3.0                                              RNS 0148
        F=FSAV                                                          RNS 0149
        R=RSAV                                                          RNS 0150
        TCEL=TSAV                                                       RNS 0151
      RETURN                                                            RNS 0152
      END IF                                                            RNS 0153
      IF(INT.EQ.1) THEN                                                 RNS 0154
        INT=3                                                           RNS 0155
        F=FM                                                            RNS 0156
        KI=2                                                            RNS 0157
      END IF                                                            RNS 0158
      IF(INT.EQ.2) THEN                                                 RNS 0159
        INT=4                                                           RNS 0160
        F=FL                                                            RNS 0161
        KI=3                                                            RNS 0162
      END IF                                                            RNS 0163
      GO TO 10                                                          RNS 0164
   20 CONTINUE                                                          RNS 0165
      FDIF=FM-F                                                         RNS 0166
      FTOT=FM-FL                                                        RNS 0167
      CM=SC(KI,1)/SC(KI,2)                                              RNS 0168
      IF(CM.GT.1.0) CM=1.0                                              RNS 0169
      CL=CMT0                                                           RNS 0170
      AM=SC(KI,4)/3.0                                                   RNS 0171
      AL=G0                                                             RNS 0172
      GO TO 40                                                          RNS 0173
   30 CONTINUE                                                          RNS 0174
      FDIF=FM-F                                                         RNS 0175
      FTOT=FM-FL                                                        RNS 0176
      CM=C7500                                                          RNS 0177
      CL=SC(KI,1)/SC(KI,2)                                              RNS 0178
      IF(CL.GT.1.0) CL=1.0                                              RNS 0179
      AM=G7500                                                          RNS 0180
      AL=SC(KI,4)/3.0                                                   RNS 0181
   40 CTOT=CM-CL                                                        RNS 0182
      CAMT=FDIF*CTOT/FTOT                                               RNS 0183
      CSSA=CM-CAMT                                                      RNS 0184
      ATOT=AM-AL                                                        RNS 0185
      AAMT=FDIF*ATOT/FTOT                                               RNS 0186
      ASYMR=AM-AAMT                                                     RNS 0187
      F=FSAV                                                            RNS 0188
      R=RSAV                                                            RNS 0189
      TCEL=TSAV                                                         RNS 0190
      RETURN                                                            RNS 0191
801   FORMAT(2X,'***  THE ASYMMETRY PARAMETER DUE TO RAIN IS BASED ON', RNS 0192
     1 'DATA BETWEEN 19 AND 231 GHZ',                                   RNS 0193
     2 /2X,'***  EXTRAPOLATION IS USED FOR FREQUENCIES LOWER AND',      RNS 0194
     3 'HIGHER THAN THIS RANGE')                                        RNS 0195
802   FORMAT(2X,'***  TEMPERATURE RANGE OF DATA IS -10 TO +10 ',        RNS 0196
     1'DEGREES CELSIUS',/2X,'***  BEYOND THESE VALUES IT IS ',          RNS 0197
     2'TREATED AS IF AT THE EXTREMES')                                  RNS 0198
803   FORMAT(2X,'***  RAIN RATES BETWEEN 0 AND 50 MM/HR ARE',           RNS 0199
     1'WITHIN THIS DATA RANGE',/2X,'***  ABOVE THAT THE ASYMMETRY',     RNS 0200
     2' PARAMETER IS CALCULATED FOR 50 MM/HR')                          RNS 0201
      END                                                               RNS 0202
