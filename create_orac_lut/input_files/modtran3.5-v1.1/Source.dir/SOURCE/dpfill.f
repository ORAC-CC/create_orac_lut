      SUBROUTINE DPFILL(HA,HB,JNEXT)                                    FIL 0001
C********************************************************************   FIL 0002
C     THIS SUBROUTINE DEFINES THE ATMOSPHERIC BOUNDARIES OF THE PATH    FIL 0003
C     FROM HA TO HB AND INTERPOLATES (EXTRAPOLATES) THE DENSITIES TO    FIL 0004
C     THESE BOUNDARIES ASSUMING THE DENSITIES VARY EXPONENTIALLY        FIL 0005
C     WITH HEIGHT                                                       FIL 0006
C********************************************************************   FIL 0007
      INTEGER IRD,IPR,IPU,NPR,IPR1,I,J,JNEXT,IA,IB,I2,K,I1              FIL 0008
      LOGICAL LSMALL                                                    FIL 0009
      DOUBLE PRECISION ZP,PP,TP,RFNDXP,                                 FIL 0010
     $     SP,PPSUM,TPSUM,RHOPSM,DENP,AMTP                              FIL 0011
      INCLUDE 'PARAM.LST'                                               FIL 0012
      DOUBLE PRECISION HA,HB,Z(LAYDIM),P(LAYDIM),T(LAYDIM),             FIL 0013
     1  DRFNDX(LAYDIM),DDNSTY(KMAX,LAYDIM),A                            FIL 0014
C                                                                       FIL 0015
C                                                                       FIL 0016
C                                                                       FIL 0017
C     CONVENTION                                                        FIL 0018
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         FIL 0019
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            FIL 0020
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             FIL 0021
C                                                                       FIL 0022
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              FIL 0023
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       FIL 0024
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     FIL 0025
C                                                                       FIL 0026
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           FIL 0027
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      FIL 0028
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        FIL 0029
C                                                                       FIL 0030
      DOUBLE PRECISION DI1X, DI2X, DJNXTX                               FIL 0031
C                                                                       FIL 0032
      COMMON /MODELX/ DNSTYX(MMOLX,LAYDIM)                              FIL 0033
      COMMON /RFRPTX/ DENPX(MMOLX,LAYDIM+1),AMTPX(MMOLX,LAYDIM+1)       FIL 0034
C                                                                       FIL 0035
C                                                                       FIL 0036
C                                                                       FIL 0037
C     SSI COMMENTS ON DOUBLE PRECISION VARIABLES:                       FIL 0038
C     RFRPTH IS THE OLD COMMON BLOCK IN SINGLE PRECISION.               FIL 0039
C     DPRFRP IS THE SAME COMMON BLOCK IN DOUBLE PRECISION; IT IS NEW.   FIL 0040
C     IN THIS ROUTINE SP IS USED AS A PREFIX TO DENOTE THE              FIL 0041
C     SINGLE PRECISION VARIABLES OF RFRPTH.                             FIL 0042
C     THE FOLLOWING ARE THE EXCEPTIONS:                                 FIL 0043
C     SPRFN  STANDS FOR THE OLD SINGLE PRECISION RFNDXP                 FIL 0044
C     SPPPSU STANDS FOR THE OLD SINGLE PRECISION PPSUM                  FIL 0045
C     SPTPSU STANDS FOR THE OLD SINGLE PRECISION TPSUM                  FIL 0046
C     SPRHOP STANDS FOR THE OLD SINGLE PRECISION RHOPSM                 FIL 0047
C     THE VARIABLES OF THE DOUBLE PRECISION BLOCK DPRFRP HAVE THE       FIL 0048
C     SAME NAMES AS THOSE OF THE ORIGINAL SINGLE PRECISION BLOCK;       FIL 0049
C     THAT IS, WITHOUT ANY PREFIXES.                                    FIL 0050
C                                                                       FIL 0051
C     VARIABLES IN THE COMMON BLOCK /MODEL/ ARE TEMPORARILY STORED      FIL 0052
C     IN DOUBLE PRECISION ARRAYS FOR DPFILL CALCULATIONS:               FIL 0053
C       ZM => Z, PM => P, TM => T, RFNDX => DRFNDX & DENSTY => DDNSTY   FIL 0054
C                                                                       FIL 0055
C     SOME OTHER VARIABLES WERE DECLARED IN DOUBLE PRECISION.           FIL 0056
C     THEIR SP COUNTERPARTS WERE PREFIXED WITH "SP".                    FIL 0057
C     SINCE THESE DO NOT INVOLVE COMMON BLOCKS NOTHING MORE ABOUT       FIL 0058
C     THEM IS SAID.  SEE THE DECLARATIONS ABOVE FOR THE SPECIFIC        FIL 0059
C     VARIABLES.                                                        FIL 0060
C                                                                       FIL 0061
      COMMON /IFIL/ IRD,IPR,IPU,NPR,IPR1,ISCRCH                         FIL 0062
      REAL ZM,PM,TM,RFNDX,DENSTY                                        FIL 0063
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    FIL 0064
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               FIL 0065
      REAL RE,ZMAX                                                      FIL 0066
      INTEGER IMAX,IMOD,IPATH                                           FIL 0067
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             FIL 0068
      COMMON /RFRPTH/ SPZP(LAYDIM+1),SPPP(LAYDIM+1),SPTP(LAYDIM+1),     FIL 0069
     $ SPRFN(LAYDIM+1),SPSP(LAYDIM+1),SPPPSU(LAYDIM+1),SPTPSU(LAYDIM+1) FIL 0070
     $ ,SPRHOP(LAYDIM+1),SPDENP(KMAX,LAYDIM+1),SPAMTP(KMAX,LAYDIM+1)    FIL 0071
      COMMON /DPRFRP/ ZP(LAYDIM+1),PP(LAYDIM+1),TP(LAYDIM+1),           FIL 0072
     $ RFNDXP(LAYDIM+1),SP(LAYDIM+1), PPSUM(LAYDIM+1),TPSUM(LAYDIM+1),  FIL 0073
     $ RHOPSM(LAYDIM+1),DENP(KMAX,LAYDIM+1), AMTP(KMAX,LAYDIM+1)        FIL 0074
      COMMON /SMALL2/LSMALL                                             FIL 0075
      DO 9 I=1, LAYDIM                                                  FIL 0076
         Z(I)=DBLE(ZM(I))                                               FIL 0077
         P(I)=DBLE(PM(I))                                               FIL 0078
         T(I)=DBLE(TM(I))                                               FIL 0079
         DRFNDX(I)=DBLE(RFNDX(I))                                       FIL 0080
         DO 9999 J=1, KMAX                                              FIL 0081
            DDNSTY(J,I)=DBLE(DENSTY(J,I))                               FIL 0082
 9999    CONTINUE                                                       FIL 0083
C        FOR CFC'S, WE ALWAYS USE SINGLE PRECISION VARIABLES, SO NO     FIL 0084
C        NEED FOR THE ABOVE LOOP FOR CFC'S.                             FIL 0085
C                                                                       FIL 0086
 9    CONTINUE                                                          FIL 0087
C                                                                       FIL 0088
      DO 99 I=1, LAYDIM+1                                               FIL 0089
         ZP(I)=DBLE(SPZP(I))                                            FIL 0090
         PP(I)=DBLE(SPPP(I))                                            FIL 0091
         TP(I)=DBLE(SPTP(I))                                            FIL 0092
         RFNDXP(I)=DBLE(SPRFN(I))                                       FIL 0093
         DO 999 J=1, KMAX                                               FIL 0094
            DENP(J,I)=DBLE(SPDENP(J,I))                                 FIL 0095
 999     CONTINUE                                                       FIL 0096
C        FOR CFC'S, WE ALWAYS USE SINGLE PRECISION VARIABLES, SO NO     FIL 0097
C        NEED FOR THE ABOVE LOOP FOR CFC'S.                             FIL 0098
C                                                                       FIL 0099
 99   CONTINUE                                                          FIL 0100
C                                                                       FIL 0101
C     IF(HA.GE.HB .AND. LSMALL .NE. .TRUE.) THEN                        FIL 0102
      IF(HA.GE.HB .AND.    (.NOT. LSMALL )) THEN                        FIL 0103
         WRITE(IPR,22) HA,HB,JNEXT                                      FIL 0104
 22      FORMAT('0SUBROUTINE DPFILL- ERROR, HA .GE. HB',//,             FIL 0105
     $        10X,'HA, HB, JNEXT = ',2E25.15,I6)                        FIL 0106
         STOP                                                           FIL 0107
      ENDIF                                                             FIL 0108
C***  FIND Z(IA):  THE SMALLEST Z(I).GT.HA                              FIL 0109
      DO 100 I=1,IMAX                                                   FIL 0110
         IF(HA.LT.Z(I)) THEN                                            FIL 0111
            IA=I                                                        FIL 0112
            GO TO 110                                                   FIL 0113
         ENDIF                                                          FIL 0114
 100  CONTINUE                                                          FIL 0115
      IA=IMAX+1                                                         FIL 0116
      IB=IA                                                             FIL 0117
      GO TO 130                                                         FIL 0118
C***  FIND Z(IB):  THE SMALLEST Z(I).GE.HB                              FIL 0119
 110  CONTINUE                                                          FIL 0120
      DO 120 I=IA,IMAX                                                  FIL 0121
         IF(HB-Z(I).LE. .0001) THEN                                     FIL 0122
            IB=I                                                        FIL 0123
            GO TO 130                                                   FIL 0124
         ENDIF                                                          FIL 0125
 120  CONTINUE                                                          FIL 0126
      IB=IMAX+1                                                         FIL 0127
 130  CONTINUE                                                          FIL 0128
C***  INTERPOLATE DENSITIES TO HA,HB                                    FIL 0129
      ZP(JNEXT)=HA                                                      FIL 0130
      I2=IA                                                             FIL 0131
      IF(I2.EQ.1) I2=2                                                  FIL 0132
      IF(I2.GT.IMAX) I2=IMAX                                            FIL 0133
      I1=I2-1                                                           FIL 0134
      A=(HA-Z(I1))/(Z(I2)-Z(I1))                                        FIL 0135
      CALL DPEXNT(PP(JNEXT),P(I1),P(I2),A)                              FIL 0136
      TP(JNEXT)=T(I1)+(T(I2)-T(I1))*A                                   FIL 0137
      CALL DPEXNT(RFNDXP(JNEXT),DRFNDX(I1),DRFNDX(I2),A)                FIL 0138
      DO 140 K=1,KMAX                                                   FIL 0139
          IF(K.EQ.66 .OR. K.EQ.67 .OR. K.EQ.16)THEN                     FIL 0140
C                                                                       FIL 0141
C             LINEARLY INTERPOLATE CLOUD DENSITIES                      FIL 0142
              DENP(K,JNEXT)=DDNSTY(K,I1)+A*(DDNSTY(K,I2)-DDNSTY(K,I1))  FIL 0143
          ELSE                                                          FIL 0144
C                                                                       FIL 0145
C             EXPONENTIALLY INTERPOLATE                                 FIL 0146
              CALL DPEXNT(DENP(K,JNEXT),DDNSTY(K,I1),DDNSTY(K,I2),A)    FIL 0147
          ENDIF                                                         FIL 0148
 140  CONTINUE                                                          FIL 0149
C                                                                       FIL 0150
      DO 145 KX=1, NSPECX                                               FIL 0151
C        STORE SINGLE PRECISION VARIABLES IN DOUBLE PRECISION.          FIL 0152
         DI1X=DBLE(DNSTYX(KX,I1))                                       FIL 0153
         DI2X=DBLE(DNSTYX(KX,I2))                                       FIL 0154
         CALL DPEXNT(DJNXTX,DI1X,DI2X,A)                                FIL 0155
C        BACK TO SINGLE PRECISION                                       FIL 0156
         DENPX(KX,JNEXT)=REAL(DJNXTX)                                   FIL 0157
 145  CONTINUE                                                          FIL 0158
C                                                                       FIL 0159
      IF(IA.NE.IB) THEN                                                 FIL 0160
C***     FILL IN DENSITIES BETWEEN HA AND HB                            FIL 0161
         I1=IA                                                          FIL 0162
         I2=IB-1                                                        FIL 0163
         DO 151 I=I1,I2                                                 FIL 0164
            JNEXT=JNEXT+1                                               FIL 0165
            ZP(JNEXT)=Z(I)                                              FIL 0166
            PP(JNEXT)=P(I)                                              FIL 0167
            TP(JNEXT)=T(I)                                              FIL 0168
            RFNDXP(JNEXT)=DRFNDX(I)                                     FIL 0169
            DO 150 K=1,KMAX                                             FIL 0170
               DENP(K,JNEXT)=DDNSTY(K,I)                                FIL 0171
 150        CONTINUE                                                    FIL 0172
C                                                                       FIL 0173
            DO 155 KX=1,NSPECX                                          FIL 0174
               DENPX(KX,JNEXT)=DNSTYX(KX,I)                             FIL 0175
 155        CONTINUE                                                    FIL 0176
C                                                                       FIL 0177
 151     CONTINUE                                                       FIL 0178
      ENDIF                                                             FIL 0179
C***  INTERPOLATE THE DENSITIES TO HB                                   FIL 0180
      JNEXT=JNEXT+1                                                     FIL 0181
      ZP(JNEXT)=HB                                                      FIL 0182
      I2=IB                                                             FIL 0183
      IF(I2.EQ.1) I2=2                                                  FIL 0184
      IF(I2.GT.IMAX) I2=IMAX                                            FIL 0185
      I1=I2-1                                                           FIL 0186
      A=(HB-Z(I1))/(Z(I2)-Z(I1))                                        FIL 0187
      CALL DPEXNT(PP(JNEXT),P(I1),P(I2),A)                              FIL 0188
      TP(JNEXT)=T(I1)+(T(I2)-T(I1))*A                                   FIL 0189
      CALL DPEXNT(RFNDXP(JNEXT),DRFNDX(I1),DRFNDX(I2),A)                FIL 0190
      DO 170 K=1,KMAX                                                   FIL 0191
          IF(K.EQ.66 .OR. K.EQ.67 .OR. K.EQ.16)THEN                     FIL 0192
C                                                                       FIL 0193
C             LINEARLY INTERPOLATE CLOUD DENSITIES                      FIL 0194
              DENP(K,JNEXT)=DDNSTY(K,I1)+A*(DDNSTY(K,I2)-DDNSTY(K,I1))  FIL 0195
          ELSE                                                          FIL 0196
C                                                                       FIL 0197
C             EXPONENTIALLY INTERPOLATE                                 FIL 0198
              CALL DPEXNT(DENP(K,JNEXT),DDNSTY(K,I1),DDNSTY(K,I2),A)    FIL 0199
          ENDIF                                                         FIL 0200
 170  CONTINUE                                                          FIL 0201
C                                                                       FIL 0202
      DO 175 KX=1, NSPECX                                               FIL 0203
C        STORE SINGLE PRECISION VARIABLES IN DOUBLE PRECISION.          FIL 0204
         DI1X=DBLE(DNSTYX(KX,I1))                                       FIL 0205
         DI2X=DBLE(DNSTYX(KX,I2))                                       FIL 0206
         CALL DPEXNT(DJNXTX, DI1X, DI2X, A)                             FIL 0207
C        BACK TO SINGLE PRECISION                                       FIL 0208
         DENPX(KX,JNEXT)=REAL(DJNXTX)                                   FIL 0209
 175  CONTINUE                                                          FIL 0210
C                                                                       FIL 0211
      DO 1 I=1, LAYDIM                                                  FIL 0212
         ZM(I)=REAL(Z(I))                                               FIL 0213
         PM(I)=REAL(P(I))                                               FIL 0214
         TM(I)=REAL(T(I))                                               FIL 0215
         RFNDX(I)=REAL(DRFNDX(I))                                       FIL 0216
         DO 1111 J=1, KMAX                                              FIL 0217
            DENSTY(J,I)=REAL(DDNSTY(J,I))                               FIL 0218
 1111    CONTINUE                                                       FIL 0219
C        NO NEED FOR A LOOP CORRESPONDING TO ABOVE FOR CFC'S.           FIL 0220
 1    CONTINUE                                                          FIL 0221
      DO 11 I=1, LAYDIM+1                                               FIL 0222
         SPZP(I)=REAL(ZP(I))                                            FIL 0223
         SPPP(I)=REAL(PP(I))                                            FIL 0224
         SPTP(I)=REAL(TP(I))                                            FIL 0225
         SPRFN(I)=REAL(RFNDXP(I))                                       FIL 0226
         DO 111 J=1, KMAX                                               FIL 0227
            SPDENP(J,I)=REAL(DENP(J,I))                                 FIL 0228
 111     CONTINUE                                                       FIL 0229
C        NO NEED FOR A LOOP CORRESPONDING TO ABOVE FOR CFC'S.           FIL 0230
 11   CONTINUE                                                          FIL 0231
      RETURN                                                            FIL 0232
      END                                                               FIL 0233
