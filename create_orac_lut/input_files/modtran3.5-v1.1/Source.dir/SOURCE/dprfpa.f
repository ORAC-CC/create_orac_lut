      SUBROUTINE DPRFPA(H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,BETA,RANGE,BENDNG)PTH 0001
C                                                                       PTH 0002
C     THIS ROUTINE, A DOUBLE PRECISION VERSION OF THE PREVIOUS ROUTINE  PTH 0003
C     "RFPATH", TRACES THE REFRACTED RAY FROM H1 WITH INITIAL ZENITH    PTH 0004
C     ANGLE "ANGLE" TO H2 WHERE THE ZENITH ANGLE IS PHI, AND CALCULATES PTH 0005
C     THE ABSORBER AMOUNTS (IF IAMT.EQ.1) ALONG THE PATH.  IT STARTS    PTH 0006
C     FROM THE LOWEST POINT ALONG THE PATH (THE TANGENT HEIGHT HMIN     PTH 0007
C     IF LEN=1 OR HA=MIN(H1,H2) IF LEN=0) AND PROCEEDS TO THE HIGHEST   PTH 0008
C     POINT.  BETA AND RANGE ARE THE EARTH CENTERED ANGLE AND THE TOTAL PTH 0009
C     DISTANCE RESPECTIVELY FOR THE REFRACTED PATH FROM H1 TO H2.       PTH 0010
C                                                                       PTH 0011
C     DECLARE INPUTS                                                    PTH 0012
      DOUBLE PRECISION H1,H2,ANGLE,PHI,HMIN,BETA,RANGE,BENDNG           PTH 0013
      INTEGER LEN,IAMT                                                  PTH 0014
C                                                                       PTH 0015
C     INCLUDE PARAMETERS                                                PTH 0016
      INCLUDE 'PARAM.LST'                                               PTH 0017
C                                                                       PTH 0018
C     LIST COMMONS                                                      PTH 0019
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               PTH 0020
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           PTH 0021
      REAL RE,ZMAX                                                      PTH 0022
      INTEGER IMAX,IMOD,IPATH                                           PTH 0023
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             PTH 0024
C                                                                       PTH 0025
C       PI       THE CONSTANT PI                                        PTH 0026
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       PTH 0027
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       PTH 0028
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         PTH 0029
      REAL PI,DEG,BIGNUM,BIGEXP                                         PTH 0030
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                PTH 0031
      REAL SPZP,SPPP,SPTP,SPRFN,SPSP,SPPPSU,SPTPSU,SPRHOP,SPDENP,SPAMTP PTH 0032
      COMMON/RFRPTH/SPZP(LAYDIM+1),SPPP(LAYDIM+1),                      PTH 0033
     1  SPTP(LAYDIM+1),SPRFN(LAYDIM+1),SPSP(LAYDIM+1),                  PTH 0034
     2  SPPPSU(LAYDIM+1),SPTPSU(LAYDIM+1),SPRHOP(LAYDIM+1),             PTH 0035
     3  SPDENP(KMAX,LAYDIM+1),SPAMTP(KMAX,LAYDIM+1)                     PTH 0036
      DOUBLE PRECISION ZP,PP,TP,RFNDXP,SP,PPSUM,TPSUM,RHOPSM,DENP,AMTP  PTH 0037
      COMMON/DPRFRP/ZP(LAYDIM+1),PP(LAYDIM+1),TP(LAYDIM+1),             PTH 0038
     1  RFNDXP(LAYDIM+1),SP(LAYDIM+1),PPSUM(LAYDIM+1),TPSUM(LAYDIM+1),  PTH 0039
     2  RHOPSM(LAYDIM+1),DENP(KMAX,LAYDIM+1),AMTP(KMAX,LAYDIM+1)        PTH 0040
      INTEGER JTURN,LJ                                                  PTH 0041
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       PTH 0042
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  PTH 0043
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      PTH 0044
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    PTH 0045
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     PTH 0046
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                PTH 0047
C                                                                       PTH 0048
C     /PATH/                                                            PTH 0049
C       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.               PTH 0050
C       AHT     ALTITUDES AT PATH BOUNDARIES.                           PTH 0051
C       TPH     TEMPERATURE AT PATH BOUNDARIES.                         PTH 0052
C       IMAP    MAPPING FROM PATH SEGMENTS TO LAYERS.                   PTH 0053
      INTEGER IMAP                                                      PTH 0054
      REAL QTHETA,AHT,TPH                                               PTH 0055
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),TPH(LAYTWO),IMAP(LAYTWO)   PTH 0056
      DOUBLE PRECISION DHALFR,DPRNG2                                    PTH 0057
      COMMON/SMALL1/DHALFR,DPRNG2                                       PTH 0058
      LOGICAL LSMALL                                                    PTH 0059
      COMMON/SMALL2/LSMALL                                              PTH 0060
      DOUBLE PRECISION SMH1,SMH2,TRANGE                                 PTH 0061
      COMMON/SMALL4/SMH1,SMH2,TRANGE                                    PTH 0062
      LOGICAL LPRINT                                                    PTH 0063
      COMMON/CPRINT/LPRINT                                              PTH 0064
C                                                                       PTH 0065
C     DECLARE LOCAL VARIABLES                                           PTH 0066
      INTEGER IHIGH,KA,JD,J2A,J2B,JO,IORDER,                            PTH 0067
     1  MAX,I,JNEXT,JHA,JMAX,J2,J,IHLOW,J1                              PTH 0068
      DOUBLE PRECISION DPANDX,HA,HB,SH,GAMMA,CPATH,S,SINAI,THETA,DS,    PTH 0069
     1  DBEND,DBETA,COSAI,ANGLEA,RHOBAR,TBAR,PBAR,SAVCOS,SAVSIN,RNG,DPREPTH 0070
C                                                                       PTH 0071
C     LIST DATA                                                         PTH 0072
      CHARACTER*2 HLOW(2)                                               PTH 0073
      DOUBLE PRECISION DPDEG                                            PTH 0074
      DATA HLOW/'H1','H2'/,DPDEG/57.2957795131D0/                       PTH 0075
      IF(LSMALL)THEN                                                    PTH 0076
C         IF (LSMALL) SET H1 AND H2 EQUAL TO DOUBLE PRECISION           PTH 0077
C         VALUES SAVED VIA COMMON STATEMENTS.                           PTH 0078
          H1=SMH1                                                       PTH 0079
          H2=SMH2                                                       PTH 0080
      ENDIF                                                             PTH 0081
      DPRE=DBLE(RE)                                                     PTH 0082
      DO 10 I=1,LAYDIM+1                                                PTH 0083
          ZP(I)=DBLE(SPZP(I))                                           PTH 0084
          RFNDXP(I)=DBLE(SPRFN(I))                                      PTH 0085
          PPSUM(I)=DBLE(SPPPSU(I))                                      PTH 0086
          TPSUM(I)=DBLE(SPTPSU(I))                                      PTH 0087
          RHOPSM(I)=DBLE(SPRHOP(I))                                     PTH 0088
   10 CONTINUE                                                          PTH 0089
      MAX=0                                                             PTH 0090
      IF(H1.LE.H2)THEN                                                  PTH 0091
          IORDER=1                                                      PTH 0092
          HA=H1                                                         PTH 0093
          HB=H2                                                         PTH 0094
          ANGLEA=ANGLE                                                  PTH 0095
      ELSE                                                              PTH 0096
          IORDER=-1                                                     PTH 0097
          HA=H2                                                         PTH 0098
          HB=H1                                                         PTH 0099
          ANGLEA=PHI                                                    PTH 0100
      ENDIF                                                             PTH 0101
      IF(LSMALL)THEN                                                    PTH 0102
          SAVSIN=SIN(ANGLEA/DPDEG)                                      PTH 0103
          SAVCOS=-COS(ANGLEA/DPDEG)                                     PTH 0104
          SINAI=SIN(ANGLEA/DPDEG)                                       PTH 0105
          COSAI=-COS(ANGLEA/DPDEG)                                      PTH 0106
      ENDIF                                                             PTH 0107
      JNEXT=1                                                           PTH 0108
      IF(IAMT.EQ.1 .AND. NPR.LT.1 .AND. LPRINT)WRITE(IPR,'(/A,///(3A))')PTH 0109
     1  '1CALCULATION OF THE REFRACTED PATH THROUGH THE ATMOSPHERE',    PTH 0110
     2  '  I START ALT   END ALT     THETA    DRANGE     RANGE',        PTH 0111
     3  '     DBETA      BETA       PHI     DBEND   BENDING',           PTH 0112
     4  '       PBAR    TBAR    RHOBAR',                                PTH 0113
     5  '       (KM)      (KM)       (DEG)     (KM)       (KM)',        PTH 0114
     6  '     (DEG)     (DEG)     (DEG)     (DEG)    (DEG) ',           PTH 0115
     7  '       (MB)     (K)  (GM/CM3)'                                 PTH 0116
      IF(LEN.EQ.1) THEN                                                 PTH 0117
C                                                                       PTH 0118
C         LONG PATH:  FILL IN THE SYMMETRIC PART FROM                   PTH 0119
C                     THE TANGENT HEIGHT TO HA                          PTH 0120
          CALL DPFILL(HMIN,HA,JNEXT)                                    PTH 0121
          JHA=JNEXT                                                     PTH 0122
      ENDIF                                                             PTH 0123
C                                                                       PTH 0124
C     IF LEN=0, OR IF LEN=1 TO FILL IN THE REMAINING PATH FROM HA TO HB PTH 0125
      IF(HA.NE.HB)CALL DPFILL(HA,HB,JNEXT)                              PTH 0126
      JMAX=JNEXT                                                        PTH 0127
      IPATH=JMAX                                                        PTH 0128
C                                                                       PTH 0129
C     INTEGRATE EACH SEGMENT OF THE PATH                                PTH 0130
C     CALCULATE CPATH SEPERATELY FOR LEN=0,1                            PTH 0131
      IF(LEN.EQ.0)THEN                                                  PTH 0132
          CALL DPFISH(HA,SH,GAMMA)                                      PTH 0133
          CPATH=(DPRE+HA)*DPANDX(HA,SH,GAMMA)*SIN(ANGLEA/DPDEG)         PTH 0134
      ELSE IF (LEN .EQ. 1) THEN                                         PTH 0135
          CALL DPSCHT(ZP(1),ZP(2),RFNDXP(1),RFNDXP(2),SH,GAMMA)         PTH 0136
          CPATH=(DPRE+HMIN)*DPANDX(HMIN,SH,GAMMA)                       PTH 0137
      ENDIF                                                             PTH 0138
      BETA=DBLE(0.)                                                     PTH 0139
      S=DBLE(0.)                                                        PTH 0140
      BENDNG=DBLE(0.)                                                   PTH 0141
      IF(LEN.EQ.0) THEN                                                 PTH 0142
C                                                                       PTH 0143
C         SHORT PATH                                                    PTH 0144
C           ANGLEA   ZENITH ANGLE AT HA [DEG].                          PTH 0145
C           SINAI    SINE OF THE INCIDENT ANGLE.                        PTH 0146
C           COSAI    COSINE OF THE INCIDENT ANGLE.                      PTH 0147
          JNEXT=1                                                       PTH 0148
          THETA=ANGLEA                                                  PTH 0149
          SINAI=COS((DBLE(90.)-ANGLEA)/DPDEG)                           PTH 0150
          COSAI=-SIN((DBLE(90.)-ANGLEA)/DPDEG)                          PTH 0151
          IF(SINAI.GT. DBLE(1.))SINAI= DBLE(1.)                         PTH 0152
          IF(SINAI.LT.-DBLE(1.))SINAI=-DBLE(1.)                         PTH 0153
      ELSE                                                              PTH 0154
C                                                                       PTH 0155
C         DO SYMMETRIC PART,FROM TANGENT HEIGHT(HMIN) TO HA             PTH 0156
          IHLOW=1                                                       PTH 0157
          IF(IORDER.EQ.-1)IHLOW=2                                       PTH 0158
          IF(IAMT.EQ.1 .AND. NPR.LT.1 .AND. LPRINT)WRITE(IPR,           PTH 0159
     1      '(T7,A,T17,A2,/T7,A,/)')'TANGENT',HLOW(IHLOW),'HEIGHT'      PTH 0160
          SINAI=DBLE(1.)                                                PTH 0161
          COSAI=DBLE(0.)                                                PTH 0162
          THETA=DBLE(90.)                                               PTH 0163
          J2=JHA-1                                                      PTH 0164
          RNG=DHALFR/J2                                                 PTH 0165
          DO 20 J=1,J2                                                  PTH 0166
              CALL DPSCHT(ZP(J),ZP(J+1),RFNDXP(J),RFNDXP(J+1),SH,GAMMA) PTH 0167
              CALL DPLAYR(J,SINAI,COSAI,CPATH,                          PTH 0168
     1          SH,GAMMA,IAMT,DS,DBEND,RNG)                             PTH 0169
              IF(LSMALL)THEN                                            PTH 0170
                  SINAI=SAVSIN                                          PTH 0171
                  COSAI=SAVCOS                                          PTH 0172
              ENDIF                                                     PTH 0173
              IF(SINAI.GT. DBLE(1.))SINAI= DBLE(1.)                     PTH 0174
              IF(SINAI.LT.-DBLE(1.))SINAI=-DBLE(1.)                     PTH 0175
              DBEND=DBEND*DPDEG                                         PTH 0176
              PHI=ASIN(SINAI)*DPDEG                                     PTH 0177
              DBETA=THETA-PHI+DBEND                                     PTH 0178
              IF(LSMALL)DBETA=(2*DPDEG*SAVSIN*DS)/(ZP(J)+ZP(J+1)+2*DPRE)PTH 0179
              PHI=DBLE(180.)-PHI                                        PTH 0180
              S=S+DS                                                    PTH 0181
C                                                                       PTH 0182
C             SAVE REFRACTED RAY PATH LENGTH FOR MULTIPLE SCATTERING.   PTH 0183
              BENDNG=BENDNG+DBEND                                       PTH 0184
              BETA=BETA+DBETA                                           PTH 0185
              IF(IAMT.EQ.1)THEN                                         PTH 0186
                  PBAR=PPSUM(J)/RHOPSM(J)                               PTH 0187
                  TBAR=TPSUM(J)/RHOPSM(J)                               PTH 0188
                  RHOBAR=RHOPSM(J)/DS                                   PTH 0189
                  IF(IAMT.EQ.1 .AND. NPR.LT.1 .AND. LPRINT)             PTH 0190
     1              WRITE(IPR,'(I3,10F10.5,F11.5,F8.2,1PE10.3)')        PTH 0191
     2              J,ZP(J),ZP(J+1),THETA,DS,S,DBETA,BETA,              PTH 0192
     3              PHI,DBEND,BENDNG,PBAR,TBAR,RHOBAR                   PTH 0193
                  J2A=J2-J+1                                            PTH 0194
                  J2B=J2+J                                              PTH 0195
                  QTHETA(J2A)=COS(REAL(PHI)/DEG)                        PTH 0196
                  QTHETA(J2B)=COS(REAL(THETA)/DEG)                      PTH 0197
                  MAX=J2B                                               PTH 0198
              ENDIF                                                     PTH 0199
              IF(ISSGEO.NE.1)THEN                                       PTH 0200
                  ATHETA(J)=REAL(THETA)                                 PTH 0201
                  ADBETA(J)=REAL(DBETA)                                 PTH 0202
              ENDIF                                                     PTH 0203
              THETA=DBLE(180.)-PHI                                      PTH 0204
   20     CONTINUE                                                      PTH 0205
C                                                                       PTH 0206
C         DOUBLE PATH QUANTITIES FOR OTHER PART OF THE SYMMETRIC PATH.  PTH 0207
          BENDNG=2*BENDNG                                               PTH 0208
          BETA=2*BETA                                                   PTH 0209
          S=2*S                                                         PTH 0210
          IF(IAMT.EQ.1 .AND. NPR.LT.1 .AND. LPRINT)                     PTH 0211
     1      WRITE(IPR,'(2(/T10,A),T40,F9.3,T58,F9.3,T85,F9.3,/)')       PTH 0212
     2      'DOUBLE RANGE, BETA, BENDING',                              PTH 0213
     3      'FOR SYMMETRIC PART OF PATH',S,BETA,BENDNG                  PTH 0214
          JNEXT=JHA                                                     PTH 0215
      ENDIF                                                             PTH 0216
C                                                                       PTH 0217
C     DO PATH FROM HA TO HB                                             PTH 0218
      IF(HA.NE.HB)THEN                                                  PTH 0219
          JO=MAX-JNEXT+1                                                PTH 0220
          J1=JNEXT                                                      PTH 0221
          J2=JMAX-1                                                     PTH 0222
          IF(H1.GT.H2)THEN                                              PTH 0223
              JD=J2-J1+1                                                PTH 0224
              DO 30 KA=MAX,1,-1                                         PTH 0225
                  QTHETA(KA+JD)=QTHETA(KA)                              PTH 0226
   30         CONTINUE                                                  PTH 0227
              MAX=MAX+JD                                                PTH 0228
              JO=0                                                      PTH 0229
          ENDIF                                                         PTH 0230
          IHLOW=1                                                       PTH 0231
          IF(IORDER.EQ.-1)IHLOW=2                                       PTH 0232
          IHIGH=MOD(IHLOW,2)+1                                          PTH 0233
          IF(IAMT.EQ.1 .AND. NPR.LT.1 .AND. LPRINT)                     PTH 0234
     1      WRITE(IPR,'(T11,A2,A,A2,/)') HLOW(IHLOW),' TO ',HLOW(IHIGH) PTH 0235
          DO 40 J=J1,J2                                                 PTH 0236
              CALL DPSCHT(ZP(J),ZP(J+1),RFNDXP(J),RFNDXP(J+1),SH,GAMMA) PTH 0237
              RNG=DPRNG2*((ZP(J+1)-ZP(J))/(ZP(J2+1)-ZP(J1)))            PTH 0238
              CALL DPLAYR(J,SINAI,COSAI,CPATH,                          PTH 0239
     1          SH,GAMMA,IAMT,DS,DBEND,RNG)                             PTH 0240
              IF(LSMALL)THEN                                            PTH 0241
                  SINAI=SAVSIN                                          PTH 0242
                  COSAI=SAVCOS                                          PTH 0243
              ENDIF                                                     PTH 0244
              IF(SINAI.GT. DBLE(1.))SINAI= DBLE(1.)                     PTH 0245
              IF(SINAI.LT.-DBLE(1.))SINAI=-DBLE(1.)                     PTH 0246
              DBEND=DBEND*DPDEG                                         PTH 0247
              PHI=ASIN(SINAI)*DPDEG                                     PTH 0248
              DBETA=THETA-PHI+DBEND                                     PTH 0249
              IF(LSMALL)DBETA=(DPDEG*2*SAVSIN*DS)/(ZP(J)+ZP(J+1)+2*DPRE)PTH 0250
              PHI=DBLE(180.)-PHI                                        PTH 0251
              S=S+DS                                                    PTH 0252
              BENDNG=BENDNG+DBEND                                       PTH 0253
              BETA=BETA+DBETA                                           PTH 0254
              IF(IAMT.EQ.1)THEN                                         PTH 0255
                  PBAR=PPSUM(J)/RHOPSM(J)                               PTH 0256
                  TBAR=TPSUM(J)/RHOPSM(J)                               PTH 0257
                  RHOBAR=RHOPSM(J)/DS                                   PTH 0258
                  IF(IAMT.EQ.1 .AND. NPR.LT.1 .AND. LPRINT)             PTH 0259
     1              WRITE(IPR,'(I3,10F10.5,F11.5,F8.2,1PE10.3)')        PTH 0260
     2              J,ZP(J),ZP(J+1),THETA,DS,S,DBETA,BETA,              PTH 0261
     3              PHI,DBEND,BENDNG,PBAR,TBAR,RHOBAR                   PTH 0262
              ENDIF                                                     PTH 0263
C                                                                       PTH 0264
C             SAVE LAYER REFRACTED PATH ANGLE FOR MULTIPLE SCATTERING.  PTH 0265
              IF(H2.GT.H1)THEN                                          PTH 0266
                  QTHETA(J+JO)=COS(REAL(THETA)/DEG)                     PTH 0267
                  MAX=J+JO                                              PTH 0268
              ELSE                                                      PTH 0269
                  J2B=J2-J+1                                            PTH 0270
                  QTHETA(J2B)=COS(REAL(PHI)/DEG)                        PTH 0271
              ENDIF                                                     PTH 0272
              IF(ISSGEO.NE.1)THEN                                       PTH 0273
                  ATHETA(J)=REAL(THETA)                                 PTH 0274
                  ADBETA(J)=REAL(DBETA)                                 PTH 0275
              ENDIF                                                     PTH 0276
              THETA=DBLE(180.)-PHI                                      PTH 0277
   40     CONTINUE                                                      PTH 0278
      ENDIF                                                             PTH 0279
      IF(ISSGEO.EQ.0)ATHETA(JMAX)=REAL(THETA)                           PTH 0280
      IF(IORDER.EQ.-1)PHI=ANGLEA                                        PTH 0281
      RANGE=S                                                           PTH 0282
      DO 50 I=1,LAYDIM+1                                                PTH 0283
          SPZP(I)=REAL(ZP(I))                                           PTH 0284
          SPRFN(I)=REAL(RFNDXP(I))                                      PTH 0285
          SPPPSU(I)=REAL(PPSUM(I))                                      PTH 0286
          SPTPSU(I)=REAL(TPSUM(I))                                      PTH 0287
          SPRHOP(I)=REAL(RHOPSM(I))                                     PTH 0288
   50 CONTINUE                                                          PTH 0289
      RETURN                                                            PTH 0290
      END                                                               PTH 0291
