      SUBROUTINE CRUSPC(REFWAV,REFDEP,                                  CRU 0001
     1  ICLD0,ICIR0,XNRMWD,XNRMIP,CWDCOL,CIPCOL,CTHICK)                 CRU 0002
C                                                                       CRU 0003
C     THIS ROUTINE DEFINES USER-DEFINED CLOUD SPECTRAL DATA.            CRU 0004
C     EXTC(6,IWAV)   WATER DROPLET EXTINCTION COEFS [KM-1 M3/GM]        CRU 0005
C     ABSC(6,IWAV)   WATER DROPLET ABSORPTION COEFS [KM-1 M3/GM]        CRU 0006
C     ASYM(6,IWAV)   WATER DROPLET HENYEY-GREENSTEIN ASYMMETRY FACTORS  CRU 0007
C     EXTC(7,IWAV)   ICE PARTICLE EXTINCTION COEFS [KM-1 M3/GM]         CRU 0008
C     ABSC(7,IWAV)   ICE PARTICLE ABSORPTION COEFS [KM-1 M3/GM]         CRU 0009
C     ASYM(7,IWAV)   ICE PARTICLE HENYEY-GREENSTEIN ASYMMETRY FACTORS   CRU 0010
C                                                                       CRU 0011
C     DECLARE INPUTS                                                    CRU 0012
C     REFWAV   REFERENCE WAVELENGTH [MICRONS]                           CRU 0013
C     REFDEP   REFERENCE VERTICAL OPTICAL DEPTH AT REFWAV               CRU 0014
C     ICLD0    INDEX OF CLOUD WATER DROPLET MODEL                       CRU 0015
C     ICIR0    INDEX OF CLOUD ICE PARTICLE MODEL                        CRU 0016
C     XNRMWD   DEFAULT WATER DROPLET MODEL EXTINCTION                   CRU 0017
C              COEFFICIENT AT 0.55 MICRONS [KM-1 M3/GM]                 CRU 0018
C     XNRMIP   DEFAULT ICE PARTICLE MODEL EXTINCTION                    CRU 0019
C              COEFFICIENT AT 0.55 MICRONS [KM-1 M3/GM]                 CRU 0020
C     CWDCOL   CLOUD WATER DROPLET VERTICAL COLUMN DENSITY [KM GM/M3]   CRU 0021
C     CIPCOL   CLOUD ICE PARTICLE VERTICAL COLUMN DENSITY [KM GM/M3]    CRU 0022
C     CTHICK   CLOUD VERTICAL THICKNESS [KM]                            CRU 0023
      INTEGER ICLD0,ICIR0                                               CRU 0024
      REAL REFWAV,REFDEP,XNRMWD,XNRMIP,CWDCOL,CIPCOL,CTHICK             CRU 0025
C                                                                       CRU 0026
C     LIST PARAMETERS                                                   CRU 0027
      INTEGER NCLDS,NCIRS                                               CRU 0028
      PARAMETER(NCLDS=5,NCIRS=2)                                        CRU 0029
      INCLUDE 'PARAM.LST'                                               CRU 0030
C                                                                       CRU 0031
C     LIST COMMONS                                                      CRU 0032
      INTEGER KPOINT                                                    CRU 0033
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     CRU 0034
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   CRU 0035
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   CRU 0036
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     CRU 0037
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               CRU 0038
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           CRU 0039
      REAL VX0,DUMDAT,CLDSPC                                            CRU 0040
      COMMON/EXTD/VX0(NWAVLN),DUMDAT(NWAVLN,66),CLDSPC(NWAVLN,3,NCLDS)  CRU 0041
      REAL CIRSPC                                                       CRU 0042
      COMMON/CIRR/CIRSPC(NWAVLN,3,NCIRS)                                CRU 0043
      INTEGER NCRALT,NCRSPC                                             CRU 0044
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    CRU 0045
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      CRU 0046
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       CRU 0047
C                                                                       CRU 0048
C     LIST LOCAL VARIABLES AND ARRAYS                                   CRU 0049
      INTEGER NWAVM1,I0LO,I0HI,I,N,IWAV,IWAVM1                          CRU 0050
      REAL WAVOLD,FACTOR,X,EXTWD,EXTIP,RATDEP                           CRU 0051
C                                                                       CRU 0052
C     CHECK THAT NCRSPC IS NOT TOO LARGE                                CRU 0053
      IF(NCRSPC.GT.MXWVLN)THEN                                          CRU 0054
          WRITE(IPR,'(/2A,I3,/14X,A,I3,A)')' FATAL ERROR:  THE INPUT',  CRU 0055
     1      ' NUMBER OF CLOUD SPECTRAL DATA POINTS (NCRSPC =',NCRSPC,   CRU 0056
     2      ' EXCEEDS THE MAXIMUM NUMBER (PARAMETER MXWVLN =',MXWVLN,   CRU 0057
     3      ' IN PARAMETER.LST).'                                       CRU 0058
          STOP 'INPUT NUMBER OF CLOUD SPECTRAL DATA POINTS IS TOO LARGE'CRU 0059
      ENDIF                                                             CRU 0060
C                                                                       CRU 0061
C     ANNOUNCE READING OF CLOUD SPECTRAL DATA                           CRU 0062
      WRITE(IPR,'(/A,I3,A)')' USER-DEFINED CLOUD SPECTRAL DATA AT',     CRU 0063
     1  NCRSPC,' SPECTRAL POINTS.'                                      CRU 0064
C                                                                       CRU 0065
C     LOOP OVER SPECTRAL DATA                                           CRU 0066
      NWAVM1=NWAVLN-1                                                   CRU 0067
      I0LO=1                                                            CRU 0068
      WAVOLD=-.1                                                        CRU 0069
      DO 30 I=1,NCRSPC                                                  CRU 0070
C                                                                       CRU 0071
C         READ SPECTRAL DATA                                            CRU 0072
          READ(IRD,'(7F10.5)',ERR=70)WAVLEN(I),                         CRU 0073
     1      EXTC(6,I),ABSC(6,I),ASYM(6,I),EXTC(7,I),ABSC(7,I),ASYM(7,I) CRU 0074
C                                                                       CRU 0075
C         CHECK WAVELENGTH                                              CRU 0076
          IF(WAVLEN(I).LE.WAVOLD)THEN                                   CRU 0077
              WRITE(IPR,'(/3A,/14X,2A,/(14X,I6,F10.6))')' FATAL',       CRU 0078
     1          ' ERROR:  CLOUD SPECTRAL DATA MUST BE READ IN',         CRU 0079
     2          ' INCREASING WAVELENGTH ORDER',' THE WAVELENGTHS',      CRU 0080
     3          ' ENCOUNTERED THUS FAR ARE:',(N,WAVLEN(N),N=1,I)        CRU 0081
              STOP 'NON-INCREASING WAVELENGTH IN CLOUD SPECTRAL TABLE'  CRU 0082
          ENDIF                                                         CRU 0083
          WAVOLD=WAVLEN(I)                                              CRU 0084
C                                                                       CRU 0085
C         SET UP INTERPOLATION                                          CRU 0086
          IF(WAVLEN(I).LT.VX0(1))THEN                                   CRU 0087
              I0HI=1                                                    CRU 0088
              FACTOR=0.                                                 CRU 0089
          ELSE                                                          CRU 0090
              DO 10 I0HI=I0LO+1,NWAVM1                                  CRU 0091
                  IF(WAVLEN(I).LE.VX0(I0HI))THEN                        CRU 0092
                      FACTOR=(WAVLEN(I)-VX0(I0LO))/(VX0(I0HI)-VX0(I0LO))CRU 0093
                      GOTO20                                            CRU 0094
                  ENDIF                                                 CRU 0095
   10         I0LO=I0HI                                                 CRU 0096
              I0HI=NWAVM1                                               CRU 0097
              FACTOR=0.                                                 CRU 0098
          ENDIF                                                         CRU 0099
   20     CONTINUE                                                      CRU 0100
C                                                                       CRU 0101
C         CHECK WATER DROPLET SPECTRAL DATA                             CRU 0102
C                                                                       CRU 0103
C         IF THE INPUT SPECTRAL EXTINCTION COEFFICIENT IS NEGATIVE, IT  CRU 0104
C         IS REPLACED BY THE WAVELENGTH-INTERPOLATED CLOUD MODEL VALUE. CRU 0105
C         IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS LESS THAN -1. CRU 0106
C         OR IF IT EXCEEDS THE EXTINCTION COEFFICIENT, IT IS REPLACED   CRU 0107
C         BY THE VALUE WHICH YIELDS THE CLOUD MODEL SINGLE SCATTERING   CRU 0108
C         ALBEDO.  IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS      CRU 0109
C         NEGATIVE BUT NOT LESS THAN -1., THAN THE INPUT VALUE IS TAKEN CRU 0110
C         TO BE THE SINGLE SCATTERING ALBEDO MINUS ONE, -(1-OMEGA).     CRU 0111
          IF(EXTC(6,I).LT.0.)THEN                                       CRU 0112
              EXTC(6,I)=XNRMWD*(CLDSPC(I0LO,1,ICLD0)                    CRU 0113
     1          +FACTOR*(CLDSPC(I0HI,1,ICLD0)-CLDSPC(I0LO,1,ICLD0)))    CRU 0114
              IF(ABSC(6,I).LT.-1. .OR. ABSC(6,I).GT.EXTC(6,I))THEN      CRU 0115
                  ABSC(6,I)=XNRMWD*(CLDSPC(I0LO,2,ICLD0)                CRU 0116
     1              +FACTOR*(CLDSPC(I0HI,2,ICLD0)-CLDSPC(I0LO,2,ICLD0)))CRU 0117
              ELSEIF(ABSC(6,I).LT.0.)THEN                               CRU 0118
                  ABSC(6,I)=-ABSC(6,I)*EXTC(6,I)                        CRU 0119
              ENDIF                                                     CRU 0120
          ELSEIF(ABSC(6,I).LT.-1. .OR. ABSC(6,I).GT.EXTC(6,I))THEN      CRU 0121
              X=CLDSPC(I0LO,1,ICLD0)                                    CRU 0122
     1          +FACTOR*(CLDSPC(I0HI,1,ICLD0)-CLDSPC(I0LO,1,ICLD0))     CRU 0123
              ABSC(6,I)=0.                                              CRU 0124
              IF(X.GT.0.)ABSC(6,I)=EXTC(6,I)*(CLDSPC(I0LO,2,ICLD0)      CRU 0125
     1          +FACTOR*(CLDSPC(I0HI,2,ICLD0)-CLDSPC(I0LO,2,ICLD0)))/X  CRU 0126
          ELSEIF(ABSC(6,I).LT.0.)THEN                                   CRU 0127
              ABSC(6,I)=-ABSC(6,I)*EXTC(6,I)                            CRU 0128
          ENDIF                                                         CRU 0129
C                                                                       CRU 0130
C         IF CARD2A INPUT IS AN ACCEPTABLE VALUE (-1 < ASYMWD < 1) FOR  CRU 0131
C         THE WATER DROPLET HENYEY-GREENSTEIN PHASE FUNCTION ASYMMETRY  CRU 0132
C         FACTOR, THE CONSTANT VALUE IS USED AT ALL WAVELENGTHS.        CRU 0133
C         IF THE ASYMWD IS OUT OF RANGE, THE USER-DEFINED SPECTRAL      CRU 0134
C         ASYMMETRY FACTOR IS USED UNLESS IT IS ALSO OUT OF RANGE.      CRU 0135
C         IF BOTH ASYMMETRY FACTOR INPUTS ARE OUT OF RANGE, THE         CRU 0136
C         WAVELENGTH-INTERPOLATED CLOUD MODEL VALUE IS USED.            CRU 0137
          IF(ABS(ASYMWD).LT.1.)THEN                                     CRU 0138
              ASYM(6,I)=ASYMWD                                          CRU 0139
          ELSEIF(ABS(ASYM(6,I)).GE.1.)THEN                              CRU 0140
              ASYM(6,I)=CLDSPC(I0LO,3,ICLD0)                            CRU 0141
     1          +FACTOR*(CLDSPC(I0HI,3,ICLD0)-CLDSPC(I0LO,3,ICLD0))     CRU 0142
          ENDIF                                                         CRU 0143
C                                                                       CRU 0144
C         CHECK ICE PARTICLE SPECTRAL DATA                              CRU 0145
C                                                                       CRU 0146
C         IF THE INPUT SPECTRAL EXTINCTION COEFFICIENT IS NEGATIVE, IT  CRU 0147
C         IS REPLACED BY THE WAVELENGTH-INTERPOLATED CIRRUS MODEL VALUE.CRU 0148
C         IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS NEGATIVE OR   CRU 0149
C         IF IT EXCEEDS THE EXTINCTION COEFFICIENT, IT IS REPLACED BY   CRU 0150
C         THE VALUE WHICH YIELDS THE CIRRUS MODEL SINGLE SCATTERING     CRU 0151
C         ALBEDO.  IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS      CRU 0152
C         NEGATIVE BUT NOT LESS THAN -1., THAN THE INPUT VALUE IS TAKEN CRU 0153
C         TO BE THE SINGLE SCATTERING ALBEDO MINUS ONE, -(1-OMEGA).     CRU 0154
          IF(EXTC(7,I).LT.0.)THEN                                       CRU 0155
              EXTC(7,I)=XNRMIP*(CIRSPC(I0LO,1,ICIR0)                    CRU 0156
     1          +FACTOR*(CIRSPC(I0HI,1,ICIR0)-CIRSPC(I0LO,1,ICIR0)))    CRU 0157
              IF(ABSC(7,I).LT.-1. .OR. ABSC(7,I).GT.EXTC(7,I))THEN      CRU 0158
                  ABSC(7,I)=XNRMIP*(CIRSPC(I0LO,2,ICIR0)                CRU 0159
     1              +FACTOR*(CIRSPC(I0HI,2,ICIR0)-CIRSPC(I0LO,2,ICIR0)))CRU 0160
              ELSEIF(ABSC(7,I).LT.0.)THEN                               CRU 0161
                  ABSC(7,I)=-ABSC(7,I)*EXTC(7,I)                        CRU 0162
              ENDIF                                                     CRU 0163
          ELSEIF(ABSC(7,I).LT.-1. .OR. ABSC(7,I).GT.EXTC(7,I))THEN      CRU 0164
              X=CIRSPC(I0LO,1,ICIR0)                                    CRU 0165
     1          +FACTOR*(CIRSPC(I0HI,1,ICIR0)-CIRSPC(I0LO,1,ICIR0))     CRU 0166
              ABSC(7,I)=0.                                              CRU 0167
              IF(X.GT.0.)ABSC(7,I)=EXTC(7,I)*(CIRSPC(I0LO,2,ICIR0)      CRU 0168
     1          +FACTOR*(CIRSPC(I0HI,2,ICIR0)-CIRSPC(I0LO,2,ICIR0)))/X  CRU 0169
          ELSEIF(ABSC(7,I).LT.0.)THEN                                   CRU 0170
              ABSC(7,I)=-ABSC(7,I)*EXTC(7,I)                            CRU 0171
          ENDIF                                                         CRU 0172
C                                                                       CRU 0173
C         IF CARD2A INPUT IS AN ACCEPTABLE VALUE (-1 < ASYMIP < 1) FOR  CRU 0174
C         THE ICE PARTICLE HENYEY-GREENSTEIN PHASE FUNCTION ASYMMETRY   CRU 0175
C         FACTOR, THE CONSTANT VALUE IS USED AT ALL WAVELENGTHS.        CRU 0176
C         IF THE ASYMIP IS OUT OF RANGE, THE USER-DEFINED SPECTRAL      CRU 0177
C         ASYMMETRY FACTOR IS USED UNLESS IT IS ALSO OUT OF RANGE.      CRU 0178
C         IF BOTH ASYMMETRY FACTOR INPUTS ARE OUT OF RANGE, THE         CRU 0179
C         WAVELENGTH-INTERPOLATED CIRRUS MODEL VALUE IS USED.           CRU 0180
          IF(ABS(ASYMIP).LT.1.)THEN                                     CRU 0181
              ASYM(7,I)=ASYMIP                                          CRU 0182
          ELSEIF(ABS(ASYM(7,I)).GE.1.)THEN                              CRU 0183
              ASYM(7,I)=CIRSPC(I0LO,3,ICIR0)                            CRU 0184
     1          +FACTOR*(CIRSPC(I0HI,3,ICIR0)-CIRSPC(I0LO,3,ICIR0))     CRU 0185
          ENDIF                                                         CRU 0186
   30 CONTINUE                                                          CRU 0187
C                                                                       CRU 0188
C     RENORMALIZE SPECTRAL DATA IF REFERENCE OPTICAL DEPTH WAS INPUT    CRU 0189
      IF(REFDEP.GT.0.)THEN                                              CRU 0190
C                                                                       CRU 0191
C         DETERMINE BRACKETING WAVELENGTHS FOR MODEL DATA               CRU 0192
          IF(REFWAV.LT.WAVLEN(1) .OR. REFWAV.GT.WAVLEN(NCRSPC))THEN     CRU 0193
              WRITE(IPR,'(/A,F10.6,A,/14X,2A,/14X,A,2(F10.6,A))')       CRU 0194
     1          ' FATAL ERROR:  THE WAVELENGTH (',REFWAV,               CRU 0195
     2          ' MICRONS) USED TO DEFINE THE CLOUD',                   CRU 0196
     3          ' VERTICAL OPTICAL DEPTH IS OUTSIDE THE RANGE',         CRU 0197
     4          ' OF THE USER-DEFINED',' CLOUD SPECTRAL DATA (',        CRU 0198
     5          WAVLEN(1),' TO',WAVLEN(NCRSPC),' MICRONS).'             CRU 0199
              STOP 'INPUT SPECTRAL CLOUD DEPTH OUTSIDE SPECTRAL RANGE'  CRU 0200
          ENDIF                                                         CRU 0201
          IWAVM1=1                                                      CRU 0202
          DO 40 IWAV=2,NCRSPC-1                                         CRU 0203
              IF(REFWAV.LE.WAVLEN(IWAV))GOTO50                          CRU 0204
   40     IWAVM1=IWAV                                                   CRU 0205
          IWAV=NCRSPC                                                   CRU 0206
   50     CONTINUE                                                      CRU 0207
C                                                                       CRU 0208
C         DETERMINE DEFAULT EXTINCTION COEFFICIENTS AT REFWAV           CRU 0209
          FACTOR=(REFWAV-WAVLEN(IWAVM1))/(WAVLEN(IWAV)-WAVLEN(IWAVM1))  CRU 0210
          EXTWD=EXTC(6,IWAVM1)+FACTOR*(EXTC(6,IWAV)-EXTC(6,IWAVM1))     CRU 0211
          EXTIP=EXTC(7,IWAVM1)+FACTOR*(EXTC(7,IWAV)-EXTC(7,IWAVM1))     CRU 0212
C                                                                       CRU 0213
C         DETERMINE RATIO OF INPUT TO CURRENT CLOUD DEPTH               CRU 0214
          RATDEP=REFDEP/(EXTWD*CWDCOL+EXTIP*CIPCOL)                     CRU 0215
C                                                                       CRU 0216
C         SCALE THE CLOUD PARTICLE SPECTRAL DATA                        CRU 0217
          DO 60 IWAV=1,NCRSPC                                           CRU 0218
              EXTC(6,IWAV)=RATDEP*EXTC(6,IWAV)                          CRU 0219
              ABSC(6,IWAV)=RATDEP*ABSC(6,IWAV)                          CRU 0220
              EXTC(7,IWAV)=RATDEP*EXTC(7,IWAV)                          CRU 0221
              ABSC(7,IWAV)=RATDEP*ABSC(7,IWAV)                          CRU 0222
   60     CONTINUE                                                      CRU 0223
      ENDIF                                                             CRU 0224
C                                                                       CRU 0225
C     RETURN TO CRUSPC IF SPECTRAL DATA IS NOT TO BE OUTPUT.            CRU 0226
      IF(NPR.GE.0)RETURN                                                CRU 0227
C                                                                       CRU 0228
C     WRITE SPECTRAL DATA HEADER                                        CRU 0229
      WRITE(IPR,'(A,/A,//54X,A,34X,A,/38X,A,2X,A,/(3A))')'1',           CRU 0230
     1  ' CLOUD SPECTRAL DATA','WATER DROPLETS','ICE PARTICLES',        CRU 0231
     2  '----------------------------------------------',               CRU 0232
     3  '----------------------------------------------',               CRU 0233
     4  ' IWAV   WAVLEN       FREQ   VERT EXT',                         CRU 0234
     5  '  EXT COEF  ABS COEF  SCT COEF     ASYM  SCT ALB',             CRU 0235
     6  '  EXT COEF  ABS COEF  SCT COEF     ASYM  SCT ALB',             CRU 0236
     7  '      (MICRON)     (CM-1)     (KM-1)',                         CRU 0237
     8  '  (        KM-1 M3/GM        )                  ',             CRU 0238
     9  '  (        KM-1 M3/GM        )'                                CRU 0239
C                                                                       CRU 0240
C     WRITE SPECTRAL DATA                                               CRU 0241
      WRITE(IPR,'((I4,F10.4,F11.3,F11.5,2(3F10.5,2F9.5)))')             CRU 0242
     1  (IWAV,WAVLEN(IWAV),10000./WAVLEN(IWAV),                         CRU 0243
     2  (EXTC(6,IWAV)*CWDCOL+EXTC(7,IWAV)*CIPCOL)/CTHICK,               CRU 0244
     3  EXTC(6,IWAV),ABSC(6,IWAV),EXTC(6,IWAV)-ABSC(6,IWAV),            CRU 0245
     4  ASYM(6,IWAV),1.-ABSC(6,IWAV)/EXTC(6,IWAV),                      CRU 0246
     5  EXTC(7,IWAV),ABSC(7,IWAV),EXTC(7,IWAV)-ABSC(7,IWAV),            CRU 0247
     6  ASYM(7,IWAV),1.-ABSC(7,IWAV)/EXTC(7,IWAV),IWAV=1,NCRSPC)        CRU 0248
      WRITE(IPR,'(/A,//)')' END OF CLOUD PARTICLE SPECTRAL DATA'        CRU 0249
C                                                                       CRU 0250
C     RETURN TO ROUTINE CRSPEC                                          CRU 0251
      RETURN                                                            CRU 0252
C                                                                       CRU 0253
C     FATAL ERROR READING CLOUD SPECTRAL DATA                           CRU 0254
   70 CONTINUE                                                          CRU 0255
      WRITE(IPR,'(/A,I3,A)')' FATAL ERROR:  UNABLE TO READ LINE',       CRU 0256
     1  I,' OF CLOUD SPECTRAL DATA.'                                    CRU 0257
      IF(I.GT.1)WRITE(IPR,'(/A,I3,A,/(5X,7F10.6))')                     CRU 0258
     1  ' THE FIRST',I-1,' LINES OF DATA ARE:',(WAVLEN(N),EXTC(6,N),    CRU 0259
     2  ABSC(6,N),ASYM(6,N),EXTC(7,N),ABSC(7,N),ASYM(7,N),N=1,I-1)      CRU 0260
      STOP 'PROBLEM READING CLOUD SPECTRAL DATA'                        CRU 0261
      END                                                               CRU 0262
