      SUBROUTINE CRSPEC(CTHICK,CWDCOL,CIPCOL)                           CRS 0001
C                                                                       CRS 0002
C     THIS ROUTINE DEFINES CLOUD PARTICLE SPECTRAL DATA:                CRS 0003
C     EXTC(6,IWAV)   WATER DROPLET EXTINCTION COEFS [KM-1 M3/GM]        CRS 0004
C     ABSC(6,IWAV)   WATER DROPLET ABSORPTION COEFS [KM-1 M3/GM]        CRS 0005
C     ASYM(6,IWAV)   WATER DROPLET HENYEY-GREENSTEIN ASYMMETRY FACTORS  CRS 0006
C     EXTC(7,IWAV)   ICE PARTICLE EXTINCTION COEFS [KM-1 M3/GM]         CRS 0007
C     ABSC(7,IWAV)   ICE PARTICLE ABSORPTION COEFS [KM-1 M3/GM]         CRS 0008
C     ASYM(7,IWAV)   ICE PARTICLE HENYEY-GREENSTEIN ASYMMETRY FACTORS   CRS 0009
C                                                                       CRS 0010
C     LIST PARAMETERS                                                   CRS 0011
      INTEGER NMDLS,NCLDS,NCIRS                                         CRS 0012
      PARAMETER(NMDLS=10,NCLDS=5,NCIRS=2)                               CRS 0013
      INCLUDE 'PARAM.LST'                                               CRS 0014
C                                                                       CRS 0015
C     LIST ARGUMENTS                                                    CRS 0016
C     CTHICK   INPUT CLOUD VERTICAL THICKNESS [KM]                      CRS 0017
C     CWDCOL   INPUT WATER DROPLET VERTICAL COLUMN DENSITY [KM GM/M3]   CRS 0018
C     CIPCOL   INPUT ICE PARTICLE VERTICAL COLUMN DENSITY [KM GM/M3]    CRS 0019
      REAL CTHICK,CWDCOL,CIPCOL                                         CRS 0020
C                                                                       CRS 0021
C     LIST COMMONS                                                      CRS 0022
      INTEGER KPOINT                                                    CRS 0023
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     CRS 0024
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   CRS 0025
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   CRS 0026
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     CRS 0027
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               CRS 0028
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           CRS 0029
      REAL VX0,DUMDAT,CLDSPC                                            CRS 0030
      COMMON/EXTD/VX0(NWAVLN),DUMDAT(NWAVLN,66),CLDSPC(NWAVLN,3,NCLDS)  CRS 0031
      REAL CIRSPC                                                       CRS 0032
      COMMON/CIRR/CIRSPC(NWAVLN,3,NCIRS)                                CRS 0033
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA                       CRS 0034
      REAL VIS,WSS,WHH,RAINRT                                           CRS 0035
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 CRS 0036
     1  VIS,WSS,WHH,RAINRT                                              CRS 0037
      INTEGER NCRALT,NCRSPC                                             CRS 0038
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    CRS 0039
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      CRS 0040
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       CRS 0041
C                                                                       CRS 0042
C     LIST LOCAL VARIABLES AND ARRAYS                                   CRS 0043
      INTEGER ICLD0,ICIR0,IWAV,IWAVM1                                   CRS 0044
      REAL REFWAV,REFDEP,FACTOR,EXTWD,EXTIP,XNRMWD,XNRMIP,RATDEP        CRS 0045
C                                                                       CRS 0046
C     LIST DATA                                                         CRS 0047
      INTEGER MDLCLD(NMDLS),MDLCIR(NMDLS)                               CRS 0048
      REAL CLDX55(NCLDS),CIRX55(NCIRS)                                  CRS 0049
C                                                                       CRS 0050
C     MDLCLD    MAPS CLOUD/RAIN MODEL INDICES INTO CLOUD MODEL INDICES  CRS 0051
C     MDLCIR    MAPS CLOUD/RAIN MODEL INDICES INTO CIRRUS MODEL INDICES CRS 0052
C     CLDX55    0.55 MICRON CLOUD MODEL EXTINCTION COEFS [KM-1 M3/GM]   CRS 0053
      DATA MDLCLD/1,2,3,4,5,3,5,5,1,1/,MDLCIR/NMDLS*1/                  CRS 0054
      DATA CLDX55/130.16,222.78,189.68,239.41,133.01/                   CRS 0055
      DATA CIRX55/290.19, 17.21/                                        CRS 0056
C                                                                       CRS 0057
C     RETURN IF CLOUD/RAIN MODEL 1 THROUGH NMDLS WAS NOT SELECTED       CRS 0058
      IF(ICLD.LT.1 .OR. ICLD.GT.NMDLS)RETURN                            CRS 0059
C                                                                       CRS 0060
C     CHECK INPUT WATER DROPLET AND ICE PARTICLE COLUMN DENSITIES       CRS 0061
      IF(CWDCOL.LE.0.)CWDCOL=0.                                         CRS 0062
      IF(CIPCOL.LE.0.)CIPCOL=0.                                         CRS 0063
      IF(CWDCOL.EQ.0. .AND. CIPCOL.EQ.0.)THEN                           CRS 0064
C                                                                       CRS 0065
C         FATAL ERROR:  NO CLOUD VERTICAL COLUMN DENSITY                CRS 0066
          WRITE(IPR,'(/2A)')' FATAL ERROR:  CLOUD WATER DROPLET AND',   CRS 0067
     1      ' ICE PARTICLE VERTICAL COLUMN DENSITIES ARE BOTH ZERO.'    CRS 0068
          STOP 'FATAL ERROR:  NO CLOUD VERTICAL COLUMN DENSITY'         CRS 0069
      ENDIF                                                             CRS 0070
C                                                                       CRS 0071
C     ICLD0    INDEX OF CLOUD WATER DROPLET MODEL                       CRS 0072
C     ICIR0    INDEX OF CLOUD ICE PARTICLE MODEL                        CRS 0073
      ICLD0=MDLCLD(ICLD)                                                CRS 0074
      ICIR0=MDLCIR(ICLD)                                                CRS 0075
C                                                                       CRS 0076
C     FOR THE CLOUD MODELS, THE EXTINCTION COEFFICIENT NORMALIZATION    CRS 0077
C     IS INCORPORATED INTO THE SPECTRAL DATA INSTEAD OF THE COLUMN      CRS 0078
C     DENSITY.  THEREFORE, THE AWCCON INPUTS TO GAMFOG ARE ONE.         CRS 0079
      AWCCON(6)=1.                                                      CRS 0080
      AWCCON(7)=1.                                                      CRS 0081
C                                                                       CRS 0082
C     XNRMWD   WATER DROPLET EXTINCTION COEFFICIENT NORMALIZATION       CRS 0083
C     XNRMIP   ICE PARTICLE EXTINCTION COEFFICIENT NORMALIZATION        CRS 0084
      XNRMWD=CLDX55(ICLD0)                                              CRS 0085
      XNRMIP=CIRX55(ICIR0)                                              CRS 0086
C                                                                       CRS 0087
C     REFWAV    REFERENCE WAVELENGTH (MICRONS)                          CRS 0088
      REFWAV=CWAVLN                                                     CRS 0089
      IF(CWAVLN.LT.VX0(1) .OR. CWAVLN.GT.VX0(NWAVLN-1))REFWAV=.55       CRS 0090
C                                                                       CRS 0091
C     THE REFERENCE VERTICAL CLOUD DEPTH (REFDEP) IS THE                CRS 0092
C     PRODUCT OF THE VERTICAL EXTINCTION (CEXT) [KM-1]                  CRS 0093
C     AND THE CLOUD THICKNESS (CTHICK) [KM].                            CRS 0094
C                                                                       CRS 0095
C     REFDEP    VERTICAL CLOUD DEPTH AT REFWAV                          CRS 0096
      REFDEP=0.                                                         CRS 0097
      IF(CEXT.GT.0.)THEN                                                CRS 0098
          REFDEP=CEXT*CTHICK                                            CRS 0099
C                                                                       CRS 0100
C         ANNOUNCE CLOUD SPECTRAL DATA SCALING                          CRS 0101
          WRITE(IPR,'(/2A,/A,3(F9.4,A))')' CLOUD SPECTRAL DATA IS',     CRS 0102
     1      ' BEING SCALED TO MAKE THE CLOUD VERTICAL EXTINCTION',      CRS 0103
     2      ' AT',REFWAV,' MICRONS EQUAL TO',CEXT,' KM-1 FOR THE',      CRS 0104
     3      CTHICK,' KM THICK CLOUD.'                                   CRS 0105
      ENDIF                                                             CRS 0106
C                                                                       CRS 0107
C     KEY ON NCRSPC TO DETERMINE IF USER-DEFINED                        CRS 0108
C     SPECTRAL DATA IS TO BE USED.                                      CRS 0109
      IF(NCRSPC.GE.2)THEN                                               CRS 0110
C                                                                       CRS 0111
C         USER-DEFINED SPECTRAL DATA                                    CRS 0112
          CALL CRUSPC(REFWAV,REFDEP,                                    CRS 0113
     1      ICLD0,ICIR0,XNRMWD,XNRMIP,CWDCOL,CIPCOL,CTHICK)             CRS 0114
          RETURN                                                        CRS 0115
      ENDIF                                                             CRS 0116
      IF(REFDEP.GT.0.)THEN                                              CRS 0117
C                                                                       CRS 0118
C         DETERMINE BRACKETING WAVELENGTHS FOR MODEL DATA               CRS 0119
          IWAVM1=1                                                      CRS 0120
          DO 10 IWAV=2,NWAVLN-2                                         CRS 0121
              IF(REFWAV.LE.VX0(IWAV))GOTO20                             CRS 0122
   10     IWAVM1=IWAV                                                   CRS 0123
   20     CONTINUE                                                      CRS 0124
C                                                                       CRS 0125
C         DETERMINE DEFAULT EXTINCTION COEFFICIENTS AT REFWAV           CRS 0126
          FACTOR=(REFWAV-VX0(IWAVM1))/(VX0(IWAV)-VX0(IWAVM1))           CRS 0127
          EXTWD=CLDSPC(IWAVM1,1,ICLD0)                                  CRS 0128
          EXTWD=XNRMWD*(EXTWD+FACTOR*(CLDSPC(IWAV,1,ICLD0)-EXTWD))      CRS 0129
          EXTIP=CIRSPC(IWAVM1,1,ICIR0)                                  CRS 0130
          EXTIP=XNRMIP*(EXTIP+FACTOR*(CIRSPC(IWAV,1,ICIR0)-EXTIP))      CRS 0131
C                                                                       CRS 0132
C         DETERMINE RATIO OF INPUT TO CURRENT CLOUD DEPTH               CRS 0133
          RATDEP=REFDEP/(EXTWD*CWDCOL+EXTIP*CIPCOL)                     CRS 0134
C                                                                       CRS 0135
C         SCALE EXTINCTION COEFFICIENT NORMALIZATION CONSTANTS          CRS 0136
C         SO THAT INPUT CLOUD EXTINCTION RESULTS                        CRS 0137
          XNRMWD=RATDEP*XNRMWD                                          CRS 0138
          XNRMIP=RATDEP*XNRMIP                                          CRS 0139
      ENDIF                                                             CRS 0140
C                                                                       CRS 0141
C     STORE THE CLOUD PARTICLE SPECTRAL DATA                            CRS 0142
      DO 30 IWAV=1,NWAVLN-1                                             CRS 0143
          EXTC(6,IWAV)=XNRMWD*CLDSPC(IWAV,1,ICLD0)                      CRS 0144
          ABSC(6,IWAV)=XNRMWD*CLDSPC(IWAV,2,ICLD0)                      CRS 0145
          IF(ABS(ASYMWD).LT.1.)THEN                                     CRS 0146
              ASYM(6,IWAV)=ASYMWD                                       CRS 0147
          ELSE                                                          CRS 0148
              ASYM(6,IWAV)=CLDSPC(IWAV,3,ICLD0)                         CRS 0149
          ENDIF                                                         CRS 0150
          EXTC(7,IWAV)=XNRMIP*CIRSPC(IWAV,1,ICIR0)                      CRS 0151
          ABSC(7,IWAV)=XNRMIP*CIRSPC(IWAV,2,ICIR0)                      CRS 0152
          IF(ABS(ASYMIP).LT.1.)THEN                                     CRS 0153
              ASYM(7,IWAV)=ASYMIP                                       CRS 0154
          ELSE                                                          CRS 0155
              ASYM(7,IWAV)=CIRSPC(IWAV,3,ICIR0)                         CRS 0156
          ENDIF                                                         CRS 0157
   30 CONTINUE                                                          CRS 0158
C                                                                       CRS 0159
C     RETURN IF SPECTRAL DATA IS NOT TO BE OUTPUT.                      CRS 0160
      IF(NPR.GE.0)RETURN                                                CRS 0161
C                                                                       CRS 0162
C     WRITE SPECTRAL DATA HEADER                                        CRS 0163
      WRITE(IPR,'(A,/A,//54X,A,34X,A,/38X,A,2X,A,/(3A))')'1',           CRS 0164
     1  ' CLOUD SPECTRAL DATA','WATER DROPLETS','ICE PARTICLES',        CRS 0165
     2  '----------------------------------------------',               CRS 0166
     3  '----------------------------------------------',               CRS 0167
     4  ' IWAV   WAVLEN       FREQ   VERT EXT',                         CRS 0168
     5  '  EXT COEF  ABS COEF  SCT COEF     ASYM  SCT ALB',             CRS 0169
     6  '  EXT COEF  ABS COEF  SCT COEF     ASYM  SCT ALB',             CRS 0170
     7  '      (MICRON)     (CM-1)     (KM-1)',                         CRS 0171
     8  '  (        KM-1 M3/GM        )                  ',             CRS 0172
     9  '  (        KM-1 M3/GM        )'                                CRS 0173
C                                                                       CRS 0174
C     WRITE SPECTRAL DATA                                               CRS 0175
      WRITE(IPR,'((I4,F10.4,F11.3,F11.5,2(3F10.5,2F9.5)))')             CRS 0176
     1  (IWAV,VX0(IWAV),10000./VX0(IWAV),                               CRS 0177
     2  (EXTC(6,IWAV)*CWDCOL+EXTC(7,IWAV)*CIPCOL)/CTHICK,               CRS 0178
     3  EXTC(6,IWAV),ABSC(6,IWAV),EXTC(6,IWAV)-ABSC(6,IWAV),            CRS 0179
     4  ASYM(6,IWAV),1.-ABSC(6,IWAV)/EXTC(6,IWAV),                      CRS 0180
     5  EXTC(7,IWAV),ABSC(7,IWAV),EXTC(7,IWAV)-ABSC(7,IWAV),            CRS 0181
     6  ASYM(7,IWAV),1.-ABSC(7,IWAV)/EXTC(7,IWAV),IWAV=1,NWAVLN-1)      CRS 0182
      WRITE(IPR,'(/A,//)')' END OF CLOUD PARTICLE SPECTRAL DATA'        CRS 0183
      RETURN                                                            CRS 0184
      END                                                               CRS 0185
