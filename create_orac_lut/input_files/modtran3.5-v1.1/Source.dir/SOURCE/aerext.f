      SUBROUTINE AEREXT(V,MNAER)                                        AEX 0001
C                                                                       AEX 0002
C     INTERPOLATES AEROSOL EXTINCTION, ABSORPTION, AND ASYMMETRY        AEX 0003
C     COEFFICIENTS FOR THE WAVENUMBER, V.                               AEX 0004
C     MNAER   MINIMUM AEROSOL INDEX TO BE USED                          AEX 0005
C             (=1 FOR ALL AEROSOLS, =6 FOR CLOUDS ONLY).                AEX 0006
      REAL V                                                            AEX 0007
      INTEGER MNAER                                                     AEX 0008
C                                                                       AEX 0009
C     MODIFIED FOR ASYMMETRY  - JAN 1986 (A.E.R. INC.)                  AEX 0010
      INCLUDE 'PARAM.LST'                                               AEX 0011
      INTEGER KPOINT                                                    AEX 0012
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     AEX 0013
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   AEX 0014
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   AEX 0015
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     AEX 0016
      COMMON/EXTD/VX0(NWAVLN),DUMDAT(NWAVLN,81)                         AEX 0017
      REAL EXTV,ABSV,ASYV                                               AEX 0018
      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       AEX 0019
      REAL TAER                                                         AEX 0020
      COMMON/AERTM/TAER(NAER)                                           AEX 0021
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           AEX 0022
      INTEGER NCRALT,NCRSPC                                             AEX 0023
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    AEX 0024
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      AEX 0025
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       AEX 0026
      INTEGER IPH                                                       AEX 0027
      REAL G                                                            AEX 0028
      COMMON/CARD3A/IPH,G                                               AEX 0029
C                                                                       AEX 0030
C     LIST DATA                                                         AEX 0031
      LOGICAL LWARN                                                     AEX 0032
      DATA LWARN/.TRUE./                                                AEX 0033
C                                                                       AEX 0034
C     NO EXTINCTION ASSUMED FOR V=0.                                    AEX 0035
      MXAER=7                                                           AEX 0036
      IF(V.LE.1.E-5)THEN                                                AEX 0037
          DO 10 IAER=MNAER,MXAER                                        AEX 0038
              EXTV(IAER)=0.                                             AEX 0039
              ABSV(IAER)=0.                                             AEX 0040
              ASYV(IAER)=0.                                             AEX 0041
   10     CONTINUE                                                      AEX 0042
          IF(IPH.EQ.0)GOTO90                                            AEX 0043
          RETURN                                                        AEX 0044
      ENDIF                                                             AEX 0045
C                                                                       AEX 0046
C     CHECK FOR USER-DEFINED CLOUD SPECTRAL DATA                        AEX 0047
      ALAM=10000./V                                                     AEX 0048
      IF(NCRSPC.GE.2)THEN                                               AEX 0049
C                                                                       AEX 0050
C         USER-DEFINED CLOUD SPECTRAL DATA (IAER=6 AND IAER=7)          AEX 0051
C         IS BASED ON "WAVLEN" WAVELENGTHS, NOT "VX0" WAVELENGTHS.      AEX 0052
          MXAER=5                                                       AEX 0053
          IF(ALAM.LT.WAVLEN(1) .OR. ALAM.GT.WAVLEN(NCRSPC))THEN         AEX 0054
              IF(LWARN)THEN                                             AEX 0055
                  WRITE(IPR,'(2A,F8.4,A,/10X,2(A,F8.4),A,/10X,A,F8.4,   AEX 0056
     1              2A,/10X,2A,/)')' WARNING:  CLOUD SPECTRAL DATA IS', AEX 0057
     2              ' REQUIRED AT',ALAM,' MICRONS, BUT DATA WAS ONLY',  AEX 0058
     3              ' INPUT FROM',WAVLEN(1),' TO',WAVLEN(NCRSPC),       AEX 0059
     4              ' MICRONS.  END POINT DATA WILL BE USED',           AEX 0060
     5              ' FOR',ALAM,' MICRONS.  THIS WARNING IS',           AEX 0061
     6              ' NOT REPEATED IF ADDITIONAL DATA',                 AEX 0062
     7              ' IS REQUIRED THAT FALLS OUTSIDE',                  AEX 0063
     8              ' THE INPUT SPECTRAL RANGE.'                        AEX 0064
                  LWARN=.FALSE.                                         AEX 0065
              ENDIF                                                     AEX 0066
              IWAV=1                                                    AEX 0067
              IF(ALAM.GT.WAVLEN(NCRSPC))IWAV=NCRSPC                     AEX 0068
              EXTV(6)=EXTC(6,IWAV)                                      AEX 0069
              ABSV(6)=ABSC(6,IWAV)                                      AEX 0070
              ASYV(6)=ASYM(6,IWAV)                                      AEX 0071
              EXTV(7)=EXTC(7,IWAV)                                      AEX 0072
              ABSV(7)=ABSC(7,IWAV)                                      AEX 0073
              ASYV(7)=ASYM(7,IWAV)                                      AEX 0074
          ELSE                                                          AEX 0075
              IWAVM1=1                                                  AEX 0076
              DO 20 IWAV=2,NCRSPC-1                                     AEX 0077
                  IF(ALAM.LE.WAVLEN(IWAV))GOTO30                        AEX 0078
   20         IWAVM1=IWAV                                               AEX 0079
              IWAV=NCRSPC                                               AEX 0080
   30         CONTINUE                                                  AEX 0081
              FACTOR=(WAVLEN(IWAV)-ALAM)/(WAVLEN(IWAV)-WAVLEN(IWAVM1))  AEX 0082
              EXTV(6)=EXTC(6,IWAV)+FACTOR*(EXTC(6,IWAVM1)-EXTC(6,IWAV)) AEX 0083
              ABSV(6)=ABSC(6,IWAV)+FACTOR*(ABSC(6,IWAVM1)-ABSC(6,IWAV)) AEX 0084
              ASYV(6)=ASYM(6,IWAV)+FACTOR*(ASYM(6,IWAVM1)-ASYM(6,IWAV)) AEX 0085
              EXTV(7)=EXTC(7,IWAV)+FACTOR*(EXTC(7,IWAVM1)-EXTC(7,IWAV)) AEX 0086
              ABSV(7)=ABSC(7,IWAV)+FACTOR*(ABSC(7,IWAVM1)-ABSC(7,IWAV)) AEX 0087
              ASYV(7)=ASYM(7,IWAV)+FACTOR*(ASYM(7,IWAVM1)-ASYM(7,IWAV)) AEX 0088
          ENDIF                                                         AEX 0089
          IF(MNAER.GE.6)RETURN                                          AEX 0090
      ENDIF                                                             AEX 0091
C                                                                       AEX 0092
C     COMPUTE MICROWAVE ATTENUATION COEFFICIENTS (THE AWCCON FACTOR IS  AEX 0093
C     INCLUDED IN CALLS TO GAMFOG BECAUSE EXTV AND ABSV ARE MULTIPLIED  AEX 0094
C     BY COLUMN DENSITIES WHICH HAVE BEEN DIVIDED BY THIS FACTOR).      AEX 0095
      VMN=10000./VX0(NWAVLN)                                            AEX 0096
      IF(V.LE.VMN)THEN                                                  AEX 0097
          DO 40 IAER=MNAER,MXAER                                        AEX 0098
              EXTV(IAER)=GAMFOG(IAER,V,TAER(IAER),AWCCON(IAER))         AEX 0099
              ABSV(IAER)=EXTV(IAER)                                     AEX 0100
              ASYV(IAER)=0.                                             AEX 0101
   40     CONTINUE                                                      AEX 0102
          IF(IPH.EQ.0)GOTO90                                            AEX 0103
          RETURN                                                        AEX 0104
      ENDIF                                                             AEX 0105
C                                                                       AEX 0106
C     DETERMINE BRACKETING WAVELENGTHS                                  AEX 0107
      NWAVM1=NWAVLN-1                                                   AEX 0108
      IWAVM1=1                                                          AEX 0109
      DO 50 IWAV=2,NWAVM1                                               AEX 0110
          IF(ALAM.LE.VX0(IWAV))GOTO70                                   AEX 0111
   50 IWAVM1=IWAV                                                       AEX 0112
C                                                                       AEX 0113
C     ALAM IS BETWEEN 200 AND 300 MICRONS.  DEFINE THE                  AEX 0114
C     300 MICRON SPECTRAL DATA BEFORE INTERPOLATION                     AEX 0115
      FACTOR=(VX0(NWAVLN)-ALAM)/(VX0(NWAVLN)-VX0(NWAVM1))               AEX 0116
      DO 60 IAER=MNAER,MXAER                                            AEX 0117
          GMFOG=GAMFOG(IAER,VMN,TAER(IAER),AWCCON(IAER))                AEX 0118
          EXTV(IAER)=GMFOG+FACTOR*(EXTC(IAER,NWAVM1)-GMFOG)             AEX 0119
          ABSV(IAER)=GMFOG+FACTOR*(ABSC(IAER,NWAVM1)-GMFOG)             AEX 0120
          ASYV(IAER)=FACTOR*ASYM(IAER,NWAVM1)                           AEX 0121
   60 CONTINUE                                                          AEX 0122
      IF(IPH.EQ.0)GOTO90                                                AEX 0123
      RETURN                                                            AEX 0124
C                                                                       AEX 0125
C     COMPUTE INFRARED SPECTRAL DATA                                    AEX 0126
   70 CONTINUE                                                          AEX 0127
      FACTOR=(VX0(IWAV)-ALAM)/(VX0(IWAV)-VX0(IWAVM1))                   AEX 0128
      DO 80 IAER=MNAER,MXAER                                            AEX 0129
          EXTV(IAER)=EXTC(IAER,IWAV)+                                   AEX 0130
     1      FACTOR*(EXTC(IAER,IWAVM1)-EXTC(IAER,IWAV))                  AEX 0131
          ABSV(IAER)=ABSC(IAER,IWAV)+                                   AEX 0132
     1      FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))                  AEX 0133
          ASYV(IAER)=ASYM(IAER,IWAV)+                                   AEX 0134
     1      FACTOR*(ASYM(IAER,IWAVM1)-ASYM(IAER,IWAV))                  AEX 0135
   80 CONTINUE                                                          AEX 0136
      IF(IPH.NE.0)RETURN                                                AEX 0137
C                                                                       AEX 0138
C     OVERWRITE AEROSOL ASYMMETRY FACTORS WITH USER-DEFINED VALUE       AEX 0139
   90 CONTINUE                                                          AEX 0140
      DO 100 IAER=MNAER,5                                               AEX 0141
          ASYV(IAER)=G                                                  AEX 0142
  100 CONTINUE                                                          AEX 0143
      RETURN                                                            AEX 0144
      END                                                               AEX 0145
