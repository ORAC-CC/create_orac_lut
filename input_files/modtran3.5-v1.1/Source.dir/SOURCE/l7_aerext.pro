;===============================================================================
;+
;      SUBROUTINE AEREXT(V,0)                                        AEX 0001
;                                                                       AEX 0002
;     INTERPOLATES AEROSOL EXTINCTION, ABSORPTION, AND ASYMMETRY        AEX 0003
;     COEFFICIENTS FOR THE WAVENUMBER, V.                               AEX 0004
;
; R.S. 13/11/97 Extracted from MODTRAN aerext.f
;-
;===============================================================================
pro l7_AEREXT,V
;C                                                                       ;AEX 0002
;C     INTERPOLATES AEROSOL EXTINCTION, ABSORPTION, AND ASYMMETRY        ;AEX 0003
;C     COEFFICIENTS FOR THE WAVENUMBER, V.                               ;AEX 0004
;C     0   MINIMUM AEROSOL INDEX TO BE USED                          ;AEX 0005
;C             (=1 FOR ALL AEROSOLS, =6 FOR CLOUDS ONLY).                ;AEX 0006
;      REAL V                                                            ;AEX 0007
;      INTEGER 0                                                     ;AEX 0008
;C                                                                       ;AEX 0009
;C     MODIFIED FOR ASYMMETRY  - JAN 1986 (A.E.R. INC.)                  ;AEX 0010
;      INCLUDE 'PARAM.LST'                                               ;AEX 0011
;      INTEGER KPOINT                                                    ;AEX 0012
;      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     ;AEX 0013
;      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   ;AEX 0014
;     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   ;AEX 0015
;     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     ;AEX 0016
;      COMMON/EXTD/VX0(NWAVLN),DUMDAT(NWAVLN,81)                         ;AEX 0017
;      REAL EXTV,ABSV,ASYV                                               ;AEX 0018
;      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       ;AEX 0019
;      REAL TAER                                                         ;AEX 0020
;      COMMON/AERTM/TAER(NAER)                                           ;AEX 0021
;      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           ;AEX 0022
;      INTEGER NCRALT,NCRSPC                                             ;AEX 0023
;      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    ;AEX 0024
;      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      ;AEX 0025
;     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       ;AEX 0026
;      REAL G                                                            ;AEX 0028
;C                                                                       ;AEX 0030
;C     LIST DATA                                                         ;AEX 0031
;      LOGICAL LWARN                                                     ;AEX 0032
;      DATA LWARN/.TRUE./                                                ;AEX 0033
;C                                                                       ;AEX 0034
;C     NO EXTINCTION ASSUMED FOR V=0.                                    ;AEX 0035
      MXAER=7                                                           ;AEX 0036
      x=l7_extdta()
      l7_exabin,x,extc,absx,asym
stop
      vx0=x.vx0
      ALAM=10000./V
      nwavln=n_elements(vx0)
      extv=fltarr(mxaer)
      absv=fltarr(mxaer)
      asyv=fltarr(mxaer)
      if V le 1.E-5 THEN BEGIN                                          ;AEX 0037
          for IAER=0,MXAER-1 do begin
              EXTV(IAER)=0.                                             ;AEX 0039
              ABSV(IAER)=0.                                             ;AEX 0040
              ASYV(IAER)=0.                                             ;AEX 0041
	  endfor
      endif                                                             ;AEX 0045
;C                                                                       ;AEX 0092
;C     COMPUTE MICROWAVE ATTENUATION COEFFICIENTS (THE AWCCON FACTOR IS  ;AEX 0093
;C     INCLUDED IN CALLS TO GAMFOG BECAUSE EXTV AND ABSV ARE MULTIPLIED  ;AEX 0094
;C     BY COLUMN DENSITIES WHICH HAVE BEEN DIVIDED BY THIS FACTOR).      ;AEX 0095
      VMN=10000./VX0(NWAVLN-1)                                            ;AEX 0096
      if V le VMN then begin                                                  ;AEX 0097
          for IAER=0,MXAER-1 do begin                                        ;AEX 0098
              EXTV(IAER)=GAMFOG(IAER,V,TAER(IAER),AWCCON(IAER))         ;AEX 0099
              ABSV(IAER)=EXTV(IAER)                                     ;AEX 0100
              ASYV(IAER)=0.                                             ;AEX 0101
          endfor
          return                                                        ;AEX 0104
      endif                                                             ;AEX 0105
;C                                                                       ;AEX 0106
;C     DETERMINE BRACKETING WAVELENGTHS                                  ;AEX 0107
      NWAVM1=NWAVLN-1                                                   ;AEX 0108
      IWAVM1=0                                                          ;AEX 0109
      for IWAV=1,NWAVM1-1 do begin                                               ;AEX 0110
          if ALAM le VX0(IWAV) then GOTO,aext70                                   ;AEX 0111
      endfor
      IWAVM1=IWAV                                                       ;AEX 0112
;C                                                                       ;AEX 0113
;C     ALAM IS BETWEEN 200 AND 300 MICRONS.  DEFINE THE                  ;AEX 0114
;C     300 MICRON SPECTRAL DATA BEFORE INTERPOLATION                     ;AEX 0115
      FACTOR=(VX0(NWAVLN)-ALAM)/(VX0(NWAVLN)-VX0(NWAVM1))               ;AEX 0116
      for IAER=0,MXAER-1 do begin                                            ;AEX 0117
          GMFOG=GAMFOG(IAER,VMN,TAER(IAER),AWCCON(IAER))                ;AEX 0118
          EXTV(IAER)=GMFOG+FACTOR*(EXTC(IAER,NWAVM1)-GMFOG)             ;AEX 0119
          ABSV(IAER)=GMFOG+FACTOR*(ABSC(IAER,NWAVM1)-GMFOG)             ;AEX 0120
          ASYV(IAER)=FACTOR*ASYM(IAER,NWAVM1)                           ;AEX 0121
      endfor
      return                                                            ;AEX 0124
;C                                                                       ;AEX 0125
;C     COMPUTE INFRARED SPECTRAL DATA                                    ;AEX 0126
aext70:
      FACTOR=(VX0(IWAV)-ALAM)/(VX0(IWAV)-VX0(IWAVM1))                   ;AEX 0128
      for IAER=0,MXAER-1 do begin                                            ;AEX 0129
          EXTV(IAER)=EXTC(IAER,IWAV)+$                                   ;AEX 0130
            FACTOR*(EXTC(IAER,IWAVM1)-EXTC(IAER,IWAV))                  ;AEX 0131
          ABSV(IAER)=ABSC(IAER,IWAV)+$                                   ;AEX 0132
            FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))                  ;AEX 0133
          ASYV(IAER)=ASYM(IAER,IWAV)+$                                  ;AEX 0134
            FACTOR*(ASYM(IAER,IWAVM1)-ASYM(IAER,IWAV))                  ;AEX 0135
      endfor
      END                                                               ;AEX 0145
