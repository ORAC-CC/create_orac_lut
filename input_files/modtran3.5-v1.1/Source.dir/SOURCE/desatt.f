      SUBROUTINE DESATT(WSPD,VIS)                                       DES 0001
C********************************************************************** DES 0002
C*                                                                    * DES 0003
C*    THIS SUBROUTINE CALCULATES THE ATTENUATION COEFFICIENTS AND     * DES 0004
C*    ASYMMETRY PARAMETER FOR THE DESERT AEROSOL BASED ON THE WIND    * DES 0005
C*    SPEED AND METEOROLOGICAL RANGE                                  * DES 0006
C*                                                                    * DES 0007
C*                                                                    * DES 0008
C*                                                                    * DES 0009
C*    PROGRAMMED BY:  D. R. LONGTIN         OPTIMETRICS, INC.         * DES 0010
C*                                          BURLINGTON, MASSACHUSETTS * DES 0011
C*                                          JULY 1987                 * DES 0012
C*                                                                    * DES 0013
C*                                                                    * DES 0014
C*    INPUTS:    WSPD    -  WIND SPEED (IN M/S) AT 10 M               * DES 0015
C*               VIS     -  METEOROLOGICAL RANGE (KM)                 * DES 0016
C*                                                                    * DES 0017
C*    OUTPUTS: DESEXT - EXTINCTION COEFFICIENT AT "NWAVLN" WAVELENGTHS* DES 0018
C*    *****    DESABS - ABSORPTION COEFFICIENT AT "NWAVLN" WAVELENGTHS* DES 0019
C*             DESG   - ASYMMETRY PARAMETER AT "NWAVLN" WAVELENGTHS   * DES 0020
C*                                                                    * DES 0021
C********************************************************************** DES 0022
C                                                                       DES 0023
      INCLUDE 'PARAM.LST'                                               DES 0024
      INTEGER KPOINT                                                    DES 0025
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     DES 0026
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   DES 0027
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   DES 0028
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     DES 0029
      REAL DSEXT,DSABS,DSASYM                                           DES 0030
      COMMON/DESAER/DSEXT(NWAVLN,4),DSABS(NWAVLN,4),DSASYM(NWAVLN,4)    DES 0031
      REAL DESEXT(NWAVLN),DESABS(NWAVLN),DESG(NWAVLN),WIND(4)           DES 0032
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          DES 0033
      INTEGER WAVEL                                                     DES 0034
      DATA WIND/0., 10., 20., 30./                                      DES 0035
      DATA RAYSCT / 0.01159 /                                           DES 0036
      IF(WSPD .LT. 0.) WSPD = 10.                                       DES 0037
C                                                                       DES 0038
      NWSPD = INT(WSPD/10) + 1                                          DES 0039
      IF (NWSPD.GE.4) WRITE(IPR,999)                                    DES 0040
      IF (NWSPD.GE.4) NWSPD = 3                                         DES 0041
C                                                                       DES 0042
C     INTERPOLATE THE RADIATIVE PROPERTIES AT WIND SPEED WSPD           DES 0043
C                                                                       DES 0044
      DO 100 WAVEL=1,NWAVLN                                             DES 0045
C                                                                       DES 0046
C     EXTINCTION COEFFICIENT                                            DES 0047
C                                                                       DES 0048
         SLOPE = LOG(DSEXT(WAVEL,NWSPD+1)/DSEXT(WAVEL,NWSPD))/          DES 0049
     *           (WIND(NWSPD+1)-WIND(NWSPD))                            DES 0050
         B = LOG(DSEXT(WAVEL,NWSPD+1)) - SLOPE*WIND(NWSPD+1)            DES 0051
         DESEXT(WAVEL) = EXP(SLOPE*WSPD + B)                            DES 0052
C                                                                       DES 0053
C     ABSORPTION COEFFICIENT                                            DES 0054
C                                                                       DES 0055
         SLOPE = LOG(DSABS(WAVEL,NWSPD+1)/DSABS(WAVEL,NWSPD))/          DES 0056
     *           (WIND(NWSPD+1)-WIND(NWSPD))                            DES 0057
         B = LOG(DSABS(WAVEL,NWSPD+1)) - SLOPE*WIND(NWSPD+1)            DES 0058
         DESABS(WAVEL) = EXP(SLOPE*WSPD + B)                            DES 0059
C                                                                       DES 0060
C     ASYMMETRY PARAMETER                                               DES 0061
C                                                                       DES 0062
         SLOPE = (DSASYM(WAVEL,NWSPD+1)-DSASYM(WAVEL,NWSPD))            DES 0063
     1           /(WIND(NWSPD+1)-WIND(NWSPD))                           DES 0064
         B = DSASYM(WAVEL,NWSPD+1) - SLOPE*WIND(NWSPD+1)                DES 0065
         DESG(WAVEL) = SLOPE*WSPD + B                                   DES 0066
100   CONTINUE                                                          DES 0067
C                                                                       DES 0068
          EXT55 = DESEXT(4)                                             DES 0069
C                                                                       DES 0070
C         DETERMINE METEROLOGICAL RANGE FROM 0.55 EXTINCTION            DES 0071
C          AND KOSCHMIEDER FORMULA                                      DES 0072
C                                                                       DES 0073
      IF (VIS. LE .0.) THEN                                             DES 0074
               VIS = 3.912/(DESEXT(4) + RAYSCT )                        DES 0075
      ENDIF                                                             DES 0076
C                                                                       DES 0077
C        RENORMALIZE ATTENUATION COEFFICIENTS TO 1.0 KM-1 AT            DES 0078
C        0.55 MICRONS FOR CAPABILTY WITH LOWTRAN                        DES 0079
C                                                                       DES 0080
          DO 200 WAVEL=1,NWAVLN                                         DES 0081
            EXTC(1,WAVEL) = DESEXT(WAVEL)       /EXT55                  DES 0082
            ABSC(1,WAVEL) = DESABS(WAVEL)       /EXT55                  DES 0083
            ASYM(1,WAVEL) = DESG(WAVEL)                                 DES 0084
200      CONTINUE                                                       DES 0085
       WRITE(IPR,900) VIS,WSPD                                          DES 0086
900    FORMAT(//,'  VIS = ',F10.3,' WIND = ',F10.3)                     DES 0087
       RETURN                                                           DES 0088
C                                                                       DES 0089
999    FORMAT(' WARNING: WIND SPEED IS BEYOND 30 M/S; RADIATIVE',       DES 0090
     *'PROPERTIES',/,'OF THE DESERT AEROSOL HAVE BEEN EXTRAPOLATED')    DES 0091
      END                                                               DES 0092
