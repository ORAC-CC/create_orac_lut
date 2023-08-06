      SUBROUTINE XINTRP(LAYX,IXMOLS)                                    XIN 0001
C                                                                       XIN 0002
C     THIS SUBROUTINE EXPONENTIALLY INTERPOLATES THE PROFILE            XIN 0003
C     DENX ON THE ALTITUDE GRID ZX INTO DENM ON THE GRID ZM.            XIN 0004
C                                                                       XIN 0005
C     DECLARE ARGUMENTS.                                                XIN 0006
      INTEGER LAYX,IXMOLS                                               XIN 0007
C                                                                       XIN 0008
C     DECLARE PARAMETERS.                                               XIN 0009
      INCLUDE 'PARAM.LST'                                               XIN 0010
C                                                                       XIN 0011
C     LIST COMMONS.                                                     XIN 0012
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               XIN 0013
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           XIN 0014
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     XIN 0015
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                XIN 0016
      REAL PZERO,TZERO,AVOGAD,ALOSMT,GASCON,PLANK,BOLTZ,                XIN 0017
     1  CLIGHT,ADCON,ALZERO,AVMWT,AIRMWT,AMWT                           XIN 0018
      COMMON/CONSTN/PZERO,TZERO,AVOGAD,ALOSMT,GASCON,PLANK,BOLTZ,       XIN 0019
     1  CLIGHT,ADCON,ALZERO,AVMWT,AIRMWT,AMWT(35)                       XIN 0020
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2                             XIN 0021
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         XIN 0022
     1  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)    XIN 0023
      REAL ZM,PM,TM,RFNDX,DENSTY                                        XIN 0024
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    XIN 0025
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               XIN 0026
      REAL DENM,DENP                                                    XIN 0027
      COMMON/DEAMT/DENM(35,LAYTWO),DENP(35,LAYTHR+1)                    XIN 0028
      REAL ZX,DTMP,DENX                                                 XIN 0029
      COMMON/XPDEM/ZX(LAYDIM),DTMP(35),DENX(35,LAYDIM)                  XIN 0030
C                                                                       XIN 0031
C     DECLARE FUNCTIONS.                                                XIN 0032
      REAL EXPINT                                                       XIN 0033
C                                                                       XIN 0034
C     DECLARE LOCAL VARIABLES.                                          XIN 0035
      INTEGER LX,L,K                                                    XIN 0036
      REAL A,RHOAIR                                                     XIN 0037
C                                                                       XIN 0038
C     BEGIN CALCULATION.                                                XIN 0039
      LX=2                                                              XIN 0040
      DO 30 L=1,ML                                                      XIN 0041
C                                                                       XIN 0042
C         FIND THE SMALLEST ZX GE ZM(L)                                 XIN 0043
   10     CONTINUE                                                      XIN 0044
          IF(ZM(L).LE.ZX(LX) .OR. LX.EQ.LAYX)THEN                       XIN 0045
              A=(ZM(L)-ZX(LX-1))/(ZX(LX)-ZX(LX-1))                      XIN 0046
              IF(A.LT.0.)WRITE(IPR,'(/2A)')' XINTRP warning:  ',        XIN 0047
     1          'Extrapolating MODEL atmosphere below minimum altitude.'XIN 0048
              IF(A.GT.1.)WRITE(IPR,'(/2A)')' XINTRP warning:  ',        XIN 0049
     1          'Extrapolating MODEL atmosphere above maximum altitude.'XIN 0050
C                                                                       XIN 0051
C             CALCULATE THE NUMBER DENSITY OF AIR                       XIN 0052
              RHOAIR=ALOSMT*(P(L)/PZERO)/(T(L)/TZERO)                   XIN 0053
              DO 20 K=1,IXMOLS                                          XIN 0054
                  DENM(K,L)=EXPINT(DENX(K,LX-1),DENX(K,LX),A)           XIN 0055
C                                                                       XIN 0056
C                 CONVERT MIXING RATIO (PPMV) TO NUMBER DENSITY         XIN 0057
                  DENM(K,L)=RHOAIR*DENM(K,L)*1.E-6                      XIN 0058
   20         CONTINUE                                                  XIN 0059
          ELSE                                                          XIN 0060
              LX=LX+1                                                   XIN 0061
              GOTO10                                                    XIN 0062
          ENDIF                                                         XIN 0063
   30 CONTINUE                                                          XIN 0064
      RETURN                                                            XIN 0065
      END                                                               XIN 0066
