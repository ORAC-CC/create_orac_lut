      SUBROUTINE WATVAP(P,T)                                            WAT 0001
C*************************************************************          WAT 0002
C                                                                       WAT 0003
C        WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR       WAT 0004
C        UNIFORM DATA INPUT -                                           WAT 0005
C        RE-WRITTEN JULY, 1996 TO CORRECT COMMENTS, CORRECT RELATIVE    WAT 0006
C        HUMIDITY CHECK, SIMPLIFY                                       WAT 0007
C                                                                       WAT 0008
C     JCHAR    JUNIT                                                    WAT 0009
C                                                                       WAT 0010
C    " ",A       10    VOLUME MIXING RATIO (PPMV)                       WAT 0011
C        B       11    NUMBER DENSITY (CM-3)                            WAT 0012
C        C       12    MASS MIXING RATIO (GM(K)/KG(AIR))                WAT 0013
C        D       13    MASS DENSITY (GM M-3)                            WAT 0014
C        E       14    PARTIAL PRESSURE (MB)                            WAT 0015
C        F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY           WAT 0016
C        G       16     "    "     "  (TD IN T(C)) - H2O ONLY           WAT 0017
C        H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY     WAT 0018
C        I       18    AVAILABLE FOR USER DEFINITION                    WAT 0019
C        J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE    WAT 0020
C                                                                       WAT 0021
C     THIS SUBROUTINE COMPUTES THE WATER VAPOR MASS DENSITY (GM M-3)    WAT 0022
C     GIVEN HUMIDITY  # TD = DEW POINT TEMP(K,C), RH = RELATIVE         WAT 0023
C     (PERCENT), PPH2O = WATER VAPOR PARTIAL PRESSURE (MB), DENH2O =    WAT 0024
C     WATER VAPOR MASS DENSITY (GM M-3),AMSMIX = MASS MIXING RATIO      WAT 0025
C     (GM/KG).                                                          WAT 0026
C                     THE FUNCTION DENSAT FOR THE SATURATION            WAT 0027
C     WATER VAPOR DENSITY OVER WATER IS ACCURATE TO BETTER THAN 1       WAT 0028
C     PERCENT FROM -50 TO +50 DEG C. (SEE THE LOWTRAN3 OR 5 REPORT)     WAT 0029
C                                                                       WAT 0030
C       'JUNIT' GOVERNS CHOICE OF UNITS -                               WAT 0031
C                                                                       WAT 0032
C***********************************************************************WAT 0033
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          WAT 0034
      COMMON /CARD1B/ JUNITP,JUNITT,JUNIT1(13),WMOL1(12),WAIR,JLOW      WAT 0035
      COMMON /CONSTN/ PZERO,TZERO,AVOGAD,ALOSMT,GASCON,PLANK,BOLTZ,     WAT 0036
     1    CLIGHT,ADCON,ALZERO,AVMWT,AIRMWT,AMWT(35)                     WAT 0037
      DATA C1/18.9766/,C2/-14.9595/,C3/-2.43882/                        WAT 0038
      DENSAT(ATEMP) = ATEMP*B*EXP(C1+C2*ATEMP+C3*ATEMP**2)*1.0E-6       WAT 0039
C*****                                                                  WAT 0040
      PSS = P / PZERO                                                   WAT 0041
      A = TZERO / T                                                     WAT 0042
      B = AVOGAD / AMWT(1)                                              WAT 0043
      WAIR = ALOSMT * PSS * A                                           WAT 0044
      JUNIT = JUNIT1(1)                                                 WAT 0045
      WMOL  = WMOL1(1)                                                  WAT 0046
      IF (JUNIT .EQ. 10) THEN                                           WAT 0047
C         given volume mixing ratio (ppmv)                              WAT 0048
          WMOL1(1) = WMOL * WAIR / B                                    WAT 0049
      ELSE IF (JUNIT .EQ. 11) THEN                                      WAT 0050
C         given number density                                          WAT 0051
          WMOL1(1) = 1.0E6 * WMOL / B                                   WAT 0052
      ELSE IF (JUNIT .EQ. 12) THEN                                      WAT 0053
C         given mass mixing ratio (gm kg-1)                             WAT 0054
          WMOL1(1) = 1.0E3 * WMOL * WAIR * AIRMWT / AVOGAD              WAT 0055
      ELSE IF (JUNIT .EQ. 13) THEN                                      WAT 0056
C         given mass density (gm m-3)                                   WAT 0057
          WMOL1(1) = WMOL                                               WAT 0058
      ELSE IF (JUNIT .EQ. 14) THEN                                      WAT 0059
C         given water vapor partial pressure (mb)                       WAT 0060
          WMOL1(1) = 1.0E6 * WAIR * WMOL / (B * P)                      WAT 0061
      ELSE IF (JUNIT .EQ. 15) THEN                                      WAT 0062
C         given dewpoint (deg K)                                        WAT 0063
          WMOL1(1) = 1.0E6 * DENSAT(TZERO/WMOL) * WMOL / (B * T)        WAT 0064
      ELSE IF (JUNIT .EQ. 16) THEN                                      WAT 0065
C         given dewpoint (deg C)                                        WAT 0066
          WMOL1(1)=1.E6*DENSAT(TZERO/(TZERO+WMOL))*(TZERO+WMOL)/(B*T)   WAT 0067
      ELSE IF (JUNIT .EQ. 17) THEN                                      WAT 0068
C         given relative humidity (percent)                             WAT 0069
          WMOL1(1) = 1.0E6 * DENSAT(A) * (WMOL/100.0) / B               WAT 0070
      ELSE                                                              WAT 0071
C         given invalid JUNIT                                           WAT 0072
          WRITE (IPR,100) JUNIT                                         WAT 0073
  100     FORMAT(/,'  **** ERROR IN WATVAP ****, JUNIT = ',I5)          WAT 0074
          STOP 'JUNIT'                                                  WAT 0075
      END IF                                                            WAT 0076
C     insure that relative humidity is valid                            WAT 0077
      RHP = 100.0 * WMOL1(1) * B / (1.0E6 * DENSAT(A))                  WAT 0078
      IF(RHP.GT.100.0) WRITE(IPR,200) RHP                               WAT 0079
  200 FORMAT(/,' ********WARNING (FROM WATVAP) # RELATIVE HUMIDTY = ',  WAT 0080
     1    G10.3,' IS GREATER THAN 100 PERCENT')                         WAT 0081
C                                                                       WAT 0082
      RETURN                                                            WAT 0083
      END                                                               WAT 0084
