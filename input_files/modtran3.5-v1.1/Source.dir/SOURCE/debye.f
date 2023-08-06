      SUBROUTINE DEBYE(WAVL,TC,KEY,RE,AI)                               DBY 0001
CCC                                                                     DBY 0002
CCC    CALCULATES WAVENUMBER DEPENDENCE OF DIELECTRIC CONSTANT          DBY 0003
CCC    OF WATER                                                         DBY 0004
C                                                                       DBY 0005
C     LIST COMMONS:                                                     DBY 0006
C       PI       THE CONSTANT PI                                        DBY 0007
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       DBY 0008
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       DBY 0009
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         DBY 0010
      REAL PI,DEG,BIGNUM,BIGEXP                                         DBY 0011
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                DBY 0012
      T=TC+273.15                                                       DBY 0013
      IF(KEY.NE.0) GO TO 1                                              DBY 0014
      GO TO 2                                                           DBY 0015
    1 EFIN=5.27137+.0216474*TC-.00131198*TC*TC                          DBY 0016
      ALPHA=-16.8129/T+.0609265                                         DBY 0017
      TAU=.00033836*EXP(2513.98/T)                                      DBY 0018
      SIG=12.5664E+08                                                   DBY 0019
      ES=78.54*(1.-.004579*(TC-25.)+.0000119*(TC-25.)**2-.000000028*    DBY 0020
     1(TC-25.)**3)                                                      DBY 0021
      GO TO 3                                                           DBY 0022
    2 EFIN=3.168                                                        DBY 0023
      ALPHA=.00023*TC*TC+.0052*TC+.288                                  DBY 0024
      T125 = 12500./(T*1.9869)                                          DBY 0025
      IF(T125. LE.BIGEXP) THEN                                          DBY 0026
           SIG=1.26*EXP(-T125)                                          DBY 0027
      ELSE                                                              DBY 0028
           SIG = 0.                                                     DBY 0029
      ENDIF                                                             DBY 0030
      TAU=9.990288E-05*EXP(13200./(T*1.9869))                           DBY 0031
      ES=3.168+.15*TC*TC+2.5*TC+200.                                    DBY 0032
    3 C1=TAU/WAVL                                                       DBY 0033
CCC                                                                     DBY 0034
CCC    TEMPORARY FIX TO CLASSICAL DEBYE EQUATION                        DBY 0035
CCC    TO HANDLE ZERO CM-1 PROBLEM                                      DBY 0036
CCC                                                                     DBY 0037
      ALPHA=0.0                                                         DBY 0038
      SIG=0.0                                                           DBY 0039
CCC                                                                     DBY 0040
      C2=1.5708*ALPHA                                                   DBY 0041
      DEM=1.+2.*C1**(1.-ALPHA)*SIN(C2)+C1**(2.*(1.-ALPHA))              DBY 0042
      E1=EFIN+(ES-EFIN)*(1.+(C1**(1.-ALPHA)*SIN(C2)))/DEM               DBY 0043
      IF(KEY.NE.0.AND.WAVL.GE.300.) E1=87.53-0.3956*TC                  DBY 0044
      IF(KEY.NE.0 .AND. WAVL.GE.300.) E1=ES                             DBY 0045
      E2=(ES-EFIN)*C1**(1.-ALPHA)*COS(C2)/DEM+SIG*WAVL/18.8496E+10      DBY 0046
CCC                                                                     DBY 0047
CCC    PERMANENT FIX TO CLSSICAL DEBYE EQUATION                         DBY 0048
CCC    TO HANDLE ZERO CM-1 PROBLEM                                      DBY 0049
CCC                                                                     DBY 0050
      E1=EFIN+(ES-EFIN)/(1.0+C1**2)                                     DBY 0051
CCC                                                                     DBY 0052
      E2=((ES-EFIN)*C1)/(1.0+C1**2)                                     DBY 0053
CCC                                                                     DBY 0054
      RE=SQRT((E1+SQRT(E1*E1+E2*E2))/2.)                                DBY 0055
      AI=-E2/(2.*RE)                                                    DBY 0056
      RETURN                                                            DBY 0057
      END                                                               DBY 0058
