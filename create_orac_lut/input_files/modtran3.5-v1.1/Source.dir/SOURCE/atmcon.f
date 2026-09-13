      BLOCK DATA ATMCON                                                 BAT 0001
C***********************************************************************BAT 0002
C     THIS SUBROUTINE INITIALIZES THE CONSTANTS USED IN THE             BAT 0003
C     PROGRAM. CONSTANTS RELATING TO THE ATMOSPHERIC PROFILES ARE STOREDBAT 0004
C     IN BLOCK DATA MLATMB.                                             BAT 0005
      COMMON /CONSTN/ PZERO,TZERO,AVOGAD,ALOSMT,GASCON,PLANK,BOLTZ,     BAT 0006
     1    CLIGHT,ADCON,ALZERO,AVMWT,AIRMWT,AMWT(35)                     BAT 0007
      COMMON/HMOLC/HMOLC(35)                                            BAT 0008
      CHARACTER*8 HMOLC                                                 BAT 0009
      DATA PZERO/1013.25/,TZERO/273.15/                                 BAT 0010
      DATA AVOGAD/6.022045E+23/,ALOSMT/2.68675E+19/,                    BAT 0011
     1    GASCON/8.31441E+7/,PLANK/6.626176E-27/,BOLTZ/1.380662E-16/,   BAT 0012
     2    CLIGHT/2.9979246E10/                                          BAT 0013
C                                                                       BAT 0014
C**   ALZERO IS THE MEAN LORENTZ HALFWIDTH AT PZERO AND 296.0 K.        BAT 0015
C**   AVMWT IS THE MEAN MOLECULAR WEIGHT USED TO AUTOMATICALLY          BAT 0016
C**   GENERATE THE FASCODE BOUNDARIES IN AUTLAY                         BAT 0017
C                                                                       BAT 0018
      DATA ALZERO/0.1/,AVMWT/36.0/                                      BAT 0019
C                                                                       BAT 0020
C**   ORDER OF MOLECULES H2O(1), CO2(2), O3(3), N2O(4), CO(5), CH4(6),  BAT 0021
C**       O2(7), NO(8), SO2(9), NO2(10), NH3(11), HNO3(12), OH(13),     BAT 0022
C**       HF(14 ), HCL(15), HBR(16), HI(17), CLO(18), OCS(19), H2CO(20) BAT 0023
C**       HOCL(21), N2(22), HCN(23), CH3CL(24), H2O2(25), C2H2(26),     BAT 0024
C**       C2H6(27), PH3(28),COF2(29),   SF6(30)                         BAT 0025
C                                                                       BAT 0026
       DATA HMOLC   /'  H2O   ','  CO2   ','   O3   ',                  BAT 0027
     1    '  N2O   ','   CO   ','  CH4   ','   O2   ',                  BAT 0028
     2    '   NO   ','  SO2   ','  NO2   ','  NH3   ',                  BAT 0029
     3    ' HNO3   ','   OH   ','   HF   ','  HCL   ',                  BAT 0030
     4    '  HBR   ','   HI   ','  CLO   ','  OCS   ',                  BAT 0031
     5    ' H2CO   ',' HOCL   ','   N2   ','  HCN   ',                  BAT 0032
     6    ' CH3CL  ',' H2O2   ',' C2H2   ',' C2H6   ',                  BAT 0033
     7    '  PH3   ',' COF2   ','  SF6   ',' H2S    ',                  BAT 0034
     8    ' HCOOH  ',3*'        '/                                      BAT 0035
C                                                                       BAT 0036
C**   MOLECULAR WEIGHTS                                                 BAT 0037
C                                                                       BAT 0038
      DATA AIRMWT/28.964/,AMWT/18.015,44.010,47.998,44.01,28.011,       BAT 0039
     1    16.043,31.999,30.01,64.06,46.01,17.03,63.01,17.00,20.01,      BAT 0040
     2    36.46,80.92,127.91,51.45,60.08,30.03,52.46,28.014,            BAT 0041
     3    27.03, 50.49, 34.01, 26.03, 30.07, 34.00,66.0,146.,           BAT 0042
     4    34.08,46.016,3*0./                                            BAT 0043
      END                                                               BAT 0044
