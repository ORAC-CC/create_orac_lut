      SUBROUTINE SUBSOL(THETAS,PHIS,TIME,IDAY)                          SUB 0001
C                                                                       SUB 0002
C     SUBROUTINE SUBSOL CALCULATES THE SUBSOLAR POINT ANGLES            SUB 0003
C     THETAS AND PHIS BASED UPON IDAY AND TIME. SINCE EACH              SUB 0004
C     YEAR IS 365.25 DAYS LONG THE EXACT VALUE OF THE DECLINATION       SUB 0005
C     ANGLE CHANGES FROM YEAR TO YEAR.  FOR PRECISE VALUES CONSULT      SUB 0006
C     ' THE AMERICAN EPHEMERIS AND NAUTICAL ALMANAC' PUBLISHED YEARLY   SUB 0007
C     BY THE U.S. GOVT. PRINTING OFFICE.  ALSO, THE SOLAR POSITION      SUB 0008
C     IS CHARACTERIZED BY 25 POINTS BELOW; THIS SHOULD PREDICT THE SUBSOSUB 0009
C     ANGLES WITHIN ONE DEGREE.  FOR INCREASED ACCURACY ADD MORE DATA   SUB 0010
C     POINTS                                                            SUB 0011
C                                                                       SUB 0012
C     THE EQUATION OF TIME, EQT, IS IN MINUTES                          SUB 0013
C     THE DECLINATION ANGLE, DEC IS IN DEGREES                          SUB 0014
C                                                                       SUB 0015
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          SUB 0016
      DIMENSION NDAY(25),EQT(25),DEC(25)                                SUB 0017
      DATA NDAY /1,9,21,32,44,60,91,121,141,152,160,172,182,            SUB 0018
     1 190,202,213,244,274,305,309,325,335,343,355,366/                 SUB 0019
      DATA DEC /-23.07,-22.22,-20.08,-17.32,-13.62,-7.88,4.23,          SUB 0020
     1 14.83, 20.03,21.95,22.87,23.45,23.17,22.47,20.63,18.23,8.58,     SUB 0021
     2 -2.88,-14.18,-15.45,-19.75,-21.68,-22.75,-23.43,-23.07/          SUB 0022
      DATA EQT /-3.23,-6.83,-11.17,-13.57,-14.33,-12.63,-4.2,           SUB 0023
     1 2.83,3.57,2.45,1.10,-1.42,-3.52,-4.93,-6.25,-6.28,-0.25,         SUB 0024
     2 10.02,16.35,16.38,14.3,11.27,8.02,2.32,-3.23/                    SUB 0025
      IF(IDAY.LT.1.OR.IDAY.GT.366)  GO TO 900                           SUB 0026
      IF(TIME.LT.0.0.OR.TIME.GT.24.0)  GO TO 910                        SUB 0027
      DO 10 I=1,25                                                      SUB 0028
      IF(NDAY(I).EQ.IDAY) GO TO 30                                      SUB 0029
10    IF(NDAY(I).GT.IDAY) GO TO 20                                      SUB 0030
20    I=I-1                                                             SUB 0031
      EQTIME=EQT(I)+(EQT(I+1)-EQT(I))*(IDAY-NDAY(I))/(NDAY(I+1)-NDAY(I))SUB 0032
      DECANG=DEC(I)+(DEC(I+1)-DEC(I))*(IDAY-NDAY(I))/(NDAY(I+1)-NDAY(I))SUB 0033
      GO TO 40                                                          SUB 0034
30    EQTIME=EQT(I)                                                     SUB 0035
      DECANG=DEC(I)                                                     SUB 0036
40    THETAS=DECANG                                                     SUB 0037
      EQTIME=EQTIME/60.0                                                SUB 0038
      PHIS=15.0*(TIME+EQTIME)-180.0                                     SUB 0039
      IF(PHIS.LT.0.0) PHIS=PHIS+360.0                                   SUB 0040
      RETURN                                                            SUB 0041
900   WRITE(IPR,901) IDAY                                               SUB 0042
901   FORMAT(' FROM SUBSOL - IDAY OUT OF RANGE, IDAY=',I6)              SUB 0043
      STOP                                                              SUB 0044
910   WRITE(IPR,902) TIME                                               SUB 0045
902   FORMAT(' FROM SUBSOL - TIME OUT OF RANGE, TIME=',E12.5)           SUB 0046
      STOP                                                              SUB 0047
      END                                                               SUB 0048
