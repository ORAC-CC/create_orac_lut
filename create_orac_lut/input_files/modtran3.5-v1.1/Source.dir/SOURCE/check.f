      SUBROUTINE CHECK(V,IV,KEY)                                        CHK 0001
C                                                                       CHK 0002
C      UNITS CONVERSION FOR P AND T                                     CHK 0003
C                                                                       CHK 0004
C     V = P OR T     AND  IV =JUNITP(I.E. MB,ATM,TORR)                  CHK 0005
C                            =JUNITT(I.E. DEG K OR C)                   CHK 0006
C                            =JUNITR(I.E. KM,M,OR CM)                   CHK 0007
C                                                                       CHK 0008
      DATA PMB/1013.25/,PTORR/760./,DEGK/273.15/                        CHK 0009
      IF(IV.LE.10) RETURN                                               CHK 0010
      GO TO(100,200,300) KEY                                            CHK 0011
C                                                                       CHK 0012
C      PRESSURE CONVERSIONS                                             CHK 0013
C                                                                       CHK 0014
  100 IF(IV.EQ.11)GO TO 110                                             CHK 0015
      IF(IV.EQ.12)GO TO 120                                             CHK 0016
      STOP'CHECK(P)'                                                    CHK 0017
  110 V=V*PMB                                                           CHK 0018
      RETURN                                                            CHK 0019
  120 V=V*PMB/PTORR                                                     CHK 0020
      RETURN                                                            CHK 0021
C                                                                       CHK 0022
C      TEMPERATURE COMVERSIONS                                          CHK 0023
C                                                                       CHK 0024
  200 IF(IV.GT.11)STOP'CHECK(T)'                                        CHK 0025
      V=V+DEGK                                                          CHK 0026
      RETURN                                                            CHK 0027
C                                                                       CHK 0028
C      RANGE CONVERSIONS                                                CHK 0029
C                                                                       CHK 0030
  300 IF(IV.EQ.11)GO TO 310                                             CHK 0031
      IF(IV.EQ.12)GO TO 320                                             CHK 0032
      STOP'CHECK(R)'                                                    CHK 0033
  310 V=V/1.E3                                                          CHK 0034
      RETURN                                                            CHK 0035
  320 V=V/1.E5                                                          CHK 0036
      RETURN                                                            CHK 0037
      END                                                               CHK 0038
