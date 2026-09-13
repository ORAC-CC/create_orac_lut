      SUBROUTINE INDX (WAVL,TC,KEY,REIL,ZIMAG)                          INX 0001
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * INX 0002
C * *                                                                   INX 0003
C * * WAVELENGTH IS IN CENTIMETERS.  TEMPERATURE IS IN DEG. C.      * * INX 0004
C * *                                                                   INX 0005
C * * KEY IS SET TO 1 IN SUBROUTINE GAMFOG                          * * INX 0006
C * *                                                                   INX 0007
C * * KEY IS SET TO 0 IN SUBROUTINE GAMFOG    FOR CIRRUS            * * INX 0008
C * *                                                                   INX 0009
C * * REAL IS THE REAL PART OF THE REFRACTIVE INDEX.                * * INX 0010
C * *                                                                   INX 0011
C * * ZIMAG IS THE IMAGINARY PART OF THE REFRACTIVE INDEX IT IS     * * INX 0012
C * *                                                                   INX 0013
C * * RETURNED NEG. I.E.  M= REAL - I*ZIMAG  .                      * * INX 0014
C * *                                                                   INX 0015
C * * A SERIES OF CHECKS ARE MADE AND WARNINGS GIVEN.               * * INX 0016
C * *                                                                   INX 0017
C * * RAY APPLIED OPTICS VOL 11,NO.8,AUG 72, PG. 1836-1844          * * INX 0018
C * *                                                                   INX 0019
C * * CORRECTIONS HAVE BEEN MADE TO RAYS ORIGINAL PAPER             * * INX 0020
C * *                                                                   INX 0021
C * *                                                                   INX 0022
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * INX 0023
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          INX 0024
      R1=0.0                                                            INX 0025
      R2=0.0                                                            INX 0026
      IF(WAVL.LT..0001) WRITE(IPR,1)                                    INX 0027
      IF(TC.LT.-20.) WRITE(IPR,2)                                       INX 0028
    1 FORMAT(///,30X,'ATTEMPTING TO EVALUATE FOR A WAVELENGTH LESS THAN INX 0029
     1ONE MICRON',//)                                                   INX 0030
    2 FORMAT(///,30X,'ATTEMPTING TO EVALUATE FOR A TEMPERATURE LESS THANINX 0031
     2 -20. DEGREES CENTIGRADE',//)                                     INX 0032
      CALL DEBYE(WAVL,TC,KEY,REIL,ZIMAG)                                INX 0033
C * *  TABLE 3 WATER PG. 1840                                           INX 0034
    4 IF(WAVL.GT..034) GO TO 5                                          INX 0035
      GO TO 7                                                           INX 0036
    5 IF(WAVL.GT..1) GO TO 11                                           INX 0037
    6  R2 =DOP(WAVL,1.83899,1639.,52340.4,10399.2,588.24,345005.,       INX 0038
     3259913.,161.29,43319.7,27661.2)                                   INX 0039
      R2=R2+R2*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)                  INX 0040
      REIL=REIL*(WAVL-.034)/.066+R2*(.1-WAVL)/.066                      INX 0041
      GO TO 11                                                          INX 0042
    7 IF(WAVL.GT..0006) GO TO 8                                         INX 0043
      GO TO 10                                                          INX 0044
    8 REIL=DOP(WAVL,1.83899,1639.,52340.4,10399.2,588.24,345005.,       INX 0045
     4259913.,161.29,43319.7,27661.2)                                   INX 0046
      REIL=REIL+REIL*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)            INX 0047
      IF(WAVL.GT..0007) GO TO 11                                        INX 0048
    9 R1=DOP(WAVL,1.79907,3352.27,99.914E+04,15.1963E+04,1639.,50483.5, INX 0049
     59246.27,588.24,84.4697E+04,10.7615E+05)                           INX 0050
      R1=R1+R1*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)                  INX 0051
      REIL=R1*(.0007-WAVL)/.0001+REIL*(WAVL-.0006)/.0001                INX 0052
      GO TO 11                                                          INX 0053
   10 REIL=DOP(WAVL,1.79907,3352.27,99.914E+04,15.1963E+04,1639.,       INX 0054
     650483.5,9246.27,588.24,84.4697E+04,10.7615E+05)                   INX 0055
      REIL=REIL+REIL*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)            INX 0056
C * *  TABLE 2 WATER PG. 1840                                           INX 0057
   11 IF(WAVL.GE..3) GO TO 57                                           INX 0058
      IF(WAVL.GE..03) GO TO 12                                          INX 0059
      GO TO 13                                                          INX 0060
   12 ZIMAG=ZIMAG+AB(WAVL,.25,300.,.47,3.)+AB(WAVL,.39,17.,.45,1.3)+    INX 0061
     7AB(WAVL,.41,62.,.35,1.7)                                          INX 0062
      GO TO 57                                                          INX 0063
   13 IF(WAVL.GE..0062) GO TO 14                                        INX 0064
      GO TO 15                                                          INX 0065
   14 ZIMAG=ZIMAG+AB(WAVL,.41,62.,.35,1.7)+AB(WAVL,.39,17.,.45,1.3)+    INX 0066
     8AB(WAVL,.25,300.,.4,2.)                                           INX 0067
      GO TO 57                                                          INX 0068
   15 IF(WAVL.GE..0017) GO TO 16                                        INX 0069
      GO TO 17                                                          INX 0070
   16 ZIMAG=ZIMAG+AB(WAVL,.39,17.,.45,1.3)+AB(WAVL,.41,62.,.22,1.8)+    INX 0071
     9AB(WAVL,.25,300.,.4,2.)                                           INX 0072
      GO TO 57                                                          INX 0073
   17 IF(WAVL.GE..00061) GO TO 18                                       INX 0074
      GO TO 19                                                          INX 0075
   18 ZIMAG=ZIMAG+AB(WAVL,.12,6.1,.042,.6)+AB(WAVL,.39,17.,.165,2.4)+   INX 0076
     1AB(WAVL,.41,62.,.22,1.8)                                          INX 0077
      GO TO 57                                                          INX 0078
   19 IF(WAVL.GE..000495) GO TO 20                                      INX 0079
      GO TO 21                                                          INX 0080
   20 ZIMAG=ZIMAG+AB(WAVL,.01,4.95,.05,1.)+AB(WAVL,.12,6.1,.009,2.)     INX 0081
      GO TO 57                                                          INX 0082
   21 IF(WAVL.GE..000297) GO TO 22                                      INX 0083
      GO TO 23                                                          INX 0084
   22 ZIMAG=ZIMAG+AB(WAVL,.27,2.97,.04,2.)+AB(WAVL,.01,4.95,.06,1.)     INX 0085
      GO TO 57                                                          INX 0086
   23 ZIMAG=ZIMAG+AB(WAVL,.27,2.97,.025,2.)+AB(WAVL,.01,4.95,.06,1.)    INX 0087
   57 CONTINUE                                                          INX 0088
      RETURN                                                            INX 0089
      END                                                               INX 0090
