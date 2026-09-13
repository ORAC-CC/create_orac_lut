      INTEGER FUNCTION NUNIT()                                          UNT 0001
C                                                                       UNT 0002
C     THIS INTEGER FUNCTION DETERMINES AN UNUSED                        UNT 0003
C     FILE UNIT NUMBER BETWEEN 10 AND 99.                               UNT 0004
C                                                                       UNT 0005
C     COMMON BLOCKS                                                     UNT 0006
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               UNT 0007
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           UNT 0008
C                                                                       UNT 0009
C     LOCAL VARIABLES                                                   UNT 0010
      LOGICAL LOPEN                                                     UNT 0011
C                                                                       UNT 0012
C     FIND UNUSED UNIT NUMBER AND RETURN.                               UNT 0013
      DO 10 NUNIT=99,10,-1                                              UNT 0014
          INQUIRE(UNIT=NUNIT,OPENED=LOPEN)                              UNT 0015
          IF(.NOT.LOPEN)RETURN                                          UNT 0016
   10 CONTINUE                                                          UNT 0017
C                                                                       UNT 0018
C     NO FILE UNIT NUMBER AVAILABLE                                     UNT 0019
      WRITE(IPR ,'(/A)')                                                UNT 0020
     1  ' Error in routine NUNIT:  No file unit numbers available?'     UNT 0021
      STOP ' Error in routine NUNIT:  No file unit numbers available?'  UNT 0022
      END                                                               UNT 0023
