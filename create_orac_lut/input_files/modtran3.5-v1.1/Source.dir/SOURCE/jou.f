      FUNCTION   JOU(CHAR)                                              JOU 0001
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,ISCRCH                      JOU 0002
C                                                                       JOU 0003
      CHARACTER*1 CHAR,HOLVEC(22)                                       JOU 0004
      DIMENSION INDX1(22)                                               JOU 0005
      DATA  HOLVEC                                                      JOU 0006
     X /'1','2','3','4','5','6','0','0','0','0',' ','A',                JOU 0007
     X  'B','C','D','E','F','G','H','I','J','K'/                        JOU 0008
      DATA  INDX1                                                       JOU 0009
     X /  1,  2,  3,  4,  5,  6,  0,  0,  0,  0, 10, 10,                JOU 0010
     X   11, 12, 13, 14, 15, 16, 17, 18, 19, 20/                        JOU 0011
C                                                                       JOU 0012
       INDX=0                                                           JOU 0013
      DO 100 I=1,22                                                     JOU 0014
       IF (HOLVEC(I) .NE. CHAR) GO TO 100                               JOU 0015
       INDX=INDX1(I)                                                    JOU 0016
       GO TO 110                                                        JOU 0017
100   CONTINUE                                                          JOU 0018
110   IF (INDX .EQ. 0) THEN                                             JOU 0019
        WRITE(IPR,910) CHAR                                             JOU 0020
910     FORMAT('0 INVALID PARAMETER :',2X,A1)                           JOU 0021
        STOP ' JOU: BAD PARAM '                                         JOU 0022
      END IF                                                            JOU 0023
920   FORMAT(5X,A1,I5)                                                  JOU 0024
      JOU=INDX                                                          JOU 0025
                    RETURN                                              JOU 0026
      END                                                               JOU 0027
