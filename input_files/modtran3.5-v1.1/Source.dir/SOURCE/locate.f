      SUBROUTINE LOCATE(OLAT,OLONG,OPSI,DBETA,RLAT,RLONG)               LOC 0001
C                                                                       LOC 0002
C     THIS ROUTINE CALCULATES THE LATITUDE AND LONGITUDE AT A POINT "R" LOC 0003
C     ALONG A LINE OF SIGHT FROM AN OBSERVER AT "O"                     LOC 0004
C                                                                       LOC 0005
C     INPUT                                                             LOC 0006
C       OLAT    THE LATITUDE AT "O" [DEG NORTH, -90 TO 90]              LOC 0007
C       OLONG   THE LONGITUDE AT "O" [DEG WEST, 0 TO 360]               LOC 0008
C       OPSI    THE PATH AZIMUTH AT "O" [DEG EAST OF NORTH, 0 TO 360]   LOC 0009
C       DBETA   THE EARTH CENTER ANGLE FROM "O" TO "R" [DEG, 0 TO 180]  LOC 0010
C                                                                       LOC 0011
C     OUTPUT                                                            LOC 0012
C       RLAT    THE LATITUDE AT "R" [DEG, -90 TO 90]                    LOC 0013
C       RLONG   THE LONGITUDE AT "R" [DEG WEST, 0 TO 360]               LOC 0014
      REAL OLAT,OLONG,OPSI,DBETA,RLAT,RLONG                             LOC 0015
C                                                                       LOC 0016
C     CAUTION:  OLONG AND OPSI ARE NOT WELL DEFINED AT THE NORTH        LOC 0017
C               (OLAT=90) AND SOUTH (OLAT=-90) POLES.  THIS ROUTINE     LOC 0018
C               RETURNS THE VALUES OF LATITUDE AND LONGITUDE OBTAINED   LOC 0019
C               IN THE LIMIT AS OLAT APPROACHES + OR - 90.              LOC 0020
C                                                                       LOC 0021
C     DECLARE LOCAL VARIABLES                                           LOC 0022
      REAL DEG,SOLAT,COLAT,SOLONG,COLONG,SDBETA,CDBETA,                 LOC 0023
     1  SOPSI,COPSI,SRLAT,XP,YP,X,Y                                     LOC 0024
      DATA DEG/57.2957795/                                              LOC 0025
      SOLAT=SIN(OLAT/DEG)                                               LOC 0026
      COLAT=COS(OLAT/DEG)                                               LOC 0027
      SOLONG=SIN(OLONG/DEG)                                             LOC 0028
      COLONG=COS(OLONG/DEG)                                             LOC 0029
      SDBETA=SIN(DBETA/DEG)                                             LOC 0030
      CDBETA=COS(DBETA/DEG)                                             LOC 0031
      SOPSI=SIN(OPSI/DEG)                                               LOC 0032
      COPSI=COS(OPSI/DEG)                                               LOC 0033
      SRLAT=SOLAT*CDBETA+COLAT*SDBETA*COPSI                             LOC 0034
      IF(SRLAT.GE.1.)THEN                                               LOC 0035
          RLAT=90.                                                      LOC 0036
      ELSEIF(SRLAT.GT.-1.)THEN                                          LOC 0037
          RLAT=DEG*ASIN(SRLAT)                                          LOC 0038
      ELSE                                                              LOC 0039
          RLAT=-90.                                                     LOC 0040
      ENDIF                                                             LOC 0041
      XP=COLAT*CDBETA-SOLAT*SDBETA*COPSI                                LOC 0042
      YP=SDBETA*SOPSI                                                   LOC 0043
      X=COLONG*XP+SOLONG*YP                                             LOC 0044
      Y=SOLONG*XP-COLONG*YP                                             LOC 0045
      IF(X.EQ.0.)THEN                                                   LOC 0046
          RLONG=90.                                                     LOC 0047
          IF(Y.LT.0.)RLONG=270.                                         LOC 0048
      ELSE                                                              LOC 0049
          RLONG=DEG*ATAN(Y/X)                                           LOC 0050
          IF(X.LT.0)RLONG=RLONG+180.                                    LOC 0051
          IF(RLONG.LT.0.)RLONG=RLONG+360.                               LOC 0052
      ENDIF                                                             LOC 0053
      RETURN                                                            LOC 0054
      END                                                               LOC 0055
