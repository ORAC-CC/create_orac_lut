      SUBROUTINE PSIECA(OLAT,OLONG,SLAT,SLONG,PSI,ECA)                  ECA 0001
C                                                                       ECA 0002
C     THIS ROUTINE CALCULATES THE EARTH CENTER ANGLE AND THE PATH       ECA 0003
C     AZIMUTH FROM AN OBSERVER TO A SOURCE POINT GIVEN LATITUDE AND     ECA 0004
C     LONGITUDE INFORMATION.                                            ECA 0005
C                                                                       ECA 0006
C     DECLARE INPUTS                                                    ECA 0007
C       OLAT   OBSERVER LATITUDE [DEG NORTH, -90 TO 90]                 ECA 0008
C       OLONG  OBSERVER LONGITUDE [DEG WEST, 0 TO 360]                  ECA 0009
C       SLAT   SOURCE POINT LATITUDE [DEG NORTH, -90 TO 90]             ECA 0010
C       SLONG  SOURCE POINT LONGITUDE [DEG WEST, 0 TO 360]              ECA 0011
      REAL OLAT,OLONG,SLAT,SLONG                                        ECA 0012
C                                                                       ECA 0013
C     DECLARE OUTPUTS                                                   ECA 0014
C       PSI    PATH AZIMUTH AT OBSERVER [DEG EAST OF NORTH, -180 TO 180]ECA 0015
C       ECA    EARTH CENTER ANGLE [DEG, 0 TO 180]                       ECA 0016
      REAL PSI,ECA                                                      ECA 0017
C                                                                       ECA 0018
C     NOTE:  OLONG AND PSI ARE NOT WELL DEFINED AT THE NORTH (OLAT=90)  ECA 0019
C            AND SOUTH (OLAT=-90) POLES.  THIS ROUTINE RETURNS THE      ECA 0020
C            VALUE OF PATH AZIMUTH OBTAINED IN THE LIMIT AS |OLAT|      ECA 0021
C            APPROACHES 90 FROM BELOW.                                  ECA 0022
C                                                                       ECA 0023
C     DECLARE LOCAL VARIABLES                                           ECA 0024
      REAL DEG,SOLAT,COLAT,SSLAT,CSLAT,DLONG,SDLONG,CDLONG,COSECA,X,Y   ECA 0025
C                                                                       ECA 0026
C     LIST DATA                                                         ECA 0027
      DATA DEG/57.2957795/                                              ECA 0028
C                                                                       ECA 0029
C     TREAT SAME LONGITUDE AS SPECIAL CASE                              ECA 0030
      IF(OLONG.EQ.SLONG)THEN                                            ECA 0031
          IF(OLAT.LE.SLAT)THEN                                          ECA 0032
              ECA=SLAT-OLAT                                             ECA 0033
              PSI=0.                                                    ECA 0034
          ELSE                                                          ECA 0035
              ECA=OLAT-SLAT                                             ECA 0036
              PSI=180.                                                  ECA 0037
          ENDIF                                                         ECA 0038
          RETURN                                                        ECA 0039
      ENDIF                                                             ECA 0040
      SOLAT=SIN(OLAT/DEG)                                               ECA 0041
      COLAT=COS(OLAT/DEG)                                               ECA 0042
      SSLAT=SIN(SLAT/DEG)                                               ECA 0043
      CSLAT=COS(SLAT/DEG)                                               ECA 0044
      DLONG=(OLONG-SLONG)/DEG                                           ECA 0045
      SDLONG=SIN(DLONG)                                                 ECA 0046
      CDLONG=COS(DLONG)                                                 ECA 0047
      COSECA=COLAT*CSLAT*CDLONG+SOLAT*SSLAT                             ECA 0048
      IF(COSECA.GE.1.)THEN                                              ECA 0049
          ECA=0.                                                        ECA 0050
      ELSEIF(COSECA.GT.-1.)THEN                                         ECA 0051
          ECA=DEG*ACOS(COSECA)                                          ECA 0052
      ELSE                                                              ECA 0053
          ECA=180.                                                      ECA 0054
      ENDIF                                                             ECA 0055
      X=COLAT*SSLAT-SOLAT*CSLAT*CDLONG                                  ECA 0056
      Y=CSLAT*SDLONG                                                    ECA 0057
      IF(X.EQ.0.)THEN                                                   ECA 0058
          PSI=90.                                                       ECA 0059
          IF(Y.LT.0.)PSI=-90.                                           ECA 0060
      ELSE                                                              ECA 0061
          PSI=DEG*ATAN(Y/X)                                             ECA 0062
          IF(X.LT.0)PSI=PSI+180.                                        ECA 0063
          IF(PSI.GT.180.)PSI=PSI-360.                                   ECA 0064
      ENDIF                                                             ECA 0065
      RETURN                                                            ECA 0066
      END                                                               ECA 0067
