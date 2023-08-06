      REAL FUNCTION SOURCE(IV,ISOURC,IDAY,ANGLE)                        SRC 0001
C                                                                       SRC 0002
C     THIS ROUTINE RETURNS SOURCE IRRADIANCE [W CM-2 / CM-1].           SRC 0003
C     CORRECTIONS ARE MADE FOR THE SUN'S ELLIPTIC ORBIT.  IF THE        SRC 0004
C     SOURCE IS THE MOON RATHER THAN THE EARTH,  THE PHASE ANGLE        SRC 0005
C     BETWEEN THE SUN, MOON AND EARTH IS TAKEN IN ACCOUNT.              SRC 0006
C                                                                       SRC 0007
C     DECLARE INPUTS.                                                   SRC 0008
C       IV       FREQUENCY [CM-1].                                      SRC 0009
C       ISOURC   SOURCE FLAG [0=SUN AND 1=MOON].                        SRC 0010
C       IDAY     DAY OF YEAR [0-366, DEFAULT (0) IS DAY 91].            SRC 0011
C       ANGLE    LUNAR PHASE ANGLE [0 TO 180 DEG].                      SRC 0012
      INTEGER IV,ISOURC,IDAY                                            SRC 0013
      REAL ANGLE                                                        SRC 0014
C                                                                       SRC 0015
C     DECLARE PARAMETERS                                                SRC 0016
      INTEGER MAXSUN                                                    SRC 0017
      PARAMETER(MAXSUN=50000)                                           SRC 0018
C                                                                       SRC 0019
C     LIST COMMON BLOCKS                                                SRC 0020
      INTEGER ICALL                                                     SRC 0021
      REAL FPHS,FALB,FORBIT                                             SRC 0022
      COMMON/ICLL/ICALL,FPHS,FALB,FORBIT                                SRC 0023
      REAL SUN                                                          SRC 0024
      COMMON/SOLAR1/SUN(0:MAXSUN)                                       SRC 0025
C                                                                       SRC 0026
C     DECLARE LOCAL VARIABLE                                            SRC 0027
      INTEGER I,IM1                                                     SRC 0028
      REAL V                                                            SRC 0029
C                                                                       SRC 0030
C     DECLARE DATA                                                      SRC 0031
      INTEGER NDAY(13)                                                  SRC 0032
      REAL RAT(13),PHS(0:16),ALB(29)                                    SRC 0033
      DATA NDAY/1,32,60,91,121,152,182,213,244,274,305,335,366/         SRC 0034
      DATA RAT/1.034,1.030,1.019,1.001,.985,.972,.967,.971,.982,        SRC 0035
     1    .998,1.015,1.029,1.034/                                       SRC 0036
      DATA PHS/100.,73.2,57.8,42.3,32.0,23.3,16.7,12.4,8.7,6.7,         SRC 0037
     1    4.7,3.6,2.4,1.2,0.9,0.4,.002/                                 SRC 0038
      DATA ALB/.001,.01,.03,.075,.1,.13,.155,.17,.178,.185,.2,.211,     SRC 0039
     1    .231,.25,.275,.289,.285,.287,.3,.29,.3,.31,.313,.319,.329     SRC 0040
     1    ,.339,.345,.350,.4/                                           SRC 0041
C                                                                       SRC 0042
C     CHECK FOR LUNAR SOURCE.                                           SRC 0043
      IF(ISOURC.EQ.1)THEN                                               SRC 0044
C                                                                       SRC 0045
C         CHECK FOR INITIAL CALL.                                       SRC 0046
          IF(ICALL.NE.1)THEN                                            SRC 0047
C                                                                       SRC 0048
C             LUNAR PHASE ANGLE FACTOR                                  SRC 0049
              FPHS=0.                                                   SRC 0050
              I=ANGLE/10                                                SRC 0051
              IF(I.LT.16)FPHS=PHS(I)+(ANGLE/10-I)*(PHS(I+1)-PHS(I))     SRC 0052
          ENDIF                                                         SRC 0053
C                                                                       SRC 0054
C         GEOMETRICAL ALBEDO OF THE MOON                                SRC 0055
          FALB=.4                                                       SRC 0056
          IF(IV.GT.3571)THEN                                            SRC 0057
              V=10000./IV                                               SRC 0058
              I=10*V                                                    SRC 0059
              FALB=ALB(I)+(ALB(I+1)-ALB(I))*(10*V-I)                    SRC 0060
          ELSEIF(IV.GT.2000)THEN                                        SRC 0061
              V=10000./IV                                               SRC 0062
              FALB=ALB(28)+(ALB(29)-ALB(28))*(V-2.8)/2.2                SRC 0063
          ENDIF                                                         SRC 0064
      ENDIF                                                             SRC 0065
C                                                                       SRC 0066
C     CHECK FOR INITIAL CALL.                                           SRC 0067
      IF(ICALL.NE.1)THEN                                                SRC 0068
          ICALL=1                                                       SRC 0069
C                                                                       SRC 0070
C         SUN ELLIPTIC ORBIT FACTOR                                     SRC 0071
          IF(IDAY.LT.NDAY(1) .OR. IDAY.GT.NDAY(13))THEN                 SRC 0072
              FORBIT=1.                                                 SRC 0073
          ELSEIF(IDAY.EQ.1)THEN                                         SRC 0074
              FORBIT=RAT(1)                                             SRC 0075
          ELSE                                                          SRC 0076
              IM1=1                                                     SRC 0077
              DO 10 I=2,12                                              SRC 0078
                  IF(IDAY.LE.NDAY(I))GOTO20                             SRC 0079
   10         IM1=I                                                     SRC 0080
              I=13                                                      SRC 0081
   20         FORBIT=RAT(IM1)+(IDAY-NDAY(IM1))*                         SRC 0082
     1          (RAT(I)-RAT(IM1))/(NDAY(I)-NDAY(I-1))                   SRC 0083
          ENDIF                                                         SRC 0084
      ENDIF                                                             SRC 0085
C                                                                       SRC 0086
C     SOLAR INTENSITY [W CM-2 / CM-1]                                   SRC 0087
      SOURCE=FORBIT*SUN(IV)                                             SRC 0088
C                                                                       SRC 0089
C     LUNAR INTENSITY [W CM-2 / CM-1]                                   SRC 0090
      IF(ISOURC.EQ.1)SOURCE=2.04472E-7*FPHS*FALB*SOURCE                 SRC 0091
      RETURN                                                            SRC 0092
      END                                                               SRC 0093
