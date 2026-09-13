      SUBROUTINE CIRR18                                                 C18 0001
C*********************************************************************  C18 0002
C*  ROUTINE TO SET CTHIK CALT CEXT  FOR  CIRRUS CLOUDS 18 19        **  C18 0003
C*  INPUTS]                                                         **  C18 0004
C*           CTHIK    -  CIRRUS THICKNESS (KM)                      **  C18 0005
C*                       0 = USE THICKNESS STATISTICS               **  C18 0006
C*                       .NE. 0 = USER DEFINES THICKNESS            **  C18 0007
C*                                                                  **  C18 0008
C*           CALT     -  CIRRUS BASE ALTITUDE (KM)                  **  C18 0009
C*                       0 = USE CALCULATED VALUE                   **  C18 0010
C*                       .NE. 0 = USER DEFINES BASE ALTITUDE        **  C18 0011
C*                                                                  **  C18 0012
C*           ICLD     -  CIRRUS PRESENCE FLAG                       **  C18 0013
C*                       0 = NO CIRRUS                              **  C18 0014
C*                       18  19 = USE CIRRUS PROFILE                **  C18 0015
C*                                                                  **  C18 0016
C*           MODEL    -  ATMOSPHERIC MODEL                          **  C18 0017
C*                       1-5  AS IN MAIN PROGRAM                    **  C18 0018
C*                       MODEL = 0,6,7 NOT USED SET TO 2            **  C18 0019
C*                                                                  **  C18 0020
C*  OUTPUTS]                                                        **  C18 0021
C*         CTHIK        -  CIRRUS THICKNESS (KM)                    **  C18 0022
C*         CALT         -  CIRRUS BASE ALTITUDE (KM)                **  C18 0023
C          CEXT IS THE EXTINCTION COEFFIENT(KM-1) AT 0.55               C18 0024
C               DEFAULT VALUE 0.14*CTHIK                                C18 0025
C*                                                                  **  C18 0026
C*********************************************************************  C18 0027
C                                                                       C18 0028
      COMMON /CARD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB  C18 0029
     1  ,MODTRN                                                         C18 0030
      LOGICAL MODTRN                                                    C18 0031
      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   C18 0032
     1    RAINRT                                                        C18 0033
      INTEGER NCRALT,NCRSPC                                             C18 0034
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    C18 0035
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      C18 0036
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       C18 0037
      INCLUDE 'PARAM.LST'                                               C18 0038
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           C18 0039
      REAL CAMEAN(5)                                                    C18 0040
      DATA CAMEAN/11.,10.,8.,7.,5./                                     C18 0041
      MDL = MODEL                                                       C18 0042
C                                                                       C18 0043
C  CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES             C18 0044
C  DEFAULTED MEAN CIRRUS THICKNESS IS 1.0KM  OR 0.2 KM.                 C18 0045
C                                                                       C18 0046
      IF ( CTHIK .GT. 0.0 ) GO TO 25                                    C18 0047
      IF(ICLD.EQ.18) CTHIK=1.0                                          C18 0048
      IF(ICLD.EQ.19) CTHIK=0.2                                          C18 0049
25    IF(CEXT .EQ. 0.) CEXT = 0.14 * CTHIK                              C18 0050
C                                                                       C18 0051
C  BASE HEIGHT CALCULATIONS                                             C18 0052
C                                                                       C18 0053
      IF ( MODEL .LT. 1  .OR.  MODEL .GT. 5 ) MDL = 2                   C18 0054
C                                                                       C18 0055
      IF(CALT.LE.0.)CALT=CAMEAN(MDL)                                    C18 0056
C                                                                       C18 0057
      IF(ICLD. EQ. 18) WRITE(IPR,1219)                                  C18 0058
1219  FORMAT(15X,'CIRRUS ATTENUATION INCLUDED   (STANDARD CIRRUS)')     C18 0059
      IF(ICLD. EQ. 19) WRITE(IPR,1220)                                  C18 0060
1220  FORMAT(15X,'CIRRUS ATTENUATION INCLUDED   (THIN     CIRRUS)')     C18 0061
      WRITE(IPR,1221) CTHIK                                             C18 0062
1221  FORMAT(15X,'CIRRUS THICKNESS ',                                   C18 0063
     X F10.3,'KM')                                                      C18 0064
      WRITE(IPR,1224)CALT                                               C18 0065
1224  FORMAT(15X,'CIRRUS BASE ALTITUDE ',                               C18 0066
     X F10.3,' KM')                                                     C18 0067
       WRITE(IPR,1226) CEXT                                             C18 0068
1226    FORMAT(15X,'CIRRUS PROFILE EXTINCT ',F10.3)                     C18 0069
C                                                                       C18 0070
C       END OF CIRRUS MODEL SET UP                                      C18 0071
      RETURN                                                            C18 0072
      END                                                               C18 0073
