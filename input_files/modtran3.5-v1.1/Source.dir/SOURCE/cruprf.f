      SUBROUTINE CRUPRF(NCRALT)                                         CRU 0001
C                                                                       CRU 0002
C     THIS ROUTINE READS IN USER-DEFINED CLOUD/RAIN MODEL PROFILES.     CRU 0003
C                                                                       CRU 0004
C     LIST PARAMETERS                                                   CRU 0005
      INCLUDE 'PARAM.LST'                                               CRU 0006
C                                                                       CRU 0007
C     DECLARE INPUTS                                                    CRU 0008
C       NCRALT    NUMBER OF CLOUD/RAIN PROFILES BOUNDARY ALTITUDES      CRU 0009
      INTEGER NCRALT                                                    CRU 0010
C                                                                       CRU 0011
C     LIST COMMONS                                                      CRU 0012
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               CRU 0013
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           CRU 0014
      REAL ZCLD,CLD,CLDICE,RR                                           CRU 0015
      COMMON/CLDRR/ZCLD(1:NZCLD,0:1),CLD(1:NZCLD,0:5),                  CRU 0016
     1  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)                             CRU 0017
C                                                                       CRU 0018
C     DECLARE LOCAL VARIABLES                                           CRU 0019
      INTEGER IALTM1,ICRALT                                             CRU 0020
C                                                                       CRU 0021
C     CHECK NUMBER OF BOUNDARY ALTITUDES (NCRALT)                       CRU 0022
      IF(NCRALT.GT.NZCLD)THEN                                           CRU 0023
          WRITE(IPR,'(/2A,I3,A,/8X,A,I3,A,/8X,2A)')' ERROR:  ',         CRU 0024
     1      'THE INPUT NUMBER OF CLOUD/RAIN BOUNDARY ALTITUDES (',      CRU 0025
     2      NCRALT,')',' EXCEEDS PARAMETER "NZCLD" (=',                 CRU 0026
     3      NZCLD,').  INCREASE "NZCLD"',' IN "PARAM.LST"',             CRU 0027
     4      ' AND MAKE NECESSARY CHANGES IN BLOCK DATA "MDTA".'         CRU 0028
          STOP 'TOO MANY CLOUD/RAIN BOUNDARY ALTITUDES'                 CRU 0029
      ENDIF                                                             CRU 0030
C                                                                       CRU 0031
C     WRITE OUT USER-DEFINED CLOUD/RAIN MODEL PROFILES MESSAGE          CRU 0032
      WRITE(IPR,'(/2A,I3,A)')' USER-DEFINED CLOUD/RAIN MODEL',          CRU 0033
     1  ' PROFILES WITH',NCRALT,' BOUNDARY ALTITUDES.'                  CRU 0034
C                                                                       CRU 0035
C     READ IN DATA AT FIRST ALTITUDE                                    CRU 0036
      ICRALT=1                                                          CRU 0037
      READ(IRD,'(4F10.5)',ERR=20)ZCLD(1,0),CLD(1,0),CLDICE(1,0),RR(1,0) CRU 0038
C                                                                       CRU 0039
C     THE BASE ALTITUDE (DEFINED RELATIVE TO THE GROUND), THE           CRU 0040
C     DENSITIES AND THE RAIN RATE MUST BE ALL BE NON-NEGATIVE           CRU 0041
      IF(ZCLD(1,0).LT.0.)THEN                                           CRU 0042
          WRITE(IPR,'(/2A,F10.5,A)')' ERROR:  THE CLOUD/RAIN',          CRU 0043
     1      ' PROFILE BASE ALTITUDE IS',ZCLD(1,0),' KM BELOW GROUND.'   CRU 0044
          STOP 'CLOUD/RAIN PROFILE BASE ALTITUDE IS BELOW THE GROUND'   CRU 0045
      ENDIF                                                             CRU 0046
      IF(CLD(1,0).LT.0. .OR. CLDICE(1,0).LT.0. .OR. RR(1,0).LT.0.)THEN  CRU 0047
          WRITE(IPR,'(/A,(A,F10.5,A))')' ERROR:  NEGATIVE',             CRU 0048
     1      ' CLOUD/RAIN DENSITY/RATE INPUT AT',ZCLD(1,0),' KM.',       CRU 0049
     2      '         WATER DROPLET DENSITY',CLD(1,0),' GM/M3',         CRU 0050
     3      '         ICE PARTICLE DENSITY ',CLDICE(1,0),' GM/M3',      CRU 0051
     4      '         RAIN RATE            ',RR(1,0),' MM/HR'           CRU 0052
          STOP 'NEGATIVE CLOUD/RAIN DENSITY/RATE INPUT'                 CRU 0053
      ENDIF                                                             CRU 0054
C                                                                       CRU 0055
C     LOOP OVER PROFILE ALTITUDE                                        CRU 0056
      IALTM1=1                                                          CRU 0057
      DO 10 ICRALT=2,NCRALT                                             CRU 0058
C                                                                       CRU 0059
C         READ IN DATA.                                                 CRU 0060
          READ(IRD,'(4F10.5)',ERR=20)                                   CRU 0061
     1      ZCLD(ICRALT,0),CLD(ICRALT,0),CLDICE(ICRALT,0),RR(ICRALT,0)  CRU 0062
C                                                                       CRU 0063
C         CHECK FOR INCREASING ALTITUDES.                               CRU 0064
          IF(ZCLD(ICRALT,0).LE.ZCLD(IALTM1,0))THEN                      CRU 0065
              WRITE(IPR,'(/3A,/8X,A,/(8X,I5,F10.5))')' ERROR: ',        CRU 0066
     1          ' THE CLOUD/RAIN PROFILE ALTITUDES MUST BE INPUT IN',   CRU 0067
     2          ' INCREASING ORDER.',' VALUES READ IN THUS FAR ARE:',   CRU 0068
     3          (IALTM1,ZCLD(IALTM1,0),IALTM1=1,ICRALT)                 CRU 0069
              STOP 'CLOUD/RAIN ALTITUDES NOT MONOTONICALLY INCREASING'  CRU 0070
          ENDIF                                                         CRU 0071
C                                                                       CRU 0072
C         CHECK FOR NEGATIVE CLOUD DENSITIES OR RAIN RATES.             CRU 0073
          IF(CLD(ICRALT,0).LT.0. .OR. CLDICE(ICRALT,0).LT.0.            CRU 0074
     1      .OR. RR(ICRALT,0).LT.0.)THEN                                CRU 0075
              WRITE(IPR,'(/2A,F10.5,A,/(8X,A,F10.5,A))')                CRU 0076
     1          ' ERROR:  NEGATIVE CLOUD/RAIN DENSITY',                 CRU 0077
     2          '/RATE INPUT AT',ZCLD(ICRALT,0),' KM.',                 CRU 0078
     3          ' WATER DROPLET DENSITY',CLD(ICRALT,0),' GM/M3',        CRU 0079
     4          ' ICE PARTICLE DENSITY ',CLDICE(ICRALT,0),' GM/M3',     CRU 0080
     5          ' RAIN RATE            ',RR(ICRALT,0),' MM/HR'          CRU 0081
              STOP 'NEGATIVE CLOUD/RAIN DENSITY/RATE INPUT'             CRU 0082
          ENDIF                                                         CRU 0083
   10 IALTM1=ICRALT                                                     CRU 0084
C                                                                       CRU 0085
C     RETURN TO CRPROF                                                  CRU 0086
      RETURN                                                            CRU 0087
C                                                                       CRU 0088
C     ERROR READING CLOUD/RAIN DATA                                     CRU 0089
   20 CONTINUE                                                          CRU 0090
      WRITE(IPR,'(/2A,I3)')' ERROR:  UNABLE TO READ',                   CRU 0091
     1  ' CLOUD/RAIN PROFILE DATA FOR LAYER',ICRALT                     CRU 0092
      STOP 'ERROR READING CLOUD/RAIN PROFILE DATA'                      CRU 0093
      END                                                               CRU 0094
