      SUBROUTINE H2SRC(BETAH2,IPH,IPARM,PARM1,PARM2,                    H2S 0001
     1  PARM3,PARM4,PSIPO,G,MSOFF,ICH1,KNTRVL)                          H2S 0002
C                                                                       H2S 0003
C     THIS ROUTINE CALLS SSGEO FOR A VERTICAL PATH AT H2.               H2S 0004
C                                                                       H2S 0005
C     DECLARE INPUTS/OUTPUTS                                            H2S 0006
C       BETAH2   EARTH CENTER ANGLE BETWEEN H1 AND H2                   H2S 0007
C       IPH      PHASE FUNCTION FLAG                                    H2S 0008
C       IPARM    SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG                H2S 0009
C       PARM1    OBSERVER LATITUDE [DEG NORTH OF EQUATOR] (IPARM<2)     H2S 0010
C                RELATIVE SOLAR AZIMUTH [DEG EAST OF NORTH] (IPARM=2)   H2S 0011
C       PARM2    OBSERVER LONGITUDE [DEG WEST OF GREENWICH] (IPARM<2)   H2S 0012
C                SOLAR ZENITH ANGLE [DEG] (IPARM=2)                     H2S 0013
C       PARM3    SOLAR LATITUDE [DEG NORTH OF EQUATOR] (IPARM<2)        H2S 0014
C       PARM4    SOLAR LONGITUDE [DEG WEST OF GREENWICH] (IPARM<2)      H2S 0015
C       PSIPO    PATH TRUE AZIMUTH ANGLE [DEG EAST OF NORTH]            H2S 0016
C       G        HENYEY-GREENSTEIN ASYMMETRY FACTOR (IPH=0)             H2S 0017
C       MSOFF    OFFSET FOR MULTIPLE SCATTERING PATH ARRAYS             H2S 0018
C       ICH1     HAZE MODEL FLAG                                        H2S 0019
C       KNTRVL   NUMBER OF WEIGHTS IN CORRELATED-K APPROACH             H2S 0020
      REAL BETAH2,PARM1,PARM2,PARM3,PARM4,PSIPO,G                       H2S 0021
      INTEGER IPARM,IPH,MSOFF,ICH1,KNTRVL                               H2S 0022
C                                                                       H2S 0023
C     INCLUDE COMMONS                                                   H2S 0024
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               H2S 0025
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           H2S 0026
C                                                                       H2S 0027
C     DECLARE LOCAL VARIABLES                                           H2S 0028
      REAL PRM1SV,PRM2SV,PSISAV,OLAT,OLONG,RLAT,RLONG,PSISO             H2S 0029
      INTEGER IERROR                                                    H2S 0030
C                                                                       H2S 0031
C     SAVE SOLAR PATH GEOMETRY VARIABLES                                H2S 0032
      PRM1SV=PARM1                                                      H2S 0033
      PRM2SV=PARM2                                                      H2S 0034
      PSISAV=PSIPO                                                      H2S 0035
C                                                                       H2S 0036
C     WRITE OUT HEADER                                                  H2S 0037
      WRITE(IPR,'(/2A)')' THE MULTIPLE SCATTERING FLUX',                H2S 0038
     1  ' CALCULATIONS ARE BASED ON A VERTICAL PATH AT H2:'             H2S 0039
C                                                                       H2S 0040
C     DETERMINE SOLAR ANGLES AT H2                                      H2S 0041
      IF(IPARM.LT.2)THEN                                                H2S 0042
C                                                                       H2S 0043
C         DETERMINE LATITUDE, LONGITUDE AND PATH AZIMUTH AT H2.         H2S 0044
          CALL LOCATE(PRM1SV,PRM2SV,PSISAV,BETAH2,PARM1,PARM2)          H2S 0045
          CALL PSIECA(PARM1,PARM2,PRM1SV,PRM2SV,PSIPO,BETAH2)           H2S 0046
          PSIPO=PSIPO+180.                                              H2S 0047
          IF(PSIPO.GE.360.)PSIPO=PSIPO-360.                             H2S 0048
          WRITE(IPR,'((5X,A,F12.5))')                                   H2S 0049
     1      ' OBSERVER LATITUDE AT H2 [DEG NORTH OF EQUATOR]  ',PARM1,  H2S 0050
     2      ' OBSERVER LONGITUDE AT H2 [DEG WEST OF GREENWICH]',PARM2,  H2S 0051
     3      ' OBSERVER PATH AZIMUTH AT H2 [DEG EAST OF NORTH] ',PSIPO   H2S 0052
      ELSE                                                              H2S 0053
C                                                                       H2S 0054
C         RELATIVE SOLAR AZIMUTH AND SOLAR ZENITH ANGLES WERE INPUTS.   H2S 0055
C         ASSUME THAT THE OBSERVER AND SUN WERE ORIGINALLY BOTH ON THE  H2S 0056
C         EQUATOR, THAT THE OBSERVER WAS AT 0 DEG LONGITUDE, AND THAT   H2S 0057
C         THE SUN WAS WEST OF THE OBSERVER (I.E. AT PARM2 DEG LONGITUDE)H2S 0058
C         STEP 1:  DETERMINE THE ORIGINAL PATH AZIMUTH, (PSIPO = THE    H2S 0059
C                  TRUE SOLAR AZIMUTH - THE RELATIVE SOLAR AZIMUTH).    H2S 0060
          PSIPO=270.-PRM1SV                                             H2S 0061
C                                                                       H2S 0062
C         STEP 2:  DETERMINE THE NEW OBSERVER LATITUDE AND LONGITUDE.   H2S 0063
          OLAT=0.                                                       H2S 0064
          OLONG=0.                                                      H2S 0065
          CALL LOCATE(OLAT,OLONG,PSIPO,BETAH2,RLAT,RLONG)               H2S 0066
C                                                                       H2S 0067
C         STEP 3:  DETERMINE THE NEW PATH AZIMUTH.                      H2S 0068
          CALL PSIECA(RLAT,RLONG,OLAT,OLONG,PSIPO,BETAH2)               H2S 0069
          PSIPO=PSIPO+180.                                              H2S 0070
C                                                                       H2S 0071
C         STEP 4:  DETERMINE THE NEW TRUE SOLAR AZIMUTH AND SOLAR ZENITHH2S 0072
          CALL PSIECA(RLAT,RLONG,OLAT,PRM2SV,PSISO,PARM2)               H2S 0073
C                                                                       H2S 0074
C         STEP 5:  DETERMINE THE NEW RELATIVE SOLAR AZIMUTH             H2S 0075
          PARM1=PSISO-PSIPO                                             H2S 0076
          IF(PARM1.GT.180.)PARM1=PARM1-360.                             H2S 0077
          IF(PARM1.LE.-180.)PARM1=PARM1+360.                            H2S 0078
          WRITE(IPR,'((5X,2A,F12.5))')' RELATIVE SOLAR/LUNAR',          H2S 0079
     1      ' AZIMUTH AT H2 [DEG EAST OF NORTH]',PARM1,' SOLAR/LUNAR',  H2S 0080
     2      ' ZENITH ANGLE AT H2 [DEGREES]              ',PARM2         H2S 0081
      ENDIF                                                             H2S 0082
C                                                                       H2S 0083
C     CALL SSGEO                                                        H2S 0084
      IERROR=0                                                          H2S 0085
      CALL SSGEO(IERROR,IPH,IPARM,PARM1,PARM2,                          H2S 0086
     1  PARM3,PARM4,PSIPO,G,MSOFF,ICH1,KNTRVL)                          H2S 0087
C                                                                       H2S 0088
C     RETURN SAVED GEOMETRY VARIABLES                                   H2S 0089
      PARM1=PRM1SV                                                      H2S 0090
      PARM2=PRM2SV                                                      H2S 0091
      PSIPO=PSISAV                                                      H2S 0092
C                                                                       H2S 0093
C     RETURN TO DRIVER                                                  H2S 0094
      RETURN                                                            H2S 0095
      END                                                               H2S 0096
