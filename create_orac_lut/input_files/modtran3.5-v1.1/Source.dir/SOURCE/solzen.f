      SUBROUTINE SOLZEN(L,IKMAX1,REE,DEG,ANGLE)                         SZN 0001
C                                                                       SZN 0002
C     ROUTINE SOLZEN STORES THE COSINE OF THE SOLAR/LUNAR ZENITH        SZN 0003
C     ANGLE DATA FOR THE MULTIPLE SCATTERING CALCULATIONS:              SZN 0004
C       CSZEN0  COSINE OF SOLAR/LUNAR ZENITH AT THE GROUND.             SZN 0005
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH FOR VERTICAL SZN 0006
C               PATH (USED IN CALCULATION OF BACKSCATTER FRACTION).     SZN 0007
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               SZN 0008
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    SZN 0009
C                                                                       SZN 0010
C     DEFINE INPUTS                                                     SZN 0011
C       L       CURRENT LAYER BOUNDARY.                                 SZN 0012
C       IKMAX1  NUMBER OF LAYER BOUNDARIES (NUMBER OF LAYERS PLUS ONE). SZN 0013
C       REE     EARTH RADIUS [KM]                                       SZN 0014
C       DEG     NUMBER OF DEGREES IN ONE RADIAN                         SZN 0015
C       ANGLE   SOLAR ZENITH ANGLE AT BOUNDARY L [DEG]                  SZN 0016
      INTEGER L,IKMAX1                                                  SZN 0017
      REAL REE,DEG,ANGLE                                                SZN 0018
C                                                                       SZN 0019
C     LIST PARAMETERS                                                   SZN 0020
      INCLUDE 'PARAM.LST'                                               SZN 0021
C                                                                       SZN 0022
C     LIST COMMONS                                                      SZN 0023
      REAL ZM,PM,TM,RFNDX,DENSTY                                        SZN 0024
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    SZN 0025
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               SZN 0026
C                                                                       SZN 0027
C     COMMON /MSRD/                                                     SZN 0028
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            SZN 0029
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             SZN 0030
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               SZN 0031
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    SZN 0032
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  SZN 0033
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. SZN 0034
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               SZN 0035
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         SZN 0036
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          SZN 0037
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              SZN 0038
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  SZN 0039
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               SZN 0040
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           SZN 0041
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     SZN 0042
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   SZN 0043
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       SZN 0044
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     SZN 0045
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              SZN 0046
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 SZN 0047
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          SZN 0048
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             SZN 0049
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                SZN 0050
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      SZN 0051
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         SZN 0052
C                                                                       SZN 0053
C     DECLARE LOCAL VARIABLES                                           SZN 0054
      INTEGER LP1                                                       SZN 0055
      REAL RIDIF,DSTDIF,DSTRAT,PROD                                     SZN 0056
C                                                                       SZN 0057
C     COSINE OF THE SOLAR ZENITH AT UPPER LAYER BOUNDARY                SZN 0058
      CSZEN0(L)=COS(ANGLE/DEG)                                          SZN 0059
C                                                                       SZN 0060
C     COSINE OF THE SOLAR ZENITH AVERAGED FOR LAYER.                    SZN 0061
      IF(L.GT.1)CSZEN(L-1)=.5*(CSZEN0(L-1)+CSZEN0(L))                   SZN 0062
C                                                                       SZN 0063
C     COSINE OF THE SOLAR ZENITH EXITING THE LAYER ABOVE.               SZN 0064
      IF(L.GE.IKMAX1)RETURN                                             SZN 0065
      LP1=L+1                                                           SZN 0066
      RIDIF=(RFNDX(L)-RFNDX(LP1))/(1.+RFNDX(LP1))                       SZN 0067
      DSTDIF=(ZM(LP1)-ZM(L))/(REE+ZM(LP1))                              SZN 0068
      DSTRAT=1.-DSTDIF                                                  SZN 0069
      PROD=(1.+RIDIF)*DSTRAT                                            SZN 0070
      CSZENX(L)=.5*(ABS(CSZEN0(L))+                                     SZN 0071
     1  SQRT((CSZEN0(L)*PROD)**2+(1.+PROD)*(DSTDIF-RIDIF*DSTRAT)))      SZN 0072
      RETURN                                                            SZN 0073
      END                                                               SZN 0074
