      SUBROUTINE SHADE(IPH,IK,MSOFF,IPATH,KNTRVL,V,TSNOBS,TSNREF,SUMSSS)SHD 0001
C                                                                       SHD 0002
C     ROUTINE SHADE SETS THE SOLAR TRANSMISSION TERMS TO ZERO FOR       SHD 0003
C     SHADED SCATTERING POINTS                                          SHD 0004
C                                                                       SHD 0005
C     DECLARE ARGUMENTS                                                 SHD 0006
      INTEGER IPH,IK,MSOFF,IPATH,KNTRVL                                 SHD 0007
      REAL V,TSNOBS,TSNREF,SUMSSS                                       SHD 0008
      INCLUDE 'PARAM.LST'                                               SHD 0009
C                                                                       SHD 0010
C     COMMON /CORKDT/                                                   SHD 0011
C       WTKSUB   SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.           SHD 0012
C       DEPLAY   INCREMENTAL EXTINCTION OPTICAL DEPTHS                  SHD 0013
C       TRNLAY   INCREMENTAL TRANSMITTANCES                             SHD 0014
C       TRNCUM   CUMULATIVE TRANSMITTANCES                              SHD 0015
      REAL WTKSUB,DEPLAY,TRNLAY,TRNCUM                                  SHD 0016
      COMMON/CORKDT/WTKSUB(NKSUB),DEPLAY(NKSUB),                        SHD 0017
     1  TRNLAY(NKSUB),TRNCUM(NKSUB)                                     SHD 0018
      SAVE /CORKDT/                                                     SHD 0019
C                                                                       SHD 0020
C       PI       THE CONSTANT PI                                        SHD 0021
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       SHD 0022
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       SHD 0023
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         SHD 0024
      REAL PI,DEG,BIGNUM,BIGEXP                                         SHD 0025
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                SHD 0026
C                                                                       SHD 0027
C     COMMON /MSRD/                                                     SHD 0028
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            SHD 0029
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             SHD 0030
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               SHD 0031
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    SHD 0032
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  SHD 0033
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. SHD 0034
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               SHD 0035
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         SHD 0036
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          SHD 0037
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              SHD 0038
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  SHD 0039
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               SHD 0040
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           SHD 0041
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     SHD 0042
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   SHD 0043
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       SHD 0044
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     SHD 0045
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              SHD 0046
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 SHD 0047
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          SHD 0048
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             SHD 0049
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                SHD 0050
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      SHD 0051
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         SHD 0052
C                                                                       SHD 0053
C     DECLARE LOCAL VARIABLES                                           SHD 0054
      INTEGER INTRVL                                                    SHD 0055
      REAL S0,TMOL,BIGDEP                                               SHD 0056
      S0=0.                                                             SHD 0057
      IF(KNTRVL.EQ.1)THEN                                               SHD 0058
          TMOL=0.                                                       SHD 0059
          CALL SSRAD(IPH,IK,MSOFF,IPATH,V,S0,TMOL,TSNOBS,TSNREF,SUMSSS) SHD 0060
      ELSE                                                              SHD 0061
C                                                                       SHD 0062
C         DIVIDE BIGNUM BY THE LAYER NUMBER TO ARTIFICIALLY REQUIRE     SHD 0063
C         THAT OPTICAL DEPTH DECREASES WITH INCREASING LAYER NUMBER.    SHD 0064
          BIGDEP=BIGNUM/IK                                              SHD 0065
          DO 10 INTRVL=1,KNTRVL                                         SHD 0066
              TRNLAY(INTRVL)=0.                                         SHD 0067
              DEPLAY(INTRVL)=BIGDEP                                     SHD 0068
   10     CONTINUE                                                      SHD 0069
          IF(MSOFF.GT.0)THEN                                            SHD 0070
C                                                                       SHD 0071
C             FOR MULTIPLE SCATTERING, STORE THE DECREASE IN WEAK-LINE  SHD 0072
C             OPTICAL DEPTH TO THE SUN ACROSS THE CURRENT LAYER.        SHD 0073
              IF(IPATH.EQ.1)THEN                                        SHD 0074
                  DEPRAT(IK)=BIGDEP                                     SHD 0075
              ELSE                                                      SHD 0076
                  DEPRAT(IK+1)=BIGDEP                                   SHD 0077
                  DEPRAT(IK)=1.-BIGDEP/DEPRAT(IK)                       SHD 0078
              ENDIF                                                     SHD 0079
          ENDIF                                                         SHD 0080
CORK      CALL SSCORK(IPH,IK,MSOFF,IPATH,KNTRVL,V,S0,TSNOBS,SUMSSS)     SHD 0081
      ENDIF                                                             SHD 0082
      RETURN                                                            SHD 0083
      END                                                               SHD 0084
