      SUBROUTINE MAPMS(ML,IKMAX)                                        MMS 0001
C                                                                       MMS 0002
C     ROUTINE MAPMS DETERMINES WHICH ATMOSPHERIC LAYER                  MMS 0003
C     CONTAINS EACH LINE-OF-SIGHT PATH SEGMENT.                         MMS 0004
C                                                                       MMS 0005
C     LIST PARAMETERS                                                   MMS 0006
      INCLUDE 'PARAM.LST'                                               MMS 0007
C                                                                       MMS 0008
C     DECLARE INPUTS                                                    MMS 0009
C       ML      NUMBER OF LAYER BOUNDARIES                              MMS 0010
C       IKMAX   NUMBER OF LINE-OF-SIGHT PATHS SEGMENTS                  MMS 0011
      INTEGER ML,IKMAX                                                  MMS 0012
      REAL ZM,PM,TM,RFNDX,DENSTY                                        MMS 0013
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    MMS 0014
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               MMS 0015
C                                                                       MMS 0016
C     /PATH/                                                            MMS 0017
C       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.               MMS 0018
C       AHT     ALTITUDES AT PATH BOUNDARIES.                           MMS 0019
C       IMAP    MAPPING FROM PATH SEGMENTS TO LAYERS.                   MMS 0020
      INTEGER IMAP                                                      MMS 0021
      REAL QTHETA,AHT,TPH                                               MMS 0022
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),TPH(LAYTWO),IMAP(LAYTWO)   MMS 0023
C                                                                       MMS 0024
C     DECLARE LOCAL VARIABLES                                           MMS 0025
      INTEGER I,J,JP1                                                   MMS 0026
      REAL AVGALT                                                       MMS 0027
C                                                                       MMS 0028
C     LOOP OVER PATH SEGMENTS                                           MMS 0029
      J=1                                                               MMS 0030
      DO 20 JP1=2,IKMAX+1                                               MMS 0031
C                                                                       MMS 0032
C         DEFINE LINE SEGMENT MIDDLE ALTITUDE                           MMS 0033
          AVGALT=.5*(AHT(J)+AHT(JP1))                                   MMS 0034
          IMAP(J)=1                                                     MMS 0035
          DO 10 I=2,ML                                                  MMS 0036
              IF(ZM(I).GE.AVGALT)GOTO20                                 MMS 0037
   10     IMAP(J)=I                                                     MMS 0038
          STOP 'ERROR in MAPMS:  Path segment above top of atmosphere.' MMS 0039
   20 J=JP1                                                             MMS 0040
      RETURN                                                            MMS 0041
      END                                                               MMS 0042
