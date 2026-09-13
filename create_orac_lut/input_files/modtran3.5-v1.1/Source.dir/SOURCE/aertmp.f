      SUBROUTINE AERTMP                                                 AET 0001
C                                                                       AET 0002
C     THIS ROUTINE DETERMINES AEROSOL DENSITY                           AET 0003
C     WEIGHTED PATH AVERAGED TEMPERATURES                               AET 0004
      INCLUDE 'PARAM.LST'                                               AET 0005
      INTEGER KPOINT                                                    AET 0006
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     AET 0007
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   AET 0008
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   AET 0009
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     AET 0010
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     AET 0011
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                AET 0012
      REAL TAER                                                         AET 0013
      COMMON/AERTM/TAER(NAER)                                           AET 0014
C                                                                       AET 0015
C     DECLARE LOCAL VARIABLES                                           AET 0016
      INTEGER IAER,M,IK                                                 AET 0017
C                                                                       AET 0018
C     LIST DATA                                                         AET 0019
      INTEGER MAP(NAER)                                                 AET 0020
      DATA MAP/7,12,13,14,16,66,67/                                     AET 0021
      DO 20 IAER=1,NAER                                                 AET 0022
          M=MAP(IAER)                                                   AET 0023
          TAER(IAER)=TBBY(1)                                            AET 0024
          IF(W(M).LE.0)GOTO20                                           AET 0025
          TAER(IAER)=TBBY(1)*WPATH(1,M)                                 AET 0026
          DO 10 IK=2,IKMAX                                              AET 0027
              TAER(IAER)=TAER(IAER)+TBBY(IK)*WPATH(IK,M)                AET 0028
   10     CONTINUE                                                      AET 0029
          TAER(IAER)=TAER(IAER)/W(M)                                    AET 0030
   20 CONTINUE                                                          AET 0031
      RETURN                                                            AET 0032
      END                                                               AET 0033
