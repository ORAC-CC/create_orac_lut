      DOUBLE PRECISION FUNCTION TANHT(CPATH,H1)                         THT 0001
C                                                                       THT 0002
C     THIS ROUTINE FINDS THE TANGENT HEIGHT GIVEN                       THT 0003
C     CPATH (THE PATH CONSTANT). SEE ATMOSPHERIC                        THT 0004
C     TRANSMITTANCE/RADIANCE:  LOWTRAN 6, AFGL-TR-83-018                THT 0005
C                                                                       THT 0006
C     INCLUDE PARAMETERS                                                THT 0007
      INCLUDE 'PARAM.LST'                                               THT 0008
C                                                                       THT 0009
C     LIST COMMONS                                                      THT 0010
      REAL RE,ZMAX                                                      THT 0011
      INTEGER IMAX,IMOD,IPATH                                           THT 0012
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             THT 0013
      REAL ZM,PM,TM,RFNDX,DENSTY                                        THT 0014
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    THT 0015
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               THT 0016
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     THT 0017
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                THT 0018
      REAL GNDALT                                                       THT 0019
      COMMON/GRAUND/GNDALT                                              THT 0020
C                                                                       THT 0021
C     DECLARE LOCAL VARIABLES                                           THT 0022
      INTEGER J,JMAX                                                    THT 0023
      DOUBLE PRECISION CPATH,H1,X,RX,RATIO,X1,X2                        THT 0024
C                                                                       THT 0025
C     DECLARE LOCAL ARRAYS                                              THT 0026
      DOUBLE PRECISION CP(LAYDIM),R(LAYDIM)                             THT 0027
C                                                                       THT 0028
C     DECLARE FUNCTIONS                                                 THT 0029
      DOUBLE PRECISION RTBIS                                            THT 0030
      DO 10 J=1,ML                                                      THT 0031
          R(J)=DBLE(RE+ZM(J))                                           THT 0032
          X=DBLE(ZM(J))                                                 THT 0033
          CALL IRFXN(X,RX,RATIO)                                        THT 0034
          CP(J)=R(J)*RX                                                 THT 0035
          JMAX=J+1                                                      THT 0036
          IF(H1.LE.DBLE(ZM(J)))GOTO20                                   THT 0037
   10 CONTINUE                                                          THT 0038
   20 CONTINUE                                                          THT 0039
      DO 30 J=JMAX,2,-1                                                 THT 0040
          IF(CPATH.LE.CP(J) .AND. CPATH.GE.CP(J-1))THEN                 THT 0041
              X1=DBLE(ZM(J-1))                                          THT 0042
              X2=DBLE(ZM(J))                                            THT 0043
              TANHT=RTBIS(X1,X2,CPATH)                                  THT 0044
              RETURN                                                    THT 0045
          ENDIF                                                         THT 0046
   30 CONTINUE                                                          THT 0047
      TANHT=DBLE(GNDALT)                                                THT 0048
      RETURN                                                            THT 0049
      END                                                               THT 0050
