      SUBROUTINE NEWH2(H1,H2,ANGLE,RANGE,BETA,LENN,HTAN,PHI)            NH2 0001
C                                                                       NH2 0002
C     THIS ROUTINE DETERMINES                                           NH2 0003
C       HTAN     TANGENT HEIGHT [KM],                                   NH2 0004
C       H2       FINAL ALTITUDE [KM],                                   NH2 0005
C       BETA     EARTH CENTER ANGLE [DEG], AND                          NH2 0006
C       LENN     PATH SWITCH (=1 IF PATH PASSES THROUGH HTAN)           NH2 0007
C     GIVEN INPUTS                                                      NH2 0008
C       H1       INITIAL (OBSERVER) ALTITUDE [KM],                      NH2 0009
C       ANGLE    PATH ZENITH ANGLE AT H1 [DEG], AND                     NH2 0010
C       RANGE    PATH SLANT RANGE.                                      NH2 0011
C                                                                       NH2 0012
C     DECLARE ARGUMENTS:                                                NH2 0013
      DOUBLE PRECISION H1,H2,ANGLE,RANGE,BETA,HTAN,PHI                  NH2 0014
      INTEGER LENN                                                      NH2 0015
C                                                                       NH2 0016
C     LIST COMMONS:                                                     NH2 0017
      REAL RE,ZMAX                                                      NH2 0018
      INTEGER IMAX,IMOD,IPATH                                           NH2 0019
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             NH2 0020
C                                                                       NH2 0021
C     DECLARE LOCAL VARIABLES:                                          NH2 0022
      DOUBLE PRECISION RX,CPATH,RATIO                                   NH2 0023
C                                                                       NH2 0024
C     DECLARE FUNCTIONS:                                                NH2 0025
      DOUBLE PRECISION TANHT                                            NH2 0026
C                                                                       NH2 0027
C     LIST DATA:                                                        NH2 0028
      DOUBLE PRECISION DEG                                              NH2 0029
      DATA DEG/57.2957795131D0/                                         NH2 0030
C                                                                       NH2 0031
C     COMPUTE CPATH, THE PATH CONSTANT,                                 NH2 0032
      CALL IRFXN(H1,RX,RATIO)                                           NH2 0033
      CPATH=RX*SIN(ANGLE/DEG)*(H1+DBLE(RE))                             NH2 0034
C                                                                       NH2 0035
C     FIND HTAN, H2, BETA AND LENN.                                     NH2 0036
      HTAN=TANHT(CPATH,H1)                                              NH2 0037
      CALL FNDPTH(CPATH,H1,HTAN,H2,RANGE,BETA,LENN,ANGLE,PHI)           NH2 0038
      IF(ANGLE.LE.90.)HTAN=MIN(H1,H2)                                   NH2 0039
C                                                                       NH2 0040
C     ROUTINE OPERATIONS COMPLETED.                                     NH2 0041
      RETURN                                                            NH2 0042
      END                                                               NH2 0043
