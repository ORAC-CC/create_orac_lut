      SUBROUTINE DPSCHT(Z1,Z2,RFNDX1,RFNDX2,SH,GAMMA)                   SHT 0001
C***********************************************************************SHT 0002
C     THE DOUBLE PRECISION VERSION OF THE PREVIOUS ROUTINE SCALHT.      SHT 0003
C                                                                       SHT 0004
C     THIS SUBROUTINE CALCULATES THE SCALE HEIGHT SH OF THE (INDEX OF   SHT 0005
C     REFRACTION-1.0) FROM THE VALUES OF THE INDEX AT THE ALTITUDES Z1  SHT 0006
C     AND Z2 ( Z1 < Z2). IT ALSO CALCULATES THE EXTRAPOLATED VALUE      SHT 0007
C     GAMMA OF THE (INDEX-1.0) AT Z = 0.0                               SHT 0008
C***********************************************************************SHT 0009
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SHT 0010
      RF1 = RFNDX1+1.0E-20                                              SHT 0011
      RF2 = RFNDX2+1.0E-20                                              SHT 0012
      RATIO = RF1/RF2                                                   SHT 0013
      IF(ABS(RATIO-1.0).LT.1.0E-05)  GO TO 100                          SHT 0014
C     SH = (Z2-Z1)/ALOG(RATIO)                                          SHT 0015
      SH = (Z2-Z1)/LOG(RATIO)                                           SHT 0016
      GAMMA = RF1*(RF2/RF1)**(-Z1/(Z2-Z1))                              SHT 0017
      GO TO 110                                                         SHT 0018
 100  CONTINUE                                                          SHT 0019
C*****THEVARIATION IN THE INDEX OF REFRACTION WITH HEIGHT IS            SHT 0020
C*****INSIGNIFICANTOR ZERO                                              SHT 0021
      SH = 0.0                                                          SHT 0022
      GAMMA = RFNDX1                                                    SHT 0023
 110  CONTINUE                                                          SHT 0024
      RETURN                                                            SHT 0025
      END                                                               SHT 0026
