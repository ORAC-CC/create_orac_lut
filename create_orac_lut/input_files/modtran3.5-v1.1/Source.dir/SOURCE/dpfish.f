      SUBROUTINE DPFISH(H,SH,GAMMA)                                     FSH 0001
C                                                                       FSH 0002
C     DOUBLE PRECISION VERSION OF THE PREVIOUS ROUTINE FINDSH.          FSH 0003
C                                                                       FSH 0004
C*****GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARIES   FSH 0005
C*****ZM(I1) AND ZM(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE     FSH 0006
C*****HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE         FSH 0007
C*****INDEX OF REFRACTION                                               FSH 0008
C                                                                       FSH 0009
      INCLUDE 'PARAM.LST'                                               FSH 0010
      DOUBLE PRECISION H, SH, GAMMA, Z1,Z2,RFNDX1,RFNDX2                FSH 0011
      REAL RE,ZMAX                                                      FSH 0012
      INTEGER IMAX,IMOD,IPATH                                           FSH 0013
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             FSH 0014
      REAL ZM,PM,TM,RFNDX,DENSTY                                        FSH 0015
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    FSH 0016
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               FSH 0017
      DO 100 IM=2,IMOD                                                  FSH 0018
      I2 = IM                                                           FSH 0019
      IF(ZM(IM).GE.REAL(H))GOTO110                                      FSH 0020
  100 CONTINUE                                                          FSH 0021
      I2 = IMOD                                                         FSH 0022
  110 CONTINUE                                                          FSH 0023
      I1 = I2-1                                                         FSH 0024
      Z1=DBLE(ZM(I1))                                                   FSH 0025
      Z2=DBLE(ZM(I2))                                                   FSH 0026
      RFNDX1=DBLE(RFNDX(I1))                                            FSH 0027
      RFNDX2=DBLE(RFNDX(I2))                                            FSH 0028
      CALL DPSCHT(Z1,Z2,RFNDX1,RFNDX2,SH,GAMMA)                         FSH 0029
      RETURN                                                            FSH 0030
      END                                                               FSH 0031
