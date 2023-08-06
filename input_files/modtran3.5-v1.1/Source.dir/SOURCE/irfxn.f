      SUBROUTINE IRFXN(X,RX,RATIO)                                      IRF 0001
      INCLUDE 'PARAM.LST'                                               IRF 0002
C                                                                       IRF 0003
C     THIS ROUTINE FINDS INDEX OF REFRACTION AND ITS DERIVATIVE AT X.   IRF 0004
C                                                                       IRF 0005
C     X = HEIGHT FROM THE SURFACE OF THE EARTH                          IRF 0006
C     RX = INDEX OF REFRACTION                                          IRF 0007
C     DRX = DERIVATIVE                                                  IRF 0008
C     RATIO = RX/DRX                                                    IRF 0009
C                                                                       IRF 0010
      DOUBLE PRECISION  X, RX, RATIO, XN, H                             IRF 0011
      INTEGER I, J                                                      IRF 0012
C                                                                       IRF 0013
      REAL ZM,PM,TM,RFNDX,DENSTY                                        IRF 0014
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    IRF 0015
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               IRF 0016
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     IRF 0017
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                IRF 0018
C                                                                       IRF 0019
      IF (X .GT. ZM(ML)) THEN                                           IRF 0020
         J = ML                                                         IRF 0021
         GO TO 200                                                      IRF 0022
      ENDIF                                                             IRF 0023
      DO 300 I = 2, ML                                                  IRF 0024
         IF (X .LE. ZM(I)) THEN                                         IRF 0025
            J = I                                                       IRF 0026
            GO TO 200                                                   IRF 0027
         ENDIF                                                          IRF 0028
 300  CONTINUE                                                          IRF 0029
C                                                                       IRF 0030
 200  I = J -1                                                          IRF 0031
C                                                                       IRF 0032
      H = -(ZM(J)-ZM(I))/LOG(RFNDX(J)/RFNDX(I))                         IRF 0033
      XN = RFNDX(I)*EXP(-(X-ZM(I))/H)                                   IRF 0034
      RX = 1 + XN                                                       IRF 0035
      RATIO = -(RX*H)/XN                                                IRF 0036
      RETURN                                                            IRF 0037
      END                                                               IRF 0038
