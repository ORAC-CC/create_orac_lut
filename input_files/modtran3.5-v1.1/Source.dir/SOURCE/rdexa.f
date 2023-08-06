      SUBROUTINE RDEXA                                                  RXA 0001
C                                                                       RXA 0002
C     READ IN USER DEFINED EXTINCTION ABSORPTION COEFFICIENTS AND       RXA 0003
C     ASYMMETRY PARAMETERS                                              RXA 0004
C                                                                       RXA 0005
      INCLUDE 'PARAM.LST'                                               RXA 0006
      INTEGER KPOINT                                                    RXA 0007
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     RXA 0008
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   RXA 0009
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   RXA 0010
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     RXA 0011
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          RXA 0012
      COMMON /CARD2D/ IREG(4),ALTB(4),IREGC(4)                          RXA 0013
      CHARACTER*4 TITLE(18)                                             RXA 0014
      REAL VX(NWAVLN)                                                   RXA 0015
C                                                                       RXA 0016
      READ (IRD,1200) (IREG(IK),IK=1,4)                                 RXA 0017
1200  FORMAT(4I5)                                                       RXA 0018
      WRITE(IPR,1210) (IREG(IK),IK=1,4)                                 RXA 0019
1210  FORMAT('0 CARD 2D *****',4I5)                                     RXA 0020
C                                                                       RXA 0021
      DO 1300 IHC = 1,4                                                 RXA 0022
C                                                                       RXA 0023
      IF(IREG(IHC) .EQ. 0) GO TO 1300                                   RXA 0024
      READ(IRD,1220) AWCCON(IHC),TITLE                                  RXA 0025
1220  FORMAT(E10.3,18A4)                                                RXA 0026
      WRITE(IPR,1230) AWCCON(IHC),TITLE                                 RXA 0027
1230  FORMAT('0 CARD 2D1 **** EQUIVALENT WATER = ',1PE10.3,18A4)        RXA 0028
      WRITE(IPR,1250)                                                   RXA 0029
1250  FORMAT('0 CARD 2D2 ****')                                         RXA 0030
C                                                                       RXA 0031
      READ(IRD,'((3(F6.2,2F7.5,F6.4)))')                                RXA 0032
     1  (VX(I),EXTC(IHC,I),ABSC(IHC,I),ASYM(IHC,I),I=1,NWAVLN)          RXA 0033
      WRITE(IPR,'((2X,3(F6.2,2F7.5,F6.4)))')                            RXA 0034
     1  (VX(I),EXTC(IHC,I),ABSC(IHC,I),ASYM(IHC,I),I=1,NWAVLN)          RXA 0035
1300  CONTINUE                                                          RXA 0036
      RETURN                                                            RXA 0037
      END                                                               RXA 0038
