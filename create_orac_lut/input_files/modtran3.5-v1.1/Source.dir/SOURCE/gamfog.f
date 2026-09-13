      FUNCTION   GAMFOG(MR,FREQ,T,RHO)                                  GAM 0001
C                                                                       GAM 0002
C        COMPUTES ATTENUATION OF EQUIVALENT LIQUID WATER CONTENT        GAM 0003
C       IN CLOUDS OR FOG IN DB/KM                                       GAM 0004
C       CONVERTED TO NEPERS BY NEW CONSTANT 1.885                       GAM 0005
C                                                                       GAM 0006
C        FREQ = WAVENUMBER (INVERSE CM)                                 GAM 0007
C        T    = TEMPERATURE (DEGREES KELVIN)                            GAM 0008
C        RHO  = EQUIVALENT LIQUID CONTENT  (G/CUBIC METER)              GAM 0009
C      CINDEX=COMPLEX DIELECTRIC CONSTANT M  FROM INDEX                 GAM 0010
C      WAVL = WAVELENGTH IN CM                                          GAM 0011
C                                                                       GAM 0012
      COMPLEX CINDEX                                                    GAM 0013
      IF(RHO.GT.0.) GO TO 2                                             GAM 0014
      GAMFOG=0.                                                         GAM 0015
      RETURN                                                            GAM 0016
    2 CONTINUE                                                          GAM 0017
      KEY=1                                                             GAM 0018
      IF(MR. GE. 5) KEY = 0                                             GAM 0019
      WAVL=1.0/FREQ                                                     GAM 0020
      TC=T-273.2                                                        GAM 0021
CCC                                                                     GAM 0022
CCC    CHANGE TEMP SO THAT MINIMUM IS -20.0 CENT.                       GAM 0023
CCC                                                                     GAM 0024
      IF(TC.LT.-20.0) TC=-20.0                                          GAM 0025
      CALL INDX (WAVL,TC,KEY,REIL,AIMAK)                                GAM 0026
      CINDEX=CMPLX(REIL,AIMAK)                                          GAM 0027
CCC                                                                     GAM 0028
CCC   ATTENUATION = 6.0*PI*FREQ*RHO*IMAG(-K)                            GAM 0029
CCC    6.0*PI/10. = 1.885 (THE FACTOR OF 10 IS FOR UNITS CONVERSION)    GAM 0030
CCC                                                                     GAM 0031
C     GAMFOG=8.1888*FREQ*RHO*AIMAG( -  (CINDEX**2-1)/(CINDEX**2+2))     GAM 0032
      GAMFOG=1.885 *FREQ*RHO*AIMAG( -  (CINDEX**2-1)/(CINDEX**2+2))     GAM 0033
      RETURN                                                            GAM 0034
      END                                                               GAM 0035
