      SUBROUTINE HERTDA(HERZ,V)                                         HER 0001
C                                                                       HER 0002
C     HERZBERG O2 ABSORPTION                                            HER 0003
C     HALL,1987 PRIVATE COMMUNICATION, BASED ON:                        HER 0004
C                                                                       HER 0005
C     REF. JOHNSTON ET.AL, JGR,89,11661-11665,1984                      HER 0006
C         NICOLET, 1987 (RECENT STUDIES IN ATOMIC & MOLECULAR PROCESSES,HER 0007
C                        PLEMUN PUBLISHING CORP, NY 1987)               HER 0008
C     AND YOSHINO, ET.AL., 1988 (PREPRINT OF "IMPROVED ABSORPTION       HER 0009
C         CROSS SECTIONS OF OXYGEN IN THE WAVELENGTH REGION 205-240NM   HER 0010
C         OF THE HERZBERG CONTINUUM")                                   HER 0011
C                                                                       HER 0012
      HERZ=0.0                                                          HER 0013
      IF(V.LE.36000.00) RETURN                                          HER 0014
C                                                                       HER 0015
C     EXTRAPOLATE SMOOTHLY THROUGH THE HERZBERG BAND REGION             HER 0016
C     NOTE: HERZBERG BANDS ARE NOT CORRECTLY INCLUDED                   HER 0017
C                                                                       HER 0018
      CORR=0.                                                           HER 0019
      IF(V.LE.40000.)CORR=((40000.-V)/4000.)*7.917E-27                  HER 0020
C                                                                       HER 0021
C     CONVERSION TO ATM-CM /KM                                          HER 0022
C                                                                       HER 0023
      RLOSCH = 2.6868 E24 * 1.0E-5                                      HER 0024
C                                                                       HER 0025
C     HALL'S NEW HERZBERG  (LEAST SQRS FIT, LN(P))                      HER 0026
C                                                                       HER 0027
C     YRATIO=2048.7/WL(I)  ****IN ANGSTOMS****                          HER 0028
C           =.20487/WN(I)     IN MICRONS                                HER 0029
C           =WCM(I)/48811.0   IN CM-1                                   HER 0030
C                                                                       HER 0031
      YRATIO= V    /48811.0                                             HER 0032
      HERZ=6.884E-24*(YRATIO)*EXP(-69.738*(ALOG(YRATIO))**2)-CORR       HER 0033
      HERZ = HERZ * RLOSCH                                              HER 0034
      RETURN                                                            HER 0035
      END                                                               HER 0036
