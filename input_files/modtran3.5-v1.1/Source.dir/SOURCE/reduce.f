      SUBROUTINE REDUCE(H1,H2,ANGLE,PHI)                                RED 0001
C                                                                       RED 0002
C     ZMAX IS THE HIGHEST LEVEL IN THE ATMOSPHERIC PROFILE STORED IN    RED 0003
C     COMMON /MODEL/.  IF H1 AND/OR H2 ARE GREATER THAN ZMAX, THIS      RED 0004
C     SUBROUTINE REDUCES THEM TO ZMAX AND RESETS ANGLE AND/OR PHI       RED 0005
C     AS NECESSARY.  THIS REDUCTION IS NECESSARY, FOR EXAMPLE, FOR      RED 0006
C     SATELLITE ALTITUDES, BECAUSE (1) THE DENSITY PROFILES ARE POORLY  RED 0007
C     DEFINED ABOVE ZMAX AND (2) THE CALCULATION TIME FOR PATHS ABOVE   RED 0008
C     ZMAX CAN BE EXCESSIVE (E.G. FOR GEOSYNCHRONOUS ALTITUDES).        RED 0009
C                                                                       RED 0010
C     DECLARE INPUTS                                                    RED 0011
      DOUBLE PRECISION H1,H2,ANGLE,PHI                                  RED 0012
C                                                                       RED 0013
C     LIST COMMONS                                                      RED 0014
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               RED 0015
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           RED 0016
      REAL RE,ZMAX                                                      RED 0017
      INTEGER IMAX,IMOD,IPATH                                           RED 0018
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             RED 0019
C                                                                       RED 0020
C     DECLARE FUNCTIONS                                                 RED 0021
      DOUBLE PRECISION DPANDX                                           RED 0022
C                                                                       RED 0023
C     DECLARE LOCAL VARIABLES                                           RED 0024
      DOUBLE PRECISION CPATH,CZMAX,ANGMAX,SH,GAMMA,DPZMAX,DPRE          RED 0025
C                                                                       RED 0026
C     LIST DATA:                                                        RED 0027
      DOUBLE PRECISION DPDEG                                            RED 0028
      DATA DPDEG/57.2957795131D0/                                       RED 0029
      DPRE=DBLE(RE)                                                     RED 0030
      DPZMAX=DBLE(ZMAX)                                                 RED 0031
      IF(H1.LE.DPZMAX .AND. H2.LE.DPZMAX)RETURN                         RED 0032
      CALL DPFISH(H1,SH,GAMMA)                                          RED 0033
      CPATH=DPANDX(H1,SH,GAMMA)*(DPRE+H1)*SIN(ANGLE/DPDEG)              RED 0034
      CALL DPFISH(DPZMAX,SH,GAMMA)                                      RED 0035
      CZMAX=DPANDX(DPZMAX,SH,GAMMA)*(DPRE+DPZMAX)                       RED 0036
      ANGMAX=DBLE(180.)-ASIN(CPATH/CZMAX)*DPDEG                         RED 0037
      IF(H1.GT.DPZMAX)THEN                                              RED 0038
          H1=DPZMAX                                                     RED 0039
          ANGLE=ANGMAX                                                  RED 0040
      ENDIF                                                             RED 0041
      IF(H2.GT.DPZMAX)THEN                                              RED 0042
          H2=DPZMAX                                                     RED 0043
          PHI=ANGMAX                                                    RED 0044
      ENDIF                                                             RED 0045
      IF(NPR.LE.1)WRITE(IPR,'(//A,/(4X,2A,F10.3,A))')                   RED 0046
     1  ' FROM SUBROUTINE REDUCE: ','ONE OR BOTH OF H1 AND H2 ARE',     RED 0047
     2  ' ABOVE THE TOP OF THE ATMOSPHERIC PROFILE ZMAX = ',ZMAX,       RED 0048
     3  ' KM AND HAVE BEEN RESET TO ZMAX','ANGLE AND/OR PHI HAVE ALSO', RED 0049
     4  ' BEEN RESET TO THE ZENITH ANGLE AT ZMAX = ',ANGMAX,' DEG'      RED 0050
      RETURN                                                            RED 0051
      END                                                               RED 0052
