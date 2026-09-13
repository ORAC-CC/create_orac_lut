      SUBROUTINE SMGEO(SPANGL,SPBETA,SPPHI,                             SMG 0001
     1  DHALFR,DPRNG2,BENDNG,LENN,SMMIN)                                SMG 0002
      REAL SPANGL,SPPHI,BENDNG,SPBETA,AORIG                             SMG 0003
      DOUBLE PRECISION DHALFR,DPRNG2,X,Y,BETA,RANGE,H1,H2,              SMG 0004
     1  DPDEG,PHI,ANGLE,SMMIN,COSPHI,COSANG,DPRE                        SMG 0005
      INTEGER LENN                                                      SMG 0006
      REAL RE,ZMAX                                                      SMG 0007
      INTEGER IMAX,IMOD,IPATH                                           SMG 0008
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             SMG 0009
      COMMON /SMALL4/H1, H2, RANGE                                      SMG 0010
C                                                                       SMG 0011
C     AORIG IS THE ORIGINAL INPUT ANGLE IN SINGLE PRECISION             SMG 0012
      COMMON /SMALL5/AORIG                                              SMG 0013
      DATA DPDEG/57.2957795131D0/                                       SMG 0014
      DPRE=DBLE(RE)                                                     SMG 0015
      X=H2-H1                                                           SMG 0016
      Y=(RANGE+X)*(RANGE-X)                                             SMG 0017
      IF(Y.LT.0.)Y=DBLE(0.)                                             SMG 0018
      BETA=2*DPDEG*ASIN(SQRT(Y/((H1+DPRE)*(H2+DPRE)))/2)                SMG 0019
C                                                                       SMG 0020
C     IF H2 > H1, CALCULATE PHI (>90 DEG) FROM THE COSINE LAW           SMG 0021
C     AND ANGLE FROM SUBTRACTION;  IF H1 > H2, CALCULATE ANGLE          SMG 0022
C     (>90 DEG) FROM THE COSINE LAW AND PHI FROM SUBTRACTION.           SMG 0023
      IF(X.GE.0.)THEN                                                   SMG 0024
         COSPHI=-(RANGE**2+X*(H1+H2+2*DPRE))/(2*RANGE*(H2+DPRE))        SMG 0025
         PHI=DBLE(180.)                                                 SMG 0026
         IF(COSPHI.GT.-1.)PHI=DPDEG*ACOS(COSPHI)                        SMG 0027
         ANGLE=DBLE(180.)+BETA-PHI                                      SMG 0028
      ELSE                                                              SMG 0029
         COSANG=-(RANGE**2-X*(H1+H2+2*DPRE))/(2*RANGE*(H1+DPRE))        SMG 0030
         ANGLE=DBLE(180.)                                               SMG 0031
         IF(COSANG.GT.-1.)ANGLE=DPDEG*ACOS(COSANG)                      SMG 0032
         PHI=DBLE(180.)+BETA-ANGLE                                      SMG 0033
      ENDIF                                                             SMG 0034
      BENDNG=0.                                                         SMG 0035
      IF(ANGLE.GT.90.D0 .AND. PHI.GT.90.D0 .AND.                        SMG 0036
     1  (AORIG.GT.90.D0 .OR. AORIG.EQ.0.))THEN                          SMG 0037
C        THE CONDITION ON AORIG IS TO MAKE SURE THAT PATHS WHOSE        SMG 0038
C        ORIGINAL INPUT ZENITH IS 90 DEGREES OR LESS ARE TREATED AS     SMG 0039
C        LENN=0 PATHS.  CALCULATED ANGLE WILL LIKELY                    SMG 0040
C        DIFFER FROM AORIG.                                             SMG 0041
         LENN=1                                                         SMG 0042
         IF(H1.LE.H2)THEN                                               SMG 0043
             DHALFR=(H1+DPRE)*COS((DBLE(180.)-ANGLE)/DPDEG)             SMG 0044
         ELSE                                                           SMG 0045
             DHALFR=(H2+DPRE)*COS((DBLE(180.)-PHI)/DPDEG)               SMG 0046
         ENDIF                                                          SMG 0047
         DPRNG2=RANGE-2*DHALFR                                          SMG 0048
         IF (DPRNG2 .LE. 1.0D-06) THEN                                  SMG 0049
            DPRNG2=DBLE(0.)                                             SMG 0050
            DHALFR=RANGE/2                                              SMG 0051
         ENDIF                                                          SMG 0052
         SMMIN=(DPRE+H1)*SIN((DBLE(180.)-ANGLE)/DPDEG)-DPRE             SMG 0053
      ELSE                                                              SMG 0054
         LENN=0                                                         SMG 0055
         DHALFR=DBLE(0.)                                                SMG 0056
         DPRNG2=RANGE                                                   SMG 0057
         SMMIN=MIN(H1,H2)                                               SMG 0058
      ENDIF                                                             SMG 0059
      SPANGL=REAL(ANGLE)                                                SMG 0060
      SPBETA=REAL(BETA)                                                 SMG 0061
      SPPHI=REAL(PHI)                                                   SMG 0062
      RETURN                                                            SMG 0063
      END                                                               SMG 0064
