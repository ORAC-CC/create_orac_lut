      REAL FUNCTION HENGNS(G,COSANG)                                    HEN 0001
C                                                                       HEN 0002
C     CALCULATES THE HENYEY-GREENSTEIN SCATTERING PHASE FUNCTION        HEN 0003
C                                                                       HEN 0004
C     INPUTS                                                            HEN 0005
C       G        HENYEY-GREENSTEIN ASYMMETRY FACTOR                     HEN 0006
C       COSANG   COSINE OF THE SCATTERING ANGLE                         HEN 0007
      REAL G,COSANG                                                     HEN 0008
C                                                                       HEN 0009
C     OUTPUTS                                                           HEN 0010
C       HENGNS   HENYEY-GREENSTEIN SCATTERING PHASE FUNCTION [STER-1]   HEN 0011
C                                                                       HEN 0012
C     COMMONS                                                           HEN 0013
C       PI       THE CONSTANT PI                                        HEN 0014
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       HEN 0015
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       HEN 0016
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         HEN 0017
      REAL PI,DEG,BIGNUM,BIGEXP                                         HEN 0018
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                HEN 0019
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               HEN 0020
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           HEN 0021
C                                                                       HEN 0022
C     LOCAL VARIABLES                                                   HEN 0023
      REAL ONEMG,DENOM                                                  HEN 0024
      ONEMG=1.-G                                                        HEN 0025
      DENOM=ONEMG**2+2*G*(1.-COSANG)                                    HEN 0026
      DENOM=4.*PI*DENOM*SQRT(DENOM)                                     HEN 0027
      IF(DENOM.GE.0.)THEN                                               HEN 0028
          HENGNS=ONEMG*(1.+G)/DENOM                                     HEN 0029
          RETURN                                                        HEN 0030
      ENDIF                                                             HEN 0031
      WRITE(IPR,'(2A,/14X,2(A,F9.5))')' FATAL ERROR: ',                 HEN 0032
     1  ' HENYEY-GREENSTEIN PHASE FUNCTION IS NOT WELL DEFINED FOR',    HEN 0033
     2  ' ASYMMETRY FACTOR (G) =',G,' AND COS(ANGLE) =',COSANG          HEN 0034
      RETURN                                                            HEN 0035
      END                                                               HEN 0036
