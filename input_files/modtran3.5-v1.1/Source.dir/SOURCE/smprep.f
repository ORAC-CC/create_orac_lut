      SUBROUTINE SMPREP(SPH1,SPH2,SPANGL,SPRANG,SPBETA,ISLCT)           PRP 0001
C                                                                       PRP 0002
C     THIS SUBROUTINE PREPS OR PREPROCESSES PATH GEOMETRY               PRP 0003
C     PARAMETERS TO SEE IS IF THE RANGE IS SMALL.  ALL SMALL PATH       PRP 0004
C     CASES ARE CAST INTO THE EQUIVALENT CASE 2C (H1,H2,RANGE).         PRP 0005
C                                                                       PRP 0006
C     DECLARE INPUTS                                                    PRP 0007
      REAL SPH1,SPH2,SPANGL,SPRANG,SPBETA                               PRP 0008
      INTEGER ISLCT                                                     PRP 0009
C                                                                       PRP 0010
C     LIST COMMONS                                                      PRP 0011
      REAL SMALL                                                        PRP 0012
      COMMON/SMALL3/SMALL                                               PRP 0013
      REAL RE,ZMAX                                                      PRP 0014
      INTEGER IMAX,IMOD,IPATH                                           PRP 0015
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             PRP 0016
      DOUBLE PRECISION H1,H2,RANGE                                      PRP 0017
      COMMON/SMALL4/H1,H2,RANGE                                         PRP 0018
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               PRP 0019
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           PRP 0020
      REAL AORIG                                                        PRP 0021
      COMMON/SMALL5/AORIG                                               PRP 0022
C                                                                       PRP 0023
C     DECLARE LOCAL VARIABLES                                           PRP 0024
      DOUBLE PRECISION ANGLE,BETA,DPDEG,R1,DRE,STORE                    PRP 0025
C                                                                       PRP 0026
C     LIST DATA                                                         PRP 0027
      DATA DPDEG/57.2957795131D0/                                       PRP 0028
C                                                                       PRP 0029
C     AORIG IS THE ORIGINAL INPUT ANGLE IN SINGLE PRECISION             PRP 0030
      AORIG=SPANGL                                                      PRP 0031
C                                                                       PRP 0032
C     DEFINE DOUBLE PRECISION VARIABLES                                 PRP 0033
      H1=DBLE(SPH1)                                                     PRP 0034
      H2=DBLE(SPH2)                                                     PRP 0035
      RANGE=DBLE(SPRANG)                                                PRP 0036
C                                                                       PRP 0037
C     IF ALREADY IN CASE 2C (H1,H2,RANGE) RETURN.                       PRP 0038
      IF(ISLCT.EQ.23)RETURN                                             PRP 0039
C                                                                       PRP 0040
C     DETERMINE RANGE WITHOUT REFRACTION                                PRP 0041
C     (REFRACTION IS NOT IMPORTANT FOR SHORT PATHS).                    PRP 0042
      DRE=DBLE(RE)                                                      PRP 0043
      R1=H1+DRE                                                         PRP 0044
      IF(ISLCT.EQ.24)THEN                                               PRP 0045
C                                                                       PRP 0046
C         CASE 2D (H1,H2,BETA)                                          PRP 0047
          BETA=DBLE(SPBETA)/DPDEG                                       PRP 0048
          RANGE=SQRT((H1-H2)**2+4*R1*(H2+DRE)*SIN(BETA/2)**2)           PRP 0049
      ELSEIF(ISLCT.EQ.22)THEN                                           PRP 0050
C                                                                       PRP 0051
C         CASE 2B (H1,ANGLE,RANGE)                                      PRP 0052
          ANGLE=DBLE(SPANGL)                                            PRP 0053
          STORE=(RANGE-H1)*(RANGE+H1)+2*R1*(H1+RANGE*COS(ANGLE/DPDEG))  PRP 0054
          H2=STORE/(SQRT(DRE**2+STORE)+DRE)                             PRP 0055
      ENDIF                                                             PRP 0056
C                                                                       PRP 0057
C     CHECK FOR SMALL RANGE                                             PRP 0058
      IF(RANGE.LE.SMALL .AND. RANGE.GT.0.)THEN                          PRP 0059
          SPANGL=0.                                                     PRP 0060
          SPBETA=0.                                                     PRP 0061
          SPRANG=REAL(RANGE)                                            PRP 0062
          SPH1=REAL(H1)                                                 PRP 0063
          SPH2=REAL(H2)                                                 PRP 0064
          WRITE(IPR,'(/2A,F12.7,A)')' FROM SMPREP:  RANGE IS',          PRP 0065
     1      ' BELOW THE SMALL RANGE CUTOFF (',SMALL,'KM).'              PRP 0066
          WRITE(IPR,'(/14X,A,/14X,4(A,F12.7))')                         PRP 0067
     1      ' THE INPUT GEOMETRY HAS BEEN CONVERTED TO CASE 2C WITH',   PRP 0068
     2      ' H1 =',SPH1,'KM, H2 =',SPH2,'KM, AND RANGE =', RANGE,'KM.' PRP 0069
      ELSEIF(RANGE.LT.0)THEN                                            PRP 0070
          WRITE(IPR,'(/A)')' FROM SMPREP:  RANGE IS LESS THAN 0.'       PRP 0071
          STOP                                                          PRP 0072
      ENDIF                                                             PRP 0073
      RETURN                                                            PRP 0074
      END                                                               PRP 0075
