      SUBROUTINE FDBETA(H1,H2,BETA,ANGLE,PHI,LENN,HMIN,IERROR)          FBT 0001
C                                                                       FBT 0002
C     GIVEN H1,H2,AND BETA (THE EARTH CENTERED ANGLE), THIS SUBROUTINE  FBT 0003
C     CALCULATES THE ZENITH ANGLE AT H1 (ANGLE) AND AT H2 (PHI).        FBT 0004
C     BASED ON A NEWTON-RAPHSON METHOD.                                 FBT 0005
C                                                                       FBT 0006
C     DECLARE ROUTINE ARGUMENTS:                                        FBT 0007
      DOUBLE PRECISION H1,H2,BETA,ANGLE,PHI,HMIN                        FBT 0008
      INTEGER LENN,IERROR                                               FBT 0009
C                                                                       FBT 0010
C     LIST COMMONS:                                                     FBT 0011
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               FBT 0012
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          FBT 0013
      REAL RE,ZMAX                                                      FBT 0014
      INTEGER IMAX,IMOD,IPATH                                           FBT 0015
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             FBT 0016
      REAL GNDALT                                                       FBT 0017
      COMMON/GRAUND/GNDALT                                              FBT 0018
C                                                                       FBT 0019
C     DECLARE LOCAL VARIABLES:                                          FBT 0020
      INTEGER IAMTB,ITER                                                FBT 0021
      DOUBLE PRECISION HA,HB,BETA1,RANGE,BENDNG,ANGLE1,RATIOA,RATIOB,   FBT 0022
     1  STORE,DENOM,DERIV,APREV,AWGHTD,DBPREV,DBETA,DPRE                FBT 0023
C                                                                       FBT 0024
C     LIST DATA:                                                        FBT 0025
      DOUBLE PRECISION TOLRNC,DPDEG                                     FBT 0026
      INTEGER ITERMX                                                    FBT 0027
      DATA TOLRNC/2.D-4/,DPDEG/57.2957795131D0/,ITERMX/35/              FBT 0028
      DPRE=DBLE(RE)                                                     FBT 0029
      IF(BETA.EQ.0.)THEN                                                FBT 0030
         LENN=0                                                         FBT 0031
         ANGLE=DBLE(0.)                                                 FBT 0032
         HMIN=H1                                                        FBT 0033
         IF(H1.GT.H2)THEN                                               FBT 0034
            ANGLE=DBLE(180.)                                            FBT 0035
            HMIN=H2                                                     FBT 0036
         ENDIF                                                          FBT 0037
         PHI=DBLE(180.)-ANGLE                                           FBT 0038
      ENDIF                                                             FBT 0039
      IF(H1.GT.H2)THEN                                                  FBT 0040
         HA=H2                                                          FBT 0041
         HB=H1                                                          FBT 0042
      ELSE                                                              FBT 0043
         HA=H1                                                          FBT 0044
         HB=H2                                                          FBT 0045
      ENDIF                                                             FBT 0046
C                                                                       FBT 0047
C***  SET PARAMETER TO SUPPRESS CALCULATION OF AMOUNTS IN ROUTINE LAYER FBT 0048
      IAMTB = 2                                                         FBT 0049
C                                                                       FBT 0050
C***  GUESS AT ANGLE,INTEGRATETO FIND BETA, TEST FOR CONVERGENCE, AND   FBT 0051
C***  ITERATE FIRST GUESS AT ANGLE: USE GEOMETRIC SOLN (NO REFRACTION)  FBT 0052
      WRITE(IPR,'(///30H CASE 2D: GIVEN H1, H2,  BETA:,//               FBT 0053
     $     42H ITERATE AROUND ANGLE UNTIL BETA CONVERGES,//             FBT 0054
     $     14H ITER    ANGLE,T21,4HBETA,T30,5HDBETA,T40,5HRANGE,T51,    FBT 0055
     $     4HHMIN,T61,3HPHI,T70,7HBENDING,/T10,5H(DEG),T21,5H(DEG),T30, FBT 0056
     $     5H(DEG),T41,4H(KM),T51,4H(KM),T60,5H(DEG),T71,5H(DEG),/)')   FBT 0057
C                                                                       FBT 0058
C     CALCULATE ANGLE1, A GUESS VALUE FOR ANGLE                         FBT 0059
      RATIOA=(HB-HA)/(DPRE+HA)                                          FBT 0060
      RATIOB=(HB-HA)/(DPRE+HB)                                          FBT 0061
      STORE=2*SIN(DBLE(.5)*BETA/DPDEG)**2                               FBT 0062
      DENOM=RATIOB-STORE                                                FBT 0063
      ANGLE1=DBLE(90.)                                                  FBT 0064
      IF(DENOM.NE.0.)ANGLE1=DPDEG*ATAN(SIN(BETA/DPDEG)/DENOM)           FBT 0065
      IF(ANGLE1.LT.0.)ANGLE1=ANGLE1+DBLE(180.)                          FBT 0066
C                                                                       FBT 0067
C     CALCULATE THE DERIVATIVE D(ANGLE)/D(BETA)                         FBT 0068
      DERIV=(RATIOA+STORE)/(RATIOA*RATIOB+2*STORE)                      FBT 0069
C                                                                       FBT 0070
C     DERIV TENDS TO OVERSHOOT VALUE.  FUDGE (=.6) SPEEDS UP CONVERGENCEFBT 0071
      DERIV=DBLE(.6)*DERIV                                              FBT 0072
C                                                                       FBT 0073
C     BEGIN ITERATIVE PROCEDURE                                         FBT 0074
      ITER=0                                                            FBT 0075
 10   ITER=ITER+1                                                       FBT 0076
      IF(ITER.GT.ITERMX)THEN                                            FBT 0077
         WRITE(IPR,'(40H0FDBETA, CASE 2D (H1,H2,BETA): SOLUTION ,       FBT 0078
     $        16HDID NOT CONVERGE,//10X,4HH1 =,F13.6,4X,4HH2 =,F13.6,   FBT 0079
     $        4X,6HBETA =,F13.6,4X,12HITERATIONS =,I5,//10X,4HLAST,     FBT 0080
     $        10H ITERATION,//10X,7HANGLE =,F16.9,4X,6HBETA =,F16.9)')  FBT 0081
     $        H1,H2,BETA,ITER,ANGLE1,BETA1                              FBT 0082
         IERROR=1                                                       FBT 0083
         RETURN                                                         FBT 0084
      ENDIF                                                             FBT 0085
C                                                                       FBT 0086
C     DETERMINE BETA1, THE BETA CORRESPONDING TO ANGLE1                 FBT 0087
      CALL DPFNMN(HA,ANGLE1,HB,LENN,ITER,HMIN,PHI,IERROR)               FBT 0088
      CALL DPRFPA(HA,HB,ANGLE1,PHI,LENN,HMIN,                           FBT 0089
     $     IAMTB,BETA1,RANGE,BENDNG)                                    FBT 0090
      WRITE(IPR,'(I5,3F10.4,2F10.3,2F10.4)')                            FBT 0091
     $     ITER,ANGLE1,BETA1,BETA-BETA1,RANGE,HMIN,PHI,BENDNG           FBT 0092
      DBETA = BETA1-BETA                                                FBT 0093
      AWGHTD=(ANGLE1*ABS(DBETA)+APREV*ABS(DBPREV))                      FBT 0094
     $     /(ABS(DBETA)+ABS(DBPREV))                                    FBT 0095
C                                                                       FBT 0096
C     CHECK FOR CONVERGENCE                                             FBT 0097
      IF(ABS(BETA-BETA1).LT.TOLRNC)THEN                                 FBT 0098
         IF(HMIN.LT.GNDALT)THEN                                         FBT 0099
            WRITE(IPR,'(3A,//9X,A)')'0FDBETA,',                         FBT 0100
     $           ' CASE 2D(H1,H2,BETA): REFRACTED TANGENT HEIGHT',      FBT 0101
     $           ' IS LESS THAN ZERO-PATH INTERSECTS THE EARTH',        FBT 0102
     $           ' BETA IS TOO LARGE FOR THIS H1 AND H2'                FBT 0103
            IERROR=1                                                    FBT 0104
         ELSEIF(H1.LE.H2)THEN                                           FBT 0105
            BETA=BETA1                                                  FBT 0106
            ANGLE=ANGLE1                                                FBT 0107
         ELSE                                                           FBT 0108
            BETA=BETA1                                                  FBT 0109
            ANGLE=PHI                                                   FBT 0110
            PHI=ANGLE1                                                  FBT 0111
         ENDIF                                                          FBT 0112
         RETURN                                                         FBT 0113
      ENDIF                                                             FBT 0114
C                                                                       FBT 0115
      APREV =  ANGLE1                                                   FBT 0116
      IF (DBPREV*DBETA .LT. 0.0D00 .AND. ITER .GT. 5) THEN              FBT 0117
         ANGLE1=AWGHTD                                                  FBT 0118
         DBPREV = BETA1 - BETA                                          FBT 0119
      ELSE                                                              FBT 0120
         DBPREV = BETA1 - BETA                                          FBT 0121
         ANGLE1=ANGLE1-DERIV*DBPREV                                     FBT 0122
      ENDIF                                                             FBT 0123
      GOTO10                                                            FBT 0124
C                                                                       FBT 0125
      END                                                               FBT 0126
