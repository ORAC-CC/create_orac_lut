      SUBROUTINE CONVRT (P,T)                                           CNV 0001
C*************************************************************          CNV 0002
C                                                                       CNV 0003
C     WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR          CNV 0004
C     UNIFORM DATA INPUT -                                              CNV 0005
C                                                                       CNV 0006
C     JCHAR    JUNIT                                                    CNV 0007
C                                                                       CNV 0008
C     " ",A       10    VOLUME MIXING RATIO (PPMV)                      CNV 0009
C     B       11    NUMBER DENSITY (CM-3)                               CNV 0010
C     C       12    MASS MIXING RATIO (GM(K)/KG(AIR))                   CNV 0011
C     D       13    MASS DENSITY (GM M-3)                               CNV 0012
C     E       14    PARTIAL PRESSURE (MB)                               CNV 0013
C     F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY              CNV 0014
C     G       16     "    "     "  (TD IN T(C)) - H2O ONLY              CNV 0015
C     H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY        CNV 0016
C     I       18    AVAILABLE FOR USER DEFINITION                       CNV 0017
C     J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE       CNV 0018
C                                                                       CNV 0019
C***************************************************************        CNV 0020
C                                                                       CNV 0021
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          CNV 0022
      COMMON /CONSTN/ PZERO,TZERO,AVOGAD,ALOSMT,GASCON,PLANK,BOLTZ,     CNV 0023
     1     CLIGHT,ADCON,ALZERO,AVMWT,AIRMWT,AMWT(35)                    CNV 0024
      COMMON /CARD1B/ JUNITP,JUNITT,JUNIT1(13),WMOL(12),WAIR1,JLOW      CNV 0025
C                                                                       CNV 0026
C                                                                       CNV 0027
C                                                                       CNV 0028
C     CONVENTION                                                        CNV 0029
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         CNV 0030
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            CNV 0031
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             CNV 0032
C                                                                       CNV 0033
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              CNV 0034
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       CNV 0035
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     CNV 0036
C                                                                       CNV 0037
      INCLUDE 'PARAM.LST'                                               CNV 0038
      REAL AMWTX(MMOLX)                                                 CNV 0039
      COMMON /ATMWTX/AMWTX                                              CNV 0040
      COMMON /CRD1BX/JUNITX, WMOLX(MMOLX)                               CNV 0041
C*****                                                                  CNV 0042
C                                                                       CNV 0043
C     VARIABLES ENDING WITH "X" ARE INTRODUCED TO DEAL WITH THE         CNV 0044
C     NSPECX EXTRA MOLECULES.                                           CNV 0045
C                                                                       CNV 0046
      RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)                               CNV 0047
C                                                                       CNV 0048
      DO 200 K = 1,NSPC+NSPECX                                          CNV 0049
         IF (K .LE. NSPC) THEN                                          CNV 0050
            B = AVOGAD/AMWT(K)                                          CNV 0051
            R = AIRMWT/AMWT(K)                                          CNV 0052
            JUNIT = JUNIT1(K)                                           CNV 0053
            WHOLD  = WMOL(K)                                            CNV 0054
         ELSE IF ( K .GT. NSPC) THEN                                    CNV 0055
            J = K-NSPC                                                  CNV 0056
            B = AVOGAD/AMWTX(J)                                         CNV 0057
            R = AIRMWT/AMWTX(J)                                         CNV 0058
            JUNIT = JUNITX                                              CNV 0059
            WHOLD  = WMOLX(J)                                           CNV 0060
         ENDIF                                                          CNV 0061
C                                                                       CNV 0062
         IF(K.EQ.1) THEN                                                CNV 0063
            CALL WATVAP(P,T)                                            CNV 0064
            GO TO 200                                                   CNV 0065
         ENDIF                                                          CNV 0066
C                                                                       CNV 0067
         IF (JUNIT .LE. 10) THEN                                        CNV 0068
C           ***** GIVEN VOL. MIXING RATIO                               CNV 0069
            GO TO 200                                                   CNV 0070
         ELSE IF (JUNIT .EQ. 11) THEN                                   CNV 0071
C           *****GIVEN NUMBER DENSITY (CM-3)                            CNV 0072
            IF (K .LE. NSPC) WMOL(K)=WHOLD/(RHOAIR*1.E-6)               CNV 0073
            IF (K .GT. NSPC) WMOLX(J)=WHOLD/(RHOAIR*1.E-6)              CNV 0074
            GO TO 200                                                   CNV 0075
         ELSE IF (JUNIT .EQ. 12) THEN                                   CNV 0076
C           *****GIVEN MASS MIXING RATIO (GM KG-1)                      CNV 0077
            IF (K .LE. NSPC) WMOL(K)= R*WHOLD*1.0E+3                    CNV 0078
            IF (K .GT. NSPC) WMOLX(J)= R*WHOLD*1.0E+3                   CNV 0079
            GO TO 200                                                   CNV 0080
         ELSE IF (JUNIT .EQ. 13) THEN                                   CNV 0081
C           *****GIVEN MASS DENSITY (GM M-3)                            CNV 0082
            IF (K .LE. NSPC) WMOL(K) = B*WHOLD/RHOAIR                   CNV 0083
            IF (K .GT. NSPC) WMOLX(J)= B*WHOLD/RHOAIR                   CNV 0084
            GO TO 200                                                   CNV 0085
         ELSE IF (JUNIT .EQ. 14) THEN                                   CNV 0086
C           *****GIVEN  PARTIAL PRESSURE (MB)                           CNV 0087
            WTEM    = ALOSMT*(WHOLD/PZERO)*(TZERO/T)                    CNV 0088
            IF (K . LE. NSPC) WMOL(K)= WTEM/(RHOAIR*1.E-6)              CNV 0089
            IF (K . GT. NSPC) WMOLX(J)=WTEM/(RHOAIR*1.E-6)              CNV 0090
            GO TO 200                                                   CNV 0091
         ELSE                                                           CNV 0092
            WRITE (IPR,951)JUNIT                                        CNV 0093
 951        FORMAT(/,'   **** ERROR IN CONVERT ****, JUNIT = ',I5)      CNV 0094
            GO TO 200                                                   CNV 0095
         ENDIF                                                          CNV 0096
 200  CONTINUE                                                          CNV 0097
      END                                                               CNV 0098
