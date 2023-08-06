      SUBROUTINE RDNSM                                                  NSM 0001
C                                                                       NSM 0002
C     THIS SUBROUTINE READS MODEL 7 DATA WHEN ISVA EQ 1                 NSM 0003
C                                                                       NSM 0004
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          NSM 0005
      COMMON /CARD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB  NSM 0006
     1  ,MODTRN                                                         NSM 0007
      LOGICAL MODTRN                                                    NSM 0008
      COMMON /CARD1A/ M4,M5,M6,MDEF,IRD1,IRD2                           NSM 0009
      COMMON /CARD1B/ JUNIT(15),WMOL(12),WAIR1,JLOW                     NSM 0010
      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICIR,IVSA,VIS,WSS,WHH,   NSM 0011
     1    RAINRT                                                        NSM 0012
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     NSM 0013
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                NSM 0014
      COMMON /NSINP/ ZMDL(40),P(40),T(40),WMDL(40,13)                   NSM 0015
      CHARACTER*1 JCHAR                                                 NSM 0016
      DIMENSION  JCHAR(15)                                              NSM 0017
      IF(ML. GT. 24) THEN                                               NSM 0018
         WRITE(IPR,900) ML                                              NSM 0019
900      FORMAT('  ML = ',I5,'  GT 24 ML RESET TO 24')                  NSM 0020
         ML = 24                                                        NSM 0021
      ENDIF                                                             NSM 0022
      JLOW = 1                                                          NSM 0023
      DO 200   K=1,ML                                                   NSM 0024
      DO 10 KM = 1,15                                                   NSM 0025
      JCHAR(KM) = ' '                                                   NSM 0026
      IF(KM. GT. 12) GO TO 10                                           NSM 0027
      WMOL(KM) = 0.                                                     NSM 0028
10    CONTINUE                                                          NSM 0029
      IRD0 = 1                                                          NSM 0030
      ICONV = 1                                                         NSM 0031
      IF((MODEL .GT. 0.) .AND. (MODEL .LT. 7)) IRD0 = 0                 NSM 0032
      IF((IRD0  .EQ. 1)  .AND. (IVSA.EQ.1)   ) THEN                     NSM 0033
           IRD0 = 0                                                     NSM 0034
           ICONV =0                                                     NSM 0035
      ENDIF                                                             NSM 0036
C                                                                       NSM 0037
C        PARAMETERS - JCHAR = INPUT KEY (SEE BELOW)                     NSM 0038
C                                                                       NSM 0039
C                                                                       NSM 0040
C     ***  ROUTINE ALSO ACCEPTS VARIABLE UNITS ON PRESS AND TEMP        NSM 0041
C                                                                       NSM 0042
C          SEE INPUT KEY BELOW                                          NSM 0043
C                                                                       NSM 0044
C                                                                       NSM 0045
C                                                                       NSM 0046
C     FOR MOLECULAR SPECIES ONLY                                        NSM 0047
C                                                                       NSM 0048
C       JCHAR   JUNIT                                                   NSM 0049
C                                                                       NSM 0050
C     " ",A      10    VOLUME MIXING RATIO (PPMV)                       NSM 0051
C         B      11    NUMBER DENSITY (CM-3)                            NSM 0052
C         C      12    MASS MIXING RATIO (GM(K)/KG(AIR))                NSM 0053
C         D      13    MASS DENSITY (GM M-3)                            NSM 0054
C         E      14    PARTIAL PRESSURE (MB)                            NSM 0055
C         F      15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY           NSM 0056
C         G      16     "    "     "  (TD IN T(C)) - H2O ONLY           NSM 0057
C         H      17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY     NSM 0058
C         I      18    AVAILABLE FOR USER DEFINITION                    NSM 0059
C        1-6    1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            NSM 0060
C                                                                       NSM 0061
C     ****************************************************************  NSM 0062
C     ****************************************************************  NSM 0063
C                                                                       NSM 0064
C     ***** OTHER 'JCHAR' SPECIFICATIONS -                              NSM 0065
C                                                                       NSM 0066
C       JCHAR   JUNIT                                                   NSM 0067
C                                                                       NSM 0068
C      " ",A     10    PRESSURE IN (MB)                                 NSM 0069
C          B     11       "     "  (ATM)                                NSM 0070
C          C     12       "     "  (TORR)                               NSM 0071
C         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            NSM 0072
C                                                                       NSM 0073
C      " ",A     10    AMBIENT TEMPERATURE IN DEG(K)                    NSM 0074
C          B     11       "         "       "  " (C)                    NSM 0075
C          C     12       "         "       "  " (F)                    NSM 0076
C         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            NSM 0077
C                                                                       NSM 0078
C     ***** DEFINITION OF "DEFAULT" CHOICES FOR PROFILE SELECTION ***** NSM 0079
C                                                                       NSM 0080
C      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL          NSM 0081
C      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS   NSM 0082
C      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE                NSM 0083
C                                                                       NSM 0084
C     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)                         NSM 0085
C                                                                       NSM 0086
C      FOR MOLECULES 8-35, ONLY US STD PROFILES ARE AVIALABLE           NSM 0087
C      THEREFORE, WHEN  'JCHAR(K) = 1-5', JCHAR(K) WILL BE RESET TO 6   NSM 0088
C                                                                       NSM 0089
C                                                                       NSM 0090
      READ(IRD,80)ZMDL(K),P(K),T(K),WMOL(1),WMOL(2),WMOL(3),            NSM 0091
     X (JCHAR(KM),KM=1,15)                                              NSM 0092
80    FORMAT ( F10.3,5E10.3,15A1)                                       NSM 0093
       WRITE(IPR,81)ZMDL(K),P(K),T(K),WMOL(1),WMOL(2),WMOL(3),          NSM 0094
     X (JCHAR(KM),KM=1,15)                                              NSM 0095
81    FORMAT ( F10.3,1P5E10.3,10X,15A1)                                 NSM 0096
      IF(ZMDL(K) .LE. 2.0)JLOW = K                                      NSM 0097
      IF(IRD1 .EQ. 1) THEN                                              NSM 0098
           READ(IRD,83)(WMOL(KM),KM=4,12)                               NSM 0099
83         FORMAT(8E10.3)                                               NSM 0100
           WRITE(IPR,84)(WMOL(KM),KM=4,12)                              NSM 0101
84         FORMAT(1P8E10.3)                                             NSM 0102
      ENDIF                                                             NSM 0103
C                                                                       NSM 0104
C                                                                       NSM 0105
C     AHAZE =  AEROSOL VISIBLE EXTINCTION COFF (KM-1)                   NSM 0106
C     AT A WAVELENGTH OF 0.55 MICROMETERS                               NSM 0107
C                                                                       NSM 0108
C     EQLWCZ=LIQUID WATER CONTENT (PPMV) AT ALT Z                       NSM 0109
C            FOR AEROSOL, CLOUD OR FOG MODELS                           NSM 0110
C                                                                       NSM 0111
C     RRATZ=RAIN RATE (MM/HR) AT ALT Z                                  NSM 0112
C                                                                       NSM 0113
C     IHA1 AEROSOL MODEL USED FOR SPECTRAL DEPENDENCE OF EXTINCTION     NSM 0114
C                                                                       NSM 0115
C     IVUL1 STRATOSPHERIC AEROSOL MODEL USED FOR SPECTRAL DEPENDENCE    NSM 0116
C     OF EXT AT Z                                                       NSM 0117
C                                                                       NSM 0118
C     ICLD1 CLOUD MODEL USED FOR SPECTRAL DEPENDENCE OF EXT AT Z        NSM 0119
C                                                                       NSM 0120
C     ONLY ONE OF IHA1,ICLD1  OR IVUL1 IS ALLOWED                       NSM 0121
C     IHA1 NE 0 OTHERS IGNORED                                          NSM 0122
C     IHA1 EQ 0 AND ICLD1 NE 0 USE ICLD1                                NSM 0123
C                                                                       NSM 0124
C     IF AHAZE AND EQLWCZ ARE BOUTH ZERO                                NSM 0125
C        DEFAULT PROFILE ARE LOADED FROM IHA1,ICLD1,IVUL1               NSM 0126
C     ISEA1 = AEROSOL SEASON CONTROL FOR ALTITUDE Z                     NSM 0127
C                                                                       NSM 0128
      IF(IRD2 .EQ. 1) THEN                                              NSM 0129
           READ(IRD,82)    AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1    NSM 0130
           WRITE(IPR,82)    AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1   NSM 0131
82         FORMAT(10X,3F10.3,4I5)                                       NSM 0132
      ENDIF                                                             NSM 0133
      DO 12 KM = 1,15                                                   NSM 0134
12    JUNIT(KM) = JOU(JCHAR(KM))                                        NSM 0135
      IF(M1 .NE. 0) JUNIT(1) = M1                                       NSM 0136
      IF(M1 .NE. 0) JUNIT(2) = M1                                       NSM 0137
      CALL CHECK(P(K),JUNIT(1),1)                                       NSM 0138
      CALL CHECK(T(K),JUNIT(2),2)                                       NSM 0139
      CALL DEFALT(ZMDL(K),P(K),T(K))                                    NSM 0140
      CALL CONVRT (P(K),T(K) )                                          NSM 0141
      DO 20 KM = 1,12                                                   NSM 0142
20    WMDL(K,KM) = WMOL(KM)                                             NSM 0143
      WMDL(K,13) = WAIR1                                                NSM 0144
200   CONTINUE                                                          NSM 0145
      RETURN                                                            NSM 0146
      END                                                               NSM 0147
