      SUBROUTINE FRQ5DT(LOOP0,IV)                                       5DT 0001
C                                                                       5DT 0002
C     THIS ROUTINE DEFINES THE LAYER INDEPENDENT 5 CM-1 DATA            5DT 0003
      LOGICAL LOOP0,MODTRN                                              5DT 0004
C                                                                       5DT 0005
C       PI       THE CONSTANT PI                                        5DT 0006
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       5DT 0007
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       5DT 0008
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         5DT 0009
      REAL PI,DEG,BIGNUM,BIGEXP                                         5DT 0010
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                5DT 0011
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   5DT 0012
     1  MODTRN                                                          5DT 0013
      COMMON/AABBCC/A1(11),B1(11),C1(11),IBND(11),QA(11),CPS(11)        5DT 0014
      COMMON/WNLOHI/IWLH2O(15),IWLO3(6),IWLCO2(11),IWLCO(4),IWLCH4(5),  5DT 0015
     1  IWLN2O(12),IWLO2(7),IWLNH3(3),IWLNO(2),IWLNO2(4),IWLSO2(5),     5DT 0016
     2              IWHH2O(15),IWHO3(6),IWHCO2(11),IWHCO(4),IWHCH4(5),  5DT 0017
     3  IWHN2O(12),IWHO2(7),IWHNH3(3),IWHNO(2),IWHNO2(4),IWHSO2(5)      5DT 0018
      COMMON/H2O/CPH2O(3515)                                            5DT 0019
      COMMON/O3/CPO3(447)                                               5DT 0020
      COMMON/UFMIX1/CPCO2(1219)                                         5DT 0021
      COMMON/UFMIX2/CPCO(173),CPCH4(493),CPN2O(704),CPO2(382)           5DT 0022
      COMMON/TRACEG/CPNH3(431),CPNO(62),CPNO2(142),CPSO2(226)           5DT 0023
      COMMON/FRQ5/SIGO20,SIGO2A,SIGO2B,C0,CT1,CT2,ABB(19),CRSNO2,CRSSO2 5DT 0024
      SAVE INDH2O,INDO3,INDCO2,INDCO,INDCH4,INDN2O,                     5DT 0025
     1  INDO2,INDNH3,INDNO,INDNO2,INDSO2                                5DT 0026
C     INITIALIZE DATA                                                   5DT 0027
      DATA MNAER/1/                                                     5DT 0028
      IF(LOOP0)THEN                                                     5DT 0029
C                                                                       5DT 0030
C         INITIALIZE FREQUENCY REGION INDEX FOR EACH MOLECULE           5DT 0031
          INDH2O=1                                                      5DT 0032
          INDO3 =1                                                      5DT 0033
          INDCO2=1                                                      5DT 0034
          INDCO =1                                                      5DT 0035
          INDCH4=1                                                      5DT 0036
          INDN2O=1                                                      5DT 0037
          INDO2 =1                                                      5DT 0038
          INDNH3=1                                                      5DT 0039
          INDNO =1                                                      5DT 0040
          INDNO2=1                                                      5DT 0041
          INDSO2=1                                                      5DT 0042
          RETURN                                                        5DT 0043
      ENDIF                                                             5DT 0044
      V=FLOAT(IV)                                                       5DT 0045
      IF(.NOT.MODTRN)THEN                                               5DT 0046
          CPS( 1)=CXDTA(V,IWLH2O,IWHH2O,CPH2O,INDH2O)                   5DT 0047
          CPS( 2)=CXDTA(V,IWLCO2,IWHCO2,CPCO2,INDCO2)                   5DT 0048
          CPS( 3)=CXDTA(V,IWLO3 ,IWHO3 ,CPO3 ,INDO3 )                   5DT 0049
          CPS( 4)=CXDTA(V,IWLN2O,IWHN2O,CPN2O,INDN2O)                   5DT 0050
          CPS( 5)=CXDTA(V,IWLCO, IWHCO, CPCO, INDCO )                   5DT 0051
          CPS( 6)=CXDTA(V,IWLCH4,IWHCH4,CPCH4,INDCH4)                   5DT 0052
          CPS( 7)=CXDTA(V,IWLO2, IWHO2, CPO2, INDO2 )                   5DT 0053
          CPS( 8)=CXDTA(V,IWLNO, IWHNO, CPNO, INDNO )                   5DT 0054
          CPS( 9)=CXDTA(V,IWLSO2,IWHSO2,CPSO2,INDSO2)                   5DT 0055
          CPS(10)=CXDTA(V,IWLNO2,IWHNO2,CPNO2,INDNO2)                   5DT 0056
          CPS(11)=CXDTA(V,IWLNH3,IWHNH3,CPNH3,INDNH3)                   5DT 0057
      ENDIF                                                             5DT 0058
C                                                                       5DT 0059
C     N2 CONTINUUM                                                      5DT 0060
      CALL C4DTA(ABB(4),V)                                              5DT 0061
      CALL ABCDTA(IV)                                                   5DT 0062
C                                                                       5DT 0063
      CALL SLF296(V,SH2OT0)                                             5DT 0064
      CALL SLF260(V,SH2OT1)                                             5DT 0065
      CALL FRN296(V,FH2O)                                               5DT 0066
C                                                                       5DT 0067
C                                                                       5DT 0068
C     RADIATION FIELD                                                   5DT 0069
C                                                                       5DT 0070
                                                                        5DT 0071
      VTEMP=1.438786*V                                                  5DT 0072
      IF(VTEMP/260..GT.BIGEXP)THEN                                      5DT 0073
          RADFN0=V                                                      5DT 0074
          RADFN1=V                                                      5DT 0075
      ELSE                                                              5DT 0076
          STORE=EXP(-VTEMP/296.)                                        5DT 0077
          RADFN0=V*(1.-STORE)/(1.+STORE)                                5DT 0078
          STORE=EXP(-VTEMP/260.)                                        5DT 0079
          RADFN1=V*(1.-STORE)/(1.+STORE)                                5DT 0080
      ENDIF                                                             5DT 0081
      ABB(5)=SH2OT0*RADFN0                                              5DT 0082
C                                                                       5DT 0083
C   CO2 CONTINUM                                                        5DT 0084
C                                                                       5DT 0085
CC    CALL CCO2  (V,FCO2)                                               5DT 0086
C                                                                       5DT 0087
CC    FCO2 = FCO2 * RADFN0                                              5DT 0088
      CALL C6DTA(ABB(6),V)                                              5DT 0089
CC    ABB(9)=SH2OT1*RADFN1-SH2OT0*RADFN0                                5DT 0090
      ABB(9)=SH2OT1*RADFN1                                              5DT 0091
      ABB(10)=FH2O*RADFN0                                               5DT 0092
C     HNO3 ABSORPTION CALCULATION                                       5DT 0093
      ABB(11)=0.                                                        5DT 0094
      IF(.NOT.MODTRN)CALL HNO3(V,ABB(11))                               5DT 0095
      CALL AEREXT(V,MNAER)                                              5DT 0096
C                                                                       5DT 0097
C     O2 CONTINUUM CONTRIBUTIONS                                        5DT 0098
          CALL HERTDA(ABB(17),V)                                        5DT 0099
      CALL O2CONT(V,SIGO20,SIGO2A,SIGO2B)                               5DT 0100
      IF(V.GE.49600.)CALL SCHRUN(V,CPS(7))                              5DT 0101
C                                                                       5DT 0102
C     DIFFUSE OZONE                                                     5DT 0103
      C0=0.                                                             5DT 0104
      IF(V.GE.40800.)THEN                                               5DT 0105
          CALL O3UV(V,C0)                                               5DT 0106
          CT1=0.                                                        5DT 0107
          CT2=0.                                                        5DT 0108
          ABB(8)=0.                                                     5DT 0109
      ELSEIF(V.GT.24565.)THEN                                           5DT 0110
          CALL O3HHT0(V,C0)                                             5DT 0111
          CALL O3HHT1(V,CT1)                                            5DT 0112
          CALL O3HHT2(V,CT2)                                            5DT 0113
          ABB(8)=0.                                                     5DT 0114
      ELSE                                                              5DT 0115
          CALL O3CHAP(V,C0,CT1,CT2)                                     5DT 0116
          ABB(8)=0.                                                     5DT 0117
      ENDIF                                                             5DT 0118
C                                                                       5DT 0119
C     NO2 CROSS-SECTIONS (14000 CM-1 TO 50000 CM-1)                     5DT 0120
      CALL NO2XS(V,CRSNO2)                                              5DT 0121
C                                                                       5DT 0122
C     SO2 CROSS-SECTIONS (14000 CM-1 TO 50000 CM-1)                     5DT 0123
      CALL SO2XS(V,CRSSO2)                                              5DT 0124
      RETURN                                                            5DT 0125
      END                                                               5DT 0126
