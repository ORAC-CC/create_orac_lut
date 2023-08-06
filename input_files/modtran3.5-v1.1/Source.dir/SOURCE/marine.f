      SUBROUTINE MARINE(VIS,MODEL,WS,WH,ICSTL,BEXT,BABS,NL)             MAR 0001
      INCLUDE 'PARAM.LST'                                               MAR 0002
C                                                                       MAR 0003
C        THIS SUBROUTINE DETERMINES AEROSOL EXT + ABS COEFFICIENTS      MAR 0004
C          FOR THE NAVY MARITIME MODEL                                  MAR 0005
C            CODED BY STU GATHMAN                  -  NRL               MAR 0006
C                                                                       MAR 0007
C        INPUTS-                                                        MAR 0008
C        WSS = CURRENT WIND SPEED (M/S)                                 MAR 0009
C        WHH = 24 HOUR AVERAGE WIND SPEED (M/S)                         MAR 0010
C        RHH = RELATIVE HUMIDITY (PERCENTAGE)                           MAR 0011
C        VIS = METEOROLOGICAL RANGE (KM)                                MAR 0012
C        ICTL = AIR MASS CHARACTER  1 = OPEN OCEAN                      MAR 0013
C                      10 = STRONG CONTINENTAL INFLUENCE                MAR 0014
C        MODEL = MODEL ATMOSPHERE                                       MAR 0015
C                                                                       MAR 0016
C        OUTPUTS-                                                       MAR 0017
C        BEXT = EXTINCTION COEFFICIENT (KM-1)                           MAR 0018
C        BABS = ABSORPTION COEFFICIENT (KM-1)                           MAR 0019
C                                                                       MAR 0020
      COMMON /MART/ RHH                                                 MAR 0021
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          MAR 0022
C                                                                       MAR 0023
C       PI       THE CONSTANT PI                                        MAR 0024
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       MAR 0025
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       MAR 0026
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         MAR 0027
      REAL PI,DEG,BIGNUM,BIGEXP                                         MAR 0028
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                MAR 0029
      COMMON/A/T1QEXT(40,4),T2QEXT(40,4),T3QEXT(40,4),                  MAR 0030
     1T1QABS(40,4),T2QABS(40,4),T3QABS(40,4),ALAM(40),AREL(4)           MAR 0031
      REAL WSPD(8),RHD(8),BEXT(NAER,MXWVLN),BABS(NAER,MXWVLN)           MAR 0032
      DATA WSPD/6.9, 4.1, 4.1, 10.29, 6.69, 12.35, 7.2, 6.9/            MAR 0033
      DATA RHD/80., 75.63, 76.2, 77.13, 75.24, 80.53, 45.89, 80./       MAR 0034
      PISC = PI/1000.0                                                  MAR 0035
      WRITE(IPR,890)                                                    MAR 0036
C                                                                       MAR 0037
C     CHECK LIMITS OF MODEL VALIDITY                                    MAR 0038
C                                                                       MAR 0039
      RH = RHH                                                          MAR 0040
      IF(RHH.GT.0.) GO TO 10                                            MAR 0041
      RH=RHD(MODEL+1)                                                   MAR 0042
10    IF(WS.GT.20.0) WS=20.                                             MAR 0043
      IF(WH.GT.20.0) WH = 20.                                           MAR 0044
      IF(RH.GT.98.0) RH = 98.                                           MAR 0045
      IF(RH.LT.50.0.AND.RH.GE.0.0) RH = 50.                             MAR 0046
      IF(ICSTL.LT.1.OR.ICSTL.GT.10) ICSTL = 3                           MAR 0047
C                                                                       MAR 0048
C     FIND SIZE DISTRIBUTION PARAMETERS FROM METEOROLOGY INPUT          MAR 0049
C                                                                       MAR 0050
      IF(WH.LE.0.) WRITE(IPR,920)                                       MAR 0051
      IF(WH .LE. 0.0) WH = WSPD(MODEL + 1)                              MAR 0052
      IF(WS.LE.0.) WRITE(IPR,930)                                       MAR 0053
      IF(WS.LE.0.0)WS=WH                                                MAR 0054
      WRITE(IPR,910)WS,WH,RH,ICSTL                                      MAR 0055
C                                                                       MAR 0056
C        F IS A RELATIVE HUMIDITY DEPENDENT GROWTH CORRECTION           MAR 0057
C        TO THE ATTENUATION COEFFICIENT.                                MAR 0058
C                                                                       MAR 0059
      F=((2.-RH/100.)/(6.*(1.-RH/100.)))**0.33333                       MAR 0060
      A1=2000.0*ICSTL*ICSTL                                             MAR 0061
      A2 = AMAX1(5.866*(WH-2.2), 0.5)                                   MAR 0062
CCC   A3 = AMAX1(0.01527*(WS-2.2), 1.14E-5)                             MAR 0063
      A3 = 10**(0.06*WS-2.8)                                            MAR 0064
C                                                                       MAR 0065
C     FIND EXTINCTION AT 0.55 MICRONS AND NORMALIZE TO 1.               MAR 0066
C                                                                       MAR 0067
C     INTERPOLATE FOR RELATIVE HUMIDITY                                 MAR 0068
C                                                                       MAR 0069
      DO 40   J=2,4                                                     MAR 0070
      IF(RH.LE.AREL(J)) GO TO 42                                        MAR 0071
 40   CONTINUE                                                          MAR 0072
 42   DELRH=AREL(J)-AREL(J-1)                                           MAR 0073
      DELRHV=RH-AREL(J-1)                                               MAR 0074
      RATIO=DELRHV/DELRH                                                MAR 0075
      QE1=T1QEXT(4,J-1)+(T1QEXT(4,J)-T1QEXT(4,J-1))*RATIO               MAR 0076
      QE2=T2QEXT(4,J-1)+(T2QEXT(4,J)-T2QEXT(4,J-1))*RATIO               MAR 0077
      QE3=T3QEXT(4,J-1)+(T3QEXT(4,J)-T3QEXT(4,J-1))*RATIO               MAR 0078
      TOTAL = A1*10.**QE1 + A2*10.**QE2 + A3*10.**QE3                   MAR 0079
      EXT55=PISC*TOTAL/F                                                MAR 0080
C                                                                       MAR 0081
C     IF METEOROLOLICAL RANGE NOT SPECIFIED,FIND FROM METEOR DATA       MAR 0082
C                                                                       MAR 0083
      IF(VIS.LE.0.) VIS=3.912/(EXT55+0.01159)                           MAR 0084
      C=(1./EXT55)*(PISC/F)                                             MAR 0085
      A1=C*A1                                                           MAR 0086
      A2=C*A2                                                           MAR 0087
      A3=C*A3                                                           MAR 0088
C                                                                       MAR 0089
C     CALCULATE NORMALIZED ATTENUATION COEFICIENTS                      MAR 0090
C                                                                       MAR 0091
      DO 45   I=1,40                                                    MAR 0092
      T1XV = T1QEXT(I,J-1) + (T1QEXT(I,J) - T1QEXT(I,J-1))*RATIO        MAR 0093
      T2XV = T2QEXT(I,J-1) + (T2QEXT(I,J) - T2QEXT(I,J-1))*RATIO        MAR 0094
      T3XV = T3QEXT(I,J-1) + (T3QEXT(I,J) - T3QEXT(I,J-1))*RATIO        MAR 0095
      T1AV = T1QABS(I,J-1) + (T1QABS(I,J) - T1QABS(I,J-1))*RATIO        MAR 0096
      T2AV = T2QABS(I,J-1) + (T2QABS(I,J) - T2QABS(I,J-1))*RATIO        MAR 0097
      T3AV = T3QABS(I,J-1) + (T3QABS(I,J) - T3QABS(I,J-1))*RATIO        MAR 0098
      BEXT(NL,I)=A1*10**T1XV+A2*10**T2XV+A3*10**T3XV                    MAR 0099
      BABS(NL,I)=A1*10**T1AV+A2*10**T2AV+A3*10**T3AV                    MAR 0100
 45   CONTINUE                                                          MAR 0101
      WRITE(IPR,900) VIS                                                MAR 0102
      RETURN                                                            MAR 0103
890   FORMAT('0MARINE AEROSOL MODEL USED')                              MAR 0104
900   FORMAT('0',T10,'VIS = ',F10.2,' KM')                              MAR 0105
910   FORMAT(T10,'WIND SPEED = ',F8.2,' M/SEC',/,T10,                   MAR 0106
     1 'WIND SPEED (24 HR AVERAGE) = ',F8.2,' M/SEC',/,                 MAR 0107
     2 T10,'RELATIVE HUMIDITY = ',F8.2,' PERCENT',/,                    MAR 0108
     3 T10,'AIRMASS CHARACTER =' ,I3)                                   MAR 0109
920   FORMAT('0  WS NOT SPECIFIED, A DEFAULT VALUE IS USED')            MAR 0110
930   FORMAT('0  WH NOT SPECIFIED, A DEFAULT VALUE IS USED')            MAR 0111
      END                                                               MAR 0112
