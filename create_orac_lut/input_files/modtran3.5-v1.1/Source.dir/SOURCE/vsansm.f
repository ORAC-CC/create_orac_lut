      SUBROUTINE VSANSM(K,AHAZE,IHA1,ZNEW)                              VNM 0001
      INCLUDE 'PARAM.LST'                                               VNM 0002
      INTEGER KPOINT                                                    VNM 0003
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     VNM 0004
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   VNM 0005
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   VNM 0006
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     VNM 0007
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          VNM 0008
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     VNM 0009
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                VNM 0010
      COMMON /CARD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB  VNM 0011
     1  ,MODTRN                                                         VNM 0012
      LOGICAL MODTRN                                                    VNM 0013
      COMMON /CARD1B/ JUNIT(15),WMOLI(12),WAIR1,JLOW                    VNM 0014
      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICIR,IVSA,VIS,WSS,WHH,   VNM 0015
     1    RAINRT                                                        VNM 0016
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)            VNM 0017
      COMMON /NSINP/ ZMDL(40),PM(40),TM(40),WMDL(40,13)                 VNM 0018
C                                                                       VNM 0019
      COMMON /MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),        VNM 0020
     X WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)     VNM 0021
      COMMON /MDATA1/ WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),            VNM 0022
     X WNH3(LAYDIM),WAIR(LAYDIM),WHNO3(LAYDIM)                          VNM 0023
      DIMENSION WMOL(13)                                                VNM 0024
C                                                                       VNM 0025
C     OUTPUT COMMON MDATA AND MDATA1                                    VNM 0026
C                                                                       VNM 0027
C                                                                       VNM 0028
C     MODEL 7 CODING                                                    VNM 0029
C     OLD LAYERS  AEROSOL RETURNED                                      VNM 0030
C     NEW LAYERS P,T,DP,AEROSOL                                         VNM 0031
C                                                                       VNM 0032
C                                                                       VNM 0033
C                                                                       VNM 0034
      JML=ML                                                            VNM 0035
      J=1                                                               VNM 0036
      KN=K                                                              VNM 0037
110   IF(KN.GT.10)GO TO 140                                             VNM 0038
      JL=J-1                                                            VNM 0039
      IF(JL.LT.1)JL=1                                                   VNM 0040
      JP=JL+1                                                           VNM 0041
      JLS = JL                                                          VNM 0042
      IF(ZVSA(KN).EQ.ZMDL  (JL))GO TO 140                               VNM 0043
      JLS = JP                                                          VNM 0044
      IF(ZVSA(KN).EQ.ZMDL  (JP))GO TO 140                               VNM 0045
      IF(ZVSA(KN).GT.ZMDL  (JL).AND.ZVSA(KN).LT.ZMDL  (JP))GO TO 115    VNM 0046
      IF(J. GE. JML) GO TO 115                                          VNM 0047
      J = J + 1                                                         VNM 0048
      GO TO 110                                                         VNM 0049
115   ZNEW=ZVSA(KN)                                                     VNM 0050
      DIF=ZMDL  (JP)-ZMDL  (JL)                                         VNM 0051
      DZ=ZVSA(KN)-ZMDL  (JL)                                            VNM 0052
      DLIN=DZ/DIF                                                       VNM 0053
      P(K)  = (PM(JP)-PM(JL))*DLIN+PM(JL)                               VNM 0054
      T(K)   =(TM(JP)-TM(JL))*DLIN+TM(JL)                               VNM 0055
      DO 120 KM = 1,13                                                  VNM 0056
      WMOL(KM)=(WMDL(JP,KM)-WMDL(JL,KM))*DLIN+WMDL(JL,KM)               VNM 0057
120   CONTINUE                                                          VNM 0058
      IHA1  =IHVSA(KN)                                                  VNM 0059
      AHAZE  =AHVSA(KN)                                                 VNM 0060
      FAC=(ZVSA(KN)-ZMDL  (JL))/DIF                                     VNM 0061
      IF(PM(JP).GT.0.0.AND.PM(JL).GT.0.) THEN                           VNM 0062
           P(K)  =PM(JL)*(PM(JP)/PM(JL))**FAC                           VNM 0063
      ENDIF                                                             VNM 0064
      IF(TM(JP).GT.0.0.AND.TM(JL).GT.0.) THEN                           VNM 0065
           T(K)   =TM(JL)*(TM(JP)/TM(JL))**FAC                          VNM 0066
      ENDIF                                                             VNM 0067
      DO 130 KM = 1,13                                                  VNM 0068
      IF(WMDL(JP,KM) .GT.0.0.AND.WMDL(JL,KM).GT.0.0) THEN               VNM 0069
           WMOL(KM)=(WMDL(JL,KM)*(WMDL(JP,KM))/WMDL(JL,KM))**FAC        VNM 0070
      ENDIF                                                             VNM 0071
130   CONTINUE                                                          VNM 0072
       WH(K)    = WMOL(1)                                               VNM 0073
       WCO2(K)  = WMOL(2)                                               VNM 0074
       WO(K)    = WMOL(3)                                               VNM 0075
       WN2O(K)  = WMOL(4)                                               VNM 0076
       WCO(K)   = WMOL(5)                                               VNM 0077
       WCH4(K)  = WMOL(6)                                               VNM 0078
       WO2(K)   = WMOL(7)                                               VNM 0079
       WNO(K)   = WMOL(8)                                               VNM 0080
       WSO2(K)  = WMOL(9)                                               VNM 0081
       WNO2(K)  = WMOL(10)                                              VNM 0082
       WNH3(K)  = WMOL(11)                                              VNM 0083
       WHNO3(K) = WMOL(12)                                              VNM 0084
       WAIR(K)  = WMOL(13)                                              VNM 0085
      RETURN                                                            VNM 0086
140   CONTINUE                                                          VNM 0087
      J = JLS                                                           VNM 0088
      IF(K.GT.10) THEN                                                  VNM 0089
         J = K - 10 + JLOW                                              VNM 0090
         IHA1  =0                                                       VNM 0091
         AHAZE  =0.                                                     VNM 0092
      ENDIF                                                             VNM 0093
      ZNEW = ZMDL(J)                                                    VNM 0094
      P(K)  =PM(J)                                                      VNM 0095
      T(K) = TM(J)                                                      VNM 0096
      DO 135 KM = 1,13                                                  VNM 0097
      WMOL(KM)= WMDL(J,KM)                                              VNM 0098
135   CONTINUE                                                          VNM 0099
       WH(K)    = WMOL(1)                                               VNM 0100
       WCO2(K)  = WMOL(2)                                               VNM 0101
       WO(K)    = WMOL(3)                                               VNM 0102
       WN2O(K)  = WMOL(4)                                               VNM 0103
       WCO(K)   = WMOL(5)                                               VNM 0104
       WCH4(K)  = WMOL(6)                                               VNM 0105
       WO2(K)   = WMOL(7)                                               VNM 0106
       WNO(K)   = WMOL(8)                                               VNM 0107
       WSO2(K)  = WMOL(9)                                               VNM 0108
       WNO2(K)  = WMOL(10)                                              VNM 0109
       WNH3(K)  = WMOL(11)                                              VNM 0110
       WHNO3(K) = WMOL(12)                                              VNM 0111
       WAIR(K)  = WMOL(13)                                              VNM 0112
      IF(KN.LE.9) IHA1  =IHVSA(KN)                                      VNM 0113
      IF(KN.LE.9)AHAZE  =AHVSA(KN)                                      VNM 0114
      RETURN                                                            VNM 0115
      END                                                               VNM 0116
