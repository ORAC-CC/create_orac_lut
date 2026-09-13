      SUBROUTINE CLDPRF(K,ICLD,IHA1,IC1,ICHIC1,HAZECZ)                  CLD 0001
C***********************************************************************CLD 0002
C     WILL COMPUTE DENSITY    PROFILES FOR CLOUDS                       CLD 0003
C***********************************************************************CLD 0004
      REAL MDLWC                                                        CLD 0005
      INCLUDE 'PARAM.LST'                                               CLD 0006
      INTEGER KPOINT                                                    CLD 0007
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     CLD 0008
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   CLD 0009
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   CLD 0010
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     CLD 0011
      REAL ZM,PM,TM,RFNDX,DENSTY                                        CLD 0012
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    CLD 0013
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               CLD 0014
      DIMENSION RHZONE(4)                                               CLD 0015
      DIMENSION ELWCR(4),ELWCU(4),ELWCM(4),ELWCT(4)                     CLD 0016
      DATA RHZONE/0.,70.,80.,99./                                       CLD 0017
      DATA ELWCR/3.517E-04,3.740E-04,4.439E-04,9.529E-04/               CLD 0018
      DATA ELWCM/4.675E-04,6.543E-04,1.166E-03,3.154E-03/               CLD 0019
      DATA ELWCU/3.102E-04,3.802E-04,4.463E-04,9.745E-04/               CLD 0020
      DATA ELWCT/1.735E-04,1.820E-04,2.020E-04,2.408E-04/               CLD 0021
      DATA AFLWC/1.295E-02/,RFLWC/1.804E-03/,CULWC/7.683E-03/           CLD 0022
      DATA ASLWC/4.509E-03/,STLWC/5.272E-03/,SCLWC/4.177E-03/           CLD 0023
      DATA SNLWC/7.518E-03/,BSLWC/1.567E-04/,FVLWC/5.922E-04/           CLD 0024
      DATA AVLWC/1.675E-04/,MDLWC/4.775E-04/                            CLD 0025
      DATA TNLWC/3.446E-3/ ,TKLWC/5.811E-2/                             CLD 0026
      HAZECZ=0.                                                         CLD 0027
      IF(DENSTY(66,K).LE.0.)GOTO15                                      CLD 0028
      IH=ICLD                                                           CLD 0029
      IF(IH .EQ. 0) GO TO 200                                           CLD 0030
      IF(ICLD .EQ. 18) THEN                                             CLD 0031
           HAZECZ=DENSTY(66,K)/TNLWC                                    CLD 0032
           RETURN                                                       CLD 0033
      ENDIF                                                             CLD 0034
      IF(ICLD .EQ. 19) THEN                                             CLD 0035
           HAZECZ=DENSTY(66,K)/TKLWC                                    CLD 0036
           RETURN                                                       CLD 0037
      ENDIF                                                             CLD 0038
      IF(ICLD .EQ. 20) THEN                                             CLD 0039
           RETURN                                                       CLD 0040
      ENDIF                                                             CLD 0041
      GO TO(114,115,116,117,118,116,118,118,114,114,114),IH             CLD 0042
114   HAZECZ=DENSTY(66,K)/CULWC                                         CLD 0043
      RETURN                                                            CLD 0044
115   HAZECZ=DENSTY(66,K)/ASLWC                                         CLD 0045
      RETURN                                                            CLD 0046
116   HAZECZ=DENSTY(66,K)/STLWC                                         CLD 0047
      RETURN                                                            CLD 0048
117   HAZECZ=DENSTY(66,K)/SCLWC                                         CLD 0049
      RETURN                                                            CLD 0050
118   HAZECZ=DENSTY(66,K)/SNLWC                                         CLD 0051
 15   RETURN                                                            CLD 0052
200   IF(IHA1 .GT. 0) GO TO 205                                         CLD 0053
      PRINT*,' WARNING ICLD NOT SET '                                   CLD 0054
      RETURN                                                            CLD 0055
205   CONTINUE                                                          CLD 0056
      WRH = RELHUM(K)                                                   CLD 0057
C                                                                       CLD 0058
      M = IC1                                                           CLD 0059
      IF(ICHIC1.EQ.6 .AND. M.NE.1)WRH=70.                               CLD 0060
C     THIS CODING  DOES NOT ALLOW TROP RH DEPENDENT  ABOVE EH(7,I)      CLD 0061
C     DEFAULTS TO TROPOSPHERIC AT 70. PERCENT                           CLD 0062
      DO 210 I=2,4                                                      CLD 0063
      IF (WRH.LT.RHZONE(I)) GO TO 215                                   CLD 0064
  210 CONTINUE                                                          CLD 0065
      I=4                                                               CLD 0066
  215 II=I-1                                                            CLD 0067
      IF(WRH.GT.0.0.AND.WRH.LT.99.)X=ALOG(100.0-WRH)                    CLD 0068
      X1=ALOG(100.0-RHZONE(II))                                         CLD 0069
      X2=ALOG(100.0-RHZONE(I))                                          CLD 0070
      IF (WRH.GE.99.0) X=X2                                             CLD 0071
      IF (WRH.LE.0.0) X=X1                                              CLD 0072
      ITA=ICHIC1                                                        CLD 0073
      IF(ITA.EQ.3. AND. M.EQ.1) GO TO 218                               CLD 0074
      IF(ITA.GT.6) GO TO 245                                            CLD 0075
CC                                                                      CLD 0076
CC    MICROWAVE                                                         CLD 0077
      N = 41                                                            CLD 0078
CC                                                                      CLD 0079
218   IF(N.GE.41. AND. ITA.EQ.3) ITA = 4                                CLD 0080
C     RH DEPENDENT AEROSOLS                                             CLD 0081
      GO TO (220,220,222,225,230,235), ITA                              CLD 0082
 220  E2=ALOG(ELWCR(I))                                                 CLD 0083
      E1=ALOG(ELWCR(II))                                                CLD 0084
      GO TO 240                                                         CLD 0085
 222  IF(M.GT.1) GO TO 225                                              CLD 0086
      E2=ALOG(ELWCM(I))                                                 CLD 0087
      E1=ALOG(ELWCM(II))                                                CLD 0088
      GO TO 240                                                         CLD 0089
  225 E2=ALOG(ELWCM(I))                                                 CLD 0090
      E1=ALOG(ELWCM(II))                                                CLD 0091
      GO TO 240                                                         CLD 0092
  230 E2=ALOG(ELWCU(I))                                                 CLD 0093
      E1=ALOG(ELWCU(II))                                                CLD 0094
      GO TO 240                                                         CLD 0095
  235 E2=ALOG(ELWCT(I))                                                 CLD 0096
      E1=ALOG(ELWCT(II))                                                CLD 0097
  240 EC=E1+(E2-E1)*(X-X1)/(X2-X1)                                      CLD 0098
      CON=EXP(EC)                                                       CLD 0099
      HAZECZ=DENSTY(66,K)/CON                                           CLD 0100
      RETURN                                                            CLD 0101
  245 IF (ITA.GT.19) GO TO 275                                          CLD 0102
      ITC=ICHIC1-7                                                      CLD 0103
      IF (ITC.LT.1) RETURN                                              CLD 0104
      GO TO (250,255,280,260,265,270,265,270,260,260,270,275), ITC      CLD 0105
250   CON=AFLWC                                                         CLD 0106
      GO TO 280                                                         CLD 0107
255   CON=RFLWC                                                         CLD 0108
      GO TO 280                                                         CLD 0109
260   CON=BSLWC                                                         CLD 0110
      GO TO 280                                                         CLD 0111
265   CON=AVLWC                                                         CLD 0112
      GO TO 280                                                         CLD 0113
270   CON=FVLWC                                                         CLD 0114
      GO TO 280                                                         CLD 0115
275   CON=MDLWC                                                         CLD 0116
280   CONTINUE                                                          CLD 0117
      HAZECZ=DENSTY(66,K)/CON                                           CLD 0118
      RETURN                                                            CLD 0119
      END                                                               CLD 0120
