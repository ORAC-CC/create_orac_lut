      SUBROUTINE LAYCLD(K,EQLWCZ,RRATZ,ICLD1,GNDALT)                    LCL 0001
C                                                                       LCL 0002
C     ROUTINE TO DETERMINE CLOUD DENSITY AND RAIN RATE FOR LAYER K      LCL 0003
C                                                                       LCL 0004
C     ZM      COMMON /MODEL/ FINAL ALTITUDES FOR LOWTRAN                LCL 0005
C     ZK      EFFECTIVE CLOUD ALTITUDES                                 LCL 0006
C     ZCLD    CLOUD ALTITUDE ARRAY                                      LCL 0007
      INCLUDE 'PARAM.LST'                                               LCL 0008
      REAL ZCLD,CLD,CLDICE,RR                                           LCL 0009
      COMMON/CLDRR/ZCLD(1:NZCLD,0:1),CLD(1:NZCLD,0:5),                  LCL 0010
     1  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)                             LCL 0011
      REAL ZM,PM,TM,RFNDX,DENSTY                                        LCL 0012
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    LCL 0013
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               LCL 0014
      IF(ICLD1.LE.0 .OR. ICLD1.GT.11)RETURN                             LCL 0015
      ZK=ZM(K)-GNDALT                                                   LCL 0016
      IF(ZK.LT.0.)ZK=0.                                                 LCL 0017
      IF(ICLD1.LE.5)THEN                                                LCL 0018
C                                                                       LCL 0019
C         ICLD1 IS 1- 5 ONE OF 5 SPECIFIC CLOUD MODELS IS CHOSEN        LCL 0020
          MC=ICLD1                                                      LCL 0021
          MR=6                                                          LCL 0022
      ELSE                                                              LCL 0023
C                                                                       LCL 0024
C         ICLD1 IS 6-11 ONE OF 5 SPECIFIC CLOUD/RAIN MODELS CHOSEN      LCL 0025
          MC=1                                                          LCL 0026
          IF(ICLD1.EQ.6)MC=3                                            LCL 0027
          IF(ICLD1.EQ.7 .OR. ICLD1.EQ.8)MC=5                            LCL 0028
          MR=ICLD1-5                                                    LCL 0029
      ENDIF                                                             LCL 0030
      EQLWCZ=0.                                                         LCL 0031
      RRATZ=0.                                                          LCL 0032
      IF(ZK.GE.ZCLD(1,1))THEN                                           LCL 0033
          MKM1=1                                                        LCL 0034
          DO 10 MK=2,16                                                 LCL 0035
              IF(ZK.LE.ZCLD(MK,1))THEN                                  LCL 0036
                  FAC=(ZCLD(MK,1)-ZK)/(ZCLD(MK,1)-ZCLD(MKM1,1))         LCL 0037
                  EQLWCZ=CLD(MK,MC)+FAC*(CLD(MKM1,MC)-CLD(MK,MC))       LCL 0038
                  IF(MR.LE.5)RRATZ=RR(MK,MR)+FAC*(RR(MKM1,MR)-RR(MK,MR))LCL 0039
                  RETURN                                                LCL 0040
              ENDIF                                                     LCL 0041
   10     MKM1=MK                                                       LCL 0042
      ENDIF                                                             LCL 0043
      RETURN                                                            LCL 0044
      END                                                               LCL 0045
