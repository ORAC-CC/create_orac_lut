      SUBROUTINE LAYVSA(K,RH,AHAZE,IHA1,ZNEW)                           LVS 0001
C                                                                       LVS 0002
C     RETURNS HAZE FOR VSA OPTION                                       LVS 0003
C                                                                       LVS 0004
      INCLUDE 'PARAM.LST'                                               LVS 0005
      REAL ZM,PM,TM,RFNDX,DENSTY                                        LVS 0006
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    LVS 0007
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               LVS 0008
      COMMON /CARD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB  LVS 0009
     1  ,MODTRN                                                         LVS 0010
      LOGICAL MODTRN                                                    LVS 0011
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)            LVS 0012
C                                                                       LVS 0013
      DIMENSION ZNEW(LAYDIM)                                            LVS 0014
      RH=0.                                                             LVS 0015
      AHAZE=0                                                           LVS 0016
      IHA1=0                                                            LVS 0017
C     HMXVSA=ZVSA(9)                                                    LVS 0018
      IF(MODEL.EQ.0.OR.MODEL.EQ.7) RETURN                               LVS 0019
      IF(K.GT.9) RETURN                                                 LVS 0020
      ZM(K)=ZVSA(K)                                                     LVS 0021
      RH=RHVSA(K)                                                       LVS 0022
      AHAZE=AHVSA(K)                                                    LVS 0023
      IHA1=IHVSA(K)                                                     LVS 0024
      RETURN                                                            LVS 0025
      END                                                               LVS 0026
