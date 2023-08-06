      SUBROUTINE AERPRF (I,  VIS,HAZE,IHAZE,     ISEASN,IVULCN,N)       AER 0001
C***********************************************************************AER 0002
C     WILL COMPUTE DENSITY    PROFILES FOR AEROSOLS                     AER 0003
C***********************************************************************AER 0004
      COMMON/PRFD  /ZHT(34),HZ2K(34,5),FAWI50(34),                      AER 0005
     1FAWI23(34),SPSU50(34),SPSU23(34),BASTFW(34),                      AER 0006
     2VUMOFW(34),HIVUFW(34),EXVUFW(34),BASTSS(34),                      AER 0007
     3VUMOSS(34),HIVUSS(34),EXVUSS(34),UPNATM(34),                      AER 0008
     3VUTONO(34),VUTOEX(34),EXUPAT(34)                                  AER 0009
      DIMENSION VS(5)                                                   AER 0010
      DATA VS/50.,23.,10.,5.,2./                                        AER 0011
      HAZE=0.                                                           AER 0012
      N=7                                                               AER 0013
      IF (IHAZE.EQ.0) RETURN                                            AER 0014
      IF (ZHT(I).GT.2.0) GO TO 15                                       AER 0015
      DO 5 J=2,5                                                        AER 0016
      IF (VIS.GE.VS(J)) GO TO 10                                        AER 0017
    5 CONTINUE                                                          AER 0018
      J=5                                                               AER 0019
   10 CONST=1./(1./VS(J)-1./VS(J-1))                                    AER 0020
      HAZE=CONST*( (HZ2K(I,J)-HZ2K(I,J-1))/VIS +                        AER 0021
     1 HZ2K(I,J-1)/VS(J) - HZ2K(I,J )/VS(J-1) )                         AER 0022
      IF(ZHT(I).GT.2.0) GO TO 15                                        AER 0023
      RETURN                                                            AER 0024
   15 IF (ZHT(I).GT.10.) GO TO 35                                       AER 0025
      CONST=1./(1./23.-1./50.)                                          AER 0026
      IF (ISEASN.GT.1) GO TO 25                                         AER 0027
      IF (VIS.LE.23.) HAZI=SPSU23(I)                                    AER 0028
      IF (VIS.LE.23.) GO TO 200                                         AER 0029
      IF (ZHT(I).GT.4.0) GO TO 20                                       AER 0030
      HAZI=CONST*((SPSU23(I)-SPSU50(I))/VIS+SPSU50(I)/23.-SPSU23(I)/50.)AER 0031
      GO TO 200                                                         AER 0032
   20 HAZI=SPSU50(I)                                                    AER 0033
      GO TO 200                                                         AER 0034
   25 IF (VIS.LE.23.) HAZI=FAWI23(I)                                    AER 0035
      IF (VIS.LE.23.) GO TO 200                                         AER 0036
      IF (ZHT(I).GT.4.0) GO TO 30                                       AER 0037
      HAZI=CONST*((FAWI23(I)-FAWI50(I))/VIS+FAWI50(I)/23.-FAWI23(I)/50.)AER 0038
      GO TO 200                                                         AER 0039
   30 HAZI=FAWI50(I)                                                    AER 0040
      GO TO 200                                                         AER 0041
   35 IF (ZHT(I).GT.30.0) GO TO 75                                      AER 0042
      HAZI=BASTSS(I)                                                    AER 0043
      IF (ISEASN.GT.1) GO TO 55                                         AER 0044
      IF (IVULCN.EQ.0) HAZI=BASTSS(I)                                   AER 0045
      IF (IVULCN.EQ.0) GO TO 200                                        AER 0046
      GO TO (40,45,50,50,45,45,50,52), IVULCN                           AER 0047
   40 HAZI=BASTSS(I)                                                    AER 0048
      GO TO 200                                                         AER 0049
   45 HAZI=VUMOSS(I)                                                    AER 0050
      GO TO 200                                                         AER 0051
   50 HAZI=HIVUSS(I)                                                    AER 0052
      GO TO 200                                                         AER 0053
52    HAZI=EXVUSS(I)                                                    AER 0054
      GO TO 200                                                         AER 0055
   55 IF (IVULCN.EQ.0) HAZI=BASTFW(I)                                   AER 0056
      IF (IVULCN.EQ.0) GO TO 200                                        AER 0057
      GO TO (60,65,70,70,65,65,70,72), IVULCN                           AER 0058
   60 HAZI=BASTFW(I)                                                    AER 0059
      GO TO 200                                                         AER 0060
   65 HAZI=VUMOFW(I)                                                    AER 0061
      GO TO 200                                                         AER 0062
   70 HAZI=HIVUFW(I)                                                    AER 0063
      GO TO 200                                                         AER 0064
72    HAZI=EXVUFW(I)                                                    AER 0065
      GO TO 200                                                         AER 0066
   75 N=14                                                              AER 0067
      IF (IVULCN.GT.1) GO TO 80                                         AER 0068
      HAZI=UPNATM(I)                                                    AER 0069
      GO TO 200                                                         AER 0070
   80 HAZI=VUTONO(I)                                                    AER 0071
200   IF(HAZI.GT.0) HAZE=HAZI                                           AER 0072
      RETURN                                                            AER 0073
      END                                                               AER 0074
