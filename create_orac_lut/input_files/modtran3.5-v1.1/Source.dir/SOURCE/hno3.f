      SUBROUTINE HNO3 (V,HABS)                                          NO3 0001
C                                                                       NO3 0002
C     HNO3  STATISTICAL BAND PARAMETERS                                 NO3 0003
C                                                                       NO3 0004
      DIMENSION H1(15), H2(16), H3(13)                                  NO3 0005
C     ARRAY H1 CONTAINS HNO3 ABS, COEF(CM-1ATM-1) FROM  850 TO 920 CM-1 NO3 0006
      DATA H1/2.197,3.911,6.154,8.150,9.217,9.461,11.56,11.10,11.17,12.4NO3 0007
     10,10.49,7.509,6.136,4.899,2.866/                                  NO3 0008
C     ARRAY H2 CONTAINS HNO3 ABS, COEF(CM-1ATM-1) FROM 1275 TO1350 CM-1 NO3 0009
      DATA H2/2.828,4.611,6.755,8.759,10.51,13.74,18.00,21.51,23.09,21.6NO3 0010
     18,21.32,16.82,16.42,17.87,14.86,8.716/                            NO3 0011
C     ARRAY H3 CONTAINS HNO3 ABS, COEF(CM-1ATM-1) FROM 1675 TO1735 CM-1 NO3 0012
      DATA H3/5.003,8.803,14.12,19.83,23.31,23.58,23.22,21.09,26.99,25.8NO3 0013
     14,24.79,17.68,9.420/                                              NO3 0014
      HABS=0.                                                           NO3 0015
      IF (V.GE.850.0.AND.V.LE.920.0) GO TO 5                            NO3 0016
      IF (V.GE.1275.0.AND.V.LE.1350.0) GO TO 10                         NO3 0017
      IF (V.GE.1675.0.AND.V.LE.1735.0) GO TO 15                         NO3 0018
      RETURN                                                            NO3 0019
    5 I=(V-845.)/5.                                                     NO3 0020
      HABS=H1(I)                                                        NO3 0021
      RETURN                                                            NO3 0022
   10 I=(V-1270.)/5.                                                    NO3 0023
      HABS=H2(I)                                                        NO3 0024
      RETURN                                                            NO3 0025
   15 I=(V-1670.)/5.                                                    NO3 0026
      HABS=H3(I)                                                        NO3 0027
      RETURN                                                            NO3 0028
      END                                                               NO3 0029
