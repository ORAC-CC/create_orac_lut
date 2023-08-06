      BLOCK DATA MDTA                                                   MDT 0001
C                                                                       MDT 0002
C     CLOUD AND RAIN DATA                                               MDT 0003
      INCLUDE 'PARAM.LST'                                               MDT 0004
      REAL ZCLD,CLD,CLDICE,RR                                           MDT 0005
      COMMON/CLDRR/ZCLD(1:NZCLD,0:1),CLD(1:NZCLD,0:5),                  MDT 0006
     1  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)                             MDT 0007
C                                                                       MDT 0008
C     PROFILE ALTITUDES [KM]                                            MDT 0009
      DATA ZCLD/NZCLD*0.,                                               MDT 0010
     1    .00,  .16,  .33,  .66, 1.00, 1.50, 2.00,2.4,2.7,3.00,3.5,     MDT 0011
     2   4.00, 4.50, 5.00, 5.50, 6.00/                                  MDT 0012
C                                                                       MDT 0013
C     WATER DROPLET DENSITIES [GM/M3]                                   MDT 0014
      DATA CLD/NZCLD*0.,                                                MDT 0015
     1    .00,  .00,  .00,  .20,  .35, 1.00, 1.00,1.0, .3, .15, .0,5*.0,MDT 0016
     2    .00,  .00,  .00,  .00,  .00,  .00,  .00, .3, .4, .30, .0,5*.0,MDT 0017
     3    .00,  .00,  .15,  .30,  .15,  .00,  .00, .0, .0, .00, .0,5*.0,MDT 0018
     4    .00,  .00,  .00,  .10,  .15,  .15,  .10, .0, .0, .00, .0,5*.0,MDT 0019
     5    .00,  .30,  .65,  .40,  .00,  .00,  .00, .0, .0, .00, .0,5*.0/MDT 0020
C                                                                       MDT 0021
C     ICE PARTICLE DENSITIES [GM/M3]                                    MDT 0022
      DATA CLDICE/NZCLD*0.,NZCLD*0./                                    MDT 0023
C                                                                       MDT 0024
C     RAIN RATES [MM/HR]                                                MDT 0025
      DATA RR/NZCLD*0.,                                                 MDT 0026
     1   2.00, 1.78, 1.43, 1.22,  .86,  .22,  .00, .0, .0, .00, .0,5*.0,MDT 0027
     2   5.00, 4.00, 3.40, 2.60,  .80,  .20,  .00, .0, .0, .00, .0,5*.0,MDT 0028
     3  12.50,10.50, 8.00, 6.00, 2.50,  .80,  .20, .0, .0, .00, .0,5*.0,MDT 0029
     4  25.00,21.50,17.50,12.00, 7.50, 4.20, 2.50,1.0, .7, .20, .0,5*.0,MDT 0030
     5  75.00,70.00,65.00,60.00,45.00,20.00,12.50,7.0,3.5,1.00, .2,5*.0/MDT 0031
      END                                                               MDT 0032
