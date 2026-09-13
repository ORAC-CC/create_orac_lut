      BLOCK DATA WVBNRG                                                 WVB 0001
C>    BLOCK DATA                                                        WVB 0002
C     WAVENUMBER-LOW AND WAVENUMBER-HIGH SPECIFY A BAND REGION          WVB 0003
C     FOR A MOLECULAR ABSORBER.                                         WVB 0004
C     THE UNIT FOR WAVENUMBER IS 1/CM.                                  WVB 0005
C     -999 IS AN INDICATOR TO INDICATE THE END OF ABSORPTION BANDS      WVB 0006
C     FOR ANY SPECIFIC ABSORBER.                                        WVB 0007
      COMMON /WNLOHI/                                                   WVB 0008
     L   IWLH2O(15),IWLO3 ( 6),IWLCO2(11),IWLCO ( 4),IWLCH4( 5),        WVB 0009
     L   IWLN2O(12),IWLO2 ( 7),IWLNH3( 3),IWLNO ( 2),IWLNO2( 4),        WVB 0010
     L   IWLSO2( 5),                                                    WVB 0011
     H   IWHH2O(15),IWHO3 ( 6),IWHCO2(11),IWHCO ( 4),IWHCH4( 5),        WVB 0012
     H   IWHN2O(12),IWHO2 ( 7),IWHNH3( 3),IWHNO ( 2),IWHNO2( 4),        WVB 0013
     H   IWHSO2( 5)                                                     WVB 0014
C                                                                       WVB 0015
      DATA IWLH2O/   0,    350,   1005,   1645,   2535,   3425,   4315, WVB 0016
     L    6155,   8005,   9620,  11545,  13075,  14865,  16340,   -999/ WVB 0017
      DATA IWHH2O/ 345,   1000,   1640,   2530,   3420,   4310,   6150, WVB 0018
     H    8000,   9615,  11540,  13070,  14860,  16045,  17860,   -999/ WVB 0019
C                                                                       WVB 0020
      DATA IWLO3 /   0,    515,   1630,   2670,   2850,   -999/         WVB 0021
      DATA IWHO3 / 200,   1275,   2295,   2845,   3260,   -999/         WVB 0022
C                                                                       WVB 0023
      DATA IWLCO2/ 425,    840,   1805,   3070,   3760,   4530,   5905, WVB 0024
     L    7395,   8030,   9340,   -999/                                 WVB 0025
      DATA IWHCO2/ 835,   1440,   2855,   3755,   4065,   5380,   7025, WVB 0026
     H    7785,   8335,   9670,   -999/                                 WVB 0027
C                                                                       WVB 0028
      DATA IWLCO /   0,   1940,   4040,   -999/                         WVB 0029
      DATA IWHCO / 175,   2285,   4370,   -999/                         WVB 0030
C                                                                       WVB 0031
      DATA IWLCH4/1065,   2345,   4110,   5865,   -999/                 WVB 0032
      DATA IWHCH4/1775,   3230,   4690,   6135,   -999/                 WVB 0033
C                                                                       WVB 0034
      DATA IWLN2O/   0,    490,    865,   1065,   1545,   2090,   2705, WVB 0035
     L    3245,   4260,   4540,   4910,   -999/                         WVB 0036
      DATA IWHN2O/ 120,    775,    995,   1385,   2040,   2655,   2865, WVB 0037
     H    3925,   4470,   4785,   5165,   -999/                         WVB 0038
C                                                                       WVB 0039
      DATA IWLO2 /   0,   7650,   9235,  12850,  14300,  15695,   -999/ WVB 0040
      DATA IWHO2 / 265,   8080,   9490,  13220,  14600,  15955,   -999/ WVB 0041
C                                                                       WVB 0042
      DATA IWLNH3/   0,    390,   -999/                                 WVB 0043
      DATA IWHNH3/ 385,   2150,   -999/                                 WVB 0044
C                                                                       WVB 0045
      DATA IWLNO /1700,   -999/                                         WVB 0046
      DATA IWHNO /2005,   -999/                                         WVB 0047
C                                                                       WVB 0048
      DATA IWLNO2/ 580,   1515,   2800,   -999/                         WVB 0049
      DATA IWHNO2/ 925,   1695,   2970,   -999/                         WVB 0050
C                                                                       WVB 0051
      DATA IWLSO2/   0,    400,    950,   2415,   -999/                 WVB 0052
      DATA IWHSO2/ 185,    650,   1460,   2580,   -999/                 WVB 0053
      END                                                               WVB 0054
