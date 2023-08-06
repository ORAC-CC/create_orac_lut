      BLOCK DATA PRFDTA                                                 PRF 0001
C>    BLOCK DATA                                                        PRF 0002
C                                                                       PRF 0003
C        AEROSOL PROFILE DATA                                           PRF 0004
C                                                                       PRF 0005
CCC         0-2KM                                                       PRF 0006
CCC           HZ2K=5 VIS PROFILES- 50KM,23KM,10KM,5KM,2KM               PRF 0007
CCC         >2-10KM                                                     PRF 0008
CCC           FAWI50=FALL/WINTER   50KM VIS                             PRF 0009
CCC           FAWI23=FALL/WINTER    23KM VIS                            PRF 0010
CCC           SPSU50=SPRING/SUMMER  50KM VIS                            PRF 0011
CCC           SPSU23=SPRING/SUMMER  23KM VIS                            PRF 0012
CCC         >10-30KM                                                    PRF 0013
CCC           BASTFW=BACKGROUND STRATOSPHERIC   FALL/WINTER             PRF 0014
CCC           VUMOFW=MODERATE VOLCANIC          FALL/WINTER             PRF 0015
CCC           HIVUFW=HIGH VOLCANIC              FALL/WINTER             PRF 0016
CCC           EXVUFW=EXTREME VOLCANIC           FALL/WINTER             PRF 0017
CCC           BASTSS,VUMOSS,HIVUSS,EXVUSS=      SPRING/SUMMER           PRF 0018
CCC         >30-100KM                                                   PRF 0019
CCC           UPNATM=NORMAL UPPER ATMOSPHERIC                           PRF 0020
CCC           VUTONO=TRANSITION FROM VOLCANIC TO NORMAL                 PRF 0021
CCC           VUTOEX=TRANSITION FROM VOLCANIC TO EXTREME                PRF 0022
CCC           EXUPAT=EXTREME UPPER ATMOSPHERIC                          PRF 0023
      COMMON/PRFD  /ZHT(34),HZ2K(34,5),FAWI50(34),FAWI23(34),SPSU50(34),PRF 0024
     1SPSU23(34),BASTFW(34),VUMOFW(34),HIVUFW(34),EXVUFW(34),BASTSS(34),PRF 0025
     2VUMOSS(34),HIVUSS(34),EXVUSS(34),UPNATM(34),VUTONO(34),           PRF 0026
     3VUTOEX(34),EXUPAT(34)                                             PRF 0027
      DATA ZHT/                                                         PRF 0028
     *    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   PRF 0029
     *    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   PRF 0030
     *   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   PRF 0031
     *   35.,   40.,   45.,   50.,   70.,  100.,99999./                 PRF 0032
       DATA HZ2K(1,1),HZ2K(1,2),HZ2K(1,3),HZ2K(1,4),HZ2K(1,5)/          PRF 0033
     1 6.62E-02, 1.58E-01, 3.79E-01, 7.70E-01, 1.94E+00/                PRF 0034
       DATA HZ2K(2,1),HZ2K(2,2),HZ2K(2,3),HZ2K(2,4),HZ2K(2,5)/          PRF 0035
     1 4.15E-02, 9.91E-02, 3.79E-01, 7.70E-01, 1.94E+00/                PRF 0036
       DATA HZ2K(3,1),HZ2K(3,2),HZ2K(3,3),HZ2K(3,4),HZ2K(3,5)/          PRF 0037
     1 2.60E-02, 6.21E-02, 6.21E-02, 6.21E-02, 6.21E-02/                PRF 0038
      DATA FAWI50  /3*0.,                                               PRF 0039
     1 1.14E-02, 6.43E-03, 4.85E-03, 3.54E-03, 2.31E-03, 1.41E-03,      PRF 0040
     2 9.80E-04,7.87E-04,23*0./                                         PRF 0041
      DATA FAWI23              /3*0.,                                   PRF 0042
     1 2.72E-02, 1.20E-02, 4.85E-03, 3.54E-03, 2.31E-03, 1.41E-03,      PRF 0043
     2 9.80E-04,7.87E-04,23*0./                                         PRF 0044
      DATA  SPSU50              / 3*0.,                                 PRF 0045
     1 1.46E-02, 1.02E-02, 9.31E-03, 7.71E-03, 6.23E-03, 3.37E-03,      PRF 0046
     2 1.82E-03  ,1.14E-03,23*0./                                       PRF 0047
      DATA  SPSU23              / 3*0.,                                 PRF 0048
     1 3.46E-02, 1.85E-02, 9.31E-03, 7.71E-03, 6.23E-03, 3.37E-03,      PRF 0049
     2 1.82E-03  ,1.14E-03,23*0./                                       PRF 0050
      DATA BASTFW       /11*0.,                                         PRF 0051
     1           7.14E-04, 6.64E-04, 6.23E-04, 6.45E-04, 6.43E-04,      PRF 0052
     2 6.41E-04, 6.00E-04, 5.62E-04, 4.91E-04, 4.23E-04, 3.52E-04,      PRF 0053
     3 2.95E-04, 2.42E-04, 1.90E-04, 1.50E-04, 3.32E-05 ,7*0./          PRF 0054
      DATA    VUMOFW       /11*0.,                                      PRF 0055
     1           1.79E-03, 2.21E-03, 2.75E-03, 2.89E-03, 2.92E-03,      PRF 0056
     2 2.73E-03, 2.46E-03, 2.10E-03, 1.71E-03, 1.35E-03, 1.09E-03,      PRF 0057
     3 8.60E-04, 6.60E-04, 5.15E-04, 4.09E-04, 7.60E-05 ,7*0./          PRF 0058
      DATA    HIVUFW       /11*0.,                                      PRF 0059
     1           2.31E-03, 3.25E-03, 4.52E-03, 6.40E-03, 7.81E-03,      PRF 0060
     2 9.42E-03, 1.07E-02, 1.10E-02, 8.60E-03, 5.10E-03, 2.70E-03,      PRF 0061
     3 1.46E-03, 8.90E-04, 5.80E-04, 4.09E-04, 7.60E-05 ,7*0./          PRF 0062
      DATA    EXVUFW       /11*0.,                                      PRF 0063
     1           2.31E-03, 3.25E-03, 4.52E-03, 6.40E-03, 1.01E-02,      PRF 0064
     2 2.35E-02, 6.10E-02, 1.00E-01, 4.00E-02, 9.15E-03, 3.13E-03,      PRF 0065
     3 1.46E-03, 8.90E-04, 5.80E-04, 4.09E-04, 7.60E-05 ,7*0./          PRF 0066
      DATA    BASTSS       /11*0.,                                      PRF 0067
     1           7.99E-04, 6.41E-04, 5.17E-04, 4.42E-04, 3.95E-04,      PRF 0068
     2 3.82E-04, 4.25E-04, 5.20E-04, 5.81E-04, 5.89E-04, 5.02E-04,      PRF 0069
     3 4.20E-04, 3.00E-04, 1.98E-04, 1.31E-04, 3.32E-05 ,7*0./          PRF 0070
      DATA    VUMOSS       /11*0.,                                      PRF 0071
     1           2.12E-03, 2.45E-03, 2.80E-03, 2.89E-03, 2.92E-03,      PRF 0072
     2 2.73E-03, 2.46E-03, 2.10E-03, 1.71E-03, 1.35E-03, 1.09E-03,      PRF 0073
     3 8.60E-04, 6.60E-04, 5.15E-04, 4.09E-04, 7.60E-05 ,7*0./          PRF 0074
      DATA    HIVUSS       /11*0.,                                      PRF 0075
     1           2.12E-03, 2.45E-03, 2.80E-03, 3.60E-03, 5.23E-03,      PRF 0076
     2 8.11E-03, 1.20E-02, 1.52E-02, 1.53E-02, 1.17E-02, 7.09E-03,      PRF 0077
     3 4.50E-03, 2.40E-03, 1.28E-03, 7.76E-04, 7.60E-05 ,7*0./          PRF 0078
      DATA    EXVUSS       /11*0.,                                      PRF 0079
     1           2.12E-03, 2.45E-03, 2.80E-03, 3.60E-03, 5.23E-03,      PRF 0080
     2 8.11E-03, 1.27E-02, 2.32E-02, 4.85E-02, 1.00E-01, 5.50E-02,      PRF 0081
     3 6.10E-03, 2.40E-03, 1.28E-03, 7.76E-04, 7.60E-05 ,7*0./          PRF 0082
      DATA UPNATM       /26*0.,                                         PRF 0083
     1 3.32E-05, 1.64E-05, 7.99E-06, 4.01E-06, 2.10E-06, 1.60E-07,      PRF 0084
     2 9.31E-10, 0.      /                                              PRF 0085
      DATA VUTONO       /26*0.,                                         PRF 0086
     1 7.60E-05, 2.45E-05, 7.99E-06, 4.01E-06, 2.10E-06, 1.60E-07,      PRF 0087
     2 9.31E-10, 0.      /                                              PRF 0088
      DATA VUTOEX       /26*0.,                                         PRF 0089
     1 7.60E-05, 7.20E-05, 6.95E-05, 6.60E-05, 5.04E-05, 1.03E-05,      PRF 0090
     2 4.50E-07, 0.      /                                              PRF 0091
      DATA EXUPAT       /26*0.,                                         PRF 0092
     1 3.32E-05, 4.25E-05, 5.59E-05, 6.60E-05, 5.04E-05, 1.03E-05,      PRF 0093
     2 4.50E-07, 0.      /                                              PRF 0094
      END                                                               PRF 0095
