      SUBROUTINE BMLOAD                                                 BML 0001
C                                                                       BML 0002
C  THIS SUBROUTINE (CALLED BY BMOD) LOADS BAND MODEL DATA FOR A SINGLE  BML 0003
C  PARAMETER SET INTO THE MATRICES SD, OD AND ALF0 (FOUND IN /BMDCOM/). BML 0004
C                                                                       BML 0005
C  ORIGINALLY, THE BAND MODEL PARAMETERS ARE STORE IN COMPACT PARAMETER BML 0006
C  SETS.  THERE ARE TWO TYPES OF PARAMETER SETS CORRESPONDING TO LINE   BML 0007
C  CENTER AND LINE TAIL CONTRIBUTIONS.  IF IMOL IS 12 OR LESS, THE      BML 0008
C  PARAMETER SET CONTAINS LINE CENTER MOLECULAR ABSORPTION INFO:        BML 0009
C      IBIN = THE BIN NUMBER (IBIN*IBNDWD = BIN CENTER FREQUENCY)       BML 0010
C      IMOL =  1 FOR H2O                                                BML 0011
C           =  2 FOR CO2                                                BML 0012
C           =  3 FOR O3                                                 BML 0013
C           =  4 FOR N2O                                                BML 0014
C           =  5 FOR CO                                                 BML 0015
C           =  6 FOR CH4                                                BML 0016
C           =  7 FOR O2                                                 BML 0017
C           =  8 FOR NO                                                 BML 0018
C           =  9 FOR SO2                                                BML 0019
C           = 10 FOR NO2                                                BML 0020
C           = 11 FOR NH3                                                BML 0021
C           = 12 FOR HNO3                                               BML 0022
C    SDZ(IT)= AVERAGE MOLECULAR ABSORPTION COEFFICIENT PARAMETER        BML 0023
C                 FOR THE IT'TH TEMPERATURE (CM-1/AMAGAT)               BML 0024
C      IALF = LORENTZ LINE WIDTH PARAMETER AT STANDARD PRESSURE AND     BML 0025
C                 TEMPERATURE (1.E-04 CM-1/ATM)                         BML 0026
C    ODZ(IT)= LINE DENSITY PARAMETER FOR THE TEMPERATURE IT (CM-1)      BML 0027
C                                                                       BML 0028
C  IF IMOL IS 12 OR MORE, THE PARAMETER SET CONTAINS CONTINUUM          BML 0029
C  ABSORPTION INFORMATION AND                                           BML 0030
C              IBIN = THE BIN NUMBER                                    BML 0031
C     IMOL,IALF = 13 FOR  H2O TAIL CONTRIBUTION                         BML 0032
C                       = 14 FOR  CO2 TAIL CONTRIBUTION                 BML 0033
C                       = 15 FOR   O3 TAIL CONTRIBUTION                 BML 0034
C                       = 16 FOR  N2O TAIL CONTRIBUTION                 BML 0035
C                       = 17 FOR   CO TAIL CONTRIBUTION                 BML 0036
C                       = 18 FOR  CH4 TAIL CONTRIBUTION                 BML 0037
C                       = 19 FOR   O2 TAIL CONTRIBUTION                 BML 0038
C                       = 20 FOR   NO TAIL CONTRIBUTION                 BML 0039
C                       = 21 FOR  SO2 TAIL CONTRIBUTION                 BML 0040
C                       = 22 FOR  NO2 TAIL CONTRIBUTION                 BML 0041
C                       = 23 FOR  NH3 TAIL CONTRIBUTION                 BML 0042
C                       = 24 FOR HNO3 TAIL CONTRIBUTION                 BML 0043
C            SDZ(IT) = MOLECULAR LINE WING CONTINUUM ABSORPTION         BML 0044
C                         COEFFICIENT FOR SPECIES IMOL AT THE           BML 0045
C                         IT'TH TEMPERATURE (CM-1/AMAGAT)               BML 0046
C            ODZ(IT) = MOLECULAR LINE WING CONTINUUM ABSORPTION         BML 0047
C                         COEFFICIENT FOR SPECIES IALF AT THE           BML 0048
C                         IT'TH TEMPERATURE (CM-1/AMAGAT)               BML 0049
C                                                                       BML 0050
C                                                                       BML 0051
C                                                                       BML 0052
C                                                                       BML 0053
C     CONVENTION                                                        BML 0054
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         BML 0055
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            BML 0056
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             BML 0057
C                                                                       BML 0058
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              BML 0059
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       BML 0060
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     BML 0061
C                                                                       BML 0062
      INCLUDE 'PARAM.LST'                                               BML 0063
C                                                                       BML 0064
C     LIST COMMONS                                                      BML 0065
      INTEGER IBNDWD,IP,ITB,NTEMP,IBIN,IMOL,IALF,JJ,JJS                 BML 0066
      REAL TBAND,SDZ,ODZ,SD,OD,ALF0,T5,PTM75,FF,T5S,                    BML 0067
     1  PTM75S,FFS,DOPFAC,DOP0,DOPSUM,COLSUM,ODSUM,SDSUM                BML 0068
      COMMON/BMDCOM/IBNDWD,IP,ITB,NTEMP,TBAND(MXTEMP),SDZ(MXTEMP),      BML 0069
     1  ODZ(MXTEMP),IBIN,IMOL,IALF,SD(MXTEMP,MMOLT2),OD(MXTEMP,MMOLT),  BML 0070
     2  ALF0(MMOLT),T5(LAYTHR),PTM75(LAYTHR),JJ(LAYTHR),FF(LAYTHR),     BML 0071
     3  T5S(LAYTHR,NSPC),PTM75S(LAYTHR,NSPC),JJS(LAYTHR,MMOLT),         BML 0072
     4  FFS(LAYTHR,MMOLT),DOPFAC(MMOLT),DOP0(MMOLT),SDSUM(MMOLT2),      BML 0073
     5  ODSUM(MMOLT),DOPSUM(MMOLT),COLSUM(MMOLT)                        BML 0074
C                                                                       BML 0075
C  FILL THE SD, OD AND ALF0 MATRICES.                                   BML 0076
C                                                                       BML 0077
      IM=IMOL                                                           BML 0078
      IF(IM.LE.NSPC)THEN                                                BML 0079
C                                                                       BML 0080
C         LINE CENTER CONTRIBUTIONS                                     BML 0081
C                                                                       BML 0082
          IMTAIL=IM+NSPC                                                BML 0083
          DO 10 IT=1,NTEMP                                              BML 0084
              SD(IT,IM)=SDZ(IT)                                         BML 0085
              OD(IT,IM)=ODZ(IT)                                         BML 0086
              SD(IT,IMTAIL)=0.                                          BML 0087
   10     CONTINUE                                                      BML 0088
          ALF0(IM)=1.E-04*IALF                                          BML 0089
      ELSE                                                              BML 0090
C                                                                       BML 0091
C         LINE TAIL CONTRIBUTIONS                                       BML 0092
C                                                                       BML 0093
          JM=IALF                                                       BML 0094
          IF(JM.EQ.0)THEN                                               BML 0095
C                                                                       BML 0096
C             ONE TAIL                                                  BML 0097
C                                                                       BML 0098
              DO 20 IT=1,NTEMP                                          BML 0099
                  SD(IT,IM)=SDZ(IT)                                     BML 0100
   20         CONTINUE                                                  BML 0101
          ELSE                                                          BML 0102
C                                                                       BML 0103
C             TWO TAILS                                                 BML 0104
C                                                                       BML 0105
              DO 30 IT=1,NTEMP                                          BML 0106
                  SD(IT,JM)=ODZ(IT)                                     BML 0107
                  SD(IT,IM)=SDZ(IT)                                     BML 0108
   30         CONTINUE                                                  BML 0109
          ENDIF                                                         BML 0110
      ENDIF                                                             BML 0111
      RETURN                                                            BML 0112
      END                                                               BML 0113
