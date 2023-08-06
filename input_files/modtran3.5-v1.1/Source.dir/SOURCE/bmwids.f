      BLOCK DATA BMWIDS                                                 BMW 0001
C                                                                       BMW 0002
C     THIS BLOCK DATA DEFINED DOPPLER HALF-WIDTH                        BMW 0003
C     CONSTANTS FOR THE MOLECULAR SPECIES.                              BMW 0004
C                                                                       BMW 0005
C     INCLUDE PARAMETERS                                                BMW 0006
      INCLUDE 'PARAM.LST'                                               BMW 0007
C                                                                       BMW 0008
C     LIST COMMONS                                                      BMW 0009
      INTEGER IBNDWD,IP,ITB,NTEMP,IBIN,IMOL,IALF,JJ,JJS                 BMW 0010
      REAL TBAND,SDZ,ODZ,SD,OD,ALF0,T5,PTM75,FF,T5S,                    BMW 0011
     1  PTM75S,FFS,DOPFAC,DOP0,DOPSUM,COLSUM,ODSUM,SDSUM                BMW 0012
      COMMON/BMDCOM/IBNDWD,IP,ITB,NTEMP,TBAND(MXTEMP),SDZ(MXTEMP),      BMW 0013
     1  ODZ(MXTEMP),IBIN,IMOL,IALF,SD(MXTEMP,MMOLT2),OD(MXTEMP,MMOLT),  BMW 0014
     2  ALF0(MMOLT),T5(LAYTHR),PTM75(LAYTHR),JJ(LAYTHR),FF(LAYTHR),     BMW 0015
     3  T5S(LAYTHR,NSPC),PTM75S(LAYTHR,NSPC),JJS(LAYTHR,MMOLT),         BMW 0016
     4  FFS(LAYTHR,MMOLT),DOPFAC(MMOLT),DOP0(MMOLT),SDSUM(MMOLT2),      BMW 0017
     5  ODSUM(MMOLT),DOPSUM(MMOLT),COLSUM(MMOLT)                        BMW 0018
C                                                                       BMW 0019
C     LIST DATA                                                         BMW 0020
      DATA (DOPFAC(I),I=1,NSPECT)/                                      BMW 0021
     1  1.3945E-6,  0.8922E-6,  0.8543E-6,  0.8921E-6,  1.1183E-6,      BMW 0022
     2  1.4777E-6,  1.0463E-6,  1.0805E-6,  0.7395E-6,  0.8726E-6,      BMW 0023
     3  1.4342E-6,  0.7456E-6,  5.04981E-7, 5.38245E-7, 5.79088E-7,     BMW 0024
     4  6.30907E-7, 6.36486E-7, 4.32375E-7, 4.52709E-7, 4.76212E-7,     BMW 0025
     5  5.99528E-7, 6.65841E-7, 5.83391E-7, 4.7721E-7,  5.69488E-7/     BMW 0026
      END                                                               BMW 0027
