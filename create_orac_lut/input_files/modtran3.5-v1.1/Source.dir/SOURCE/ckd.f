      SUBROUTINE  CKD(SH2OT0,SH2OT1,TAVE,SH2O,FH2O,V)                   CKD 0001
C                                                                       CKD 0002
C                                                                       CKD 0003
C     REVISION:  3.3                                                    CKD 0004
C     CREATED:   28 APR 1994  10:14:28                                  CKD 0005
C     PRESENTLY: 28 APR 1994  10:15:53                                  CKD 0006
C     CLOUGH CKD_2.1                                                    CKD 0007
C     HVRCNT = '3.3'                                                    CKD 0008
C                                                                       CKD 0009
      REAL XFAC(0:50)                                                   CKD 0010
C                                                                       CKD 0011
      SAVE S260, S296, SFAC                                             CKD 0012
      DATA P0 / 1013. /,T0 / 296. /                                     CKD 0013
      DATA XLOSMT / 2.68675E+19 /                                       CKD 0014
C                                                                       CKD 0015
C     THESE ARE SELF-CONTINUUM MODIFICATION FACTORS FROM 700-1200 CM-1  CKD 0016
C                                                                       CKD 0017
C                                                                       CKD 0018
      DATA (XFAC(I),I=0,50)/                                            CKD 0019
     1    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,              CKD 0020
     2    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,              CKD 0021
     3    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,              CKD 0022
     4    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,              CKD 0023
     5    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,              CKD 0024
     6    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,              CKD 0025
     7    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,              CKD 0026
     8    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,              CKD 0027
     9    1.00000,1.00000,1.00000/                                      CKD 0028
C                                                                       CKD 0029
C                                                                       CKD 0030
      DATA S260/-12345678./, S296/-12345678./                           CKD 0031
C--------------------------------------------------------------------   CKD 0032
C                             SELF                                      CKD 0033
                                                                        CKD 0034
      ALPHA2 = 200.**2                                                  CKD 0035
C                                                                       CKD 0036
      ALPHS2= 120.**2                                                   CKD 0037
      BETAS = 5.E-06                                                    CKD 0038
      V0S=1310.                                                         CKD 0039
      FACTRS= 0.15                                                      CKD 0040
C                                                                       CKD 0041
C--------------------------------------------------------------------   CKD 0042
C                                                                       CKD 0043
C                             FOREIGN                                   CKD 0044
      HWSQF= 330.**2                                                    CKD 0045
      BETAF = 8.  E-11                                                  CKD 0046
      V0F =1130.                                                        CKD 0047
      FACTRF = 0.97                                                     CKD 0048
C                                                                       CKD 0049
      V0F2 =1900.                                                       CKD 0050
      HWSQF2 = 150.**2                                                  CKD 0051
      BETA2 = 3.E-06                                                    CKD 0052
C                                                                       CKD 0053
C--------------------------------------------------------------------   CKD 0054
         VJ = V                                                         CKD 0055
         VS2 = (VJ-V0S)**2                                              CKD 0056
         SH2O = 0.                                                      CKD 0057
         IF(SH2OT0.GT.0.)THEN                                           CKD 0058
         TFAC = (TAVE-T0)/(260.-T0)                                     CKD 0059
            SH2O = SH2OT0*(SH2OT1/SH2OT0)**TFAC                         CKD 0060
         SFAC = 1.                                                      CKD 0061
         IF (VJ.GE.700. .AND.  VJ.LE.1200.) THEN                        CKD 0062
            JFAC = (VJ-700.)/10. + 0.00001                              CKD 0063
            SFAC = XFAC(JFAC)                                           CKD 0064
         ENDIF                                                          CKD 0065
C                                                                       CKD 0066
C     CORRECTION TO SELF CONTINUUM (1 SEPT 85); FACTOR OF 0.78 AT 1000  CKD 0067
C                             AND  .......                              CKD 0068
C                                                                       CKD 0069
      SH2O = SFAC * SH2O*(1.-0.2333*(ALPHA2/((VJ-1050.)**2+ALPHA2))) *  CKD 0070
     C                 (1.-FACTRS*(ALPHS2/(VS2+(BETAS*VS2**2)+ALPHS2))) CKD 0071
         ENDIF                                                          CKD 0072
C-------------------------------------------                            CKD 0073
C                                                                       CKD 0074
C                                                                       CKD 0075
C     CORRECTION TO FOREIGN CONTINUUM                                   CKD 0076
C                                                                       CKD 0077
        VF2 = (VJ-V0F)**2                                               CKD 0078
        VF6 = VF2 * VF2 * VF2                                           CKD 0079
        FSCAL  = (1.-FACTRF*(HWSQF/(VF2+(BETAF*VF6)+HWSQF)))            CKD 0080
        VF2 = (VJ-V0F2)**2                                              CKD 0081
        VF4 = VF2*VF2                                                   CKD 0082
        FSCAL = FSCAL* (1.- 0.6*(HWSQF2/(VF2 + BETA2*VF4 + HWSQF2)))    CKD 0083
C                                                                       CKD 0084
        FH2O=FH2O*FSCAL                                                 CKD 0085
        RETURN                                                          CKD 0086
        END                                                             CKD 0087
