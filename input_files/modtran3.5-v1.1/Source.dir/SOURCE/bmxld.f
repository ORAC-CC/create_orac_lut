      SUBROUTINE BMXLD                                                  BMX 0001
C                                                                       BMX 0002
C     THIS ROUTINE IS THE EQUIVALENT OF BMLOAD FOR THE X-MOLECULES.     BMX 0003
C                                                                       BMX 0004
C     THIS SUBROUTINE (CALLED BY BMOD) LOADS BAND MODEL DATA FOR A SINGLBMX 0005
C     PARAMETER SET INTO THE MATRICES SDX, ODX AND ALF0X (FOUND IN /BMDCBMX 0006
C                                                                       BMX 0007
C     CONVENTION                                                        BMX 0008
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         BMX 0009
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            BMX 0010
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             BMX 0011
C                                                                       BMX 0012
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              BMX 0013
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       BMX 0014
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     BMX 0015
C                                                                       BMX 0016
      INCLUDE 'PARAM.LST'                                               BMX 0017
C                                                                       BMX 0018
C     TRANS VARIABLES                                                   BMX 0019
      INTEGER IBINX,IMOLX,IALFX                                         BMX 0020
      REAL SDZX,ODZX                                                    BMX 0021
      COMMON/BMDCMX/SDZX(MXTEMP),ODZX(MXTEMP),IBINX,IMOLX,IALFX         BMX 0022
      INTEGER IBNDWD,IP,ITB,NTEMP,IBIN,IMOL,IALF,JJ,JJS                 BMX 0023
      REAL TBAND,SDZ,ODZ,SD,OD,ALF0,T5,PTM75,FF,T5S,                    BMX 0024
     1  PTM75S,FFS,DOPFAC,DOP0,DOPSUM,COLSUM,ODSUM,SDSUM                BMX 0025
      COMMON/BMDCOM/IBNDWD,IP,ITB,NTEMP,TBAND(MXTEMP),SDZ(MXTEMP),      BMX 0026
     1  ODZ(MXTEMP),IBIN,IMOL,IALF,SD(MXTEMP,MMOLT2),OD(MXTEMP,MMOLT),  BMX 0027
     2  ALF0(MMOLT),T5(LAYTHR),PTM75(LAYTHR),JJ(LAYTHR),FF(LAYTHR),     BMX 0028
     3  T5S(LAYTHR,NSPC),PTM75S(LAYTHR,NSPC),JJS(LAYTHR,MMOLT),         BMX 0029
     4  FFS(LAYTHR,MMOLT),DOPFAC(MMOLT),DOP0(MMOLT),SDSUM(MMOLT2),      BMX 0030
     5  ODSUM(MMOLT),DOPSUM(MMOLT),COLSUM(MMOLT)                        BMX 0031
C                                                                       BMX 0032
C     FILL THE SDX, ODX AND ALF0X MATRICES.                             BMX 0033
      IMX=IMOLX                                                         BMX 0034
      NSPED=NSPC+NSPC                                                   BMX 0035
      IF(IMX.LE.NSPECX)THEN                                             BMX 0036
C                                                                       BMX 0037
C        LINE CENTER CONTRIBUTIONS                                      BMX 0038
C                                                                       BMX 0039
         IMTAIL=IMX+NSPECX                                              BMX 0040
         DO 10 IT=1,NTEMP                                               BMX 0041
            SD(IT,IMX+NSPED)=SDZX(IT)                                   BMX 0042
            OD(IT,IMX+NSPC)=ODZX(IT)                                    BMX 0043
            SD(IT,IMTAIL+NSPED)=0.                                      BMX 0044
 10      CONTINUE                                                       BMX 0045
         ALF0(IMX+NSPC)=1.E-04*IALFX                                    BMX 0046
      ELSE                                                              BMX 0047
C                                                                       BMX 0048
C        LINE TAIL CONTRIBUTIONS                                        BMX 0049
C                                                                       BMX 0050
         JMX=IALFX                                                      BMX 0051
         IF(JMX.EQ.0)THEN                                               BMX 0052
C                                                                       BMX 0053
C           ONE TAIL                                                    BMX 0054
C                                                                       BMX 0055
            DO 20 IT=1,NTEMP                                            BMX 0056
               SD(IT,IMX+NSPED)=SDZX(IT)                                BMX 0057
 20         CONTINUE                                                    BMX 0058
         ELSE                                                           BMX 0059
C                                                                       BMX 0060
C           TWO TAILS                                                   BMX 0061
C                                                                       BMX 0062
            DO 30 IT=1,NTEMP                                            BMX 0063
               SD(IT,JMX+NSPED)=ODZX(IT)                                BMX 0064
               SD(IT,IMX+NSPED)=SDZX(IT)                                BMX 0065
 30         CONTINUE                                                    BMX 0066
         ENDIF                                                          BMX 0067
      ENDIF                                                             BMX 0068
      RETURN                                                            BMX 0069
      END                                                               BMX 0070
