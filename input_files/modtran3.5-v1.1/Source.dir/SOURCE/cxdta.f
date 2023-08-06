      FUNCTION   CXDTA(V,IWL,IWH,CP,IND)                                CXD 0001
C                                                                       CXD 0002
C     THIS FUNCTION IS A BIT DIFFERENT FROM AND REPLACES                CXD 0003
C     THE OLD ROUTINE, A SUBROUTINE, CXDTA                              CXD 0004
C                                                                       CXD 0005
C     THIS SUBROUTINE FINDS THE C' FOR THE WAVENUMBER V.                CXD 0006
C     INPUT:         V --- WAVENUMBER                                   CXD 0007
C            (IWL,IWH) --- WAVENUMBER PAIR SPECIFIES THE ABSORPTION     CXD 0008
C                          REGION. BOTH ARE ARRAYS AND TERMINATED       CXD 0009
C                          WITH THE VALUE -999                          CXD 0010
C                   CP --- ARRAY CONTAINS THE C'S                       CXD 0011
C     I/O:         IND --- INDICATOR INDICATES THE ABSORPTION REGION    CXD 0012
C                          WHERE THE WAVENUMBER IS EXPECTED TO BE IN    CXD 0013
C                          OR NEARBY (IT SERVES FOR THE PURPOSE         CXD 0014
C                          TO SPEED UP THE SEARCHING PROCESS)           CXD 0015
      DIMENSION IWL(*),IWH(*),CP(*)                                     CXD 0016
      IV=V                                                              CXD 0017
      CXDTA=-20.0                                                       CXD 0018
      IF (IWL(IND+1) .EQ. -999 .AND. IV .GT. IWH(IND)) RETURN           CXD 0019
      IF (IV .LT. IWL(1)) RETURN                                        CXD 0020
      IC=0                                                              CXD 0021
  100 IF (IV .GE. IWL(IND) .AND. IV .LE. IWH(IND)) GO TO 200            CXD 0022
      IF (IV .GT. IWH(IND) .AND. IV .LT. IWL(IND+1)) RETURN             CXD 0023
      IND=IND+1                                                         CXD 0024
      IF (IWL(IND) .NE. -999) GO TO 100                                 CXD 0025
      IND=IND-1                                                         CXD 0026
      IF (IV .GT. IWH(IND)) RETURN                                      CXD 0027
      IND=1                                                             CXD 0028
      GO TO 100                                                         CXD 0029
  200 IF (IND .EQ. 1) GO TO 400                                         CXD 0030
      INDM1=IND-1                                                       CXD 0031
      DO 300 I=1,INDM1                                                  CXD 0032
        IC=IC+(IWH(I)-IWL(I))/5+1                                       CXD 0033
  300 CONTINUE                                                          CXD 0034
  400 IC=IC+(IV-IWL(IND))/5+1                                           CXD 0035
      CXDTA=CP(IC)                                                      CXD 0036
      RETURN                                                            CXD 0037
      END                                                               CXD 0038
