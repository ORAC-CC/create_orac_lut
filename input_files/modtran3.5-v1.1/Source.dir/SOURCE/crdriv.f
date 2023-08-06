      SUBROUTINE CRDRIV                                                 CRD 0001
C                                                                       CRD 0002
C     THIS ROUTINE IS THE DRIVER FOR CLOUD/RAIN MODELS 1 THROUGH 10.    CRD 0003
      REAL CTHICK,CWDCOL,CIPCOL                                         CRD 0004
C                                                                       CRD 0005
C     DETERMINE CLOUD/RAIN DENSITY PROFILES                             CRD 0006
      CALL CRPROF(CTHICK,CWDCOL,CIPCOL)                                 CRD 0007
C                                                                       CRD 0008
C     MERGE OLD ATMOSPHERIC AND CLOUD/RAIN LAYER BOUNDARY DATA.         CRD 0009
      CALL CRMERG                                                       CRD 0010
C                                                                       CRD 0011
C     DETERMINE CLOUD/RAIN SPECTRAL DATA                                CRD 0012
      CALL CRSPEC(CTHICK,CWDCOL,CIPCOL)                                 CRD 0013
C                                                                       CRD 0014
C     RETURN TO MODTRAN DRIVER ROUTINE                                  CRD 0015
      RETURN                                                            CRD 0016
      END                                                               CRD 0017
