      REAL FUNCTION DBLTX(W,CPRIME,QA)                                  DBL 0001
C                                                                       DBL 0002
C     DOUBLE EXPONENTIAL TRANSMITTANCE ROUTINE.                         DBL 0003
      REAL W,CPRIME,QA,QAWS                                             DBL 0004
C                                                                       DBL 0005
C     NO OPTICAL DEPTH                                                  DBL 0006
      DBLTX=1.                                                          DBL 0007
      IF(W.LT.1.E-20 .OR. CPRIME.LE.-20.)RETURN                         DBL 0008
      QAWS=QA*(ALOG10(W)+CPRIME)                                        DBL 0009
      IF(QAWS.LE.-6.)RETURN                                             DBL 0010
C                                                                       DBL 0011
C     OPTICAL DEPTH OVERFLOW                                            DBL 0012
      DBLTX=0.                                                          DBL 0013
      IF(QAWS.GE.2.)RETURN                                              DBL 0014
C                                                                       DBL 0015
C     STANDARD TRANSMITTANCE                                            DBL 0016
      DBLTX=EXP(-10.**QAWS)                                             DBL 0017
      RETURN                                                            DBL 0018
      END                                                               DBL 0019
