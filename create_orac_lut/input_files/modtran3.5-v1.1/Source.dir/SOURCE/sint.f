      SUBROUTINE SINT(V1,V1C,DV,NPT,CONTI,CONTO)                        SIN 0001
C                                                                       SIN 0002
C     INTERPOLATION  FOR CONTINUUM WITH LOWTRAN                         SIN 0003
C                                                                       SIN 0004
      DIMENSION CONTI(2003)                                             SIN 0005
      CONTO=0.                                                          SIN 0006
      I=(V1C-V1)/DV+1.00001                                             SIN 0007
      IF(I.GE.NPT)GO TO 10                                              SIN 0008
      CONTO=CONTI(I)                                                    SIN 0009
      IMOD=AMOD(V1C,10.)                                                SIN 0010
      IF(IMOD.GT.0) CONTO=(CONTI(I)+CONTI(I+1))/2.                      SIN 0011
10    CONTINUE                                                          SIN 0012
      RETURN                                                            SIN 0013
      END                                                               SIN 0014
