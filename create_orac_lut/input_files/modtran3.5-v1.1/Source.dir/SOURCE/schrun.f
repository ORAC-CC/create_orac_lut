      SUBROUTINE SCHRUN(V,CPRUN)                                        SCH 0001
      COMMON /SHUR/ SHN(430)                                            SCH 0002
      DATA V1,V2,DV,INUM /49600.,51710.,5.,425/                         SCH 0003
      CPRUN = -20.                                                      SCH 0004
      IF(V .LT. V1) GO TO 20                                            SCH 0005
      IF(V .GT. V2) GO TO 20                                            SCH 0006
      IND = (V - V1)/DV + 1.0001                                        SCH 0007
      IF(IND . GT. INUM) THEN                                           SCH 0008
            PRINT*,'  IND GT INUM  V IND ',V,IND                        SCH 0009
            GO TO 20                                                    SCH 0010
      ENDIF                                                             SCH 0011
      CPRUN = SHN(IND)                                                  SCH 0012
20    RETURN                                                            SCH 0013
      END                                                               SCH 0014
