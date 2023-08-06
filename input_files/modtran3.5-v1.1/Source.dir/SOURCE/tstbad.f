      LOGICAL FUNCTION  TSTBAD( VARNAM, RELERR )                        TBD 0001
                                                                        TBD 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            TBD 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                TBD 0004
C       WRITE NAME (-VARNAM-) OF VARIABLE FAILING SELF-TEST AND ITS     TBD 0005
C       PERCENT ERROR FROM THE CORRECT VALUE;  RETURN  'FALSE'.         TBD 0006
                                                                        TBD 0007
      CHARACTER*(*)  VARNAM                                             TBD 0008
      REAL*8           RELERR                                           TBD 0009
                                                                        TBD 0010
                                                                        TBD 0011
      TSTBAD = .FALSE.                                                  TBD 0012
      WRITE( *, '(/,3A,1P,E11.2,A)' )                                   TBD 0013
     $       ' OUTPUT VARIABLE  ', VARNAM,'  DIFFERED BY', 100.*RELERR, TBD 0014
     $       '  PER CENT FROM CORRECT VALUE.  SELF-TEST FAILED.'        TBD 0015
      RETURN                                                            TBD 0016
      END                                                               TBD 0017
