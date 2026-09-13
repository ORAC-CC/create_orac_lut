      LOGICAL FUNCTION  WRTDIM ( DIMNAM, MINVAL )                       WDM 0001
                                                                        WDM 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            WDM 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                WDM 0004
C          WRITE NAME OF TOO-SMALL SYMBOLIC DIMENSION AND               WDM 0005
C          THE VALUE IT SHOULD BE INCREASED TO;  RETURN 'TRUE'          WDM 0006
                                                                        WDM 0007
C      INPUT :  DIMNAM = NAME OF SYMBOLIC DIMENSION WHICH IS TOO SMALL  WDM 0008
C                        ( CHARACTER, ANY LENGTH )                      WDM 0009
C               MINVAL = VALUE TO WHICH THAT DIMENSION SHOULD BE        WDM 0010
C                        INCREASED (AT LEAST)                           WDM 0011
C ----------------------------------------------------------------------WDM 0012
      CHARACTER*(*)  DIMNAM                                             WDM 0013
      INTEGER        MINVAL                                             WDM 0014
                                                                        WDM 0015
                                                                        WDM 0016
      WRITE ( *, '(3A,I7)' )  ' ****  SYMBOLIC DIMENSION  ', DIMNAM,    WDM 0017
     $                     '  SHOULD BE INCREASED TO AT LEAST ', MINVAL WDM 0018
      WRTDIM = .TRUE.                                                   WDM 0019
      RETURN                                                            WDM 0020
      END                                                               WDM 0021
