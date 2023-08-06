      LOGICAL FUNCTION  WRTBAD ( VARNAM )                               WBD 0001
                                                                        WBD 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            WBD 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                WBD 0004
C          WRITE NAMES OF ERRONEOUS VARIABLES AND RETURN 'TRUE'         WBD 0005
                                                                        WBD 0006
C      INPUT :   VARNAM = NAME OF ERRONEOUS VARIABLE TO BE WRITTEN      WBD 0007
C                         ( CHARACTER, ANY LENGTH )                     WBD 0008
C ----------------------------------------------------------------------WBD 0009
      CHARACTER*(*)  VARNAM                                             WBD 0010
      INTEGER        MAXMSG, NUMMSG                                     WBD 0011
      SAVE  NUMMSG, MAXMSG                                              WBD 0012
      DATA  NUMMSG / 0 /,  MAXMSG / 50 /                                WBD 0013
                                                                        WBD 0014
                                                                        WBD 0015
      WRTBAD = .TRUE.                                                   WBD 0016
      NUMMSG = NUMMSG + 1                                               WBD 0017
      WRITE ( *, '(3A)' )  ' ****  INPUT VARIABLE  ', VARNAM,           WBD 0018
     $                     '  IN ERROR  ****'                           WBD 0019
      IF ( NUMMSG.EQ.MAXMSG )                                           WBD 0020
     $   CALL  ERRMSG ( 'TOO MANY INPUT ERRORS.  ABORTING...$', .TRUE. )WBD 0021
      RETURN                                                            WBD 0022
      END                                                               WBD 0023
