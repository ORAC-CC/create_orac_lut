      SUBROUTINE  ERRMSG( MESSAG, FATAL )                               ERR 0001
                                                                        ERR 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            ERR 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                ERR 0004
C        PRINT OUT A WARNING OR ERROR MESSAGE;  ABORT IF ERROR          ERR 0005
                                                                        ERR 0006
      LOGICAL       FATAL, ONCE                                         ERR 0007
      CHARACTER*(*) MESSAG                                              ERR 0008
      INTEGER       MAXMSG, NUMMSG                                      ERR 0009
      SAVE          MAXMSG, NUMMSG, ONCE                                ERR 0010
      DATA NUMMSG / 0 /,  MAXMSG / 100 /,  ONCE / .FALSE. /             ERR 0011
                                                                        ERR 0012
                                                                        ERR 0013
      IF ( FATAL )  THEN                                                ERR 0014
         WRITE ( *, '(/,2A)' )  ' ******* ERROR >>>>>>  ', MESSAG       ERR 0015
         STOP                                                           ERR 0016
      END IF                                                            ERR 0017
                                                                        ERR 0018
      NUMMSG = NUMMSG + 1                                               ERR 0019
      IF ( NUMMSG.GT.MAXMSG )  THEN                                     ERR 0020
         IF ( .NOT.ONCE )  WRITE ( *,99 )                               ERR 0021
         ONCE = .TRUE.                                                  ERR 0022
      ELSE                                                              ERR 0023
         WRITE ( *, '(/,2A)' )  ' ******* WARNING >>>>>>  ', MESSAG     ERR 0024
      ENDIF                                                             ERR 0025
                                                                        ERR 0026
      RETURN                                                            ERR 0027
                                                                        ERR 0028
   99 FORMAT( ///,' >>>>>>  TOO MANY WARNING MESSAGES --  ',            ERR 0029
     $   'THEY WILL NO LONGER BE PRINTED  <<<<<<<', /// )               ERR 0030
      END                                                               ERR 0031
