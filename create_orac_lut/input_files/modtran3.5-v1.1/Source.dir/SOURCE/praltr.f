      SUBROUTINE  PRALTR( UMU, NUMU, ALBMED, TRNMED )                   PAT 0001
                                                                        PAT 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            PAT 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                PAT 0004
C        PRINT PLANAR ALBEDO AND TRANSMISSIVITY OF MEDIUM               PAT 0005
C        AS A FUNCTION OF INCIDENT BEAM ANGLE                           PAT 0006
                                                                        PAT 0007
      REAL*8     UMU(*), ALBMED(*), TRNMED(*)                           PAT 0008
                                                                        PAT 0009
                                                                        PAT 0010
      WRITE( *,110 )                                                    PAT 0011
      DO 20  IU = 1, NUMU                                               PAT 0012
         ANGL = 180.0/3.14159265 * DACOS( UMU(IU) )                     PAT 0013
         WRITE(*,111)  ANGL, UMU(IU), ALBMED(IU), TRNMED(IU)            PAT 0014
 20   CONTINUE                                                          PAT 0015
                                                                        PAT 0016
      RETURN                                                            PAT 0017
                                                                        PAT 0018
110   FORMAT( ///, ' *******  FLUX ALBEDO AND/OR TRANSMISSIVITY OF ',   PAT 0019
     $ 'ENTIRE MEDIUM  ********', //,                                   PAT 0020
     $ ' BEAM ZEN ANG  dcos(BEAM ZEN ANG)      ALBEDO   TRANSMISSIVITY')PAT 0021
111   FORMAT( 0P,F13.4, F20.6, F12.5, 1P,E17.4 )                        PAT 0022
      END                                                               PAT 0023
