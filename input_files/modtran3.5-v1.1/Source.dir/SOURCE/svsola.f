      SUBROUTINE SVSOLA(IPARM,IPH,IDAY,ISOURC,PARM1,PARM2,PARM3,PARM4,  SSL 0001
     $     TIME,G,ANGLEM,                                               SSL 0002
     $     ISAVE1,ISAVE2,ISAVE3,ISAVE4,SAVE1,SAVE2,SAVE3,SAVE4,         SSL 0003
     $     SAVE5,SAVE6,SAVE7)                                           SSL 0004
C                                                                       SSL 0005
C     SAVE THE CURRENT SOLAR PARAMETERS                                 SSL 0006
C                                                                       SSL 0007
C     IMPLICIT UNDEFINED(A-Z)                                           SSL 0008
      INTEGER IPARM,IDAY,ISOURC,ISAVE1,ISAVE2,ISAVE3,IPH,ISAVE4         SSL 0009
      REAL PARM1,PARM2,PARM3,PARM4,TIME,G,ANGLEM,                       SSL 0010
     $     SAVE1,SAVE2,SAVE3,SAVE4,SAVE5,SAVE6,SAVE7                    SSL 0011
      ISAVE1 = IPARM                                                    SSL 0012
      ISAVE2 = IPH                                                      SSL 0013
      ISAVE3 = IDAY                                                     SSL 0014
      ISAVE4 = ISOURC                                                   SSL 0015
      SAVE1 = PARM1                                                     SSL 0016
      SAVE2 = PARM2                                                     SSL 0017
      SAVE3 = PARM3                                                     SSL 0018
      SAVE4 = PARM4                                                     SSL 0019
      SAVE5 = TIME                                                      SSL 0020
      SAVE6 = G                                                         SSL 0021
      SAVE7 = ANGLEM                                                    SSL 0022
      RETURN                                                            SSL 0023
      END                                                               SSL 0024
