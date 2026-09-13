      SUBROUTINE COMPAR(IPARM,IPH,IDAY,ISOURC,PARM1,PARM2,PARM3,PARM4,  COM 0001
     $     TIME,G,ANGLEM,                                               COM 0002
     $     ISAVE1,ISAVE2,ISAVE3,ISAVE4,SAVE1,SAVE2,SAVE3,SAVE4,         COM 0003
     $     SAVE5,SAVE6,SAVE7,LSAME)                                     COM 0004
C                                                                       COM 0005
C     IPARM, IPH, ETC ARE THE CURRENT SOLAR PARAMETERS.                 COM 0006
C     ISAVE1, ISAVE2, ETC. ARE THE SOLAR PARAMETERS OF THE              COM 0007
C     IMMEDIATELY PRECEEDING RUN.                                       COM 0008
C     LSAME IS A LOGICAL WHICH IS TRUE ONLY IF THE OLD AND THE NEW MATCHCOM 0009
C     LSAME IS MEANINGFUL ONLY IF BOTH THE CURRENT AND THE IMMEDIATELY  COM 0010
C     PRECEEDING RUNS USED SOLAR PARAMETERS.                            COM 0011
C                                                                       COM 0012
C     IMPLICIT UNDEFINED(A-Z)                                           COM 0013
      INTEGER IPARM,IDAY,IPH,ISOURC,ISAVE1,ISAVE2,ISAVE3,ISAVE4         COM 0014
      REAL PARM1,PARM2,PARM3,PARM4,TIME,G,ANGLEM,                       COM 0015
     $     SAVE1,SAVE2,SAVE3,SAVE4,SAVE5,SAVE6,SAVE7                    COM 0016
      LOGICAL LSAME                                                     COM 0017
C                                                                       COM 0018
      LSAME = .FALSE.                                                   COM 0019
      IF (ISAVE1 .NE. IPARM) RETURN                                     COM 0020
      IF (ISAVE2 .NE. IPH) RETURN                                       COM 0021
      IF (ISAVE3 .NE. IDAY) RETURN                                      COM 0022
      IF (ISAVE4 .NE. ISOURC) RETURN                                    COM 0023
      IF (SAVE1 .NE. PARM1 .AND. IPARM.NE.2) RETURN                     COM 0024
      IF (SAVE2 .NE. PARM2) RETURN                                      COM 0025
      IF (SAVE3 .NE. PARM3) RETURN                                      COM 0026
      IF (SAVE4 .NE. PARM4) RETURN                                      COM 0027
      IF (SAVE5 .NE. TIME) RETURN                                       COM 0028
      IF(IPH.EQ.0 .AND. SAVE6.NE.G)RETURN                               COM 0029
      IF (SAVE7 .NE. ANGLEM) RETURN                                     COM 0030
      LSAME = .TRUE.                                                    COM 0031
      RETURN                                                            COM 0032
      END                                                               COM 0033
