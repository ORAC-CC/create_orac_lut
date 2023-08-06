      SUBROUTINE GTSTRT(IV,IREC,ITB)                                    GST 0001
C                                                                       GST 0002
C     FIND THE RECORD (IREC) WHERE FREQUENCY IV STARTS IN THE           GST 0003
C     BAND MODEL PARAMETER TAPE.  ITB IS THE UNIT NUMBER.               GST 0004
C                                                                       GST 0005
C     DECLARE INPUTS/OUTPUTS                                            GST 0006
      INTEGER IV,IREC,ITB                                               GST 0007
C                                                                       GST 0008
C     DECLARE LOCAL VARIABLES                                           GST 0009
      INTEGER IV2,IRCLO,IRCHI,IRCTST,IDR                                GST 0010
      IRCLO=2                                                           GST 0011
C                                                                       GST 0012
C     FIND IRCHI SO THAT IT IS GREATER THAN OR EQUAL TO MAXIMUM RECORD. GST 0013
C     IRCHI = 2**16 = 65536 SHOULD EXCEED NUMBER OF RECORDS.            GST 0014
C     SO THIS RECORD IS OFF THE FILE.                                   GST 0015
C     DO THE FOLLOWING IN CASE IRCHI IS NOT LARGE ENOUGH.               GST 0016
C                                                                       GST 0017
      IRCHI=65536                                                       GST 0018
   10 READ(ITB,REC=IRCHI,ERR=20)IV2                                     GST 0019
      IRCHI=IRCHI+IRCHI                                                 GST 0020
      GOTO10                                                            GST 0021
C                                                                       GST 0022
C     START BISECTION LOOP.                                             GST 0023
C     FIND FREQUENCY AT THE MID-RECORD BETWEEN IRCLO AND IRCHI.         GST 0024
C     RESET IRCLO AND IRCHI.                                            GST 0025
C     TAKES ONLY 16 READS (MAXIMUM) TO LOCATE THE STARTING POINT.       GST 0026
C                                                                       GST 0027
   20 IDR=IRCHI-IRCLO                                                   GST 0028
      IF(IDR.GT.1)THEN                                                  GST 0029
         IRCTST=IRCLO+IDR/2                                             GST 0030
         READ(ITB,REC=IRCTST,ERR=30)IV2                                 GST 0031
         IF(IV2.LT.IV)THEN                                              GST 0032
            IRCLO=IRCTST                                                GST 0033
            GOTO20                                                      GST 0034
         ENDIF                                                          GST 0035
   30    IRCHI=IRCTST                                                   GST 0036
         GOTO20                                                         GST 0037
      ENDIF                                                             GST 0038
      READ(ITB,REC=IRCLO)IV2                                            GST 0039
      IF (IV2.EQ.IV)THEN                                                GST 0040
         IREC=IRCLO                                                     GST 0041
      ELSE                                                              GST 0042
         IREC=IRCHI                                                     GST 0043
      ENDIF                                                             GST 0044
      END                                                               GST 0045
