      SUBROUTINE  LEPOLY( NMU, M, MAXMU, TWONM1, MU, YLM )              PLY 0001
                                                                        PLY 0002
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            PLY 0003
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                PLY 0004
C       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,         PLY 0005
C       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL          PLY 0006
C       PLM = P-SUB-L-SUPER-M AS                                        PLY 0007
                                                                        PLY 0008
C             YLM(MU) = SQRT( (L-M)!/(L+M)! ) * PLM(MU)                 PLY 0009
                                                                        PLY 0010
C       FOR FIXED ORDER -M- AND ALL DEGREES FROM L = M TO TWONM1.       PLY 0011
C       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE    PLY 0012
C       FROM A PRIOR CALL TO THE ROUTINE.                               PLY 0013
                                                                        PLY 0014
C       REFERENCE: DAVE, J.V. AND B.H. ARMSTRONG, COMPUTATIONS OF       PLY 0015
C                  HIGH-ORDER ASSOCIATED LEGENDRE POLYNOMIALS,          PLY 0016
C                  J. QUANT. SPECTROSC. RADIAT. TRANSFER 10,            PLY 0017
C                  557-562, 1970.  (HEREAFTER D/A)                      PLY 0018
                                                                        PLY 0019
C       METHOD: VARYING DEGREE RECURRENCE RELATIONSHIP.                 PLY 0020
                                                                        PLY 0021
C       NOTE 1: THE D/A FORMULAS ARE TRANSFORMED BY                     PLY 0022
C               SETTING  M = N-1; L = K-1.                              PLY 0023
C       NOTE 2: ASSUMES THAT ROUTINE IS CALLED FIRST WITH  M = 0,       PLY 0024
C               THEN WITH  M = 1, ETC. UP TO  M = TWONM1.               PLY 0025
C       NOTE 3: LOOPS ARE WRITTEN IN SUCH A WAY AS TO VECTORIZE.        PLY 0026
                                                                        PLY 0027
C  I N P U T     V A R I A B L E S:                                     PLY 0028
                                                                        PLY 0029
C       NMU    :  NUMBER OF ARGUMENTS OF -YLM-                          PLY 0030
C       M      :  ORDER OF -YLM-                                        PLY 0031
C       MAXMU  :  FIRST DIMENSION OF -YLM-                              PLY 0032
C       TWONM1 :  MAX DEGREE OF -YLM-                                   PLY 0033
C       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM-                      PLY 0034
C       IF M.GT.0, YLM(M-1,I) FOR I = 1 TO NMU IS REQUIRED              PLY 0035
                                                                        PLY 0036
C  O U T P U T     V A R I A B L E:                                     PLY 0037
                                                                        PLY 0038
C       YLM(L,I) :  L = M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE     PLY 0039
C                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-           PLY 0040
C+---------------------------------------------------------------------+PLY 0041
      REAL*8     MU(*), YLM( 0:MAXMU,* )                                PLY 0042
      INTEGER  M, NMU, TWONM1                                           PLY 0043
      PARAMETER  ( MAXSQT = 1000 )                                      PLY 0044
      REAL*8     SQT( MAXSQT )                                          PLY 0045
      LOGICAL  PASS1                                                    PLY 0046
      SAVE  SQT, PASS1                                                  PLY 0047
      DATA  PASS1 / .TRUE. /                                            PLY 0048
                                                                        PLY 0049
                                                                        PLY 0050
      IF ( PASS1 )  THEN                                                PLY 0051
         PASS1 = .FALSE.                                                PLY 0052
         DO 1  NS = 1, MAXSQT                                           PLY 0053
            SQT( NS ) = DSQRT( DBLE(NS) )                               PLY 0054
    1    CONTINUE                                                       PLY 0055
      ENDIF                                                             PLY 0056
                                                                        PLY 0057
      IF ( 2*TWONM1 .GT. MAXSQT )                                       PLY 0058
     $   CALL ERRMSG( 'LEPOLY--NEED TO INCREASE PARAM MAXSQT', .TRUE. ) PLY 0059
                                                                        PLY 0060
      IF ( M .EQ. 0 )  THEN                                             PLY 0061
C                             ** UPWARD RECURRENCE FOR ORDINARY         PLY 0062
C                             ** LEGENDRE POLYNOMIALS                   PLY 0063
         DO  10  I = 1, NMU                                             PLY 0064
            YLM( 0,I ) = 1.                                             PLY 0065
            YLM( 1,I ) = MU( I )                                        PLY 0066
  10     CONTINUE                                                       PLY 0067
                                                                        PLY 0068
         DO  20  L = 2, TWONM1                                          PLY 0069
            DO  20  I = 1, NMU                                          PLY 0070
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )          PLY 0071
     $                      - ( L-1 ) * YLM( L-2,I ) ) / L              PLY 0072
  20     CONTINUE                                                       PLY 0073
                                                                        PLY 0074
      ELSE                                                              PLY 0075
                                                                        PLY 0076
         DO  30  I = 1, NMU                                             PLY 0077
C                               ** Y-SUB-M-SUPER-M; DERIVED FROM        PLY 0078
C                               ** D/A EQS. (11,12)                     PLY 0079
                                                                        PLY 0080
            YLM( M,I) = - SQT( 2*M-1 ) / SQT( 2*M )                     PLY 0081
     $                  * DSQRT( 1. - MU(I)**2 ) * YLM( M-1,I )         PLY 0082
                                                                        PLY 0083
C                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM     PLY 0084
C                              ** D/A EQS. (13,14) USING EQS. (11,12)   PLY 0085
                                                                        PLY 0086
            YLM( M+1,I ) = SQT( 2*M+1 ) * MU(I) * YLM( M,I )            PLY 0087
30       CONTINUE                                                       PLY 0088
C                                   ** UPWARD RECURRENCE; D/A EQ. (10)  PLY 0089
         DO  40  L = M+2, TWONM1                                        PLY 0090
            TMP1 = SQT( L-M ) * SQT( L+M )                              PLY 0091
            TMP2 = SQT( L-M-1 ) * SQT( L+M-1 )                          PLY 0092
            DO  40  I = 1, NMU                                          PLY 0093
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )          PLY 0094
     $                        - TMP2 * YLM( L-2,I ) ) / TMP1            PLY 0095
40       CONTINUE                                                       PLY 0096
                                                                        PLY 0097
      END IF                                                            PLY 0098
                                                                        PLY 0099
      RETURN                                                            PLY 0100
      END                                                               PLY 0101
