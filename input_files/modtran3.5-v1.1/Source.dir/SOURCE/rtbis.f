      DOUBLE PRECISION FUNCTION RTBIS(X1,X2,CPATH)                      RTB 0001
C                                                                       RTB 0002
C     THIS FUNCTION FINDS THE ROOT OF FUNC(X) = X*REFRACTIVE INDEX - CPARTB 0003
C     THE ROOT IS ACTUALLY THE TANGENT HEIGHT.                          RTB 0004
C     IT IS SANDWICHED BETWEEN X1 AND X2.                               RTB 0005
C     THIS ROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.          RTB 0006
C                                                                       RTB 0007
      DOUBLE PRECISION X1, X2, CPATH                                    RTB 0008
      DOUBLE PRECISION F,FMID,DX,XMID,XACC,F1,F2,RATIO,DPRE             RTB 0009
      INTEGER J, JMAX                                                   RTB 0010
      REAL RE,ZMAX                                                      RTB 0011
      INTEGER IMAX,IMOD,IPATH                                           RTB 0012
      COMMON/PARMTR/RE,ZMAX,IMAX,IMOD,IPATH                             RTB 0013
      PARAMETER (JMAX=40)                                               RTB 0014
      DATA XACC/1.D-5/                                                  RTB 0015
      DPRE=DBLE(RE)                                                     RTB 0016
      CALL IRFXN(X2,FMID, RATIO)                                        RTB 0017
      FMID = FMID*(X2+DPRE)-CPATH                                       RTB 0018
      CALL IRFXN(X1,F, RATIO)                                           RTB 0019
      F = F*(X1+DPRE)-CPATH                                             RTB 0020
      IF(F*FMID.GT.0.)                                                  RTB 0021
     1  STOP 'ERROR in RTBIS:  Root must be bracketed for bisection.'   RTB 0022
      IF(F.LT.0.)THEN                                                   RTB 0023
         RTBIS=X1                                                       RTB 0024
         DX=X2-X1                                                       RTB 0025
      ELSE                                                              RTB 0026
         RTBIS=X2                                                       RTB 0027
         DX=X1-X2                                                       RTB 0028
      ENDIF                                                             RTB 0029
      DO 11 J=1,JMAX                                                    RTB 0030
         DX=DX/2                                                        RTB 0031
         XMID=RTBIS+DX                                                  RTB 0032
         CALL IRFXN(XMID,FMID, RATIO)                                   RTB 0033
         FMID = FMID*(XMID+DPRE)-CPATH                                  RTB 0034
         IF(FMID.LE.0.)RTBIS=XMID                                       RTB 0035
         IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN                     RTB 0036
 11   CONTINUE                                                          RTB 0037
C                                                                       RTB 0038
C     COMES HERE IF UNABLE TO SOLVE.                                    RTB 0039
C                                                                       RTB 0040
      CALL IRFXN(X2,F2, RATIO)                                          RTB 0041
      F2 = F2*(X2+DPRE)-CPATH                                           RTB 0042
      CALL IRFXN(X1,F1, RATIO)                                          RTB 0043
      F1 = F1*(X1+DPRE)-CPATH                                           RTB 0044
      IF (ABS(F2) .LT. ABS(F1)) THEN                                    RTB 0045
         RTBIS = X2                                                     RTB 0046
      ELSE                                                              RTB 0047
         RTBIS = X1                                                     RTB 0048
      ENDIF                                                             RTB 0049
      END                                                               RTB 0050
