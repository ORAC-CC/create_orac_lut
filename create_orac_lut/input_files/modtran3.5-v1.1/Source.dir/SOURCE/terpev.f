      SUBROUTINE  TERPEV( CWT, EVECC, GL, GU, MAZIM, MXCMU, MXUMU,      TEV 0001
     $                    NN, NSTR, NUMU, WK, YLMC, YLMU )              TEV 0002
                                                                        TEV 0003
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            TEV 0004
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                TEV 0005
C         INTERPOLATE EIGENVECTORS TO USER ANGLES; EQ SD(8)             TEV 0006
                                                                        TEV 0007
      REAL*8  CWT(*), EVECC( MXCMU,* ), GL(0:*), GU(  MXUMU,* ), WK(*), TEV 0008
     $      YLMC(  0:MXCMU,* ), YLMU(  0:MXCMU,* )                      TEV 0009
                                                                        TEV 0010
                                                                        TEV 0011
      DO 50  IQ = 1, NSTR                                               TEV 0012
                                                                        TEV 0013
         DO 20  L = MAZIM, NSTR-1                                       TEV 0014
C                                       ** INNER SUM IN SD(8) TIMES ALL TEV 0015
C                                   ** FACTORS IN OUTER SUM BUT PLM(MU) TEV 0016
            SUM = 0.0                                                   TEV 0017
            DO 10  JQ = 1, NSTR                                         TEV 0018
               SUM = SUM + CWT(JQ) * YLMC(L,JQ) * EVECC(JQ,IQ)          TEV 0019
10          CONTINUE                                                    TEV 0020
            WK(L+1) = 0.5 * GL(L) * SUM                                 TEV 0021
20       CONTINUE                                                       TEV 0022
C                                    ** FINISH OUTER SUM IN SD(8)       TEV 0023
C                                    ** AND STORE EIGENVECTORS          TEV 0024
         DO 40  IU = 1, NUMU                                            TEV 0025
            SUM = 0.                                                    TEV 0026
            DO 30  L = MAZIM, NSTR-1                                    TEV 0027
               SUM = SUM + WK(L+1) * YLMU(L,IU)                         TEV 0028
30          CONTINUE                                                    TEV 0029
            IF ( IQ.LE.NN )  GU( IU, IQ+NN     ) = SUM                  TEV 0030
            IF ( IQ.GT.NN )  GU( IU, NSTR+1-IQ ) = SUM                  TEV 0031
40       CONTINUE                                                       TEV 0032
                                                                        TEV 0033
50    CONTINUE                                                          TEV 0034
                                                                        TEV 0035
      RETURN                                                            TEV 0036
      END                                                               TEV 0037
