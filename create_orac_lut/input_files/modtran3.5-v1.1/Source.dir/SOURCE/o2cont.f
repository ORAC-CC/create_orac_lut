      SUBROUTINE O2CONT(V,SIGMA,ALPHA,BETA)                             O2C 0001
C                                                                       O2C 0002
C     THIS ROUTINE IS DRIVEN BY FREQUENCY, RETURNING ONLY THE           O2C 0003
C     O2 COEFFICIENTS, INDEPENDENT OF TEMPERATURE.                      O2C 0004
C                                                                       O2C 0005
C  *******************************************************************  O2C 0006
C  *  THESE COMMENTS APPLY TO THE COLUME ARRAYS FOR:                 *  O2C 0007
C  *       PBAR*UBAR(O2)                                             *  O2C 0008
C  *       PBAR*UBAR(O2)*DT                                          *  O2C 0009
C  *   AND PBAR*UBAR(O2)*DT*DT    WHERE:  DT=TBAR-220.               *  O2C 0010
C  *  THAT HAVE BEEN COMPILED IN OTHER PARTS OF THE LOWTRAN CODE     *  O2C 0011
C  *                                                                 *  O2C 0012
C  *  LOWTRAN7 COMPATIBLE:                                           *  O2C 0013
C  *  O2 CONTINUUM SUBROUTINE FOR 1395-1760CM-1                      *  O2C 0014
C  *  MODIFIED BY G.P. ANDERSON, APRIL '88                           *  O2C 0015
C  *                                                                 *  O2C 0016
C  *  THE EXPONENTIAL TEMPERATURE EMPLOYED IN THE FASCOD2 ALGORITHM  *  O2C 0017
C  *  (SEE BELOW) IS NOT READILY SUITABLE FOR LOWTRAN.  THEREFORE    *  O2C 0018
C  *  THE EXPONENTIALS HAVE BEEN LINEARLY EXPANDED, KEEPING ONLY THE *  O2C 0019
C  *  LINEAR AND QUADRATIC TERMS:                                    *  O2C 0020
C  *                                                                 *  O2C 0021
C  *  EXP(A*DT)=1.+ A*DT + (A*DT)**2/2. + ....                       *  O2C 0022
C  *                                                                 *  O2C 0023
C  *     EXP(B*DT*DT)=1.+ B*DT*DT + (B*DT*DT)**2/2. + ....           *  O2C 0024
C  *                                                                 *  O2C 0025
C  *  THE PRODUCT OF THE TWO TERMS IS:                               *  O2C 0026
C  *                                                                 *  O2C 0027
C  *     (1. + A*DT + (A*A/2. + B)*DT*DT )                           *  O2C 0028
C  *                                                                 *  O2C 0029
C  *  THIS EXPANSION ONLY WORKS WELL FOR SMALL VALUES OF X IN EXP(X) *  O2C 0030
C  *                                                                 *  O2C 0031
C  *  SINCE DT = T-220., THE APPROXIMATION IS VERY GOOD UNTIL        *  O2C 0032
C  *  T.GT.260. OR DT.GT.40.   AT T=280, THE MAXIMUM ERRORS ARE STILL*  O2C 0033
C  *  LESS THAN 10% BUT AT T=300, THOSE ERRORS ARE AS LARGE AS 20%   *  O2C 0034
C  *******************************************************************  O2C 0035
C                                                                       O2C 0036
C     THE FOLLOWING COMMENTS ARE EXCERPTED DIRECTLY FROM FASCOD2        O2C 0037
C                                                                       O2C 0038
C      THIS SUBROUTINE CONTAINS THE ROGERS AND WALSHAW                  O2C 0039
C      EQUIVALENT COEFFICIENTS DERIVED FROM THE THEORETICAL             O2C 0040
C      VALUES SUPPLIED BY ROLAND DRAYSON. THESE VALUES USE              O2C 0041
C      THE SAME DATA AS TIMOFEYEV AND AGREE WITH TIMOFEYEV'S RESULTS.   O2C 0042
C      THE DATA ARE IN THE FORM OF STRENGTHS(O2SO) AND TWO              O2C 0043
C      COEFFICIENTS (O2A & O2B),  WHICH ARE USED TO CORRECT FOR         O2C 0044
C      TEMPERATURE. THE DEPENDENCY ON PRESSURE SQUARED                  O2C 0045
C      IS CONTAINED IN THE P*WO2 PART OF THE CONSTANT.                  O2C 0046
C      NOTE THAT SINCE THE COEFFICIENTS ARE FOR AIR, THE                O2C 0047
C      THE STRENGTHS ARE DIVIDED BY THE O2 MIXING RATIO FOR             O2C 0048
C      DRY AIR OF 0.20946 (THIS IS ASSUMED CONSTANT).                   O2C 0049
C      ORIGINAL FORMULATION OF THE COEFFICIENTS WAS BY LARRY GORDLEY.   O2C 0050
C      THIS VERSION WRITTEN BY EARL THOMPSON, JULY 1984.                O2C 0051
C                                                                       O2C 0052
C                                                                       O2C 0053
      COMMON/O2C/ O2DRAY(74),O2C001(74),O2S0(74),O2A(74),O2B(74),       O2C 0054
     X V1O2,V2O2,DVO2,NPTO2                                             O2C 0055
      SIGMA =0.                                                         O2C 0056
      ALPHA =0.                                                         O2C 0057
      BETA  =0.                                                         O2C 0058
      IF(V .LT. 1395) GO TO 30                                          O2C 0059
      IF(V .GT. 1760) GO TO 30                                          O2C 0060
C                                                                       O2C 0061
C                                                                       O2C 0062
      CALL O2INT(V,V1O2,DVO2,NPTO2,C,O2S0,A,O2A,B,O2B)                  O2C 0063
C                                                                       O2C 0064
C                                                                       O2C 0065
C                                                                       O2C 0066
C     OLD 'FASCOD2' TEMPERATURE DEPENDENCE USING BLOCK DATA ARRAYS      O2C 0067
C                                                                       O2C 0068
C     C(J)=O2S0(I)* EXP(O2A(I)*TD+O2B(I)*TD*TD) /(0.20946*VJ)           O2C 0069
C                                                                       O2C 0070
C     NEW COEFFICIENT DEFINITIONS FOR LOWTRAN FORMULATION               O2C 0071
C                                                                       O2C 0072
      ALPHA= A                                                          O2C 0073
      BETA=A**2/2.+B                                                    O2C 0074
      SIGMA=C/0.20946                                                   O2C 0075
C                                                                       O2C 0076
C     NEW 'LOWTRAN7' TEMPERATURE DEPENDENCE                             O2C 0077
C                                                                       O2C 0078
C     THIS WOULD BE THE CODING FOR THE LOWTRAN7 FORMULATION, BUT        O2C 0079
C       BECAUSE THE T-DEPENDENCE IS INCLUDED IN THE AMOUNTS, ONLY       O2C 0080
C       THE COEFFICIENTS (SIGMA, ALPHA & BETA) ARE BEING RETURNED       O2C 0081
C                                                                       O2C 0082
C     C(J)=SIGMA*(1.+ALPHA*TD+BETA*TD*TD)                               O2C 0083
C                                                                       O2C 0084
C     THE COEFFICIENTS FOR O2 HAVE BEEN MULTIPLIED BY A FACTOR          O2C 0085
C     OF 0.78  [RINSLAND ET AL, 1989: JGR 94; 16,303 - 16,322.].        O2C 0086
      O2FAC = 0.78                                                      O2C 0087
      SIGMA=O2FAC*SIGMA                                                 O2C 0088
30    RETURN                                                            O2C 0089
      END                                                               O2C 0090
