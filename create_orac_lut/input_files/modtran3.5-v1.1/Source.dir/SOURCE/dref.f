      REAL*8 FUNCTION  DREF( MU, HL, NSTR )                             REF 0001
C           INSERT FOR DOUBLE PRECISION - NORTH                         REF 0002
            IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                      REF 0003
C        EXACT FLUX ALBEDO FOR GIVEN ANGLE OF INCIDENCE, GIVEN          REF 0004
C        A BIDIRECTIONAL REFLECTIVITY CHARACTERIZED BY ITS              REF 0005
C        LEGENDRE COEFFICIENTS ( NOTE** THESE WILL ONLY AGREE           REF 0006
C        WITH BOTTOM-BOUNDARY ALBEDOS CALCULATED BY 'DISORT' IN         REF 0007
C        THE LIMIT AS NUMBER OF STREAMS GO TO INFINITY, BECAUSE         REF 0008
C        'DISORT' EVALUATES THE INTEGRAL 'CL' ONLY APPROXIMATELY,       REF 0009
C        BY QUADRATURE, WHILE THIS ROUTINE CALCULATES IT EXACTLY. )     REF 0010
                                                                        REF 0011
C      INPUT :   MU     COSINE OF INCIDENCE ANGLE                       REF 0012
C                HL     LEGENDRE COEFFICIENTS OF BIDIRECTIONAL REF'Y    REF 0013
C              NSTR     NUMBER OF ELEMENTS OF 'HL' TO CONSIDER          REF 0014
                                                                        REF 0015
C      INTERNAL VARIABLES (P-SUB-L IS THE L-TH LEGENDRE POLYNOMIAL) :   REF 0016
                                                                        REF 0017
C              CL    INTEGRAL FROM 0 TO 1 OF  MU * P-SUB-L(MU)          REF 0018
C                       (VANISHES FOR  L = 3, 5, 7, ... )               REF 0019
C              PL    P-SUB-L                                            REF 0020
C            PLM1    P-SUB-(L-1)                                        REF 0021
C            PLM2    P-SUB-(L-2)                                        REF 0022
                                                                        REF 0023
      PARAMETER  ( MAXTRM = 100 )                                       REF 0024
      LOGICAL      PASS1                                                REF 0025
      REAL*8         MU, HL( 0:* )                                      REF 0026
      REAL*8         C( MAXTRM )                                        REF 0027
      SAVE  PASS1, C                                                    REF 0028
      DATA  PASS1 / .TRUE. /                                            REF 0029
                                                                        REF 0030
                                                                        REF 0031
      IF ( PASS1 )  THEN                                                REF 0032
         PASS1 = .FALSE.                                                REF 0033
         CL = 0.125                                                     REF 0034
         C(2) = 10. * CL                                                REF 0035
         DO 1  L = 4, MAXTRM, 2                                         REF 0036
            CL = - CL * (L-3) / (L+2)                                   REF 0037
            C(L) = 2. * (2*L+1) * CL                                    REF 0038
    1    CONTINUE                                                       REF 0039
      END IF                                                            REF 0040
                                                                        REF 0041
      IF ( NSTR.GT.MAXTRM )  CALL                                       REF 0042
     $     ERRMSG( 'DREF--PARAMETER MAXTRM TOO SMALL', .TRUE. )         REF 0043
                                                                        REF 0044
      DREF = HL(0) - 2.*HL(1) * MU                                      REF 0045
      PLM2 = 1.0                                                        REF 0046
      PLM1 = - MU                                                       REF 0047
      DO 10  L = 2, NSTR-1                                              REF 0048
C                                ** LEGENDRE POLYNOMIAL RECURRENCE      REF 0049
                                                                        REF 0050
         PL = ( (2*L-1) * (-MU) * PLM1 - (L-1) * PLM2 ) / L             REF 0051
         IF( MOD(L,2).EQ.0 )  DREF = DREF + C(L) * HL(L) * PL           REF 0052
         PLM2 = PLM1                                                    REF 0053
         PLM1 = PL                                                      REF 0054
   10 CONTINUE                                                          REF 0055
                                                                        REF 0056
      RETURN                                                            REF 0057
      END                                                               REF 0058
