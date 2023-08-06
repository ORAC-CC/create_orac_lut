      SUBROUTINE  ZEROIT( A, LENGTH )                                   ZIT 0001
C               INSERTED LINE TO DO DOUBLE PRECISION - NORTH            ZIT 0002
                  IMPLICIT DOUBLE PRECISION ( A-H, O-Z )                ZIT 0003
                                                                        ZIT 0004
C         ZEROS A REAL ARRAY -A- HAVING -LENGTH- ELEMENTS               ZIT 0005
                                                                        ZIT 0006
      REAL*8  A(*)                                                      ZIT 0007
                                                                        ZIT 0008
      DO 10  L = 1, LENGTH                                              ZIT 0009
         A( L ) = 0.0                                                   ZIT 0010
10    CONTINUE                                                          ZIT 0011
                                                                        ZIT 0012
      RETURN                                                            ZIT 0013
      END                                                               ZIT 0014
