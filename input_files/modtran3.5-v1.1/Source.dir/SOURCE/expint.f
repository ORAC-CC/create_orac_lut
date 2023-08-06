      REAL FUNCTION EXPINT(X1,X2,A)                                     EXP 0001
C                                                                       EXP 0002
C     EXPONENTIAL INTERPOLATION FOR POSITIVE ARGUMENTS                  EXP 0003
      REAL X1,X2,A                                                      EXP 0004
      IF(X1.GT.0. .AND. X2.GT.0.)THEN                                   EXP 0005
          EXPINT=X1*(X2/X1)**A                                          EXP 0006
      ELSE                                                              EXP 0007
          EXPINT=X1+(X2-X1)*A                                           EXP 0008
      ENDIF                                                             EXP 0009
      RETURN                                                            EXP 0010
      END                                                               EXP 0011
