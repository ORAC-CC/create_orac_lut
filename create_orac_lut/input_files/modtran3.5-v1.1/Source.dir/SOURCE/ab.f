      FUNCTION   AB(WAVL,A,CEN,B,C)                                     AB  0001
CCC                                                                     AB  0002
CCC    DESCRIBES THE IMAGINARY PART OF THE DIELECTRIC CONSTANT          AB  0003
CCC                                                                     AB  0004
      AB=-A*EXP(-ABS((ALOG10(10000.*WAVL/CEN)/B))**C)                   AB  0005
      RETURN                                                            AB  0006
      END                                                               AB  0007
