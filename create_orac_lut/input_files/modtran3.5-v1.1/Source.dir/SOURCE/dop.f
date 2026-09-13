      FUNCTION   DOP(WAVL,A,CEN1,B,C,CEN2,D,E,CEN3,F,G)                 DOP 0001
CCC                                                                     DOP 0002
CCC    DESCRIBES THE REAL PART OF THE DIELECTRIC CONSTANT               DOP 0003
CCC                                                                     DOP 0004
      V=1./WAVL                                                         DOP 0005
      V2=V*V                                                            DOP 0006
      H1=CEN1**2-V2                                                     DOP 0007
      H2=CEN2**2-V2                                                     DOP 0008
      H3=CEN3**2-V2                                                     DOP 0009
      DOP=SQRT(A+B*H1/(H1*H1+C*V2)+D*H2/(H2*H2+E*V2)+F*H3/(H3*H3+G*V2)) DOP 0010
      RETURN                                                            DOP 0011
      END                                                               DOP 0012
