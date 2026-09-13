      SUBROUTINE RAIN(V,EXT,ABT,SCT,ASYMR)                              RAI 0001
      INTEGER PHASE,DIST                                                RAI 0002
      INCLUDE 'PARAM.LST'                                               RAI 0003
      INTEGER KPOINT                                                    RAI 0004
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     RAI 0005
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   RAI 0006
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   RAI 0007
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     RAI 0008
C                                                                       RAI 0009
C       PI       THE CONSTANT PI                                        RAI 0010
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       RAI 0011
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       RAI 0012
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         RAI 0013
      REAL PI,DEG,BIGNUM,BIGEXP                                         RAI 0014
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                RAI 0015
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,     RAI 0016
     1  RAINRT                                                          RAI 0017
      TMPRN=W(61)/W(3)                                                  RAI 0018
      VTEMP=1.438786*V                                                  RAI 0019
      IF(VTEMP.GE.TMPRN*BIGEXP)THEN                                     RAI 0020
          RFD=V                                                         RAI 0021
      ELSE                                                              RAI 0022
          STORE=EXP(-VTEMP/TMPRN)                                       RAI 0023
          RFD=V*(1.-STORE)/(1.+STORE)                                   RAI 0024
      ENDIF                                                             RAI 0025
      RAINAV=(W(3)/W(62))**(1./.63)                                     RAI 0026
      IF(V.LT.250.)THEN                                                 RAI 0027
          PHASE=1                                                       RAI 0028
          IF(ICLD.GT.11)PHASE=2                                         RAI 0029
          DIST=1                                                        RAI 0030
C                                                                       RAI 0031
C         CALL SCATTERING ROUTINE TO OBTAIN ASYMMETRY                   RAI 0032
C         FACTOR AND RATIO OF ABSORPTION TO EXTINCTION DUE              RAI 0033
C         TO RAIN WITHIN RANGE OF 19 TO 231 GHZ.                        RAI 0034
C         EXTRAPOLATE ABOVE AND BELOW THAT FREQ RANGE                   RAI 0035
          CALL RNSCAT(V,RAINAV,TMPRN,PHASE,DIST,CSSA,ASYMR)             RAI 0036
      ELSE                                                              RAI 0037
          CSSA=.5                                                       RAI 0038
          ASYMR=.85                                                     RAI 0039
      ENDIF                                                             RAI 0040
      RNEXPH=TNRAIN(RAINAV,V,TMPRN,RFD)*W(62)                           RAI 0041
      EXT=RNEXPH                                                        RAI 0042
      ABT=RNEXPH*CSSA                                                   RAI 0043
      SCT=RNEXPH*(1.-CSSA)                                              RAI 0044
      RETURN                                                            RAI 0045
      END                                                               RAI 0046
