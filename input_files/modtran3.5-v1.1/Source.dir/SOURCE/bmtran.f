      REAL FUNCTION BMTRAN(DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2,DV)           BMT 0001
C                                                                       BMT 0002
C     THIS FUNCTION RETURNS PATH TRANSMITTANCE AVERAGED                 BMT 0003
C     OVER A SPECTRAL INTERVAL OF WIDTH DV.                             BMT 0004
C                                                                       BMT 0005
C     OUTPUTS                                                           BMT 0006
C       BMTRAN  =   PATH TRANSMITTANCE [DIMENSIONLESS]                  BMT 0007
C                                                                       BMT 0008
C     INPUTS                                                            BMT 0009
C                    _                                                  BMT 0010
C       DEPTH   =   >_  SU/D,     THE ABSORPTION COEFFICIENT (S/D)      BMT 0011
C                                 AND COLUMN DENSITY (U) PRODUCT,       BMT 0012
C                                 SUMMED OVER THE PATH [DIMENSIONLESS]  BMT 0013
C                                                                       BMT 0014
C       ODBAR   =   <1/D>,        THE RECIPROCAL OF THE LINE            BMT 0015
C                                 SPACING, (D), PATH AVERAGED [CM]      BMT 0016
C                                                                       BMT 0017
C       ADBAR   =   <ALFDOP/D>,   THE DOPPLER LINE WIDTH (ALFDOP)       BMT 0018
C                                 OVER THE AVERAGE LINE SPACING (D),    BMT 0019
C                                 PATH AVERAGED [DIMENSIONLESS]         BMT 0020
C                                                                       BMT 0021
C       ACBAR   =   <ALFCOL/D>,   THE COLLISION (LORENTZ) LINE WIDTH    BMT 0022
C                                 (ALFCOL) OVER THE AVERAGE LINE SPACINGBMT 0023
C                                 (D), PATH AVERAGED [DIMENSIONLESS]    BMT 0024
C                                                                       BMT 0025
C                          2  2                                         BMT 0026
C       ACBAR2  =   <ALFCOL /D >, THE COLLISION (LORENTZ) LINE WIDTH    BMT 0027
C                                 (ALFCOL) SQUARED OVER THE AVERAGE     BMT 0028
C                                 LINE SPACING (D) SQUARED, PATH        BMT 0029
C                                 AVERAGED [DIMENSIONLESS] (OZONE ONLY) BMT 0030
C                                                                       BMT 0031
C       DV      =                 THE INTERVAL WIDTH [CM-1]             BMT 0032
      REAL DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2,DV                            BMT 0033
C                                                                       BMT 0034
C     DECLARE LOCAL VARIABLES                                           BMT 0035
      REAL WDW,STORE,RATIO,RHO,RHOP1,RHO2,RHO2P1,DENOM,                 BMT 0036
     1  F1,F2,F3,STR2PI,DOPTRM,WDRATD,WDRATL,WDV,ARG,WDFRAC             BMT 0037
C                                                                       BMT 0038
C     THE RODGERS AND WILLIAMS APPROXIMATION IS USED TO                 BMT 0039
C     CALCULATE THE TOTAL VOIGT EQUIVALENT WIDTH (WDV):                 BMT 0040
C                                                                       BMT 0041
C                2            2            2            2          2    BMT 0042
C       (WDV/WDW)  = (WDL/WDW)  + (WDD/WDW)  - (WDL/WDW)  (WDD/WDW)     BMT 0043
C                                                                       BMT 0044
C     WHERE WDW IS THE WEAK-LINE EQUIVALENT WIDTH,                      BMT 0045
C           WDL IS THE LORENTZ EQUIVALENT WIDTH, AND                    BMT 0046
C           WDD IS THE DOPPLER EQUIVALENT WIDTH.                        BMT 0047
C                                                                       BMT 0048
C     THE APPROXIMATIONS FOR THESE EQUIVALENT WIDTHS ARE:               BMT 0049
C                                       _                               BMT 0050
C          WDW   =   <SU>         =  [ >_ SU/D ] / <1/D>                BMT 0051
C                                                                       BMT 0052
C                             2                                         BMT 0053
C       WDRATL   =   (WDL/WDW)   =  4 / (4 + WDW/<ALFCOL>)              BMT 0054
C                                                                       BMT 0055
C                             2                                         BMT 0056
C       WDRATD   =   (WDD/WDW)   =  LN (1 + DOPTRM) / DOPTRM            BMT 0057
C                                                                       BMT 0058
C     WHERE                                                             BMT 0059
C                                                   2                   BMT 0060
C               DOPTRM  =  [LN(2)/2]  (WDW/<ALFDOP>)                    BMT 0061
C                                                                       BMT 0062
C     THE HALF-WIDTHS ARE CALCULATED FROM THE FOLLOWING EXPRESSIONS:    BMT 0063
C                                                                 _     BMT 0064
C       <ALFCOL>/WDW = [<ALFCOL/D> / <1/D>] / WDW = <ALFCOL/D> / >_ SU/DBMT 0065
C                                                                 _     BMT 0066
C       <ALFDOP>/WDW = [<ALFDOP/D> / <1/D>] / WDW = <ALFDOP/D> / >_ SU/DBMT 0067
C                                                                       BMT 0068
C     FOR OZONE AN IMPROVED APPROXIMATION FOR LORENTZ EQUIVALENT WIDTH  BMT 0069
C     IS REQUIRED BASED ON THE WORK OF R. M. GOODY, WITH REFINEMENTS    BMT 0070
C     BY L. S. BERNSTEIN TO ELIMINATE DECREASES IN THE CURVE-OF-GROWTH  BMT 0071
C                                                                       BMT 0072
C     LIST DATA                                                         BMT 0073
C        PT5LN2  =  .5 * LN(2)                                          BMT 0074
      REAL PT5LN2,PI,TWOPI                                              BMT 0075
      DATA PT5LN2,PI,TWOPI/.34657359,3.1415927,6.2831853/               BMT 0076
C                                                                       BMT 0077
C     WEAK-LINE EQUIVALENT WIDTH [CM-1]                                 BMT 0078
      WDW=DEPTH/ODBAR                                                   BMT 0079
C                                                                       BMT 0080
C     CALCULATE DENOMINATOR IN THE LORENTZ EQUIVALENT WIDTH EXPRESSION  BMT 0081
      STORE=DEPTH/ACBAR                                                 BMT 0082
      RATIO=1.                                                          BMT 0083
      IF(ACBAR2.NE.0.)THEN                                              BMT 0084
C                                                                       BMT 0085
C         REFINEMENT FOR OZONE ONLY                                     BMT 0086
          RHO=ACBAR2/ACBAR**2                                           BMT 0087
          RHOP1=RHO+1.                                                  BMT 0088
          RHO2=RHO**2                                                   BMT 0089
          RHO2P1=RHO2+1.                                                BMT 0090
          DENOM=RHOP1*(RHO2+RHOP1)                                      BMT 0091
          F1=RHO*(RHO-1)*RHO2P1/DENOM                                   BMT 0092
          DENOM=DENOM*RHOP1**2                                          BMT 0093
          F2=2.*RHO2*RHO2P1/DENOM                                       BMT 0094
C                                                                       BMT 0095
C         THE F3 EQUATION HAS BEEN IMPROVED (THE CURVE-OF-GROWTH        BMT 0096
C         WAS DECREASING FOR SOME DOWN LOOKING PATHS ON THE WING        BMT 0097
C         OF THE 9.6 MICRON O3 BAND; THIS NO LONGER OCCURS).            BMT 0098
          F3=.7200*RHO-.5314                                            BMT 0099
          STR2PI=STORE/TWOPI                                            BMT 0100
          RATIO=1.-.5*STR2PI*F1/(1.+STR2PI*(F3+STR2PI*F2))              BMT 0101
          STORE=RATIO*STORE                                             BMT 0102
      ENDIF                                                             BMT 0103
      DENOM=4.+STORE                                                    BMT 0104
C                                                                       BMT 0105
C     CALCULATE TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH         BMT 0106
      DOPTRM=PT5LN2*(DEPTH/ADBAR)**2                                    BMT 0107
C                                                                       BMT 0108
C     CALCULATE VOIGT EQUIVALENT WIDTH (RODGERS-WILLIAMS APPROXIMATION) BMT 0109
      IF(DOPTRM.GT..01)THEN                                             BMT 0110
C                                                                       BMT 0111
C         STANDARD EXPRESSIONS                                          BMT 0112
          WDRATD=LOG(1.+DOPTRM)/DOPTRM                                  BMT 0113
          WDRATL=4./DENOM*RATIO**2                                      BMT 0114
          WDV=WDW*SQRT(WDRATD+WDRATL-WDRATD*WDRATL)                     BMT 0115
      ELSEIF(DOPTRM.GT..0001)THEN                                       BMT 0116
C                                                                       BMT 0117
C         IF DOPTRM IS SMALL, WDRATD IS NEAR ONE.  REPLACE              BMT 0118
C         THE LOG AND THE SQRT WITH POWER SERIES EXPANSIONS.            BMT 0119
          STORE=DOPTRM*(.25-DOPTRM*(.16666667-.125*DOPTRM))*            BMT 0120
     1      (STORE+(1-RATIO)*(1+RATIO))/DENOM                           BMT 0121
          WDV=WDW*(1.-STORE*(1.+.5*STORE*(1.+STORE)))                   BMT 0122
      ELSE                                                              BMT 0123
C                                                                       BMT 0124
C         IF DOPTRM IS VERY SMALL, TRUNCATED EXPANSIONS SUFFICE.        BMT 0125
          WDV=WDW*(1.-.25*DOPTRM*(STORE+(1-RATIO)*(1+RATIO))/DENOM)     BMT 0126
      ENDIF                                                             BMT 0127
C                                                                       BMT 0128
C     SUBTRACT LINE TAIL CONTRIBUTIONS ASSUMING LORENTZIAN TAILS        BMT 0129
C     AND ASSUMING THE LINE IS CENTERED .2 DV FROM BIN EDGE.            BMT 0130
      ARG=SQRT(WDW*RATIO*ACBAR/(PI*ODBAR))/DV                           BMT 0131
      WDFRAC=WDV/DV-(.2*BMERFU(5.*ARG)+.8*BMERFU(1.25*ARG))             BMT 0132
C                                                                       BMT 0133
C     CALCULATE TRANSMITTANCE USING POWER LAW EXPRESSION                BMT 0134
      BMTRAN=0.                                                         BMT 0135
      IF(WDFRAC.LT.1.)BMTRAN=(1.-WDFRAC)**(DV*ODBAR)                    BMT 0136
C                                                                       BMT 0137
C     RETURN TO BMOD                                                    BMT 0138
      RETURN                                                            BMT 0139
      END                                                               BMT 0140
      REAL FUNCTION BMTRN(DEPTH,ODBAR,ADBAR,ACBAR,DV)                   BMT 0141
C                                                                       BMT 0142
C     THIS FUNCTION RETURNS PATH TRANSMITTANCE AVERAGED                 BMT 0143
C     OVER A SPECTRAL INTERVAL OF WIDTH DV.                             BMT 0144
C                                                                       BMT 0145
C     OUTPUTS                                                           BMT 0146
C       BMTRN  =   PATH TRANSMITTANCE [DIMENSIONLESS]                   BMT 0147
C                                                                       BMT 0148
C     INPUTS                                                            BMT 0149
C                    _                                                  BMT 0150
C       DEPTH   =   >_  SU/D,     THE ABSORPTION COEFFICIENT (S/D)      BMT 0151
C                                 AND COLUMN DENSITY (U) PRODUCT,       BMT 0152
C                                 SUMMED OVER THE PATH [DIMENSIONLESS]  BMT 0153
C                                                                       BMT 0154
C       ODBAR   =   <1/D>,        THE RECIPROCAL OF THE LINE            BMT 0155
C                                 SPACING, (D), PATH AVERAGED [CM]      BMT 0156
C                                                                       BMT 0157
C       ADBAR   =   <ALFDOP/D>,   THE DOPPLER LINE WIDTH (ALFDOP)       BMT 0158
C                                 OVER THE AVERAGE LINE SPACING (D),    BMT 0159
C                                 PATH AVERAGED [DIMENSIONLESS]         BMT 0160
C                                                                       BMT 0161
C       ACBAR   =   <ALFCOL/D>,   THE COLLISION (LORENTZ) LINE WIDTH    BMT 0162
C                                 (ALFCOL) OVER THE AVERAGE LINE SPACINGBMT 0163
C                                 (D), PATH AVERAGED [DIMENSIONLESS]    BMT 0164
C                                                                       BMT 0165
C       DV      =                 THE INTERVAL WIDTH [CM-1]             BMT 0166
      REAL DEPTH,ODBAR,ADBAR,ACBAR,DV                                   BMT 0167
C                                                                       BMT 0168
C     DECLARE LOCAL VARIABLES                                           BMT 0169
      REAL WDW,STORE,DENOM,DOPTRM,WDRATD,WDRATL,WDV,ARG,WDFRAC          BMT 0170
C                                                                       BMT 0171
C     THE RODGERS AND WILLIAMS APPROXIMATION IS USED TO                 BMT 0172
C     CALCULATE THE TOTAL VOIGT EQUIVALENT WIDTH (WDV):                 BMT 0173
C                                                                       BMT 0174
C                2            2            2            2          2    BMT 0175
C       (WDV/WDW)  = (WDL/WDW)  + (WDD/WDW)  - (WDL/WDW)  (WDD/WDW)     BMT 0176
C                                                                       BMT 0177
C     WHERE WDW IS THE WEAK-LINE EQUIVALENT WIDTH,                      BMT 0178
C           WDL IS THE LORENTZ EQUIVALENT WIDTH, AND                    BMT 0179
C           WDD IS THE DOPPLER EQUIVALENT WIDTH.                        BMT 0180
C                                                                       BMT 0181
C     THE APPROXIMATIONS FOR THESE EQUIVALENT WIDTHS ARE:               BMT 0182
C                                       _                               BMT 0183
C          WDW   =   <SU>         =  [ >_ SU/D ] / <1/D>                BMT 0184
C                                                                       BMT 0185
C                             2                                         BMT 0186
C       WDRATL   =   (WDL/WDW)   =  4 / (4 + WDW/<ALFCOL>)              BMT 0187
C                                                                       BMT 0188
C                             2                                         BMT 0189
C       WDRATD   =   (WDD/WDW)   =  LN (1 + DOPTRM) / DOPTRM            BMT 0190
C                                                                       BMT 0191
C     WHERE                                                             BMT 0192
C                                                   2                   BMT 0193
C               DOPTRM  =  [LN(2)/2]  (WDW/<ALFDOP>)                    BMT 0194
C                                                                       BMT 0195
C     THE HALF-WIDTHS ARE CALCULATED FROM THE FOLLOWING EXPRESSIONS:    BMT 0196
C                                                                 _     BMT 0197
C       <ALFCOL>/WDW = [<ALFCOL/D> / <1/D>] / WDW = <ALFCOL/D> / >_ SU/DBMT 0198
C                                                                 _     BMT 0199
C       <ALFDOP>/WDW = [<ALFDOP/D> / <1/D>] / WDW = <ALFDOP/D> / >_ SU/DBMT 0200
C                                                                       BMT 0201
C     FOR OZONE AN IMPROVED APPROXIMATION FOR LORENTZ EQUIVALENT WIDTH  BMT 0202
C     IS REQUIRED BASED ON THE WORK OF R. M. GOODY, WITH REFINEMENTS    BMT 0203
C     BY L. S. BERNSTEIN TO ELIMINATE DECREASES IN THE CURVE-OF-GROWTH  BMT 0204
C                                                                       BMT 0205
C     LIST DATA                                                         BMT 0206
C        PT5LN2  =  .5 * LN(2)                                          BMT 0207
      REAL PT5LN2,PI                                                    BMT 0208
      DATA PT5LN2,PI/.34657359,3.1415927/                               BMT 0209
C                                                                       BMT 0210
C     WEAK-LINE EQUIVALENT WIDTH [CM-1]                                 BMT 0211
      WDW=DEPTH/ODBAR                                                   BMT 0212
C                                                                       BMT 0213
C     CALCULATE DENOMINATOR IN THE LORENTZ EQUIVALENT WIDTH EXPRESSION  BMT 0214
      STORE=DEPTH/ACBAR                                                 BMT 0215
      DENOM=4.+STORE                                                    BMT 0216
C                                                                       BMT 0217
C     CALCULATE TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH         BMT 0218
      DOPTRM=PT5LN2*(DEPTH/ADBAR)**2                                    BMT 0219
C                                                                       BMT 0220
C     CALCULATE VOIGT EQUIVALENT WIDTH (RODGERS-WILLIAMS APPROXIMATION) BMT 0221
      IF(DOPTRM.GT..01)THEN                                             BMT 0222
C                                                                       BMT 0223
C         STANDARD EXPRESSIONS                                          BMT 0224
          WDRATD=LOG(1.+DOPTRM)/DOPTRM                                  BMT 0225
          WDRATL=4./DENOM                                               BMT 0226
          WDV=WDW*SQRT(WDRATD+WDRATL-WDRATD*WDRATL)                     BMT 0227
      ELSEIF(DOPTRM.GT..0001)THEN                                       BMT 0228
C                                                                       BMT 0229
C         IF DOPTRM IS SMALL, WDRATD IS NEAR ONE.  REPLACE              BMT 0230
C         THE LOG AND THE SQRT WITH POWER SERIES EXPANSIONS.            BMT 0231
          STORE=DOPTRM*(.25-DOPTRM*(.16666667-.125*DOPTRM))*            BMT 0232
     1      STORE/DENOM                                                 BMT 0233
          WDV=WDW*(1.-STORE*(1.+.5*STORE*(1.+STORE)))                   BMT 0234
      ELSE                                                              BMT 0235
C                                                                       BMT 0236
C         IF DOPTRM IS VERY SMALL, TRUNCATED EXPANSIONS SUFFICE.        BMT 0237
          WDV=WDW*(1.-.25*DOPTRM*STORE/DENOM)                           BMT 0238
      ENDIF                                                             BMT 0239
C                                                                       BMT 0240
C     SUBTRACT LINE TAIL CONTRIBUTIONS ASSUMING LORENTZIAN TAILS        BMT 0241
C     AND ASSUMING THE LINE IS CENTERED .2 DV FROM BIN EDGE.            BMT 0242
      ARG=SQRT(WDW*ACBAR/(PI*ODBAR))/DV                                 BMT 0243
      WDFRAC=WDV/DV-(.2*BMERFU(5.*ARG)+.8*BMERFU(1.25*ARG))             BMT 0244
C                                                                       BMT 0245
C     CALCULATE TRANSMITTANCE USING POWER LAW EXPRESSION                BMT 0246
      BMTRN=0.                                                          BMT 0247
      IF(WDFRAC.LT.1.)BMTRN=(1.-WDFRAC)**(DV*ODBAR)                     BMT 0248
C                                                                       BMT 0249
C     RETURN TO BMOD                                                    BMT 0250
      RETURN                                                            BMT 0251
      END                                                               BMT 0252
