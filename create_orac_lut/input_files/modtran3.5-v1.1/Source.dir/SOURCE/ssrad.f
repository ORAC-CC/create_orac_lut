      SUBROUTINE SSRAD(IPH,IK,MSOFF,IPATH,                              SSR 0001
     1  V,S0,TMOL,TSNOBS,TSNREF,SUMSSS)                                 SSR 0002
C                                                                       SSR 0003
C     SSRAD PERFORMS THE SINGLE SCATTERING SOLAR/LUNAR                  SSR 0004
C     RADIANCE LAYER SUM AND STORES TRANSMITTED SOLAR/LUNAR             SSR 0005
C     IRRADIANCES FOR MULTIPLE SCATTERING CALCULATIONS.                 SSR 0006
C                                                                       SSR 0007
C     DECLARE INPUTS                                                    SSR 0008
C       IPH     SCATTERING PHASE FUNCTION SWITCH                        SSR 0009
C               =2 FOR LOWTRAN MIE FUNCTIONS                            SSR 0010
C       IK      LAYER INDEX                                             SSR 0011
C       MSOFF   LAYER INDEX OFFSET (>0 FOR MULTIPLE SCATTERING PATH)    SSR 0012
C       IPATH   PATH TYPE SWITCH                                        SSR 0013
C               =1 FOR DIRECT SUN TO OBSERVER PATH ONLY                 SSR 0014
C               =2 FOR SUN TO SCATTERING POINT (MSOFF>0)                SSR 0015
C                  FOR SUN TO SCATTERING POINT TO OBSERVER (MSOFF=0)    SSR 0016
C               =3 FOR OBSERVER TO SCATTERING POINT LINE-OF-SIGHT PATH  SSR 0017
C       V       SPECTRAL FREQUENCY [CM-1]                               SSR 0018
C       S0      EXTRA-TERRESTRIAL SOURCE FUNCTION [W/(CM2-MICRON)]      SSR 0019
C       TMOL    SCATTERING POINT TO SUN (MOON) MOLECULAR TRANSMITTANCE. SSR 0020
C               [THE TOTAL TRANSMITTANCE IS TMOL*EXP(-TX(14))]          SSR 0021
      INTEGER IPH,IK,MSOFF,IPATH                                        SSR 0022
      REAL V,S0,TMOL                                                    SSR 0023
C                                                                       SSR 0024
C     DECLARE OUTPUTS                                                   SSR 0025
C       TSNOBS   SOLAR/LUNAR IRRADIANCE AT THE OBSERVER                 SSR 0026
C       TSNREF   SOLAR/LUNAR IRRADIANCE ALONG L-SHAPED PATH             SSR 0027
C       SUMSSS   SOLAR/LUNAR SINGLE SCATTERING RADIANCE SUM             SSR 0028
      REAL TSNOBS,TSNREF,SUMSSS                                         SSR 0029
      INCLUDE 'PARAM.LST'                                               SSR 0030
C                                                                       SSR 0031
C     LIST COMMONS                                                      SSR 0032
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               SSR 0033
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           SSR 0034
      INTEGER KPOINT                                                    SSR 0035
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     SSR 0036
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   SSR 0037
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   SSR 0038
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     SSR 0039
      INTEGER JTURN,LJ                                                  SSR 0040
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       SSR 0041
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  SSR 0042
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      SSR 0043
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    SSR 0044
      INTEGER NCRALT,NCRSPC                                             SSR 0045
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    SSR 0046
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      SSR 0047
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       SSR 0048
      REAL EXTV,ABSV,ASYV                                               SSR 0049
      COMMON/AER/EXTV(NAER),ABSV(NAER),ASYV(NAER)                       SSR 0050
C                                                                       SSR 0051
C     COMMON /MSRD/                                                     SSR 0052
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            SSR 0053
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             SSR 0054
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               SSR 0055
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    SSR 0056
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  SSR 0057
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. SSR 0058
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               SSR 0059
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         SSR 0060
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          SSR 0061
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              SSR 0062
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  SSR 0063
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               SSR 0064
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           SSR 0065
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     SSR 0066
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   SSR 0067
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       SSR 0068
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     SSR 0069
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              SSR 0070
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 SSR 0071
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          SSR 0072
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             SSR 0073
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                SSR 0074
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      SSR 0075
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         SSR 0076
C                                                                       SSR 0077
C       PI       THE CONSTANT PI                                        SSR 0078
C       DEG      NUMBER OF DEGREES IN ONE RADIAN.                       SSR 0079
C       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.                       SSR 0080
C       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.         SSR 0081
      REAL PI,DEG,BIGNUM,BIGEXP                                         SSR 0082
      COMMON/CNSTNS/PI,DEG,BIGNUM,BIGEXP                                SSR 0083
C                                                                       SSR 0084
C     DECLARE FUNCTIONS                                                 SSR 0085
      REAL PHASEF,HENGNS                                                SSR 0086
C                                                                       SSR 0087
C     DECLARE LOCAL VARIABLES                                           SSR 0088
      INTEGER IKP1                                                      SSR 0089
      LOGICAL LWARN                                                     SSR 0090
      REAL PMOLB,PAERB,PCLDB,PCIRB,PMOLF,PAERF,PCLDF,PCIRF,DELDEP,STRANLSSR 0091
C                                                                       SSR 0092
C     SAVE PHASE FUNCTION VARIABLES AND WARNING LOGICAL VARIABLE.       SSR 0093
      SAVE PMOLB,PAERB,PCLDB,PCIRB,PMOLF,PAERF,PCLDF,PCIRF,LWARN        SSR 0094
C                                                                       SSR 0095
C     LIST DATA                                                         SSR 0096
      DATA LWARN/.TRUE./                                                SSR 0097
C                                                                       SSR 0098
C     BRANCH TO PATH TYPE                                               SSR 0099
      GOTO(10,20,30),IPATH                                              SSR 0100
C                                                                       SSR 0101
C     IPATH=1 (DIRECT SUN TO OBSERVER PATH).                            SSR 0102
   10 CONTINUE                                                          SSR 0103
C                                                                       SSR 0104
C     STORE THE TRANSMITTED SOLAR/LUNAR IRRADIANCE TO THE               SSR 0105
C     OBSERVER AND THEN RETURN IF MULTIPLE SCATTERING PATH.             SSR 0106
      IF(TMOL.GT.0.)THEN                                                SSR 0107
          S0DEP(1,1)=TX(14)-LOG(TMOL)                                   SSR 0108
          S0TRN(1,1)=S0*TMOL*EXP(-TX(14))                               SSR 0109
      ELSE                                                              SSR 0110
          S0DEP(1,1)=BIGNUM                                             SSR 0111
          S0TRN(1,1)=0.                                                 SSR 0112
      ENDIF                                                             SSR 0113
      IF(MSOFF.GT.0)RETURN                                              SSR 0114
      TSNOBS=S0TRN(1,1)                                                 SSR 0115
C                                                                       SSR 0116
C     INITIALIZE THE SOLAR/LUNAR SINGLE SCATTERING RADIANCE LAYER SUM.  SSR 0117
C     STORE THE SUN/MOON TO OBSERVER OPTICAL DEPTH.                     SSR 0118
      SUMSSS=0.                                                         SSR 0119
C                                                                       SSR 0120
C     STORE THE RAYLEIGH (PMOLB), THE AEROSOL (PAERB),                  SSR 0121
C     THE CLOUD WATER DROPLET (PCLDB) AND THE CLOUD ICE                 SSR 0122
C     PARTICLE (PCIRB) PHASE FUNCTIONS AT THE OBSERVER.                 SSR 0123
      PMOLB=PHASFN(1,1)                                                 SSR 0124
      PAERB=PHASFN(1,2)                                                 SSR 0125
      PCLDB=PHASFN(1,3)                                                 SSR 0126
      PCIRB=PHASFN(1,4)                                                 SSR 0127
      IF(IPH.EQ.2)PAERB=PHASEF(V,AH1(1),PHASFN(1,2),ARH(1))             SSR 0128
      IF(ABS(ASYMWD).GE.1.)PCLDB=HENGNS(ASYV(6),PHASFN(1,3))            SSR 0129
      IF(ABS(ASYMIP).GE.1.)PCIRB=HENGNS(ASYV(7),PHASFN(1,4))            SSR 0130
      RETURN                                                            SSR 0131
C                                                                       SSR 0132
C     IPATH=2 (SUN TO SCATTERING POINT TO OBSERVER PATH).               SSR 0133
   20 CONTINUE                                                          SSR 0134
C                                                                       SSR 0135
C     STORE THE IPATH=2 TRANSMITTED SOLAR/LUNAR IRRADIANCE              SSR 0136
C     RETURN IF MULTIPLE SCATTERING PATH.  CALCULATE THE                SSR 0137
C     CHANGE IN AND STORE THE CURRENT IPATH=2 OPTICAL DEPTH.            SSR 0138
      IKP1=IK+1                                                         SSR 0139
      IF(TMOL.GT.0.)THEN                                                SSR 0140
          S0DEP(1,IKP1)=TX(14)-LOG(TMOL)                                SSR 0141
          S0TRN(1,IKP1)=S0*TMOL*EXP(-TX(14))                            SSR 0142
      ELSE                                                              SSR 0143
C                                                                       SSR 0144
C         DIVIDE BIGNUM BY BOUNDARY NUMBER SO THAT DEPTH                SSR 0145
C         TO THE SUN DECREASES WITH ALTITUDE IN MULTIPLE                SSR 0146
C         SCATTERING FLUX ADDING ROUTINES BMFLUX AND FLXADD.            SSR 0147
          S0DEP(1,IKP1)=BIGNUM/IKP1                                     SSR 0148
          S0TRN(1,IKP1)=0.                                              SSR 0149
      ENDIF                                                             SSR 0150
      IF(MSOFF.GT.0)THEN                                                SSR 0151
C                                                                       SSR 0152
C         STORE FRACTIONAL DECREASE IN OPTICAL DEPTH ACROSS THE LAYER.  SSR 0153
          DEPRAT(IK)=0.                                                 SSR 0154
          IF(S0DEP(1,IK).LT.S0DEP(1,IKP1))THEN                          SSR 0155
              IF(LWARN)THEN                                             SSR 0156
                  WRITE(IPR,'(/2A,2(/10X,2A,I3,A,1PE10.3),/(10X,2A))')  SSR 0157
     1              ' WARNING:  WEAK-LINE OPTICAL DEPTH TO',            SSR 0158
     2              ' THE SUN IS INCREASING WITH ALTITUDE.',            SSR 0159
     3              ' THE DEPTH TO THE SUN FROM THE BOTTOM',            SSR 0160
     4              ' OF LAYER',IK,' IS ',S0DEP(1,IK),                  SSR 0161
     5              ' THE DEPTH TO THE SUN FROM THE TOP   ',            SSR 0162
     6              ' OF LAYER',IK,' WAS',S0DEP(1,IKP1),                SSR 0163
     7              ' THIS ANOMALY CAN OCCUR BECAUSE OF',               SSR 0164
     8              ' CURVED-EARTH EFFECTS OR BECAUSE',                 SSR 0165
     9              ' OF THE CURTIS-GODSON APPROXIMATION. ',            SSR 0166
     &              ' THE DEPTH FROM LAYER TOP',                        SSR 0167
     &              ' WAS DECREASED TO MATCH THE',                      SSR 0168
     1              ' DEPTH FROM LAYER BOTTOM.',                        SSR 0169
     2              ' ***  THIS WARNING WILL NOT BE REPEATED  ***'      SSR 0170
                  LWARN=.FALSE.                                         SSR 0171
              ENDIF                                                     SSR 0172
              S0DEP(1,IKP1)=S0DEP(1,IK)                                 SSR 0173
              S0TRN(1,IKP1)=S0TRN(1,IK)                                 SSR 0174
          ELSEIF(S0DEP(1,IK).GT.0.)THEN                                 SSR 0175
              DEPRAT(IK)=1.-S0DEP(1,IKP1)/S0DEP(1,IK)                   SSR 0176
          ENDIF                                                         SSR 0177
          RETURN                                                        SSR 0178
      ENDIF                                                             SSR 0179
      TSNREF=S0TRN(1,IKP1)                                              SSR 0180
C                                                                       SSR 0181
C     STORE THE RAYLEIGH (PMOL), AEROSOL (PAER), CLOUD WATER DROPLET    SSR 0182
C     (PCLD), AND CLOUD ICE PARTICLE (PCIR) PHASE FUNCTIONS AT THE      SSR 0183
C     CURRENT (FRONT-SIDE) AND PREVIOUS (BACK-SIDE) SCATTERING POINT.   SSR 0184
      PMOLF=PMOLB                                                       SSR 0185
      PAERF=PAERB                                                       SSR 0186
      PCLDF=PCLDB                                                       SSR 0187
      PCIRF=PCIRB                                                       SSR 0188
      PMOLB=PHASFN(IKP1,1)                                              SSR 0189
      PAERB=PHASFN(IKP1,2)                                              SSR 0190
      PCLDB=PHASFN(IKP1,3)                                              SSR 0191
      PCIRB=PHASFN(IKP1,4)                                              SSR 0192
      IF(IPH.EQ.2)PAERB=PHASEF(V,AH1(IKP1),PHASFN(IKP1,2),ARH(IKP1))    SSR 0193
      IF(ABS(ASYMWD).GE.1.)PCLDB=HENGNS(ASYV(6),PHASFN(IKP1,3))         SSR 0194
      IF(ABS(ASYMIP).GE.1.)PCIRB=HENGNS(ASYV(7),PHASFN(IKP1,4))         SSR 0195
      RETURN                                                            SSR 0196
C                                                                       SSR 0197
C     IPATH=3 (OBSERVER TO SCATTERING POINT).                           SSR 0198
   30 CONTINUE                                                          SSR 0199
C                                                                       SSR 0200
C     NO PREPARATION REQUIRED FOR MULTIPLE SCATTERING PATH.             SSR 0201
      IF(MSOFF.GT.0)RETURN                                              SSR 0202
C                                                                       SSR 0203
C                           1                                           SSR 0204
C                           /                                  Z        SSR 0205
C     STRANL = S0TRN(1,IK)  |  [ S0TRN(1,IKP1) / S0TRN(1,IK) ]    DZ    SSR 0206
C                           /                                           SSR 0207
C                           0                                           SSR 0208
      IKP1=IK+1                                                         SSR 0209
      DELDEP=S0DEP(1,IK)-S0DEP(1,IKP1)                                  SSR 0210
      IF(ABS(DELDEP).LT..001)THEN                                       SSR 0211
          STRANL=.5*(S0TRN(1,IKP1)+S0TRN(1,IK))                         SSR 0212
      ELSE                                                              SSR 0213
          STRANL=(S0TRN(1,IKP1)-S0TRN(1,IK))/DELDEP                     SSR 0214
      ENDIF                                                             SSR 0215
C                                                                       SSR 0216
C     SUM THE SOLAR/LUNAR SINGLE SCATTER RADIANCE CONTRIBUTIONS.        SSR 0217
      SUMSSS=SUMSSS+STRANL*.5*(TX(2)*(PAERB+PAERF)+TX(66)*(PCLDB+PCLDF) SSR 0218
     1                       +TX(67)*(PCIRB+PCIRF)+TX(15)*(PMOLB+PMOLF))SSR 0219
      RETURN                                                            SSR 0220
      END                                                               SSR 0221
