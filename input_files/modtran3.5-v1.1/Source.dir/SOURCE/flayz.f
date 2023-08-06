      SUBROUTINE FLAYZ(ML,ICLD,GNDALT,IVSA)                             FLA 0001
C                                                                       FLA 0002
C     THIS ROUTINE DETERMINES THE LAYER BOUNDARIES WHEN                 FLA 0003
C     ONE OF THE MODEL ATMOSPHERES IS SELECTED.                         FLA 0004
C                                                                       FLA 0005
C     INCLUDE PARAMETERS.                                               FLA 0006
      INCLUDE 'PARAM.LST'                                               FLA 0007
      INTEGER NDFLT                                                     FLA 0008
      PARAMETER(NDFLT=36)                                               FLA 0009
C                                                                       FLA 0010
C     DECLARE INPUTS.                                                   FLA 0011
      INTEGER ML,ICLD,IVSA                                              FLA 0012
      REAL GNDALT                                                       FLA 0013
C                                                                       FLA 0014
C     LIST COMMONS.                                                     FLA 0015
      REAL ZM,PM,TM,RFNDX,DENSTY                                        FLA 0016
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    FLA 0017
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               FLA 0018
      INTEGER NCRALT,NCRSPC                                             FLA 0019
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP    FLA 0020
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      FLA 0021
     1  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP                       FLA 0022
      REAL ZVSA,RHVSA,AHVSA                                             FLA 0023
      INTEGER IHVSA                                                     FLA 0024
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)              FLA 0025
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               FLA 0026
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           FLA 0027
C                                                                       FLA 0028
C     DECLARE LOCAL VARIABLES.                                          FLA 0029
      INTEGER I,ILAST,INEXT,NEXT,IHI,J,IM1,ICIR(4)                      FLA 0030
      REAL DALT,CLDD,TOL,CIRRUS(4)                                      FLA 0031
C                                                                       FLA 0032
C     DECLARE DATA                                                      FLA 0033
      REAL ZDFLT(NDFLT)                                                 FLA 0034
      DATA ZDFLT/ 0., 1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,       FLA 0035
     1           12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,       FLA 0036
     2           24.,25.,30.,35.,40.,45.,50.,55.,60.,70.,80.,100./      FLA 0037
      IF(IVSA.EQ.1)THEN                                                 FLA 0038
C                                                                       FLA 0039
C         VERTICAL STRUCTURE ALGORITHM (VSA) OPTION IS ON.              FLA 0040
C         FIRST 9 BOUNDARY ALTITUDES WERE DEFINED IN ROUTINE VSA.       FLA 0041
          DO 10 I=1,9                                                   FLA 0042
              ZM(I)=ZVSA(I)                                             FLA 0043
   10     CONTINUE                                                      FLA 0044
C                                                                       FLA 0045
C         ADD A BOUNDARY ALTITUDE 10 METERS ABOVE ZVSA(9).              FLA 0046
          ML=10                                                         FLA 0047
          ZM(10)=ZM(9)+.01                                              FLA 0048
C                                                                       FLA 0049
C         FIND THE FIRST DEFAULT ALTITUDE HALF A METER ABOVE ZM(10).    FLA 0050
          ILAST=0                                                       FLA 0051
          DO 20 INEXT=1,NDFLT                                           FLA 0052
              IF(ZDFLT(INEXT).GE.ZM(10)+.0005)GOTO30                    FLA 0053
   20     ILAST=INEXT                                                   FLA 0054
C                                                                       FLA 0055
C         NO DEFAULT ALTITUDES ABOVE ZM(10)+.0005.                      FLA 0056
          RETURN                                                        FLA 0057
C                                                                       FLA 0058
C         ADD DEFAULT ALTITUDES.                                        FLA 0059
   30     CONTINUE                                                      FLA 0060
          IF(ML+NDFLT-ILAST.GT.LAYDIM)THEN                              FLA 0061
              WRITE(IPR,'(2A,I3,/8X,2A,I3,A)')                          FLA 0062
     1          ' ERROR:  VERTICAL STRUCTURE ALGORITHM INCREASES',      FLA 0063
     2          ' LAYER BOUNDARY NUMBER ABOVE',LAYDIM,' INCREASE',      FLA 0064
     3          ' PARAMETER LAYDIM TO',ML+NDFLT-ILAST,' IN PARAM.LST'   FLA 0065
               STOP 'INCREASE PARAMETER LAYDIM'                         FLA 0066
          ENDIF                                                         FLA 0067
          DO 40 I=INEXT,NDFLT                                           FLA 0068
              ML=ML+1                                                   FLA 0069
              ZM(ML)=ZM(I)                                              FLA 0070
   40     CONTINUE                                                      FLA 0071
          RETURN                                                        FLA 0072
      ENDIF                                                             FLA 0073
C                                                                       FLA 0074
C     ASSIGN DEFAULT ALTITUDES                                          FLA 0075
      ML=NDFLT                                                          FLA 0076
      DO 50 I=1,NDFLT                                                   FLA 0077
          ZM(I)=ZDFLT(I)                                                FLA 0078
   50 CONTINUE                                                          FLA 0079
C                                                                       FLA 0080
C     IF THE INPUT GROUND ALTITUDE IS LESS THAN 6 KM ABOVE              FLA 0081
C     SEA LEVEL, SCALE BOTTOM 5 LAYER BOUNDARIES.                       FLA 0082
      IF(GNDALT.LT.ZM(7))THEN                                           FLA 0083
          DALT=(ZM(7)-GNDALT)/(ZM(7)-ZM(1))                             FLA 0084
cjv 1/10/97  Lex fix. 'If a ground altitude between 0 and 6 km is used along
cjv          with one of the model atmospheres, the atmosphere was layered
cjv          erroneously.  This simple fix eliminates the problem.'
c
cjv          DO 60 I=1,6                                                   FLA 0085
          DO 60 I=6,1,-1                                                   FLA 0085
cjv ^
              ZM(I)=GNDALT+DALT*(ZM(I)-ZM(1))                           FLA 0086
   60     CONTINUE                                                      FLA 0087
      ENDIF                                                             FLA 0088
      IF(ICLD.NE.18 .AND. ICLD.NE.19)RETURN                             FLA 0089
C                                                                       FLA 0090
C     CIRRUS CLOUD.  DEFINE CLOUD BOUNDARIES.                           FLA 0091
      CLDD=.1*CTHIK                                                     FLA 0092
      CIRRUS(1)=CALT-.5*CLDD                                            FLA 0093
      IF(CIRRUS(1).LE.GNDALT)CIRRUS(1)=GNDALT                           FLA 0094
      CIRRUS(2)=CIRRUS(1)+CLDD                                          FLA 0095
      CIRRUS(3)=CIRRUS(1)+CTHIK                                         FLA 0096
      CIRRUS(4)=CIRRUS(3)+CLDD                                          FLA 0097
      IF(CIRRUS(4).GE.ZM(NDFLT))THEN                                    FLA 0098
          WRITE(IPR,'(2A,/8X,2A)')' ERROR:  CIRRUS LAYER',              FLA 0099
     1      ' IS ABOVE THE TOP OF THE MODEL ATMOSPHERE.',               FLA 0100
     2      ' IF YOU REALLY WANT A CLOUD AT THIS ALTITUDE,',            FLA 0101
     3      ' INPUT A USER-DEFINED ATMOSPHERE.'                         FLA 0102
          STOP 'CIRRUS CLOUD ABOVE TOP OF ATMOSPHERE'                   FLA 0103
      ENDIF                                                             FLA 0104
C                                                                       FLA 0105
C     CHECK ZM DIMENSION                                                FLA 0106
      ML=NDFLT+4                                                        FLA 0107
      IF(LAYDIM.LT.ML)THEN                                              FLA 0108
          WRITE(IPR,'(3A,I3,/14X,A,I3,A)')' ERROR:  CIRRUS',            FLA 0109
     1      ' LAYER INCREASES LAYER BOUNDARY NUMBER ABOVE',LAYDIM,      FLA 0110
     2      ' INCREASE PARAMETER LAYDIM TO',ML,' IN PARAM.LST'          FLA 0111
          STOP 'INCREASE PARAMETER LAYDIM'                              FLA 0112
      ENDIF                                                             FLA 0113
C                                                                       FLA 0114
C     COMBINE DEFAULT AND CIRRUS LAYERS STARTING FROM THE TOP.          FLA 0115
      NEXT=ML                                                           FLA 0116
      IHI=NDFLT                                                         FLA 0117
      DO 90 J=4,1,-1                                                    FLA 0118
          DO 70 I=IHI,1,-1                                              FLA 0119
              IF(CIRRUS(J).GE.ZM(I))GOTO80                              FLA 0120
              ZM(NEXT)=ZM(I)                                            FLA 0121
              NEXT=NEXT-1                                               FLA 0122
   70     CONTINUE                                                      FLA 0123
   80     CONTINUE                                                      FLA 0124
          ZM(NEXT)=CIRRUS(J)                                            FLA 0125
          ICIR(J)=NEXT                                                  FLA 0126
          NEXT=NEXT-1                                                   FLA 0127
          IHI=I                                                         FLA 0128
   90 CONTINUE                                                          FLA 0129
C                                                                       FLA 0130
C     SET TOLERANCE FOR MERGING LAYERS                                  FLA 0131
      TOL=MIN(.5*CLDD,.0005)                                            FLA 0132
C                                                                       FLA 0133
C     CHECK FOR MERGING                                                 FLA 0134
      DO 120 J=4,1,-1                                                   FLA 0135
          IF(ZM(ICIR(J)+1)-ZM(ICIR(J)).LT.TOL)THEN                      FLA 0136
              IM1=ICIR(J)+1                                             FLA 0137
              DO 100 I=ICIR(J)+2,ML                                     FLA 0138
                  ZM(IM1)=ZM(I)                                         FLA 0139
  100         IM1=I                                                     FLA 0140
              ML=ML-1                                                   FLA 0141
          ELSEIF(ZM(ICIR(J))-ZM(ICIR(J)-1).LT.TOL)THEN                  FLA 0142
              IM1=ICIR(J)-1                                             FLA 0143
              DO 110 I=ICIR(J),ML                                       FLA 0144
                  ZM(IM1)=ZM(I)                                         FLA 0145
  110         IM1=I                                                     FLA 0146
              ML=ML-1                                                   FLA 0147
          ENDIF                                                         FLA 0148
  120 CONTINUE                                                          FLA 0149
      RETURN                                                            FLA 0150
      END                                                               FLA 0151
