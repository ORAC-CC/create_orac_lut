      SUBROUTINE BMDATA(IV1,IFWHM,IDVX,IKMX,MXFREQ,IMSMX)               BMD 0001
C                                                                       BMD 0002
C     BMDATA (CALLED BY TRANS) MAKES THE INITIAL BAND MODEL TAPE READ   BMD 0003
C     AND CALCULATES WAVENUMBER-INDEPENDENT PARAMETERS FOR USE BY BMOD  BMD 0004
C                                                                       BMD 0005
C     CONVENTION                                                        BMD 0006
C     MMOLX = MAXIMUM NUMBER OF NEW SPECIES (IDENTIFIED BY "X")         BMD 0007
C     MMOL  = MAXIMUM NUMBER OF OLD SPECIES (W/O SUFFIX "X")            BMD 0008
C     THESE DEFINE THE MAXIMUM ARRAY SIZES.                             BMD 0009
C                                                                       BMD 0010
C     THE ACTUAL NUMBER OF PARAMETERS ARE:                              BMD 0011
C     NSPC = ACTUAL NUMBER OF OLD SPECIES (12), CAN'T EXCEED MMOL       BMD 0012
C     NSPECX = ACTUAL NUMBER OF "X" SPECIES,     CAN'T EXCEED MMOLX     BMD 0013
C                                                                       BMD 0014
C     PARAMETER KMAX DENOTES THE NUMBER OF MODTRAN "SPECIES".           BMD 0015
C     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES      BMD 0016
C     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.        BMD 0017
                                                                        BMD 0018
      INCLUDE 'PARAM.LST'                                               BMD 0019
C                                                                       BMD 0020
C     TRANS VARIABLES                                                   BMD 0021
      INTEGER IBINX,IMOLX,IALFX                                         BMD 0022
      REAL SDZX,ODZX                                                    BMD 0023
      COMMON/BMDCMX/SDZX(MXTEMP),ODZX(MXTEMP),IBINX,IMOLX,IALFX         BMD 0024
      COMMON /SOLSX/  WPTHSX(LAYTHR,MMOLX),TBBYSX(LAYTHR,MMOLX)         BMD 0025
      LOGICAL MODTRN                                                    BMD 0026
      INTEGER KPOINT                                                    BMD 0027
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     BMD 0028
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   BMD 0029
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   BMD 0030
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     BMD 0031
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               BMD 0032
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           BMD 0033
      INTEGER JTURN,LJ                                                  BMD 0034
      REAL ATHETA,ADBETA,PHASFN,AH1,ARH,ANGSUN,TBBYS,PATMS,WPATHS       BMD 0035
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),                  BMD 0036
     1  ADBETA(LAYDIM+1),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),      BMD 0037
     2  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,KMAX)    BMD 0038
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,                      BMD 0039
     1  NOPRNT,TBOUND,SALB,MODTRN                                       BMD 0040
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     BMD 0041
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                BMD 0042
      INTEGER IBNDWD,IP,ITB,NTEMP,IBIN,IMOL,IALF,JJ,JJS                 BMD 0043
      REAL TBAND,SDZ,ODZ,SD,OD,ALF0,T5,PTM75,FF,T5S,                    BMD 0044
     1  PTM75S,FFS,DOPFAC,DOP0,DOPSUM,COLSUM,ODSUM,SDSUM                BMD 0045
      COMMON/BMDCOM/IBNDWD,IP,ITB,NTEMP,TBAND(MXTEMP),SDZ(MXTEMP),      BMD 0046
     1  ODZ(MXTEMP),IBIN,IMOL,IALF,SD(MXTEMP,MMOLT2),OD(MXTEMP,MMOLT),  BMD 0047
     2  ALF0(MMOLT),T5(LAYTHR),PTM75(LAYTHR),JJ(LAYTHR),FF(LAYTHR),     BMD 0048
     3  T5S(LAYTHR,NSPC),PTM75S(LAYTHR,NSPC),JJS(LAYTHR,MMOLT),         BMD 0049
     4  FFS(LAYTHR,MMOLT),DOPFAC(MMOLT),DOP0(MMOLT),SDSUM(MMOLT2),      BMD 0050
     5  ODSUM(MMOLT),DOPSUM(MMOLT),COLSUM(MMOLT)                        BMD 0051
      DATA ITBX/31/,TZERO/273.15/                                       BMD 0052
C                                                                       BMD 0053
C     OPEN THE UNFORMATTED BAND MODEL TAPE, FILE ITB = 9.  AN ERROR     BMD 0054
C     OCCURS IF THE TAPE WAS OPENED DURING A PREVIOUS LOWTRAN RUN       BMD 0055
C     (IMPLYING IRPT WAS NOT ZERO).  IN THAT CASE, REWIND THE TAPE.     BMD 0056
      ITB=9                                                             BMD 0057
      OPEN(ITBX,FILE='DATA/CFCBMP96.ASC',STATUS='OLD',FORM='FORMATTED') BMD 0058
      REWIND(ITBX)                                                      BMD 0059
C                                                                       BMD 0060
C     READ THE BAND MODEL TAPE HEADER.                                  BMD 0061
      OPEN(ITB,FILE='DATA/MOLBMP96.BIN',                                BMD 0062
     1  STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED',RECL=60)        BMD 0063
      READ(ITB,REC=1)NTEMP,(TBAND(IT),IT=1,NTEMP),IBNDWD                BMD 0064
C                                                                       BMD 0065
C     PASS THE BAND WIDTH TO TRANS THROUGH IDVX.  ALSO DETERMINE THE    BMD 0066
C     MINIMUM FREQUENCY AT WHICH BAND MODEL DATA IS NEEDED.             BMD 0067
      IDVX=IBNDWD                                                       BMD 0068
      IV=5*((IV1-IFWHM/IBNDWD+1)/5)                                     BMD 0069
C                                                                       BMD 0070
C     PERFORM THE STANDARD LOWTRAN CALCULATION IF NO BAND MODEL DATA    BMD 0071
C     EXISTS FOR THE CHOSEN FREQUENCY RANGE.                            BMD 0072
      MXFREQ=22681                                                      BMD 0073
      IF(IV.GT.MXFREQ)THEN                                              BMD 0074
          WRITE(IPR,'(//43H (THE BAND MODEL TAPE DOES NOT HAVE DATA IN, BMD 0075
     1      32H THE REQUESTED WAVENUMBER RANGE))')                      BMD 0076
          RETURN                                                        BMD 0077
      ENDIF                                                             BMD 0078
C                                                                       BMD 0079
C     FIND 1ST RECORD (IP) IN BAND MODEL FILE WHERE FREQUENCY IV OCCURS.BMD 0080
      CALL GTSTRT(IV,IP,ITB)                                            BMD 0081
C                                                                       BMD 0082
C     READ THE FIRST RECORD                                             BMD 0083
      READ(ITB,REC=IP)                                                  BMD 0084
     1  IBIN,IMOL,(SDZ(I),I=1,MXTEMP),IALF,(ODZ(I),I=1,MXTEMP)          BMD 0085
C                                                                       BMD 0086
C     READ THE BAND MODEL PARAMETERS (ABSORPTION CROSS-SECTIONS ONLY)   BMD 0087
C     AND POSITION AT THE PROPER PLACE OF EACH FILE OF EACH X-SPECIES.  BMD 0088
   10 CONTINUE                                                          BMD 0089
      READ(ITBX,'(I6,I5,1P11E11.3)',END=20)I,IMOLX,(SDZX(IT),IT=1,NTEMP)BMD 0090
      IBINX=I                                                           BMD 0091
      READ(ITBX,'(6X,I5,1P11E11.3)')IALFX,(ODZX(IT),IT=1,NTEMP)         BMD 0092
      IF(IV.GT.IBINX)GOTO10                                             BMD 0093
   20 CONTINUE                                                          BMD 0094
C                                                                       BMD 0095
C     SET MSMAX TO LAYTWO FOR MULTIPLE SCATTERING                       BMD 0096
      MSMAX=ABS(IMULT)*LAYTWO                                           BMD 0097
C                                                                       BMD 0098
C     SET TEMPERATURE INTERPOLATION INDICES FOR EACH LAYER              BMD 0099
      IKHI=IKMX                                                         BMD 0100
      DO 60 MSOFF=0,MSMAX,LAYTWO                                        BMD 0101
          DO 50 IK=1,IKHI                                               BMD 0102
              IKOFF=IK+MSOFF                                            BMD 0103
              TT=TBBY(IKOFF)                                            BMD 0104
              IF(TT.LE.TBAND(1))THEN                                    BMD 0105
                  JJ(IKOFF)=2                                           BMD 0106
                  FF(IKOFF)=1.                                          BMD 0107
              ELSEIF(TT.GE.TBAND(NTEMP))THEN                            BMD 0108
                  JJ(IKOFF)=NTEMP                                       BMD 0109
                  FF(IKOFF)=0.                                          BMD 0110
              ELSE                                                      BMD 0111
                  DO 30 J=2,NTEMP                                       BMD 0112
                      IF(TT.LE.TBAND(J))GOTO40                          BMD 0113
   30             CONTINUE                                              BMD 0114
   40             JJ(IKOFF)=J                                           BMD 0115
                  FF(IKOFF)=(TBAND(J)-TT)/(TBAND(J)-TBAND(J-1))         BMD 0116
              ENDIF                                                     BMD 0117
C                                                                       BMD 0118
C         SET TEMPERATURE SCALING PARAMETERS                            BMD 0119
C           T5      LAYER TEMPERATURE DIVIDED BY 273.15K RAISED TO 0.5  BMD 0120
C           PTM75   LAYER TEMPERATURE DIVIDED BY 273.15K RAISED TO -0.75BMD 0121
C                   TIMES THE LAYER PRESSURE IN ATMOSPHERES.            BMD 0122
          T5(IKOFF)=SQRT(TT/TZERO)                                      BMD 0123
          PTM75(IKOFF)=PATM(IKOFF)*(TZERO/TT)**.75                      BMD 0124
   50     CONTINUE                                                      BMD 0125
          IKHI=IMSMX                                                    BMD 0126
   60  CONTINUE                                                         BMD 0127
C                                                                       BMD 0128
C     IF NO SINGLE SCATTER SOLAR, RETURN                                BMD 0129
      IF(IEMSCT.NE.2)RETURN                                             BMD 0130
C                                                                       BMD 0131
C     SET TEMPERATURE INTERPOLATION INDICES FOR SOLAR LAYERS            BMD 0132
      IKHI=IKMX+1                                                       BMD 0133
      DO 110 MSOFF=0,MSMAX,LAYTWO                                       BMD 0134
          DO 100 IK=1,IKHI                                              BMD 0135
              IKOFF=IK+MSOFF                                            BMD 0136
C                                                                       BMD 0137
C             SKIP SET UP IF THE SUN IS IN THE SHADE.                   BMD 0138
              IF(WPATHS(IKOFF,36).LT.0.)GOTO100                         BMD 0139
              DO 90 K=1,NSPECT                                          BMD 0140
                  IF(K.LE.NSPC)THEN                                     BMD 0141
                      TT=TBBYS(IKOFF,K)                                 BMD 0142
                  ELSE                                                  BMD 0143
                      TT=TBBYSX(IKOFF,K-NSPC)                           BMD 0144
                  ENDIF                                                 BMD 0145
                  IF(TT.LE.TBAND(1))THEN                                BMD 0146
                      JJS(IKOFF,K)=2                                    BMD 0147
                      FFS(IKOFF,K)=1.                                   BMD 0148
                  ELSEIF(TT.GE.TBAND(NTEMP))THEN                        BMD 0149
                      JJS(IKOFF,K)=NTEMP                                BMD 0150
                      FFS(IKOFF,K)=0.                                   BMD 0151
                  ELSE                                                  BMD 0152
                      DO 70 J=2,NTEMP                                   BMD 0153
                          IF(TT.LE.TBAND(J))GOTO80                      BMD 0154
   70                 CONTINUE                                          BMD 0155
   80                 JJS(IKOFF,K)=J                                    BMD 0156
                      FFS(IKOFF,K)=(TBAND(J)-TT)/(TBAND(J)-TBAND(J-1))  BMD 0157
                  ENDIF                                                 BMD 0158
                  IF(K.GT.NSPC)GOTO90                                   BMD 0159
C                                                                       BMD 0160
C                 SET TEMPERATURE SCALING PARAMETERS FOR SOLAR PATHS.   BMD 0161
C                   T5S     SOLAR PATH TEMPERATURE DIVIDED              BMD 0162
C                           BY 273.15K RAISED TO 0.5                    BMD 0163
C                   PTM75S  SOLAR PATH TEMPERATURE DIVIDED BY           BMD 0164
C                           273.15K RAISED TO -0.75 TIMES SOLAR         BMD 0165
C                           PATH PRESSURE IN ATMOSPHERES.               BMD 0166
                  T5S(IKOFF,K)=SQRT(TT/TZERO)                           BMD 0167
                  PTM75S(IKOFF,K)=PATMS(IKOFF,K)*(TZERO/TT)**.75        BMD 0168
   90         CONTINUE                                                  BMD 0169
  100     CONTINUE                                                      BMD 0170
          IKHI=IMSMX+1                                                  BMD 0171
  110 CONTINUE                                                          BMD 0172
      RETURN                                                            BMD 0173
      END                                                               BMD 0174
