      SUBROUTINE CIRRUS(CTHIK,CALT,ISEED,CPROB,CEXT)                    CIR 0001
C*********************************************************************  CIR 0002
C*  ROUTINE TO GENERATE ALTITUDE PROFILES OF CIRRUS DENSITY         **  CIR 0003
C*  PROGRAMMED BY   M.J. POST                                       **  CIR 0004
C*                  R.A. RICHTER        NOAA/WPL                    **  CIR 0005
C*                                      BOULDER, COLORADO           **  CIR 0006
C*                                      01/27/1981                  **  CIR 0007
C*                                                                  **  CIR 0008
C*  INPUTS!                                                         **  CIR 0009
C*           CTHIK    -  CIRRUS THICKNESS (KM)                      **  CIR 0010
C*                       0 = USE THICKNESS STATISTICS               **  CIR 0011
C*                       .NE. 0 = USER DEFINES THICKNESS            **  CIR 0012
C*                                                                  **  CIR 0013
C*           CALT     -  CIRRUS BASE ALTITUDE (KM)                  **  CIR 0014
C*                       0 = USE CALCULATED VALUE                   **  CIR 0015
C*                       .NE. 0 = USER DEFINES BASE ALTITUDE        **  CIR 0016
C*                                                                  **  CIR 0017
C*           ICIR     -  CIRRUS PRESENCE FLAG                       **  CIR 0018
C*                       0 = NO CIRRUS                              **  CIR 0019
C*                       .NE. 0 = USE CIRRUS PROFILE                **  CIR 0020
C*                                                                  **  CIR 0021
C*           MODEL    -  ATMOSPHERIC MODEL                          **  CIR 0022
C*                       1-5  AS IN MAIN PROGRAM                    **  CIR 0023
C*                       MODEL = 0,6,7 NOT USED SET TO 2            **  CIR 0024
C*                                                                  **  CIR 0025
C*           ISEED    -  RANDOM NUMBER INITIALIZATION FLAG.         **  CIR 0026
C*                       0 = USE DEFAULT MEAN VALUES FOR CIRRUS     **  CIR 0027
C*                       .NE. 0 = INITIAL VALUE OF SEED FOR RANF    **  CIR 0028
C*                       FUNCTION. CHANGE SEED VALUE EACH RUN FOR   **  CIR 0029
C*                       DIFFERENT RANDOM NUMBER SEQUENCES. THIS    **  CIR 0030
C*                       PROVIDES FOR STATISTICAL DETERMINATION     **  CIR 0031
C*                       OF CIRRUS BASE ALTITUDE AND THICKNESS.     **  CIR 0032
C*                                                                  **  CIR 0033
C*  OUTPUTS!                                                        **  CIR 0034
C*         CTHIK        -  CIRRUS THICKNESS (KM)                    **  CIR 0035
C*         CALT         -  CIRRUS BASE ALTITUDE (KM)                **  CIR 0036
C*         DENSTY(16,I) -  ARRAY, ALTITUDE PROFILE OF CIRRUS DENSITY**  CIR 0037
C*         CPROB        -  CIRRUS PROBABILITY                       **  CIR 0038
C*                                                                  **  CIR 0039
C*********************************************************************  CIR 0040
C                                                                       CIR 0041
      COMMON /CARD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB  CIR 0042
     1  ,MODTRN                                                         CIR 0043
      LOGICAL MODTRN                                                    CIR 0044
      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   CIR 0045
     1    RAINRT                                                        CIR 0046
      INCLUDE 'PARAM.LST'                                               CIR 0047
      INTEGER KPOINT                                                    CIR 0048
      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     CIR 0049
      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   CIR 0050
     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   CIR 0051
     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     CIR 0052
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     CIR 0053
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                CIR 0054
      REAL ZM,PM,TM,RFNDX,DENSTY                                        CIR 0055
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    CIR 0056
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               CIR 0057
      DIMENSION CBASE(5,2),TSTAT(11),PTAB(5),CAMEAN(5)                  CIR 0058
      DIMENSION CBASE1(5),CBASE2(5)                                     CIR 0059
      EQUIVALENCE (CBASE1(1),CBASE(1,1)),(CBASE2(1),CBASE(1,2))         CIR 0060
C                                                                       CIR 0061
      DATA  CAMEAN           / 11.0, 10.0, 8.0, 7.0, 5.0 /              CIR 0062
      DATA  PTAB           / 0.8, 0.4, 0.5, 0.45, 0.4/                  CIR 0063
      DATA  CBASE1            / 7.5, 7.3, 4.5, 4.5, 2.5 /               CIR 0064
      DATA  CBASE2            /16.5,13.5,14.0, 9.5,10.0 /               CIR 0065
      DATA  TSTAT             / 0.0,.291,.509,.655,.764,.837,.892,      CIR 0066
     + 0.928, 0.960, 0.982, 1.00 /                                      CIR 0067
C                                                                       CIR 0068
C  SET CIRRUS PROBABILITY AND PROFILE TO ALL ZEROES                     CIR 0069
C                                                                       CIR 0070
      CPROB = 0.0                                                       CIR 0071
      MDL = MODEL                                                       CIR 0072
C                                                                       CIR 0073
      DO 10 I=1,ML                                                      CIR 0074
10    DENSTY(16,I)=0.                                                   CIR 0075
C                                                                       CIR 0076
C  CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES, CALCULATE  CIR 0077
C  A STATISTICAL THICKNESS, OR USE A MEAN THICKNESS (ISEED = 0).        CIR 0078
C  DEFAULTED MEAN CIRRUS THICKNESS IS 1.0 KM.                           CIR 0079
C                                                                       CIR 0080
      IF ( CTHIK .GT. 0.0 ) GO TO 25                                    CIR 0081
      IF ( ISEED .NE. 0 ) GO TO 15                                      CIR 0082
      CTHIK = 1.0                                                       CIR 0083
      GO TO 25                                                          CIR 0084
C                                                                       CIR 0085
C  CALCULATE CLOUD THICKNESS USING LOWTRAN CIRRUS THICKNESS STATISTICS  CIR 0086
C  NOTE - THIS ROUTINE USES A UNIFORM RANDOM NUMBER GENERATOR           CIR 0087
C  FUNCTION (RANF) WHICH RETURNS A NUMBER BETWEEN 0 AND 1.              CIR 0088
C  THIS FEATURE IS MACHINE DEPENDENT!!                                  CIR 0089
C                                                                       CIR 0090
   15 CALL RANSET(ISEED)                                                CIR 0091
      URN = RANFUN(IDUM)                                                CIR 0092
      DO 20 I = 1, 10                                                   CIR 0093
         IF (URN .GE. TSTAT(I) .AND. URN .LT. TSTAT(I+1)) CTHIK = I-1   CIR 0094
   20 CONTINUE                                                          CIR 0095
      CTHIK = CTHIK / 2.0  +  RANFUN(IDUM) / 2.0                        CIR 0096
C                                                                       CIR 0097
C  DENCIR IS CIRRUS DENSITY IN KM-1                                     CIR 0098
C                                                                       CIR 0099
25    IF(CEXT .GT. 0.) THEN                                             CIR 0100
           DENCIR = CEXT / 2.                                           CIR 0101
      ELSE                                                              CIR 0102
           DENCIR = 0.07 * CTHIK                                        CIR 0103
      ENDIF                                                             CIR 0104
C                                                                       CIR 0105
C  BASE HEIGHT CALCULATIONS                                             CIR 0106
C                                                                       CIR 0107
      IF ( MODEL .LT. 1  .OR.  MODEL .GT. 5 ) MDL = 2                   CIR 0108
      CPROB = 100.0 * PTAB(MDL)                                         CIR 0109
C                                                                       CIR 0110
      HMAX = CBASE(MDL,2) - CTHIK                                       CIR 0111
      BRANGE = HMAX - CBASE(MDL,1)                                      CIR 0112
      IF ( CALT .GT. 0.0 ) GO TO 27                                     CIR 0113
      IF ( ISEED .NE. 0 ) GO TO 26                                      CIR 0114
      CALT = CAMEAN(MDL)                                                CIR 0115
      GO TO 27                                                          CIR 0116
   26 CALT = BRANGE * RANFUN(IDUM)+ CBASE(MDL,1)                        CIR 0117
C                                                                       CIR 0118
C  PUT CIRRUS DENSITY IN CORRECT ALTITUDE BINS. IF MODEL = 7,           CIR 0119
C  INTERPOLATE EH(16,I) FOR NON-STANDARD ALTITUDE BOUNDARIES.           CIR 0120
C                                                                       CIR 0121
   27 IF(MODEL .EQ. 7) GO TO 60                                         CIR 0122
      IV1=INT(CALT )                                                    CIR 0123
      IV2=INT(CALT+CTHIK )                                              CIR 0124
      DO 30 I = 2, 16                                                   CIR 0125
         IF(I .GE. IV1 .AND. I .LE. IV2) DENSTY(16,I+1) =  DENCIR       CIR 0126
   30 CONTINUE                                                          CIR 0127
C                                                                       CIR 0128
C  ADJUST FIRST AND LAST CIRRUS LEVEL IF CLOUD DOES NOT ENTIRELY        CIR 0129
C  FILL EACH LEVEL.                                                     CIR 0130
C                                                                       CIR 0131
      IHGT1 = INT( CALT )                                               CIR 0132
      IHGT2 = INT( CALT + CTHIK)                                        CIR 0133
      IF( IHGT1 . NE . IHGT2 ) GO TO 35                                 CIR 0134
      DENSTY(16,IHGT1+1) = DENSTY( 16,IHGT1+1)*CTHIK                    CIR 0135
      RETURN                                                            CIR 0136
   35 PCT1  = 1.0 - ( CALT - IHGT1 )                                    CIR 0137
      DENSTY(16,IHGT1+1) = DENSTY(16,IHGT1+1) * PCT1                    CIR 0138
      PCT2 =  ( CALT + CTHIK) - IHGT2                                   CIR 0139
      DENSTY(16,IHGT2+1) = DENSTY(16,IHGT2+1) * PCT2                    CIR 0140
      RETURN                                                            CIR 0141
C                                                                       CIR 0142
C  INTERPOLATE DENSTY(16,I) FOR USER SUPPLIED ALTITUDE BOUNDARIES       CIR 0143
C                                                                       CIR 0144
   60 TOP = CALT + CTHIK                                                CIR 0145
      BOTTOM = CALT                                                     CIR 0146
      IF(TOP.LT.ZM(1) .OR. BOTTOM.GT.ZM(ML))RETURN                      CIR 0147
      IML = ML - 1                                                      CIR 0148
      DO 70 I=1,IML                                                     CIR 0149
         ZMIN=ZM(I)                                                     CIR 0150
         ZMAX=ZM(I+1)                                                   CIR 0151
         DENOM = ZMAX - ZMIN                                            CIR 0152
         IF(BOTTOM.LE.ZMIN .AND. TOP.GE.ZMAX) DENSTY(16,I) = DENCIR     CIR 0153
         IF(BOTTOM.GE.ZMIN .AND. TOP.LT.ZMAX)                           CIR 0154
     +        DENSTY(16,I) = DENCIR * CTHIK/DENOM                       CIR 0155
         IF(BOTTOM.GE.ZMIN .AND. TOP.GE.ZMAX .AND. BOTTOM.LT.ZMAX)      CIR 0156
     +        DENSTY(16,I) = DENCIR * (ZMAX - BOTTOM)/ DENOM            CIR 0157
         IF(BOTTOM.LT.ZMIN .AND. TOP.LE.ZMAX .AND.TOP.GT.ZMIN)          CIR 0158
     +        DENSTY(16,I) = DENCIR * (TOP - ZMIN) / DENOM              CIR 0159
   70 CONTINUE                                                          CIR 0160
      RETURN                                                            CIR 0161
      END                                                               CIR 0162
