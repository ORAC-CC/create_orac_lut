      SUBROUTINE RHOEPS (IV,IBKG,SALB,BKEPS)                            RHO 0001
      INTEGER KKMAX                                                     RHO 0002
      PARAMETER(KKMAX=52)                                               RHO 0003
C                                                                       RHO 0004
C     JUN 93                                                            RHO 0005
C         THIS SUBROUTINE IS ADAPTED FROM THE SPECTRAL SCIENCES TARGET  RHO 0006
C     IR SIGNATURE CODE (SSTIRS).                                       RHO 0007
C                                                                       RHO 0008
C         IV .LT. 0    RHOEPS READS THE BACKGROUND REFLECTIVITY FILE    RHO 0009
C             GE. 0    INTERPOLATES SURFACE ALBED (SALB) TO VALUE OF V  RHO 0010
C         IBKG         LABEL FOR REQUESTED BACKGROUND DATA FILE         RHO 0011
C         SALB         RETURNED VALUE OF SURFACE ALBEDO                 RHO 0012
C         BKEPS            "      "   "     "    EMISSIVITY  (NOT USED) RHO 0013
C                                                                       RHO 0014
C        THE BACKGROUND EMISSIVITY IS NOT USED, BECAUSE MODTRAN IS ASSUMRHO 0015
C     THAT IT ALWAYS EQUALS (1 - SALB).  IF THE USER DESIRES TO DIRECTLYRHO 0016
C     BKEPS, HE MUST MAKE THE APPROPRIATE CHANGES TO BEMISS IN SUBROUTINRHO 0017
C     TRANS PLUS SEVERAL PLACES IN BMFLUX & FLXADD WHERE RUPCN=SALB AND RHO 0018
C     (1-SALB) IS USED FOR THE SURFACE EMISSIVITY.  (NOTE THAT BKEPS OR RHO 0019
C     BEMISS IS NOT PASSED TO BMFLUX & FLXADD, BUT SALB IS.)            RHO 0020
C                                                                       RHO 0021
C                                                                       RHO 0022
C     FILE  NF09   BACKGROUND DATA FILE                                 RHO 0023
C           NOUT   TAPE6 FILE                                           RHO 0024
C                                                                       RHO 0025
      SAVE          VIN,R,E,KDATA                                       RHO 0026
      CHARACTER*4   HEADER(9)                                           RHO 0027
      COMMON  /UNITS/   NOUT,NF09                                       RHO 0028
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          RHO 0029
C                                                                       RHO 0030
      REAL VIN(KKMAX),R(KKMAX,6),E(KKMAX,6),TEMP1(KKMAX),TEMP2(KKMAX)   RHO 0031
      INTEGER J                                                         RHO 0032
      LOGICAL LWARN                                                     RHO 0033
      SAVE J,LWARN                                                      RHO 0034
C                                                                       RHO 0035
      DATA NDEG/1 /                                                     RHO 0036
C                                                                       RHO 0037
C ... FILE IDENTIFIER STARTS IN COLUMN 1                                RHO 0038
C         HEADER  =  BACK   -  BACKGROUND RHO-EPS FILE                  RHO 0039
C                 =  END    -  END OF THIS DATA SET                     RHO 0040
C...THE END OF EACH REFLECTIVITY FILE IS IDENTIFIED BY                  RHO 0041
C     A BLANK CARD WITH A A NON-ZERO INTEGER IN COLUMN ONE.             RHO 0042
C                                                                       RHO 0043
C                                                                       RHO 0044
      NF09=28                                                           RHO 0045
      NOUT=IPR                                                          RHO 0046
C                                                                       RHO 0047
      V=IV                                                              RHO 0048
      IF(IV.GE.0.)GOTO80                                                RHO 0049
C                                                                       RHO 0050
C              **  **  **  **  **  **  **  **  **  **  **  **  **  **   RHO 0051
C              **  **                                          **  **   RHO 0052
C              **  **   READ  BACKGROUND  REFLECTANCE  DATA    **  **   RHO 0053
C              **  **                                          **  **   RHO 0054
C              **  **  **  **  **  **  **  **  **  **  **  **  **  **   RHO 0055
C                                                                       RHO 0056
C ... READ TITLE CARD (NOT USED)                                        RHO 0057
      READ (NF09,'(9a4)')    HEADER                                     RHO 0058
C                                                                       RHO 0059
C ... FIRST, INITIALIZE VARIOUS PARAMETERS                              RHO 0060
C               NTEST = FILE NUMBER FOR BACKGROUND SPECTRAL ALBEDOS     RHO 0061
C                   J = INDEX FOR /REFEPS/                              RHO 0062
C                   K = DATA SET NUMBER                                 RHO 0063
C               LWARN = FREQUENCY MISMATCH WARNING FLAG                 RHO 0064
      NTEST=IBKG                                                        RHO 0065
      J=1                                                               RHO 0066
      KK=1                                                              RHO 0067
      LWARN=.TRUE.                                                      RHO 0068
C                                                                       RHO 0069
C ... READ NUMBER OF BACKGROUND FILES & CHECK VALIDITY OF NTEST         RHO 0070
      READ (NF09,'(i2)') JFILE                                          RHO 0071
      IF (JFILE .LE. 0) JFILE=1                                         RHO 0072
      NLOW=1                                                            RHO 0073
C                                       ** ** **  READ DATA  ** ** **   RHO 0074
C                                                                       RHO 0075
C                             ... SKIP TO REQUESTED FILE                RHO 0076
      DO 10 L=NLOW,JFILE                                                RHO 0077
         READ (NF09,'(13x,i2)',END=16) LBKG                             RHO 0078
         READ (NF09,'(i1)'   ,END=16) NU                                RHO 0079
         IF (LBKG .EQ. NTEST) THEN                                      RHO 0080
            GO TO 20                                                    RHO 0081
         ELSE                                                           RHO 0082
   13       READ (NF09,  '(i1)'  ,END=16) N                             RHO 0083
            IF (N .EQ. 0) GO TO 13                                      RHO 0084
         END IF                                                         RHO 0085
   10 CONTINUE                                                          RHO 0086
C                                                                       RHO 0087
C                              ... PROBLEMS -- DID NOT FIND MATCHING FILRHO 0088
   16 WRITE (NOUT,17) NTEST,IBKG,L,JFILE                                RHO 0089
   17 FORMAT (5('**  '),3X,'stop  in  subroutine  rhoeps',3X,5('  **')/ RHO 0090
     &       10X,'no match was found for the requested background ',    RHO 0091
     &       'file number (=',I3,') ---- suggest checking columns '/    RHO 0092
     &       10X,'14 & 15 in the background data files (fortran ',      RHO 0093
     &       'file: nf09)'/20X,'ibkg=',I4,5X,'file',I3,                 RHO 0094
     &       ' is last background file read  (of',I3,' files')          RHO 0095
      STOP  'sub.rhoeps -- background reflectance file'                 RHO 0096
C                                                                       RHO 0097
C                                                    ... READ DATA      RHO 0098
   20 READ(NF09,22)  JDUM,VIN(KK),R(KK,J),E(KK,J)                       RHO 0099
   22 FORMAT (I1,F9.0,6(2X,2F4.2))                                      RHO 0100
C                                                                       RHO 0101
C                                        ... CHECK FOR END OF FILE      RHO 0102
      IF   (JDUM .EQ. 0)  THEN                                          RHO 0103
         IF(KK.LT.KKMAX)THEN                                            RHO 0104
            KK=KK + 1                                                   RHO 0105
         END IF                                                         RHO 0106
         GO TO 20                                                       RHO 0107
      END IF                                                            RHO 0108
      KDATA=KK - 1                                                      RHO 0109
C                                                                       RHO 0110
C       ** ** ** **  IF NEEDED, CONVERT TO WAVENUMBERS (CM-1) ** ** ** *RHO 0111
      IF (NU .EQ. 0) THEN                                               RHO 0112
         DO 30 K=1,KDATA                                                RHO 0113
           TEMP1(K)=R(K,J)                                              RHO 0114
           TEMP2(K)=E(K,J)                                              RHO 0115
   30    CONTINUE                                                       RHO 0116
         DO 40 K=1,KDATA                                                RHO 0117
           KK     =KDATA + 1 - K                                        RHO 0118
           R(KK,J)=TEMP1(K)                                             RHO 0119
           E(KK,J)=TEMP2(K)                                             RHO 0120
   40    CONTINUE                                                       RHO 0121
         DO 50 K=1,KDATA                                                RHO 0122
            IF (VIN(K) .LE. 0.)                                         RHO 0123
     1             STOP 'rhoeps:  cannot have zero wavenumber in nf09'  RHO 0124
            TEMP1(K)=10000./VIN(K)                                      RHO 0125
   50    CONTINUE                                                       RHO 0126
C                                       ... ROUND WAVENUMBERS TO 5 CM-1 RHO 0127
         DO 60 K=1,KDATA                                                RHO 0128
            KK      =KDATA + 1 - K                                      RHO 0129
            VIN(KK)=FLOAT(5*INT(.2*TEMP1(K)+.5))                        RHO 0130
   60    CONTINUE                                                       RHO 0131
      END IF                                                            RHO 0132
C                                                                       RHO 0133
C                                                                       RHO 0134
C              **  **  **  **  **  **  **  **  **  **  **  **  **  **   RHO 0135
C              **  **                                       **  **      RHO 0136
C              **  **   INTERPOLATE TO DESIRED WAVENUMBER   **  **      RHO 0137
C              **  **                                       **  **      RHO 0138
C              **  **  **  **  **  **  **  **  **  **  **  **  **  **   RHO 0139
C                                                                       RHO 0140
C                                                                       RHO 0141
      RETURN                                                            RHO 0142
   80 CONTINUE                                                          RHO 0143
      V    =FLOAT(IV)                                                   RHO 0144
      IF(V.LT.VIN(1))THEN                                               RHO 0145
          SALB =R(1,1)                                                  RHO 0146
          BKEPS=E(1,1)                                                  RHO 0147
          IF(LWARN)THEN                                                 RHO 0148
              WRITE(IPR,'(/A,I3,A,/17X,2(A,I6),A,/17X,A,/)')            RHO 0149
     1          ' RHOEPS WARNING:  SPECTRAL ALBEDO DATA SET',           RHO 0150
     2          IBKG,' COVERS THE SPECTRAL RANGE',' FROM',              RHO 0151
     3          INT(VIN(1)+.01),' TO',INT(VIN(KDATA)+.01),              RHO 0152
     4          ' CM-1.  END POINT VALUES WERE USED',                   RHO 0153
     5          ' FOR FREQUENCIES OUTSIDE OF THE FREQUENCY RANGE.'      RHO 0154
              LWARN=.FALSE.                                             RHO 0155
          ENDIF                                                         RHO 0156
      ELSEIF(V.GT.VIN(KDATA))THEN                                       RHO 0157
          SALB =R(KDATA,1)                                              RHO 0158
          BKEPS=E(KDATA,1)                                              RHO 0159
          IF(LWARN)THEN                                                 RHO 0160
              WRITE(IPR,'(/A,I3,A,/17X,2(A,I6),A,/17X,A,/)')            RHO 0161
     1          ' RHOEPS WARNING:  SPECTRAL ALBEDO DATA SET',           RHO 0162
     2          IBKG,' COVERS THE SPECTRAL RANGE',' FROM',              RHO 0163
     3          INT(VIN(1)+.01),' TO',INT(VIN(KDATA)+.01),              RHO 0164
     4          ' CM-1.  END POINT VALUES WERE USED',                   RHO 0165
     5          ' FOR FREQUENCIES OUTSIDE OF THE FREQUENCY RANGE.'      RHO 0166
              LWARN=.FALSE.                                             RHO 0167
          ENDIF                                                         RHO 0168
      ELSE                                                              RHO 0169
          SALB =XTERP(V,VIN,R(1,J),NDEG,KDATA,DERINT,IER)               RHO 0170
          BKEPS=XTERP(V,VIN,E(1,J),NDEG,KDATA,DERINT,IER)               RHO 0171
      ENDIF                                                             RHO 0172
      IF(SALB.GT.1. .OR. SALB.LT.0.)THEN                                RHO 0173
        PRINT*,' v salb ',V,SALB                                        RHO 0174
        STOP                                                            RHO 0175
      ENDIF                                                             RHO 0176
      RETURN                                                            RHO 0177
      END                                                               RHO 0178
