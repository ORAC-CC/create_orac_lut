      SUBROUTINE RDSUN(FWHM,SUNFIL)                                     RSN 0001
C                                                                       RSN 0002
C     THIS ROUTINE READS IN THE "SUNFIL" SOLAR IRRADIANCE FILE AND      RSN 0003
C     PASSES A TRIANGULAR SLIT OVER IT WITH FWHM EQUAL TO |FWHM|.       RSN 0004
C     THE RESULTING DATA IS STORED IN ARRAY SUN OF /SOLAR1/.            RSN 0005
C                                                                       RSN 0006
C     DECLARE ARGUMENTS                                                 RSN 0007
C       FWHM     FULL-WIDTH-AT-HALF-MAXIMUM USED IN TRIANGULAR SLIT     RSN 0008
C                FUNCTION [IF NEGATIVE, ABS(FWHM) IS USED IN THE        RSN 0009
C                SLIT FUNCTION AND THE DATA IS OUTPUT TO A FILE].       RSN 0010
C       SUNFIL   INPUT DATA FILE NAME.                                  RSN 0011
      INTEGER FWHM                                                      RSN 0012
      CHARACTER*(*) SUNFIL                                              RSN 0013
C                                                                       RSN 0014
C     DECLARE PARAMETERS.                                               RSN 0015
      INTEGER MAXSUN                                                    RSN 0016
      PARAMETER(MAXSUN=50000)                                           RSN 0017
C                                                                       RSN 0018
C     LIST COMMON BLOCKS.                                               RSN 0019
      INTEGER IRD,IPR,IPU,NPR,IPR1,ISCRCH                               RSN 0020
      COMMON/IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                           RSN 0021
      REAL SUN                                                          RSN 0022
      COMMON/SOLAR1/SUN(0:MAXSUN)                                       RSN 0023
C                                                                       RSN 0024
C     DECLARE FUNCTIONS.                                                RSN 0025
      INTEGER NUNIT                                                     RSN 0026
C                                                                       RSN 0027
C     DECLARE LOCAL VARIABLES.                                          RSN 0028
      INTEGER NSUN,IVLO,IVHI,IV,IRES,IRES2,IRESM1,IVOFF,I,J,IVOLD,LINE  RSN 0029
      LOGICAL LEXIST,LOPEN                                              RSN 0030
      CHARACTER*11 SUNDAT                                               RSN 0031
      REAL SUNIV,SUM                                                    RSN 0032
C                                                                       RSN 0033
C     DECLARE DATA.                                                     RSN 0034
C       A        COEFFICIENT FOR LOW FREQUENCY POWER LAW APPROXIMATION. RSN 0035
C       B        EXPONENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.    RSN 0036
C       LBLKDT   LOGICAL FLAG, TRUE TO CREATE sunbd.f FILE.             RSN 0037
      REAL A,B                                                          RSN 0038
      LOGICAL LBLKDT                                                    RSN 0039
      DATA A,B/3.50187E-13,1.93281/,LBLKDT/.FALSE./                     RSN 0040
C                                                                       RSN 0041
C     OPEN INPUT DATA FILE.                                             RSN 0042
      INQUIRE(FILE=SUNFIL,EXIST=LEXIST,OPENED=LOPEN,NUMBER=NSUN)        RSN 0043
      IF(.NOT.LEXIST)THEN                                               RSN 0044
          WRITE(IPR,'(/3A)')                                            RSN 0045
     1      ' ERROR in RDSUN:  File ',SUNFIL,' was not found.'          RSN 0046
          STOP ' ERROR in RDSUN:  Solar irradiance file was not found.' RSN 0047
      ENDIF                                                             RSN 0048
      IF(LOPEN)THEN                                                     RSN 0049
          REWIND(NSUN)                                                  RSN 0050
      ELSE                                                              RSN 0051
          NSUN=NUNIT()                                                  RSN 0052
          OPEN(NSUN,FILE=SUNFIL,STATUS='OLD')                           RSN 0053
      ENDIF                                                             RSN 0054
C                                                                       RSN 0055
C     READ FIRST LINE OF SOLAR DATA                                     RSN 0056
C       IVLO     INITIAL FREQUENCY [CM-1]                               RSN 0057
C       SUN      SOLAR IRRADIANCE AT IVLO [W CM-2 / CM-1]               RSN 0058
      READ(NSUN,'(//I7,F13.3)',ERR=100)IVLO,SUNIV                       RSN 0059
      IF(IVLO.LT.1 .OR. IVLO.GT.MAXSUN)GOTO100                          RSN 0060
      IVHI=IVLO                                                         RSN 0061
      SUN(IVHI)=SUNIV                                                   RSN 0062
C                                                                       RSN 0063
C     READ REMAINING SOLAR DATA (1 CM-1 INCREMENTS).                    RSN 0064
   10 CONTINUE                                                          RSN 0065
      READ(NSUN,'(I7,F13.3)',END=20,ERR=100)IV,SUNIV                    RSN 0066
      IF(IV.NE.IVHI+1)GOTO100                                           RSN 0067
      IVHI=IV                                                           RSN 0068
      SUN(IVHI)=SUNIV                                                   RSN 0069
      IF(IVHI.LT.MAXSUN)GOTO10                                          RSN 0070
   20 CONTINUE                                                          RSN 0071
      CLOSE(NSUN)                                                       RSN 0072
C                                                                       RSN 0073
C     PASS TRIANGULAR SLIT FUNCTION OVER DATA.                          RSN 0074
      IRES=ABS(FWHM)                                                    RSN 0075
      IF(IRES.GT.1)THEN                                                 RSN 0076
          IRES2=IRES*IRES                                               RSN 0077
          IRESM1=IRES-1                                                 RSN 0078
          DO 40 IV=IVLO+IRESM1,IVHI-IRESM1                              RSN 0079
              SUM=IRES*SUN(IV)                                          RSN 0080
              DO 30 IVOFF=1,IRESM1                                      RSN 0081
                  SUM=SUM+(IRES-IVOFF)*(SUN(IV+IVOFF)+SUN(IV-IVOFF))    RSN 0082
   30         CONTINUE                                                  RSN 0083
C                                                                       RSN 0084
C             STORE DATA OFFSET BY IRESM1 CM-1 TO                       RSN 0085
C             AVOID OVERWRITING REQUIRED DATA.                          RSN 0086
              SUN(IV-IRESM1)=SUM/IRES2                                  RSN 0087
   40     CONTINUE                                                      RSN 0088
      ELSE                                                              RSN 0089
          IRESM1=0                                                      RSN 0090
      ENDIF                                                             RSN 0091
C                                                                       RSN 0092
C     MOVE DATA TO PROPER LOCATION (FREQUENCY).                         RSN 0093
      DO 50 IV=IVHI-IRESM1,IVLO+IRESM1,-1                               RSN 0094
          SUN(IV)=SUN(IV-IRESM1)                                        RSN 0095
   50 CONTINUE                                                          RSN 0096
C                                                                       RSN 0097
C     DEFINE LOW FREQUENCY DATA [W CM-2 / CM-1] USING                   RSN 0098
C     POWER LAW APPROXIMATION.                                          RSN 0099
      SUN(0)=0.                                                         RSN 0100
      DO 60 IV=1,IVLO+IRESM1-1                                          RSN 0101
          SUN(IV)=A*FLOAT(IV)**B                                        RSN 0102
   60 CONTINUE                                                          RSN 0103
C                                                                       RSN 0104
C     INSERT A CONSTANT VALUE AT HIGH FREQUENCIES.                      RSN 0105
      DO 70 IV=IVHI-IRESM1+1,MAXSUN                                     RSN 0106
          SUN(IV)=SUN(IVHI-IRESM1)                                      RSN 0107
   70 CONTINUE                                                          RSN 0108
C                                                                       RSN 0109
C     THE FOLLOWING CODING WAS USED TO CREATE THE BLOCK DATA SUNBD.     RSN 0110
C     THE LOGICAL FLAG LBLKDT IS HARD-WIRED TO FALSE, BUT THE CODING IS RSN 0111
C     LEFT HERE TO AIDE IN INCORPORATION OF NEW SOLAR IRRADIANCE DATA.  RSN 0112
      IF(LBLKDT)THEN                                                    RSN 0113
C                                                                       RSN 0114
C         OPEN sunbd.f FILE.                                            RSN 0115
          OPEN(NSUN,FILE='sunbd.f',STATUS='unknown')                    RSN 0116
          CLOSE(NSUN,STATUS='delete')                                   RSN 0117
          OPEN(NSUN,FILE='sunbd.f',STATUS='new')                        RSN 0118
C                                                                       RSN 0119
C         WRITE HEADER.                                                 RSN 0120
          WRITE(NSUN,'((A))')                                           RSN 0121
     1      '      BLOCK DATA SUNBD',                                   RSN 0122
     2      'C',                                                        RSN 0123
     3      'C     SOLAR IRRADIANCES [W CM-2 / CM-1] TABULATED EACH',   RSN 0124
     4      'C     WAVENUMBER FROM 0 TO 50,000 CM-1 AT 5 CM-1 SPECTRAL',RSN 0125
     5      'C     RESOLUTION (FWHM TRIANGULAR SLIT).  THIS DATA IS',   RSN 0126
     6      'C     DERIVED FROM THE 1 CM-1 TABULATION (FILE SUN2) OF',  RSN 0127
     7      'C     KURUCZ HIGH SPECTRAL RESOLUTION CALCULATIONS.',      RSN 0128
     8      '      INTEGER MAXSUN',                                     RSN 0129
     9      '      PARAMETER(MAXSUN=50000)',                            RSN 0130
     &      '      REAL SUN',                                           RSN 0131
     1      '      COMMON/SOLAR1/SUN(0:MAXSUN)',                        RSN 0132
     2      '      INTEGER I'                                           RSN 0133
C                                                                       RSN 0134
C         WRITE 0 TO 100 CM-1 DATA.                                     RSN 0135
          WRITE(NSUN,'(A,2(/A),25X,1PE10.3,A)')                         RSN 0136
     1      'C','C         0 TO   100 CM-1 DATA.',                      RSN 0137
     2      '      DATA (SUN(I),I=     0,   100)/',SUN(0),','           RSN 0138
          WRITE(NSUN,'((9(I6,6(1PE10.3,A),/),A,6(1PE10.3,A)))')         RSN 0139
     1      (J,      (SUN(I),',',I=6*J-5,6*J),J=1,9),                   RSN 0140
     2      '     &',(SUN(I),',',I=55   ,60 ),                          RSN 0141
     3      (J-10,   (SUN(I),',',I=6*J-5,6*J),J=11,16),                 RSN 0142
     4       J-10,   (SUN(I),',',I=97,99),SUN(100),'/'                  RSN 0143
C                                                                       RSN 0144
C         WRITE 101 TO 100*INT(MAXSUM/100) CM-1 DATA.                   RSN 0145
          IVOLD=100                                                     RSN 0146
          DO 80 IV=200,MAXSUN,100                                       RSN 0147
              WRITE(NSUN,'(A,2(/A,2(I6,A)))')                           RSN 0148
     1          'C','C    ',IVOLD+1,' TO',IV,' CM-1 DATA.',             RSN 0149
     2          '      DATA (SUN(I),I=',IVOLD+1,',',IV,')/'             RSN 0150
              WRITE(NSUN,'((9(I6,6(1PE10.3,A),/),A,6(1PE10.3,A)))')     RSN 0151
     1          (J      ,(SUN(IVOLD+I),',',I=6*J-5,6*J),J=1,9),         RSN 0152
     2          '     &',(SUN(IVOLD+I),',',I=55,60),                    RSN 0153
     3          (J-10   ,(SUN(IVOLD+I),',',I=6*J-5,6*J),J=11,16),       RSN 0154
     4           J-10   ,(SUN(IVOLD+I),',',I=97,99),SUN(IV),'/'         RSN 0155
              IVOLD=IV                                                  RSN 0156
   80     CONTINUE                                                      RSN 0157
C                                                                       RSN 0158
C         WRITE REMAINING CM-1 DATA.                                    RSN 0159
          IF(IVOLD.LT.MAXSUN)THEN                                       RSN 0160
              WRITE(NSUN,'(A,2(/A,2(I6,A)))')                           RSN 0161
     1          'C','C    ',IVOLD+1,' TO',MAXSUN,' CM-1 DATA.',         RSN 0162
     2          '      DATA (SUN(I),I=',IVOLD+1,',',MAXSUN,')/'         RSN 0163
              LINE=1                                                    RSN 0164
   90         CONTINUE                                                  RSN 0165
              IF(IVOLD+6.LT.MAXSUN)THEN                                 RSN 0166
                  IF(LINE.LT.10)THEN                                    RSN 0167
                      WRITE(NSUN,'(I6,6(1PE10.3,A))')LINE,              RSN 0168
     1                  (SUN(IVOLD+I),',',I=1,6)                        RSN 0169
                      LINE=LINE+1                                       RSN 0170
                  ELSE                                                  RSN 0171
                      WRITE(NSUN,'(A,6(1PE10.3,A))')'     &',           RSN 0172
     1                  (SUN(IVOLD+I),',',I=1,6)                        RSN 0173
                      LINE=1                                            RSN 0174
                  ENDIF                                                 RSN 0175
                  IVOLD=IVOLD+6                                         RSN 0176
                  GOTO90                                                RSN 0177
              ENDIF                                                     RSN 0178
              IF(LINE.LT.10)THEN                                        RSN 0179
                  WRITE(NSUN,'(I6,6(1PE10.3,A))')LINE,                  RSN 0180
     1              (SUN(I),',',I=IVOLD+1,MAXSUN-1),SUN(MAXSUN),'/'     RSN 0181
              ELSE                                                      RSN 0182
                  WRITE(NSUN,'(A,6(1PE10.3,A))')'     &',               RSN 0183
     1              (SUN(I),',',I=IVOLD+1,MAXSUN-1),SUN(MAXSUN),'/'     RSN 0184
              ENDIF                                                     RSN 0185
          ENDIF                                                         RSN 0186
          WRITE(NSUN,'(A)')'      END'                                  RSN 0187
          CLOSE(NSUN)                                                   RSN 0188
      ENDIF                                                             RSN 0189
C                                                                       RSN 0190
C     RETURN UNLESS FWHM IS NEGATIVE                                    RSN 0191
      IF(FWHM.GE.0)RETURN                                               RSN 0192
C                                                                       RSN 0193
C     OPEN OUTPUT DATA FILE.                                            RSN 0194
      SUNDAT='sun0000.dat'                                              RSN 0195
      IF(IRES.LT.10)THEN                                                RSN 0196
          WRITE(SUNDAT(7:7),'(I1)')IRES                                 RSN 0197
      ELSEIF(IRES.LT.100)THEN                                           RSN 0198
          WRITE(SUNDAT(6:7),'(I2)')IRES                                 RSN 0199
      ELSEIF(IRES.LT.1000)THEN                                          RSN 0200
          WRITE(SUNDAT(5:7),'(I3)')IRES                                 RSN 0201
      ELSE                                                              RSN 0202
          WRITE(SUNDAT(4:7),'(I4)')IRES                                 RSN 0203
      ENDIF                                                             RSN 0204
      OPEN(NSUN,FILE=SUNDAT,STATUS='UNKNOWN')                           RSN 0205
      CLOSE(NSUN,STATUS='DELETE')                                       RSN 0206
      OPEN(NSUN,FILE=SUNDAT,STATUS='NEW')                               RSN 0207
C                                                                       RSN 0208
C     WRITE DATA TO SUNDAT IN ORIGINAL UNITS [W CM-2 / CM-1].           RSN 0209
      WRITE(NSUN,'(A,/A,/(I7,1PE13.3))')                                RSN 0210
     1  '  FREQ   SOLAR IRRADIANCE',                                    RSN 0211
     2  ' (CM-1)  (W CM-2 / CM-1)',(IV,SUN(IV),IV=1,MAXSUN)             RSN 0212
C                                                                       RSN 0213
C     RETURN TO DRIVER.                                                 RSN 0214
      RETURN                                                            RSN 0215
C                                                                       RSN 0216
C     WRITE OUT ERROR MESSAGE AND STOP.                                 RSN 0217
  100 CONTINUE                                                          RSN 0218
      WRITE(IPR,'(/3A)')                                                RSN 0219
     1  ' ERROR in RDSUN:  Problem reading/using file ',SUNFIL,' data.' RSN 0220
      STOP'ERROR in RDSUN:  Prolem reading/using solar irradiance data.'RSN 0221
      END                                                               RSN 0222
