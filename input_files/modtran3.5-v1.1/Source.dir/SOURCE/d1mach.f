      DOUBLE PRECISION FUNCTION D1MACH(I)                               DPM 0001
                                                                        DPM 0002
C  DOUBLE-PRECISION MACHINE CONSTANTS (SEE R1MACH FOR DOCUMENTATION)    DPM 0003
                                                                        DPM 0004
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST     DPM 0005
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.                   DPM 0006
                                                                        DPM 0007
      INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4), SC      DPM 0008
      DOUBLE PRECISION DMACH(5), EPS, EPSNEW, S                         DPM 0009
                                                                        DPM 0010
      EQUIVALENCE (DMACH(1),SMALL(1)), (DMACH(2),LARGE(1)),             DPM 0011
     $            (DMACH(3),RIGHT(1)), (DMACH(4),DIVER(1)),             DPM 0012
     $            (DMACH(5),LOG10(1))                                   DPM 0013
                                                                        DPM 0014
      LOGICAL  PASS1                                                    DPM 0015
      SAVE     PASS1                                                    DPM 0016
      DATA     PASS1/.TRUE./                                            DPM 0017
                                                                        DPM 0018
C IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T 3B SERIES AND              DPM 0019
C MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T PC 7300),          DPM 0020
C IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.                   DPM 0021
Cj HP, SUN and SIG are big-endian                                       DPM 0022
                                                                        DPM 0023
c      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2146435071,-1/,DPM 0024
c     $  (RIGHT(N),N=1,2)/1017118720,0/, (DIVER(N),N=1,2)/1018167296,0/, DPM 0025
c     $  (LOG10(N),N=1,2)/1070810131,1352628735/, SC/987/                DPM 0026
                                                                        DPM 0027
C IEEE ARITHMETIC MACHINES AND 8087-BASED MICROS, SUCH AS THE IBM PC    DPM 0028
C AND AT&T 6300, IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST.   DPM 0029
Cj ie. little-endian                                                    DPM 0030
                                                                        DPM 0031
      DATA (SMALL(N),N=1,2)/0,1048576/, (LARGE(N),N=1,2)/-1,2146435071/,DPM 0032
     $  (RIGHT(N),N=1,2)/0,1017118720/, (DIVER(N),N=1,2)/0,1018167296/, DPM 0033
     $  (LOG10(N),N=1,2)/1352628735,1070810131/, SC/987/                DPM 0034
                                                                        DPM 0035
C AMDAHL MACHINES.                                                      DPM 0036
                                                                        DPM 0037
C     DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2147483647,-1/,DPM 0038
C    $ (RIGHT(N),N=1,2)/856686592,0/, (DIVER(N),N=1,2)/ 873463808,0/,   DPM 0039
C    $ (LOG10(N),N=1,2)/1091781651,1352628735/, SC/987/                 DPM 0040
                                                                        DPM 0041
C BURROUGHS 1700 SYSTEM.                                                DPM 0042
                                                                        DPM 0043
C     DATA (SMALL(N),N=1,2)/ZC00800000,Z000000000/,                     DPM 0044
C    $ (LARGE(N),N=1,2)/ZDFFFFFFFF,ZFFFFFFFFF/,                         DPM 0045
C    $ (RIGHT(N),N=1,2)/ZCC5800000,Z000000000/,                         DPM 0046
C    $ (DIVER(N),N=1,2)/ZCC6800000,Z000000000/,                         DPM 0047
C    $ (LOG10(N),N=1,2)/ZD00E730E7,ZC77800DC0/, SC/987/                 DPM 0048
                                                                        DPM 0049
C BURROUGHS 5700 SYSTEM.                                                DPM 0050
                                                                        DPM 0051
C     DATA (SMALL(N),N=1,2)/O1771000000000000,O0000000000000000/,       DPM 0052
C    $  (LARGE(N),N=1,2)/O0777777777777777,O0007777777777777/,          DPM 0053
C    $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,          DPM 0054
C    $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,          DPM 0055
C    $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/  DPM 0056
                                                                        DPM 0057
C BURROUGHS 6700/7700 SYSTEMS.                                          DPM 0058
                                                                        DPM 0059
C     DATA (SMALL(N),N=1,2)/O1771000000000000,O7770000000000000/,       DPM 0060
C    $  (LARGE(N),N=1,2)/O0777777777777777,O7777777777777777/,          DPM 0061
C    $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,          DPM 0062
C    $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,          DPM 0063
C    $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/  DPM 0064
                                                                        DPM 0065
C FTN4 ON THE CDC 6000/7000 SERIES.                                     DPM 0066
                                                                        DPM 0067
C     DATA                                                              DPM 0068
C    $  (SMALL(N),N=1,2)/00564000000000000000B,00000000000000000000B/,  DPM 0069
C    $  (LARGE(N),N=1,2)/37757777777777777777B,37157777777777777774B/,  DPM 0070
C    $  (RIGHT(N),N=1,2)/15624000000000000000B,00000000000000000000B/,  DPM 0071
C    $  (DIVER(N),N=1,2)/15634000000000000000B,00000000000000000000B/,  DPM 0072
C    $  (LOG10(N),N=1,2)/17164642023241175717B,16367571421742254654B/,  DPM 0073
C    $  SC/987/                                                         DPM 0074
                                                                        DPM 0075
C FTN5 ON THE CDC 6000/7000 SERIES.                                     DPM 0076
                                                                        DPM 0077
C     DATA                                                              DPM 0078
C    $(SMALL(N),N=1,2)/O"00564000000000000000",O"00000000000000000000"/,DPM 0079
C    $(LARGE(N),N=1,2)/O"37757777777777777777",O"37157777777777777774"/,DPM 0080
C    $(RIGHT(N),N=1,2)/O"15624000000000000000",O"00000000000000000000"/,DPM 0081
C    $(DIVER(N),N=1,2)/O"15634000000000000000",O"00000000000000000000"/,DPM 0082
C    $(LOG10(N),N=1,2)/O"17164642023241175717",O"16367571421742254654"/,DPM 0083
C    $ SC/987/                                                          DPM 0084
                                                                        DPM 0085
C CONVEX C-1                                                            DPM 0086
                                                                        DPM 0087
C     DATA (SMALL(N),N=1,2)/'00100000'X,'00000000'X/,                   DPM 0088
C    $  (LARGE(N),N=1,2)/'7FFFFFFF'X,'FFFFFFFF'X/,                      DPM 0089
C    $  (RIGHT(N),N=1,2)/'3CC00000'X,'00000000'X/,                      DPM 0090
C    $  (DIVER(N),N=1,2)/'3CD00000'X,'00000000'X/,                      DPM 0091
C    $  (LOG10(N),N=1,2)/'3FF34413'X,'509F79FF'X/, SC/987/              DPM 0092
                                                                        DPM 0093
C CRAY 1, XMP, 2, AND 3.                                                DPM 0094
                                                                        DPM 0095
C     DATA                                                              DPM 0096
C    $ (SMALL(N),N=1,2)/201354000000000000000B,000000000000000000000B/, DPM 0097
C    $ (LARGE(N),N=1,2)/577767777777777777777B,000007777777777777776B/, DPM 0098
C    $ (RIGHT(N),N=1,2)/376434000000000000000B,000000000000000000000B/, DPM 0099
C    $ (DIVER(N),N=1,2)/376444000000000000000B,000000000000000000000B/, DPM 0100
C    $ (LOG10(N),N=1,2)/377774642023241175717B,000007571421742254654B/, DPM 0101
C    $ SC/987/                                                          DPM 0102
                                                                        DPM 0103
C DATA GENERAL ECLIPSE S/200                                            DPM 0104
C NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -          DPM 0105
C STATIC DMACH(5)                                                       DPM 0106
                                                                        DPM 0107
C     DATA SMALL/20K,3*0/, LARGE/77777K,3*177777K/,                     DPM 0108
C    $  RIGHT/31420K,3*0/, DIVER/32020K,3*0/,                           DPM 0109
C    $  LOG10/40423K,42023K,50237K,74776K/, SC/987/                     DPM 0110
                                                                        DPM 0111
C HARRIS SLASH 6 AND SLASH 7                                            DPM 0112
                                                                        DPM 0113
C     DATA (SMALL(N),N=1,2)/'20000000,'00000201/,                       DPM 0114
C    $  (LARGE(N),N=1,2)/'37777777,'37777577/,                          DPM 0115
C    $  (RIGHT(N),N=1,2)/'20000000,'00000333/,                          DPM 0116
C    $  (DIVER(N),N=1,2)/'20000000,'00000334/,                          DPM 0117
C    $  (LOG10(N),N=1,2)/'23210115,'10237777/, SC/987/                  DPM 0118
                                                                        DPM 0119
C HONEYWELL DPS 8/70 SERIES.                                            DPM 0120
                                                                        DPM 0121
C     DATA (SMALL(N),N=1,2)/O402400000000,O000000000000/,               DPM 0122
C    $  (LARGE(N),N=1,2)/O376777777777,O777777777777/,                  DPM 0123
C    $  (RIGHT(N),N=1,2)/O604400000000,O000000000000/,                  DPM 0124
C    $  (DIVER(N),N=1,2)/O606400000000,O000000000000/,                  DPM 0125
C    $  (LOG10(N),N=1,2)/O776464202324,O117571775714/, SC/987/          DPM 0126
                                                                        DPM 0127
C IBM 360/370 SERIES, XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.      DPM 0128
                                                                        DPM 0129
C     DATA (SMALL(N),N=1,2)/Z00100000,Z00000000/,                       DPM 0130
C    $  (LARGE(N),N=1,2)/Z7FFFFFFF,ZFFFFFFFF/,                          DPM 0131
C    $  (RIGHT(N),N=1,2)/Z33100000,Z00000000/,                          DPM 0132
C    $  (DIVER(N),N=1,2)/Z34100000,Z00000000/,                          DPM 0133
C    $  (LOG10(N),N=1,2)/Z41134413,Z509F79FF/, SC/987/                  DPM 0134
                                                                        DPM 0135
C INTERDATA 8/32 WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.              DPM 0136
C FOR THE INTERDATA FORTRAN VII COMPILER REPLACE                        DPM 0137
C THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.                            DPM 0138
                                                                        DPM 0139
C     DATA (SMALL(N),N=1,2)/Z'00100000',Z'00000000'/,                   DPM 0140
C    $  (LARGE(N),N=1,2)/Z'7EFFFFFF',Z'FFFFFFFF'/,                      DPM 0141
C    $  (RIGHT(N),N=1,2)/Z'33100000',Z'00000000'/,                      DPM 0142
C    $  (DIVER(N),N=1,2)/Z'34100000',Z'00000000'/,                      DPM 0143
C    $  (LOG10(N),N=1,2)/Z'41134413',Z'509F79FF'/, SC/987/              DPM 0144
                                                                        DPM 0145
C PDP-10 (KA PROCESSOR).                                                DPM 0146
                                                                        DPM 0147
C     DATA (SMALL(N),N=1,2)/"033400000000,"000000000000/,               DPM 0148
C    $  (LARGE(N),N=1,2)/"377777777777,"344777777777/,                  DPM 0149
C    $  (RIGHT(N),N=1,2)/"113400000000,"000000000000/,                  DPM 0150
C    $  (DIVER(N),N=1,2)/"114400000000,"000000000000/,                  DPM 0151
C    $  (LOG10(N),N=1,2)/"177464202324,"144117571776/, SC/987/          DPM 0152
                                                                        DPM 0153
C PDP-10 (KI PROCESSOR).                                                DPM 0154
                                                                        DPM 0155
C     DATA (SMALL(N),N=1,2)/"000400000000,"000000000000/,               DPM 0156
C    $  (LARGE(N),N=1,2)/"377777777777,"377777777777/,                  DPM 0157
C    $  (RIGHT(N),N=1,2)/"103400000000,"000000000000/,                  DPM 0158
C    $  (DIVER(N),N=1,2)/"104400000000,"000000000000/,                  DPM 0159
C    $  (LOG10(N),N=1,2)/"177464202324,"047674776746/, SC/987/          DPM 0160
                                                                        DPM 0161
C PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS                            DPM 0162
C (EXPRESSED IN INTEGER AND OCTAL).                                     DPM 0163
                                                                        DPM 0164
C     DATA (SMALL(N),N=1,2)/8388608,0/, (LARGE(N),N=1,2)/2147483647,-1/,DPM 0165
C    $  (RIGHT(N),N=1,2)/612368384,0/, (DIVER(N),N=1,2)/620756992,0/,   DPM 0166
C    $  (LOG10(N),N=1,2)/1067065498,-2063872008/, SC/987/               DPM 0167
                                                                        DPM 0168
C     DATA (SMALL(N),N=1,2)/O00040000000,O00000000000/,                 DPM 0169
C    $  (LARGE(N),N=1,2)/O17777777777,O37777777777/,                    DPM 0170
C    $  (RIGHT(N),N=1,2)/O04440000000,O00000000000/,                    DPM 0171
C    $  (DIVER(N),N=1,2)/O04500000000,O00000000000/,                    DPM 0172
C    $  (LOG10(N),N=1,2)/O07746420232,O20476747770/, SC/987/            DPM 0173
                                                                        DPM 0174
C PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS                            DPM 0175
C (EXPRESSED IN INTEGER AND OCTAL).                                     DPM 0176
                                                                        DPM 0177
C     DATA SMALL/128,3*0/, LARGE/32767,3*-1/, RIGHT/9344,3*0/,          DPM 0178
C    $  DIVER/9472,3*0/, LOG10/16282,8346,-31493,-12296/, SC/987/       DPM 0179
                                                                        DPM 0180
C     DATA SMALL/O000200,3*O000000/, LARGE/O077777,3*O177777/,          DPM 0181
C    $  RIGHT/O022200,3*O000000/, DIVER/O022400,3*O000000/,             DPM 0182
C    $  LOG10/O037632,O020232,O102373,O147770/, SC/987/                 DPM 0183
                                                                        DPM 0184
C PRIME 50 SERIES SYSTEMS WITH 32-BIT INTEGERS AND 64V MODE             DPM 0185
C INSTRUCTIONS, SUPPLIED BY IGOR BRAY.                                  DPM 0186
                                                                        DPM 0187
C     DATA (SMALL(N),N=1,2)/:10000000000,:00000100001/,                 DPM 0188
C    $  (LARGE(N),N=1,2)/:17777777777,:37777677775/,                    DPM 0189
C    $  (RIGHT(N),N=1,2)/:10000000000,:00000000122/,                    DPM 0190
C    $  (DIVER(N),N=1,2)/:10000000000,:00000000123/,                    DPM 0191
C    $  (LOG10(N),N=1,2)/:11504046501,:07674600177/, SC/987/            DPM 0192
                                                                        DPM 0193
C SEQUENT BALANCE 8000                                                  DPM 0194
                                                                        DPM 0195
C     DATA (SMALL(N),N=1,2)/$00000000, $00100000/,                      DPM 0196
C    $  (LARGE(N),N=1,2)/$FFFFFFFF, $7FEFFFFF/,                         DPM 0197
C    $  (RIGHT(N),N=1,2)/$00000000, $3CA00000/,                         DPM 0198
C    $  (DIVER(N),N=1,2)/$00000000, $3CB00000/,                         DPM 0199
C    $  (LOG10(N),N=1,2)/$509F79FF, $3FD34413/, SC/987/                 DPM 0200
                                                                        DPM 0201
C UNIVAC 1100 SERIES.                                                   DPM 0202
                                                                        DPM 0203
C     DATA (SMALL(N),N=1,2)/O000040000000,O000000000000/,               DPM 0204
C    $  (LARGE(N),N=1,2)/O377777777777,O777777777777/,                  DPM 0205
C    $  (RIGHT(N),N=1,2)/O170540000000,O000000000000/,                  DPM 0206
C    $  (DIVER(N),N=1,2)/O170640000000,O000000000000/,                  DPM 0207
C    $  (LOG10(N),N=1,2)/O177746420232,O411757177572/, SC/987/          DPM 0208
                                                                        DPM 0209
C VAX UNIX F77 COMPILER                                                 DPM 0210
                                                                        DPM 0211
C     DATA (SMALL(N),N=1,2)/128,0/, (LARGE(N),N=1,2)/-32769,-1/,        DPM 0212
C    $  (RIGHT(N),N=1,2)/9344,0/, (DIVER(N),N=1,2)/9472,0/,             DPM 0213
C    $  (LOG10(N),N=1,2)/546979738,-805796613/, SC/987/                 DPM 0214
                                                                        DPM 0215
C VAX-11 WITH FORTRAN IV-PLUS COMPILER                                  DPM 0216
                                                                        DPM 0217
C     DATA (SMALL(N),N=1,2)/Z00000080,Z00000000/,                       DPM 0218
C    $  (LARGE(N),N=1,2)/ZFFFF7FFF,ZFFFFFFFF/,                          DPM 0219
C    $  (RIGHT(N),N=1,2)/Z00002480,Z00000000/,                          DPM 0220
C    $  (DIVER(N),N=1,2)/Z00002500,Z00000000/,                          DPM 0221
C    $  (LOG10(N),N=1,2)/Z209A3F9A,ZCFF884FB/, SC/987/                  DPM 0222
                                                                        DPM 0223
C VAX/VMS VERSION 2.2                                                   DPM 0224
                                                                        DPM 0225
C     DATA (SMALL(N),N=1,2)/'80'X,'0'X/,                                DPM 0226
C    $  (LARGE(N),N=1,2)/'FFFF7FFF'X,'FFFFFFFF'X/,                      DPM 0227
C    $  (RIGHT(N),N=1,2)/'2480'X,'0'X/, (DIVER(N),N=1,2)/'2500'X,'0'X/, DPM 0228
C    $  (LOG10(N),N=1,2)/'209A3F9A'X,'CFF884FB'X/, SC/987/              DPM 0229
                                                                        DPM 0230
      IF( PASS1 )  THEN                                                 DPM 0231
                                                                        DPM 0232
         PASS1 = .FALSE.                                                DPM 0233
         IF (SC.NE.987)                                                 DPM 0234
     $       CALL ERRMSG( 'D1MACH--NO DATA STATEMENTS ACTIVE',.TRUE.)   DPM 0235
                                                                        DPM 0236
C                        ** CALCULATE MACHINE PRECISION                 DPM 0237
         EPSNEW = 0.01D0                                                DPM 0238
   10    EPS = EPSNEW                                                   DPM 0239
            EPSNEW = EPSNEW / 1.1D0                                     DPM 0240
C                               ** IMPORTANT TO STORE 'S' SINCE MAY BE  DPM 0241
C                               ** KEPT IN HIGHER PRECISION IN REGISTERSDPM 0242
            S = 1.D0 + EPSNEW                                           DPM 0243
            IF( S.GT.1.D0 ) GO TO 10                                    DPM 0244
         IF( EPS/DMACH(4).LT.0.5D0 .OR. EPS/DMACH(4).GT.2.D0 )          DPM 0245
     $       CALL ERRMSG( 'D1MACH--TABULATED PRECISION WRONG',.TRUE.)   DPM 0246
                                                                        DPM 0247
      END IF                                                            DPM 0248
                                                                        DPM 0249
      IF (I.LT.1.OR.I.GT.5)                                             DPM 0250
     $    CALL ERRMSG( 'D1MACH--ARGUMENT OUT OF BOUNDS',.TRUE.)         DPM 0251
      D1MACH = DMACH(I)                                                 DPM 0252
      RETURN                                                            DPM 0253
      END                                                               DPM 0254
