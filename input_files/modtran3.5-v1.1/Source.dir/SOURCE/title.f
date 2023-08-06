      BLOCK DATA TITLE                                                  TIT 0001
C>    BLOCK DATA                                                        TIT 0002
C     TITLE INFORMATION                                                 TIT 0003
      CHARACTER*4  HHAZE      ,HSEASN     ,HVULCN     ,BLANK,           TIT 0004
     X            HMET        ,HMODEL     ,HTRRAD                       TIT 0005
      COMMON /TITL/ HHAZE(5,16),HSEASN(5,2),HVULCN(5,8),BLANK,          TIT 0006
     X HMET(5,2),HMODEL(5,8),HTRRAD(6,4)                                TIT 0007
      COMMON /VSBD/ VSB(10)                                             TIT 0008
      DATA VSB /23.,5.,0.,23.,5.,50.,23.,0.2,0.5,0./                    TIT 0009
      DATA BLANK/'    '/                                                TIT 0010
      DATA HHAZE /                                                      TIT 0011
     1 'RURA','L   ','    ','    ','    ',                              TIT 0012
     2 'RURA','L   ','    ','    ','    ',                              TIT 0013
     3 'NAVY',' MAR','ITIM','E   ','    ',                              TIT 0014
     4 'MARI','TIME','    ','    ','    ',                              TIT 0015
     5 'URBA','N   ','    ','    ','    ',                              TIT 0016
     6 'TROP','OSPH','ERIC','    ','    ',                              TIT 0017
     7 'USER',' DEF','INED','    ','    ',                              TIT 0018
     8 'FOG1',' (AD','VECT','TION',')   ',                              TIT 0019
     9 'FOG2','(RAD','IATI','0N) ','    ',                              TIT 0020
     X 'DESE','RT A','EROS','OL  ','    ',                              TIT 0021
     A 'BACK','GROU','ND S','TRAT','O   ',                              TIT 0022
     B 'AGED',' VOL','CANI','C   ','    ',                              TIT 0023
     C 'FRES','H VO','LCAN','IC  ','    ',                              TIT 0024
     D 'AGED',' VOL','CANI','C   ','    ',                              TIT 0025
     E 'FRES','H VO','LCAN','IC  ','    ',                              TIT 0026
     F 'METE','ORIC',' DUS','T   ','    '/                              TIT 0027
      DATA HSEASN /                                                     TIT 0028
     1 'SPRI','NG-S','UMME','R   ','    ',                              TIT 0029
     2 'FALL','-WIN','TER ','    ','    ' /                             TIT 0030
      DATA HVULCN /                                                     TIT 0031
     1 'BACK','GROU','ND S','TRAT','O   ',                              TIT 0032
     2 'MODE','RATE',' VOL','CANI','C   ',                              TIT 0033
     3 'HIGH','    ',' VOL','CANI','C   ',                              TIT 0034
     4 'HIGH','    ',' VOL','CANI','C   ',                              TIT 0035
     5 'MODE','RATE',' VOL','CANI','C   ',                              TIT 0036
     6 'MODE','RATE',' VOL','CANI','C   ',                              TIT 0037
     7 'HIGH','    ',' VOL','CANI','C   ',                              TIT 0038
     8 'EXTR','EME ',' VOL','CANI','C   '/                              TIT 0039
      DATA HMET/                                                        TIT 0040
     1 'NORM','AL  ','    ','    ','    ',                              TIT 0041
     2 'TRAN','SITI','ON  ','    ','    '/                              TIT 0042
      DATA HMODEL /                                                     TIT 0043
     1 'TROP','ICAL',' MOD','EL  ','    ',                              TIT 0044
     2 'MIDL','ATIT','UDE ','SUMM','ER  ',                              TIT 0045
     3 'MIDL','ATIT','UDE ','WINT','ER  ',                              TIT 0046
     4 'SUBA','RCTI','C   ','SUMM','ER  ',                              TIT 0047
     5 'SUBA','RCTI','C   ','WINT','ER  ',                              TIT 0048
     6 '1976',' U S',' STA','NDAR','D   ',                              TIT 0049
     7 '   ','    ','    ','    ','    ',                               TIT 0050
     8 'MODE','L =0','HORI','ZONT','AL  '/                              TIT 0051
      DATA HTRRAD/                                                      TIT 0052
     1 'TRAN','SMIT','TANC','E   ','    ','    ',                       TIT 0053
     2 'RADI','ANCE','    ','    ','    ','    ',                       TIT 0054
     3 'RADI','ANCE','+SOL','AR S','CATT','ERNG',                       TIT 0055
     4 'TRAN','SMIT','TED ','SOLA','R IR','RAD.'/                       TIT 0056
      END                                                               TIT 0057
