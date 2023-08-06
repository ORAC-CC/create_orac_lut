      BLOCK DATA XSECTB                                                 XSC 0001
      IMPLICIT DOUBLE PRECISION (V)                                     XSC 0002
C                                                                       XSC 0003
C**   XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES        XSC 0004
C                                                                       XSC 0005
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                  XSC 0006
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)             XSC 0007
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),         XSC 0008
     X                NTEMPF(5,35),NSPECR(35),IXFORM(35),NUMXS          XSC 0009
C                                                                       XSC 0010
      DATA (ALIAS(1,I),I=1,35)/                                         XSC 0011
     1    'CLONO2    ', 'HNO4      ', 'CHCL2F    ', 'CCL4      ',       XSC 0012
     2    'CCL3F     ', 'CCL2F2    ', 'C2CL2F4   ', 'C2CL3F3   ',       XSC 0013
     3    'N2O5      ', 'HNO3      ', 'CF4', 'CHF2CL   ',               XSC 0014
     4    23*' ZZZZZZZZ '/                                              XSC 0015
      DATA (ALIAS(2,I),I=1,35)/                                         XSC 0016
     1    'CLNO3     ', 'ZZZZZZZZ  ', 'ZZZZZZZZ  ', 'ZZZZZZZZ  ',       XSC 0017
     2    'CFCL3     ', 'CF2CL2    ', 'C2F4CL2   ', 'C2F3CL3   ',       XSC 0018
     3    'ZZZZZZZZ  ', 'ZZZZZZZZ  ', 'ZZZZZZZZ   ', 'ZZZZZZZ  ',       XSC 0019
     4    23*' ZZZZZZZZ '/                                              XSC 0020
      DATA (ALIAS(3,I),I=1,35)/                                         XSC 0021
     1    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',       XSC 0022
     2    'CFC11     ', 'CFC12     ', 'CFC114    ', 'CFC113    ',       XSC 0023
     3    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC14     ',                     XSC 0024
     4    24*' ZZZZZZZZ ' /                                             XSC 0025
      DATA (ALIAS(4,I),I=1,35)/                                         XSC 0026
     1    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F21       ', ' ZZZZZZZZ ',       XSC 0027
     2    'F11       ', 'F12       ', 'F114      ', 'F113      ',       XSC 0028
     3    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F14       ',                     XSC 0029
     4    24*' ZZZZZZZZ ' /                                             XSC 0030
C                                                                       XSC 0031
      DATA V1FX/175*0.0/,V2FX/175*0.0/,DVFX/175*0.0/,WXM/35*0.0/        XSC 0032
      DATA NTEMPF/175*0/,NSPECR/35*0/,IXFORM/35*0/,NUMXS/0/             XSC 0033
C                                                                       XSC 0034
      END                                                               XSC 0035
