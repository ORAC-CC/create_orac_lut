      SUBROUTINE XPROFL                                                 XPR 0001
      INCLUDE 'PARAM.LST'                                               XPR 0002
C                                                                       XPR 0003
C********************************************************************** XPR 0004
C     THIS SUBROUTINE GENERATES THE DENSITY PROFILES OF THE CROSS-      XPR 0005
C     SECTION MOLECULES.  IT STORES THE PROFILES IN THE ARRAY DENM IN   XPR 0006
C     /DEAMT/ AT THE ALTITUDES ZMDL, WHICH ARE THE SAME ALTITUDES THAT  XPR 0007
C     THE PROFILES OF THE MOLECULAR AMOUNTS ARE DEFINED ON.  (NOTE: THE XPR 0008
C     ACTUAL ALTITUDES USED ARE FROM ZST WHICH IS A COPY OF ZMDL.)      XPR 0009
C     IPRFL IS A FLAG INDICATING THAT THE STANDARD PROFILES (0) OR A    XPR 0010
C     USER-INPUT PROFILE (1) IS TO BE USED.                             XPR 0011
C********************************************************************** XPR 0012
C                                                                       XPR 0013
C**   IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE   XPR 0014
C**   MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES     XPR 0015
C**   (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR   XPR 0016
C**   L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD       XPR 0017
C**   MOLECULES.                                                        XPR 0018
C                                                                       XPR 0019
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          XPR 0020
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(35),XAMNT(35,67)               XPR 0021
      COMMON/DEAMT/DENM(35,LAYTWO),DENP(35,LAYTHR+1)                    XPR 0022
C                                                                       XPR 0023
C**   COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES       XPR 0024
C**   FOR THE CROSS-SECTION MOLECULES.                                  XPR 0025
C**   XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES        XPR 0026
C                                                                       XPR 0027
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                  XPR 0028
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)             XPR 0029
C                                                                       XPR 0030
C**   AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH  XPR 0031
C**   LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL, LAYXMX LEVELS MAX     XPR 0032
C                                                                       XPR 0033
      COMMON /MLATMX/ LAYXMX,ALTX(50),AMOLX(50,35)                      XPR 0034
      COMMON /XPDEM/ZX(LAYDIM),DTMP(35),DENX(35,LAYDIM)                 XPR 0035
C                                                                       XPR 0036
C**   LOAD THE PROFILES OF ALTITUDE, PRESSURE, AND TEMPERATURE THAT     XPR 0037
C**   WERE USED TO CALCULATE THE MOLECULAR AMOUNTS BACK INTO THE        XPR 0038
C**   ARRAYS ZMDL, PM, AND TM FROM THE ARRAYS ZST, PST, AND TST         XPR 0039
C                                                                       XPR 0040
C                                                                       XPR 0041
C                                                                       XPR 0042
C**       A STANDARD PROFILE FOR X-MOLECULES DENSITY PROFILES HAS BEEN  XPR 0043
C**       SELECTED. THE PROFILES OF VOLUMNE MIXING RATIO ARE IN AMOLX   XPR 0044
C**       STORED AT THE LEVELS ALTX. LOAD THE ALTITUDES INTO ZX AND     XPR 0045
C**       DENX RESPECTIVELY.                                            XPR 0046
          IXMOLS = 3                                                    XPR 0047
          IXINDX(1) = 4                                                 XPR 0048
CC        CCL3F   IE  F11                                               XPR 0049
          IXINDX(2) = 5                                                 XPR 0050
CC        CCL2F2  IE  F12                                               XPR 0051
          IXINDX(3) = 6                                                 XPR 0052
C                                                                       XPR 0053
          LAYX=LAYXMX                                                   XPR 0054
          DO 210 L=1,LAYX                                               XPR 0055
              ZX(L) = ALTX(L)                                           XPR 0056
              DO 200 K=1,IXMOLS                                         XPR 0057
                  DENX(K,L) = AMOLX(L,IXINDX(K))                        XPR 0058
  200         CONTINUE                                                  XPR 0059
  210     CONTINUE                                                      XPR 0060
C                                                                       XPR 0061
C                                                                       XPR 0062
C**   INTERPOLATE THE DENSITY PROFILE DENX DEFINED ON ZX TO DENM        XPR 0063
C**   DEFINED ON ZMDL, THEN CONVERT MIXING RATIO TO NUMBER DENSITY.     XPR 0064
C                                                                       XPR 0065
      CALL XINTRP(LAYX,IXMOLS)                                          XPR 0066
C                                                                       XPR 0067
      RETURN                                                            XPR 0068
      END                                                               XPR 0069
