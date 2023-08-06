      SUBROUTINE MSRAD(GROUND,UANG,NSTR,V,S0,KNTRVL)                    MSR 0001
C                                                                       MSR 0002
C     ROUTINE MSRAD SETS UP OPTICAL PROPERTIES PROFILES FOR VERTICAL    MSR 0003
C     PATH THEN CALLS DISORT WHICH RETURNS MS SOURCE FUNCTIONS.         MSR 0004
C                                                                       MSR 0005
C     LIST PARAMETERS:                                                  MSR 0006
      INCLUDE 'PARAM.LST'                                               MSR 0007
      INTEGER MAXULV,MAXUMU,MAXPHI,MAXCLY,MAXCMU,MXULV                  MSR 0008
      PARAMETER(MAXULV=LAYDIM,MAXUMU=1,MAXPHI=1,MAXCLY=LAYDIM,          MSR 0009
     1  MAXCMU=16,MXULV=LAYDIM+1)                                       MSR 0010
C                                                                       MSR 0011
C     DECLARE INPUTS:                                                   MSR 0012
      LOGICAL GROUND                                                    MSR 0013
      INTEGER NSTR,KNTRVL                                               MSR 0014
      REAL V,S0                                                         MSR 0015
      DOUBLE PRECISION UANG                                             MSR 0016
C                                                                       MSR 0017
C     LIST COMMONS:                                                     MSR 0018
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT                     MSR 0019
      REAL TBOUND,SALB                                                  MSR 0020
      LOGICAL MODTRN                                                    MSR 0021
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB,   MSR 0022
     1  MODTRN                                                          MSR 0023
C                                                                       MSR 0024
C       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.       MSR 0025
C       UPFLX    LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].    MSR 0026
C       DNFLX    LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].  MSR 0027
C       UPFLXS   BOUNDARY UPWARD SCATTERED SOLAR FLUX [W CM-2 / CM-1].  MSR 0028
C       DNFLXS   BOUNDARY DOWNWARD SCATTERED SOLAR FLUX [W CM-2 / CM-1].MSR 0029
C       NTFLX    LAYER BOUNDARY NET (THERMAL PLUS SCATTERED SOLAR       MSR 0030
C                PLUS DIRECT SOLAR) UPWARD FLUX [W CM-2 / CM-1].        MSR 0031
      REAL SUBINT,UPFLX,DNFLX,UPFLXS,DNFLXS,NTFLX                       MSR 0032
      COMMON/NETFLX/SUBINT(NKSUB),UPFLX(LAYDIM),DNFLX(LAYDIM),          MSR 0033
     1  UPFLXS(LAYDIM),DNFLXS(LAYDIM),NTFLX(LAYDIM)                     MSR 0034
      INTEGER IKMAX,ML,ISSGEO,IMULT                                     MSR 0035
      COMMON/CNTRL/IKMAX,ML,ISSGEO,IMULT                                MSR 0036
C                                                                       MSR 0037
C     COMMON /MSRD/                                                     MSR 0038
C       CSZEN0  LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.            MSR 0039
C       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.             MSR 0040
C       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING               MSR 0041
C               (AWAY FROM EARTH) THE CURRENT LAYER.                    MSR 0042
C       BBGRND  THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].  MSR 0043
C       BBNDRY  LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1]. MSR 0044
C       COSBAR  LAYER HENYEY-GREENSTEIN ASYMMETRY FACTOR.               MSR 0045
C       TSCAT   LAYER SCATTERING OPTICAL DEPTH.                         MSR 0046
C       TCONT   LAYER CONTINUUM OPTICAL DEPTH.                          MSR 0047
C       TAUT    LAYER TOTAL OPTICAL DEPTH.                              MSR 0048
C       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.  MSR 0049
C       S0DEP   OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.               MSR 0050
C       S0TRN   TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]           MSR 0051
C       UPF     LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].     MSR 0052
C       DNF     LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].   MSR 0053
C       UPFS    LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].       MSR 0054
C       DNFS    LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].     MSR 0055
      REAL CSZEN0,CSZEN,CSZENX,BBGRND,BBNDRY,COSBAR,TSCAT,              MSR 0056
     1  TCONT,TAUT,DEPRAT,S0DEP,S0TRN,UPF,DNF,UPFS,DNFS                 MSR 0057
      COMMON/MSRD/CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),          MSR 0058
     1  BBGRND,BBNDRY(LAYDIM),COSBAR(LAYDIM),TSCAT(LAYDIM),             MSR 0059
     2  TCONT(LAYDIM),TAUT(NKSUB,LAYDIM),DEPRAT(LAYDIM),                MSR 0060
     3  S0DEP(NKSUB,LAYDIM),S0TRN(NKSUB,LAYDIM),UPF(NKSUB,LAYDIM),      MSR 0061
     4  DNF(NKSUB,LAYDIM),UPFS(NKSUB,LAYDIM),DNFS(NKSUB,LAYDIM)         MSR 0062
      INTEGER IV1,IV2,IDV,IFWHM                                         MSR 0063
      COMMON/CARD4/IV1,IV2,IDV,IFWHM                                    MSR 0064
      REAL ZM,PM,TM,RFNDX,DENSTY                                        MSR 0065
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    MSR 0066
     1  RFNDX(LAYDIM),DENSTY(KMAX,LAYDIM)                               MSR 0067
      DOUBLE PRECISION ALBEDO,BTEMP,FBEAM,TTEMP,WVNMLO,WVNMHI,UMU0,WN0  MSR 0068
      DOUBLE PRECISION DTAUC(MAXCLY),HL(0:MAXCMU),PHI(MAXPHI),          MSR 0069
     1  PMOM(0:MAXCMU,1:MAXCLY),SSALB(MAXCLY),TEMPER(0:MAXCLY),         MSR 0070
     2  UMU(MAXUMU),UTAU(MAXULV),RFLDIR(MAXULV),RFLDN(MAXULV),          MSR 0071
     3  FLDN(MXULV),FLUP(MAXULV),UAVG(MAXULV),DFDT(MAXULV),             MSR 0072
     4  U0U(MAXUMU,MAXULV),UU(MAXUMU,MAXULV,MAXPHI),ALBMED(MAXUMU),     MSR 0073
     5  TRNMED(MAXUMU),S0CMS(MAXUMU,MAXULV),T0CMS(MAXUMU,MAXULV)        MSR 0074
C                                                                       MSR 0075
C     DECLARE LOCAL VARIABLES                                           MSR 0076
      INTEGER IK,N,KWI,INTRVL,NEXT                                      MSR 0077
C                                                                       MSR 0078
C     DEFINE DATA                                                       MSR 0079
C       PRNT     DISORT PRINT OPTIONS                                   MSR 0080
C       ONLYFL   .FALSE. FOR USING DISORT FOR RADIATIVE                 MSR 0081
C                TRANSFER COMPUTATIONS (IT SHOULD BE LEFT               MSR 0082
C                AS ONLYFL=.TRUE. FOR THE MERGED PACKAGE).              MSR 0083
C       DELTAM   DELTA M APPROXIMATION FLAG (SET DELTAM                 MSR 0084
C                TO .FALSE. TO PERFORM CALCULATIONS                     MSR 0085
C                ANALOGOUS TO THOSE IN BMFLUX/FLXADD).                  MSR 0086
      LOGICAL USRTAU,USRANG,LAMBER,MSFLAG,PLANK,ONLYFL,DELTAM,PRNT(7)   MSR 0087
      INTEGER IBCND,NUMU,NPHI                                           MSR 0088
      DOUBLE PRECISION ACCUR,FISOT,TEMIS,PHI0,DPDEG                     MSR 0089
      CHARACTER HEADER*127                                              MSR 0090
      SAVE USRTAU,USRANG,LAMBER,MSFLAG,PLANK,ONLYFL,DELTAM,             MSR 0091
     1  PRNT,IBCND,NUMU,NPHI,ACCUR,FISOT,TEMIS,PHI0,HEADER              MSR 0092
      DATA USRTAU/.FALSE./,USRANG/.TRUE./,LAMBER/.TRUE./,MSFLAG/.TRUE./,MSR 0093
     1  PLANK/.TRUE./,ONLYFL/.TRUE./,DELTAM/.TRUE./,PRNT/7*.FALSE./     MSR 0094
      DATA IBCND/0/,NUMU/1/,NPHI/1/                                     MSR 0095
      DATA ACCUR/0./,FISOT/0./,TEMIS/0./,PHI0/0./,DPDEG/57.2957795131D0/MSR 0096
      DATA HEADER/'*'/                                                  MSR 0097
C                                                                       MSR 0098
C     DISORT INTERFACE                                                  MSR 0099
      TEMPER(0)=DBLE(TM(ML))                                            MSR 0100
      UTAU(1)=DBLE(0.)                                                  MSR 0101
      IK=ML                                                             MSR 0102
      DO 20 N=1,IKMAX                                                   MSR 0103
          IK=IK-1                                                       MSR 0104
C                                                                       MSR 0105
C         SET UP TEMPERATURE PROFILE                                    MSR 0106
          TEMPER(N)=DBLE(TM(IK))                                        MSR 0107
C                                                                       MSR 0108
C         CLEAN OUT ARRAYS CONTAINING SOURCE FUNCTIONS                  MSR 0109
          T0CMS(1,N)=DBLE(0.)                                           MSR 0110
          S0CMS(1,N)=DBLE(0.)                                           MSR 0111
C                                                                       MSR 0112
C         CALCULATE OPTICAL DEPTHS                                      MSR 0113
          DTAUC(N)=DBLE(TAUT(1,IK))                                     MSR 0114
          UTAU(N+1)=UTAU(N)+DTAUC(N)                                    MSR 0115
C                                                                       MSR 0116
C         CALCULATE SINGLE SCATTERING ALBEDO                            MSR 0117
          SSALB(N)=DBLE(0.)                                             MSR 0118
          IF(TSCAT(IK).LT.0.)TSCAT(IK)=0.                               MSR 0119
          IF(DTAUC(N).GT.0.)SSALB(N)=DBLE(TSCAT(IK))/DTAUC(N)           MSR 0120
          IF(SSALB(N).GT.1.)SSALB(N)=DBLE(1.)                           MSR 0121
C                                                                       MSR 0122
C         H-G FUNCTION MOMENTS OF LEGENDRE POLYNOMIALS                  MSR 0123
          PMOM(0,N)=DBLE(1.)                                            MSR 0124
          DO 10 KWI=1,NSTR                                              MSR 0125
               PMOM(KWI,N)=DBLE(COSBAR(IK)**KWI)                        MSR 0126
   10     CONTINUE                                                      MSR 0127
   20 CONTINUE                                                          MSR 0128
C                                                                       MSR 0129
C     OTHER DISORT VARIABLES DEFINED                                    MSR 0130
      WN0=DBLE(V)                                                       MSR 0131
      WVNMLO=DBLE(V-.5*IDV)                                             MSR 0132
      WVNMHI=DBLE(V+.5*IDV)                                             MSR 0133
      UMU(1)=-COS(UANG/DPDEG)                                           MSR 0134
      BTEMP=TEMPER(IKMAX)                                               MSR 0135
      IF(GROUND .AND. TBOUND.GT.0.)BTEMP=DBLE(TBOUND)                   MSR 0136
      TTEMP=TEMPER(0)                                                   MSR 0137
      FBEAM=DBLE(S0)                                                    MSR 0138
      UMU0=DBLE(CSZEN0(1))                                              MSR 0139
      ALBEDO=DBLE(SALB)                                                 MSR 0140
      PHI(1)=DBLE(0.)                                                   MSR 0141
C                                                                       MSR 0142
C     LOOP OVER CORRELATED-K SUB-INTERVALS                              MSR 0143
      INTRVL=1                                                          MSR 0144
      DO 40 NEXT=2,KNTRVL                                               MSR 0145
C                                                                       MSR 0146
C         CALL TO ROUTINE DISORT                                        MSR 0147
          CALL DISORT(IKMAX,DTAUC,SSALB,PMOM,TEMPER,WVNMLO,WVNMHI,      MSR 0148
     1      USRTAU,ML,UTAU,NSTR,USRANG,NUMU,UMU,NPHI,PHI,IBCND,FBEAM,   MSR 0149
     2      UMU0,PHI0,FISOT,LAMBER,ALBEDO,HL,BTEMP,TTEMP,TEMIS,DELTAM,  MSR 0150
     3      PLANK,ONLYFL,ACCUR,PRNT,HEADER,MAXCLY,MAXULV,MAXUMU,MAXCMU, MSR 0151
     4      MAXPHI,RFLDIR,RFLDN,FLDN,FLUP,DFDT,UAVG,UU,U0U,ALBMED,      MSR 0152
     5      TRNMED,MSFLAG,WN0,S0CMS,T0CMS,DNFS(INTRVL,1),DNF(INTRVL,1)) MSR 0153
C                                                                       MSR 0154
C         LOOP OVER LAYERS                                              MSR 0155
          UTAU(1)=DBLE(0.)                                              MSR 0156
          IK=ML                                                         MSR 0157
          DO 30 N=1,IKMAX                                               MSR 0158
              UPFLX(IK)=UPFLX(IK)+SUBINT(INTRVL)*REAL(FLUP(N))          MSR 0159
              DNFLX(IK)=DNFLX(IK)+SUBINT(INTRVL)*REAL(FLDN(N))          MSR 0160
              NTFLX(IK)=NTFLX(IK)+SUBINT(INTRVL)*REAL(FLUP(N)-FLDN(N))  MSR 0161
              IK=IK-1                                                   MSR 0162
C                                                                       MSR 0163
C             STORE THERMAL AND SOLAR SOURCE FUNCTIONS IN UPF AND UPFS. MSR 0164
              UPF(INTRVL,N)=REAL(T0CMS(1,N))                            MSR 0165
              UPFS(INTRVL,N)=REAL(S0CMS(1,N))                           MSR 0166
C                                                                       MSR 0167
C             CLEAN OUT ARRAYS CONTAINING SOURCE FUNCTIONS.             MSR 0168
              T0CMS(1,N)=DBLE(0.)                                       MSR 0169
              S0CMS(1,N)=DBLE(0.)                                       MSR 0170
C                                                                       MSR 0171
C             CALCULATE OPTICAL DEPTHS                                  MSR 0172
              DTAUC(N)=DBLE(TAUT(NEXT,IK))                              MSR 0173
              UTAU(N+1)=UTAU(N)+DTAUC(N)                                MSR 0174
C                                                                       MSR 0175
C             CALCULATE SINGLE SCATTERING ALBEDO                        MSR 0176
              SSALB(N)=DBLE(0.)                                         MSR 0177
              IF(TSCAT(IK).LT.0.)TSCAT(IK)=0.                           MSR 0178
              IF(DTAUC(N).GT.0.)SSALB(N)=DBLE(TSCAT(IK))/DTAUC(N)       MSR 0179
              IF(SSALB(N).GT.1.)SSALB(N)=DBLE(1.)                       MSR 0180
   30     CONTINUE                                                      MSR 0181
          UPFLX(1)=UPFLX(1)+SUBINT(INTRVL)*REAL(FLUP(ML))               MSR 0182
          DNFLX(1)=DNFLX(1)+SUBINT(INTRVL)*REAL(FLDN(ML))               MSR 0183
          NTFLX(1)=NTFLX(1)+SUBINT(INTRVL)*REAL(FLUP(ML)-FLDN(ML))      MSR 0184
          INTRVL=NEXT                                                   MSR 0185
   40 CONTINUE                                                          MSR 0186
C                                                                       MSR 0187
C     LAST CORRELATED-K SUB-INTRVL.  CALL TO ROUTINE DISORT.            MSR 0188
      CALL DISORT(IKMAX,DTAUC,SSALB,PMOM,TEMPER,WVNMLO,WVNMHI,          MSR 0189
     1  USRTAU,ML,UTAU,NSTR,USRANG,NUMU,UMU,NPHI,PHI,IBCND,FBEAM,       MSR 0190
     2  UMU0,PHI0,FISOT,LAMBER,ALBEDO,HL,BTEMP,TTEMP,TEMIS,DELTAM,      MSR 0191
     3  PLANK,ONLYFL,ACCUR,PRNT,HEADER,MAXCLY,MAXULV,MAXUMU,MAXCMU,     MSR 0192
     4  MAXPHI,RFLDIR,RFLDN,FLDN,FLUP,DFDT,UAVG,UU,U0U,ALBMED,          MSR 0193
     5  TRNMED,MSFLAG,WN0,S0CMS,T0CMS,DNFS(KNTRVL,1),DNF(KNTRVL,1))     MSR 0194
C                                                                       MSR 0195
C     LOOP OVER LAYERS                                                  MSR 0196
      IK=ML                                                             MSR 0197
      DO 50 N=1,IKMAX                                                   MSR 0198
          UPFLX(IK)=UPFLX(IK)+SUBINT(KNTRVL)*REAL(FLUP(N))              MSR 0199
          DNFLX(IK)=DNFLX(IK)+SUBINT(KNTRVL)*REAL(FLDN(N))              MSR 0200
          NTFLX(IK)=NTFLX(IK)+SUBINT(KNTRVL)*REAL(FLUP(N)-FLDN(N))      MSR 0201
          IK=IK-1                                                       MSR 0202
C                                                                       MSR 0203
C         STORE THERMAL AND SOLAR SOURCE FUNCTIONS IN UPF AND UPFS.     MSR 0204
          UPF(KNTRVL,N)=REAL(T0CMS(1,N))                                MSR 0205
          UPFS(KNTRVL,N)=REAL(S0CMS(1,N))                               MSR 0206
   50 CONTINUE                                                          MSR 0207
      UPFLX(1)=UPFLX(1)+SUBINT(KNTRVL)*REAL(FLUP(ML))                   MSR 0208
      DNFLX(1)=DNFLX(1)+SUBINT(KNTRVL)*REAL(FLDN(ML))                   MSR 0209
      NTFLX(1)=NTFLX(1)+SUBINT(KNTRVL)*REAL(FLUP(ML)-FLDN(ML))          MSR 0210
C                                                                       MSR 0211
C     RETURN TO LOOP                                                    MSR 0212
      RETURN                                                            MSR 0213
      END                                                               MSR 0214
