      FUNCTION   XTERP(XCC,X,Y,NDEG,NPTS,DINT,IER)                      XTR 0001
C                                                                       XTR 0002
C   FUNCTION PERFORMS NEWTONS INTERPOLATION FOR DISCRETE DATA           XTR 0003
C            AS A FUNCTION OF ONE VARIABLE                              XTR 0004
C                                                                       XTR 0005
C   WHERE  XC - INDEPENDENT VARIABLE AT WHICH THE INTERPOLATED VALUE    XTR 0006
C               OF THE DEPENDENT VARIABLE IS DESIRED                    XTR 0007
C           X - TABLE OF INDEPENDENT VARIABLE VALUES IN INCREASING      XTR 0008
C               ORDER                                                   XTR 0009
C           Y - CORRESPONDING TABLE OF DEPENDENT VARIABLE VALUES        XTR 0010
C        NDEG - ORDER OF THE INTERPOLATING POLYNOMIAL USED (MAX - 10)   XTR 0011
C        NPTS - NUMBER OF ENTRIES IN X AND Y                            XTR 0012
C        DINT - DERIVITIVE AT XC                                        XTR 0013
C         IER - RETURN CODE:                                            XTR 0014
C                      0 = INTERPOLATION WAS PERFORMED                  XTR 0015
C                     -1 = EXTRAPOLATION BELOW TABLE VALUES             XTR 0016
C                      1 = EXTRAPOLATION ABOVE TABLE VALUES             XTR 0017
C                                                                       XTR 0018
C  ROUTINE MODIFIED FROM 'THE COMPUTING TECHNOLOGY CENTER NUMERICAL     XTR 0019
C    ANALYSIS LIBRARY', O.R.N.L.                                        XTR 0020
C                                                                       XTR 0021
      COMMON /IFIL/IRD,IPR,IPU,NPR,IPR1,ISCRCH                          XTR 0022
      DIMENSION X(52),Y(52),Y1(11),PI(12)                               XTR 0023
      INTEGER HI                                                        XTR 0024
      DATA PI/12*1./                                                    XTR 0025
C                                                                       XTR 0026
      XC=XCC                                                            XTR 0027
      NFIT=NDEG + 1                                                     XTR 0028
      N=NFIT                                                            XTR 0029
      MFLAG=0                                                           XTR 0030
      IF(N .GT. NPTS) N=NPTS                                            XTR 0031
      IF(XC-X(1)) 50,20,10                                              XTR 0032
10    IF(XC .GT. X(NPTS)) GO TO 70                                      XTR 0033
20    IER=0                                                             XTR 0034
      DO 30 I=1,NPTS                                                    XTR 0035
      IF(XC - X(I)) 40,120,30                                           XTR 0036
120   XC=XC + .000001                                                   XTR 0037
      MFLAG=I                                                           XTR 0038
30    CONTINUE                                                          XTR 0039
40    LOW=I - (N+1)/2                                                   XTR 0040
      IF(LOW .LE. 0) GO TO 60                                           XTR 0041
      HI=LOW + N - 1                                                    XTR 0042
      IF(HI .GT. NPTS) GO TO 80                                         XTR 0043
      GO TO 90                                                          XTR 0044
C                                                                       XTR 0045
C ... XC LT X(1)                                                        XTR 0046
50    IER=-1                                                            XTR 0047
60    HI=N                                                              XTR 0048
      LOW=1                                                             XTR 0049
      GO TO 90                                                          XTR 0050
C                                                                       XTR 0051
C ... XC GT X(NPTS)                                                     XTR 0052
70    IER=1                                                             XTR 0053
80    LOW=NPTS - N + 1                                                  XTR 0054
      HI=NPTS                                                           XTR 0055
C                                                                       XTR 0056
C ... INTERPOLATE                                                       XTR 0057
90    CON=1.                                                            XTR 0058
      DINT=0.                                                           XTR 0059
      Y1(1)=Y(LOW)                                                      XTR 0060
      XTERP=Y(LOW)                                                      XTR 0061
      IM=LOW - 1                                                        XTR 0062
      IL=LOW + 1                                                        XTR 0063
      DO 110 K=IL,HI                                                    XTR 0064
      VAL=XTERP                                                         XTR 0065
      IA=K - LOW                                                        XTR 0066
      IS=IA + 1                                                         XTR 0067
      Y1(IS)=Y(K)                                                       XTR 0068
      DO 100 I=1,IA                                                     XTR 0069
      IR=IM + I                                                         XTR 0070
      IF(X(IR).EQ.X(K))GO TO 100                                        XTR 0071
      Y1(IS) = (Y1(I)-Y1(IS))/(X(IR)-X(K))                              XTR 0072
  100 CONTINUE                                                          XTR 0073
      CON=CON*(XC-X(K-1))                                               XTR 0074
      PI(IA+1)=CON                                                      XTR 0075
      CON1=PI(IA)                                                       XTR 0076
      IF(IA .LT. 2) GO TO 112                                           XTR 0077
      DO 111 I=2,IA                                                     XTR 0078
      IF(PI(I).EQ.0.0)GO TO 111                                         XTR 0079
      CON1=CON1 + CON*PI(I-1)/PI(I)                                     XTR 0080
  111 CONTINUE                                                          XTR 0081
  112 DINT=DINT + CON1*Y1(IS)                                           XTR 0082
  110 XTERP=VAL + CON*Y1(IS)                                            XTR 0083
      IF(MFLAG.NE.0)XTERP=Y(MFLAG)                                      XTR 0084
      RETURN                                                            XTR 0085
      END                                                               XTR 0086
