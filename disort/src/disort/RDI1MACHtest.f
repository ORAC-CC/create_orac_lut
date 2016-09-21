
      PROGRAM Test_RDI1MACH

c     Test f90 versions of x1MACH (old machine constant routines)

      IMPLICIT NONE
      INTEGER :: I, I1MACH
      INTEGER, DIMENSION(16) :: IMACH_temp = 0
      REAL :: R1MACH
      DOUBLE PRECISION :: D1MACH
      REAL, DIMENSION(5) :: RMACH_temp
      DOUBLE PRECISION, DIMENSION(5) :: DMACH_temp

      DO I = 1, 5
         RMACH_temp(I) = R1MACH(I)
         DMACH_temp(I) = D1MACH(I)
      END DO

      DO I = 1, 16
         IF( I==3 .OR. I==6 ) CYCLE
         IMACH_temp(I) = I1MACH(I)
      END DO

      OPEN( 10, FILE=':Test_RDI1MACH.out' )
      WRITE(10,*) ' RMACH: ', RMACH_temp(1:5)
      WRITE(10,*) ' DMACH: ', DMACH_temp(1:5)
      WRITE(10,*) ' IMACH: ', IMACH_temp(1:16)
      CLOSE(10)

      STOP
      END PROGRAM Test_RDI1MACH


c --------------------------------------------------------------------
c Output on a Unix IEEE computer:
c
c  RMACH:  1.17549E-38 3.40282E+38 5.96046E-08 1.19209E-07 .30103
c
c  DMACH:  2.225073858507202-308  1.797693134862315+308
c          1.110223024625156E-16  2.220446049250313E-16  0.3010299956639811
c
c Output on a Mac Powerbook 3400c:
c
c RMACH:   1.175494E-38  3.402823E+38  5.960464E-08  1.192093E-07  .301030
c DMACH:   2.225073858507202E-308  1.797693134862314E+308
c           1.110223024625160E-016  2.220446049250313E-016  .301029995663981
c IMACH:   5  6  0  0  32  0  2  31  2147483647  2  24 -125  128  53 -1021 1024
c --------------------------------------------------------------------
