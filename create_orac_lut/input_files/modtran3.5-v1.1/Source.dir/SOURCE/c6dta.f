      SUBROUTINE C6DTA(C6L,V)                                           C6D 0001
C     CALCULATES MOLECULAR RAYLEIGH SCATTERING COEFFICIENT              C6D 0002
C     USES APPROXIMATION OF SHETTLE ET AL., 1980 (APPL OPT, 2873-4)     C6D 0003
C        WITH THE DEPOLARIZATION = 0.0279 (INSTED OF 0.035)             C6D 0004
C     INPUT:  V = FREQUENCY IN WAVENUMBERS (CM-1)                       C6D 0005
C     OUTPUT: C6L = MOLECULAR SCATTERING COEFFICIENT (KM-1)             C6D 0006
C                    FOR TEMPERATURE = 273 K & PRESSURE 1 ATM.          C6D 0007
C                                                                       C6D 0008
      C6L=V**4/(9.38076E+18  -1.08426E+09 * V**2)                       C6D 0009
      RETURN                                                            C6D 0010
      END                                                               C6D 0011
