/* These match the corresponding Fortran functions *
 * Author:    Steven Marsh
 * Date:      August 2001
 * Modifications:
   Gareth Thomas, 08 Feb 2013.
              DISORT and GETMOM now called as integer functions, rather
              than void.
   Gareth Thomas, 10 May 2013.
              Added PLKAVG function
 */
#ifndef DISORTfunctions_h
#define DISORTfunctions_h

extern short disort_(IDL_LONG* NLYR, float* DTAUC, float* SSALB, IDL_LONG* NMOM, float* PMOM, float* TEMPER, float* WVNMLO, float* WVNMHI, IDL_LONG* USRTAU, IDL_LONG* NTAU, float* UTAU, IDL_LONG* NSTR, IDL_LONG* USRANG, IDL_LONG* NUMU, float* UMU, IDL_LONG* NPHI, float* PHI, IDL_LONG* IBCND, float* FBEAM, float* UMU0, float* PHI0, float* FISOT, IDL_LONG* LAMBER, float* ALBEDO, float* BTEMP, float* TTEMP, float* TEMIS, IDL_LONG* PLANK, IDL_LONG* ONLYFL, float* ACCUR, IDL_LONG* PRNT, IDL_LONG* MAXCLY, IDL_LONG* MAXULV, IDL_LONG* MAXUMU, IDL_LONG* MAXPHI, IDL_LONG* MAXMOM, float* RFLDIR, float* RFLDN,float*  FLUP, float* DFDT, float* UAVG, float* UU, float* ALBMED, float* TRNMED);

extern short getmom_(IDL_LONG* iphas, float* gg, IDL_LONG* nmom, float* pmom);

extern float plkavg_(float* wnlo, float* wnhi, float* temp);

#endif
