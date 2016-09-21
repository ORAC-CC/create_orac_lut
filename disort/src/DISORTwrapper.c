/* File:    DISORTwrapper.c
 * Purpose: Wrappers for the Fortran DISORT routines in
 *          DISORTfunctions.f
 * Author:  Steven Marsh
 * Date:    August 2001
 * Modifications:
   Gareth Thomas, 08 Feb 2013.
           DISORT and GETMOM now called as integer functions, rather
           than void. The return value is checked and a "1" is taken to
           represent a fatal error detected within the Fortran. If an
           error has occured, execution is halted using the
	   IDL_Message(); mechanism. This goes hand-in-hand with the
           alterations to error handling made in DISORTfunctions.F
   Gareth Thomas, 10 May 2013.
           Added PLKAVG to the list of callable routines.
 */

/* ANSI */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* IDL */
#include "export.h"

/* Local */
#include "DISORTtoIDL.h"
#include "DISORTfunctions.h"

/* IDL procedure protos */
extern void IDL_CDECL disort(int argc, IDL_VPTR argv[], char *argk);
extern void IDL_CDECL getmom(int argc, IDL_VPTR argv[], char *argk);
extern IDL_VPTR IDL_CDECL plkavg(int argc, IDL_VPTR argv[], char *argk);

/* IDL version 5.3 and greater use IDL_SYSFUN_DEF2 while earlier
 * versions use IDL_SYSFUN_DEF.  Note the addition of the final '0' in
 * each line for IDL_SYSFUN_DEF2. */

#ifdef __IDLPRE53__
/* FOR IDL < 5.3 */
/* Define the procedures */
static IDL_SYSFUN_DEF DISORTtoIDL_procedures[] = {
  {(IDL_FUN_RET) disort,"DISORT",44,44,0},
  {(IDL_FUN_RET) getmom,"GETMOM",4,4,0},
};
static IDL_SYSFUN_DEF DISORTtoIDL_functions[] = {
  {plkavg,"PLKAVG",3,3,0},
};
#else
/* FOR IDL >= 5.3 */
/* Define the procedures */
static IDL_SYSFUN_DEF2 DISORTtoIDL_procedures[] = {
  {(IDL_FUN_RET) disort,"DISORT",44,44,0,0},
  {(IDL_FUN_RET) getmom,"GETMOM",4,4,0,0},
};
static IDL_SYSFUN_DEF2 DISORTtoIDL_functions[] = {
  {plkavg,"PLKAVG",3,3,0,0},
};

#endif

/* Startup call when DLM is loaded */
int DISORTtoIDL_Startup(void)
{

  /* Register the error handler */
  IDL_ExitRegister(DISORTtoIDL_exit_handler);

  /* If IDL is pre-5.3 then change IDL_SysRtnAdd to IDL_AddSystemRoutine,
   * (NB: the parameters stay the same) */

#ifdef __IDLPRE53__
  /* FOR IDL < 5.3 */
  /* Add procedures */
  if (!IDL_AddSystemRoutine(DISORTtoIDL_procedures, FALSE,
			    ARRLEN(DISORTtoIDL_procedures))
      && !IDL_AddSystemRoutine(DISORTtoIDL_functions, TRUE,
			       ARRLEN(DISORTtoIDL_functions))) {
    return IDL_FALSE;
  }

#else
  /* FOR IDL >= 5.3 */
  /* Add procedures */
  return IDL_SysRtnAdd(DISORTtoIDL_procedures, FALSE,
		       ARRLEN(DISORTtoIDL_procedures))
    && IDL_SysRtnAdd(DISORTtoIDL_functions, TRUE,
		     ARRLEN(DISORTtoIDL_functions));

#endif

}

/* Called when IDL is shutdown */
void DISORTtoIDL_exit_handler(void) {
    /* Nothing special to do in this case */
}

/* Wrapped Fortran routines below this point */



/* =======================================================================
 * IDL Function:   disort
 * Description:    Discrete Ordinance Radiative Transfer Calculations)
 * =======================================================================
 */
/* C does not have a logical type so define a comparable boolian type */
#ifdef BOOL
#undef BOOL
#endif
#define BOOL int

/* Define the boolian false */
#ifdef FALSE
#undef FALSE
#endif
#define FALSE   (0)

/* Define the boolian true */
#ifdef TRUE
#undef TRUE
#endif
#define TRUE    (1)

void IDL_CDECL disort(int argc, IDL_VPTR argv[], char *argk) {

    /* Local */
  IDL_MEMINT RFLDIRdim[2],RFLDNdim[2],FLUPdim[2],DFDTdim[2],UAVGdim[2],UUdim[3],ALBMEDdim[2],TRNMEDdim[2];

  IDL_VPTR ivRFLDIRArray,ivRFLDNArray,ivFLUPArray,ivDFDTArray,ivUAVGArray,ivUUArray,ivALBMEDArray,ivTRNMEDArray;

  IDL_LONG NLYR,NMOM,NTAU,NSTR,NUMU,NPHI,IBCND,MAXCLY,MAXULV,MAXUMU,MAXPHI,MAXMOM;

  /* IDLBool_t USRTAU,USRANG,LAMBER,PLANK,ONLYFL;
  IDLBool_t *PRNT; */

  IDL_LONG USRTAU,USRANG,LAMBER,PLANK,ONLYFL;
  IDL_LONG *PRNT;

  float WVNMLO,WVNMHI,FBEAM,UMU0,PHI0,FISOT,ALBEDO,BTEMP,TTEMP,TEMIS,ACCUR;

  float *DTAUC,*SSALB,*PMOM,*TEMPER,*UTAU,*UMU,*PHI,*RFLDIR,*RFLDN,*FLUP,*DFDT,*UAVG,*UU,*ALBMED,*TRNMED;

  int ndims, numElmnts,i, disort_error;

  disort_error = 0;

  /*-----------------------------------------------NLYR */
  IDL_ENSURE_SCALAR(argv[0]);
  NLYR = IDL_LongScalar(argv[0]);
  /* To change maximum number of layers edit MAXCLY in
   DISORTfunctions.F  */
  if (NLYR > 100){
     IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"NLYR must be 100 or less");
  }

  /*-----------------------------------------------DTAUC*/
  IDL_ENSURE_ARRAY(argv[1]);
  ndims = argv[1]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[1]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"DTAUC array must only contain ONE dimension");

  /* Check the array is a FLOAT type */
  if (argv[1]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"DTAUC must be of FLOAT precision");

  /* Get the passed array */
  DTAUC = (float *) argv[1]->value.arr->data;

  /*-----------------------------------------------SSALB*/
  IDL_ENSURE_ARRAY(argv[2]);
  ndims = argv[2]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[2]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"SSALB array must only contain ONE dimension");

  /* Check the array is a FLOAT type */
  if (argv[2]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"SSALB must be of FLOAT precision");
  /* Get the passed array */
  SSALB = (float *) argv[2]->value.arr->data;

 /* Check all values between 0.0 and 0.9*/
  for (i=0;i<numElmnts;i++){
    if (SSALB[i] > 1.0)
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		  "SSALB contains values greater than 1.0");
  };

  /*-----------------------------------------------NMOM*/
  IDL_ENSURE_SCALAR(argv[3]);
  NMOM = IDL_LongScalar(argv[3]);

  /*-----------------------------------------------PMOM*/
  IDL_ENSURE_ARRAY(argv[4]);
  ndims = argv[4]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[4]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  /*  if (ndims != NLYR)
     IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
  		"PMOM array must contain NLYR dimensions"); */

  /* Check the array is a FLOAT type */
  if (argv[4]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"PMOM must be of FLOAT precision");
  /* Get the passed array */
  PMOM = (float *) argv[4]->value.arr->data;

/*-----------------------------------------------TEMPER*/
  IDL_ENSURE_ARRAY(argv[5]);
  ndims = argv[5]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[5]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"TEMPER array must only contain ONE dimension");

  /* Check the array is a FLOAT type */
  if (argv[5]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"TEMPER must be of FLOAT precision");
  /* Get the passed array */
  TEMPER = (float *) argv[5]->value.arr->data;

  /*-----------------------------------------------WVNMLO*/
  IDL_ENSURE_SCALAR(argv[6]);
  WVNMLO = argv[6]->value.f;
  /* printf("%f",WVNMLO); */

  /*-----------------------------------------------WVNMHI*/
  IDL_ENSURE_SCALAR(argv[7]);
  WVNMHI = argv[7]->value.f;

  /*-----------------------------------------------USRTAU*/
  IDL_ENSURE_SCALAR(argv[8]);
  USRTAU =  IDL_LongScalar(argv[8]);

  /*-----------------------------------------------NTAU*/
  IDL_ENSURE_SCALAR(argv[9]);
  NTAU = IDL_LongScalar(argv[9]);

  /*-----------------------------------------------UTAU*/
  IDL_ENSURE_ARRAY(argv[10]);
  ndims = argv[10]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[10]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"UTAU array must only contain ONE dimension");

  /* Check the array is a FLOAT type */
  if (argv[10]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"UTAU must be of FLOAT precision");
  /* Get the passed array */
  UTAU = (float *) argv[10]->value.arr->data;

  /*-----------------------------------------------NSTR*/
  IDL_ENSURE_SCALAR(argv[11]);
  NSTR = IDL_LongScalar(argv[11]);

  /*-----------------------------------------------USRANG*/
  IDL_ENSURE_SCALAR(argv[12]);
  USRANG = IDL_LongScalar(argv[12]);

  /*-----------------------------------------------NUMU*/
  IDL_ENSURE_SCALAR(argv[13]);
  NUMU = IDL_LongScalar(argv[13]);

  /*-----------------------------------------------UMU*/
  IDL_ENSURE_ARRAY(argv[14]);
  ndims = argv[14]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[14]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"UMU array must only contain ONE dimension");

  /* Check the array is a FLOAT type */
  if (argv[14]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"UMU must be of FLOAT precision");
  /* Get the passed array */
  UMU = (float *) argv[14]->value.arr->data;

  /*-----------------------------------------------NPHI*/
  IDL_ENSURE_SCALAR(argv[15]);
  NPHI = IDL_LongScalar(argv[15]);

  /*-----------------------------------------------PHI*/
  IDL_ENSURE_ARRAY(argv[16]);
  ndims = argv[16]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[16]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"PHI array must only contain ONE dimension");

  /* Check the array is a FLOAT type */
  if (argv[16]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"PHI must be of FLOAT precision");
  /* Get the passed array */
  PHI = (float *) argv[16]->value.arr->data;

  /*-----------------------------------------------IBCND*/
  IDL_ENSURE_SCALAR(argv[17]);
  IBCND = IDL_LongScalar(argv[17]);

  /*-----------------------------------------------FBEAM*/
  IDL_ENSURE_SCALAR(argv[18]);
  FBEAM = argv[18]->value.f;

  /*-----------------------------------------------UMU0*/
  IDL_ENSURE_SCALAR(argv[19]);/*UMU0*/
  UMU0 = argv[19]->value.f;

  /*-----------------------------------------------PHI0*/
  IDL_ENSURE_SCALAR(argv[20]);
  PHI0 = argv[20]->value.f;

  /*-----------------------------------------------FISOT*/
  IDL_ENSURE_SCALAR(argv[21]);
  FISOT = argv[21]->value.f;

  /*-----------------------------------------------LAMBER*/
  IDL_ENSURE_SCALAR(argv[22]);
  LAMBER = IDL_LongScalar(argv[22]);

  /*-----------------------------------------------ALBEDO*/
  IDL_ENSURE_SCALAR(argv[23]);
  ALBEDO = argv[23]->value.f;

  /*-----------------------------------------------BTEMP*/
  IDL_ENSURE_SCALAR(argv[24]);
  BTEMP = argv[24]->value.f;

  /*-----------------------------------------------TTEMP*/
  IDL_ENSURE_SCALAR(argv[25]);
  TTEMP = argv[25]->value.f;

  /*-----------------------------------------------TEMIS*/
  IDL_ENSURE_SCALAR(argv[26]);
  TEMIS = argv[26]->value.f;

  /*-----------------------------------------------PLANK*/
  IDL_ENSURE_SCALAR(argv[27]);
  PLANK = IDL_LongScalar(argv[27]);

  /*-----------------------------------------------ONLYFL*/
  IDL_ENSURE_SCALAR(argv[28]);
  ONLYFL = IDL_LongScalar(argv[28]);

  /*-----------------------------------------------ACCUR*/
  IDL_ENSURE_SCALAR(argv[29]);
  ACCUR = argv[29]->value.f;

  /*-----------------------------------------------PRNT*/
  IDL_ENSURE_ARRAY(argv[30]);
  ndims = argv[30]->value.arr->n_dim;      /* dimension */
  numElmnts = argv[30]->value.arr->n_elts; /* total number of elements */

  /* Ensure only a single dimension */
  if (ndims != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"PRNT array must only contain ONE dimension");

  /* Get the passed array */
  PRNT = (IDL_LONG *) argv[30]->value.arr->data;

  /* I forced the following five variables to ensure there was
   no additional space within matrices which was not use.
   Unless I did this the IDL call would crash with strange
   errors. */

  /*-----------------------------------------------MAXCLY
  IDL_ENSURE_SCALAR(argv[31]);
  MAXCLY = IDL_LongScalar(argv[31]);*/
  MAXCLY = NLYR; /*ensure matrices exact right size*/

  /*-----------------------------------------------MAXULV
  IDL_ENSURE_SCALAR(argv[32]);
    MAXULV = IDL_LongScalar(argv[32]);*/
  MAXULV = NTAU; /*ensure matrices exact right size*/

  /*-----------------------------------------------MAXUMU
  IDL_ENSURE_SCALAR(argv[33]);
   MAXUMU = IDL_LongScalar(argv[33]);*/
  MAXUMU = NUMU; /*ensure matrices exact right size*/

  /*-----------------------------------------------MAXPHI
  IDL_ENSURE_SCALAR(argv[34]);
  MAXPHI = IDL_LongScalar(argv[34]);*/
  MAXPHI = NPHI; /*ensure matrices exact right size*/
  /*-----------------------------------------------MAXMOM
  IDL_ENSURE_SCALAR(argv[35]);
   MAXMOM = IDL_LongScalar(argv[35]);*/
  MAXMOM = NMOM;/*ensure matrices exact right size*/

  /*-----------------------------------------------RFLDIR
    IDL_ENSURE_ARRAY(argv[36]);*/
  RFLDIRdim[0] = NTAU;
  RFLDIRdim[1] = 1; /*single dimension*/
  RFLDIR = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, RFLDIRdim,
  				       IDL_ARR_INI_ZERO, &ivRFLDIRArray);

  /*-----------------------------------------------/RFLDN
    IDL_ENSURE_ARRAY(argv[37]);*/
  RFLDNdim[0] = NTAU;
  RFLDNdim[1] = 1; /*single dimension*/
  RFLDN = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, RFLDNdim,
				       IDL_ARR_INI_ZERO, &ivRFLDNArray);

  /*-----------------------------------------------FLUP
    IDL_ENSURE_ARRAY(argv[38]);*/
  FLUPdim[0] = NTAU;
  FLUPdim[1] = 1; /*single dimension*/
  FLUP = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, FLUPdim,
				       IDL_ARR_INI_ZERO, &ivFLUPArray);

  /*----------------------------------------------- DFDT
    IDL_ENSURE_ARRAY(argv[39]);*/
  DFDTdim[0] = NTAU;
  DFDTdim[1] = 1; /*single dimension*/
  DFDT = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, DFDTdim,
				       IDL_ARR_INI_ZERO, &ivDFDTArray);

  /*-----------------------------------------------UAVG
    IDL_ENSURE_ARRAY(argv[40]);*/
  UAVGdim[0] = NTAU;
  UAVGdim[1] = 1; /*single dimension*/
  UAVG = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, UAVGdim,
				       IDL_ARR_INI_ZERO, &ivUAVGArray);

  /*-----------------------------------------------UU
    IDL_ENSURE_ARRAY(argv[41]);*/
  UUdim[0] = NUMU;
  UUdim[1] = NTAU;
  UUdim[2] = NPHI;  /*three dimensions*/
  UU = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 3, UUdim,
				       IDL_ARR_INI_ZERO, &ivUUArray);

  /*----------------------------------------------- ALBMED
    IDL_ENSURE_ARRAY(argv[42]);*/
  ALBMEDdim[0] = NUMU;
  ALBMEDdim[1] = 1; /*single dimension*/
  ALBMED = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, ALBMEDdim,
				       IDL_ARR_INI_ZERO, &ivALBMEDArray);

  /*-----------------------------------------------/TRNMED
    IDL_ENSURE_ARRAY(argv[43]);*/
  TRNMEDdim[0] = NUMU;
  TRNMEDdim[1] = 1; /*single dimension*/
  TRNMED = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, TRNMEDdim,
				       IDL_ARR_INI_ZERO, &ivTRNMEDArray);
  /*-----------------------------------------------*/

  /* Call Fortran */
  disort_error = disort_(&NLYR, DTAUC, SSALB, &NMOM, PMOM, TEMPER,
	  &WVNMLO, &WVNMHI, &USRTAU, &NTAU, UTAU, &NSTR, &USRANG,
          &NUMU, UMU, &NPHI, PHI, &IBCND, &FBEAM, &UMU0, &PHI0,
  	  &FISOT, &LAMBER, &ALBEDO, &BTEMP, &TTEMP, &TEMIS,
  	  &PLANK, &ONLYFL, &ACCUR, PRNT, &MAXCLY,
  	  &MAXULV, &MAXUMU, &MAXPHI, &MAXMOM, RFLDIR, RFLDN,
  	  FLUP, DFDT, UAVG, UU, ALBMED, TRNMED );



/* Copy the IDL_VPTRs to argv */
  IDL_VarCopy(ivRFLDIRArray, argv[36]);/*RFLDIR*/
  IDL_VarCopy(ivRFLDNArray, argv[37]);/*RFLDN*/
  IDL_VarCopy(ivFLUPArray, argv[38]);/*FLUP*/
  IDL_VarCopy(ivDFDTArray, argv[39]);/*DFDT*/
  IDL_VarCopy(ivUAVGArray, argv[40]);/*UAVG*/
  IDL_VarCopy(ivUUArray, argv[41]);/*UU*/
  IDL_VarCopy(ivALBMEDArray, argv[42]);/*ALBMED*/
  IDL_VarCopy(ivTRNMEDArray, argv[43]);/*TRNMED*/

  if (disort_error == 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"DISORT encountered an error and has aborted!");


 }


/* =======================================================================
 * IDL Function:   getmom
 * Description:    Calculates phase moments of aerosol for use in DISORT
 * =======================================================================
 */
void IDL_CDECL getmom(int argc, IDL_VPTR argv[], char *argk) {

  /* Local */
  IDL_LONG IPHAS,NMOM;
  float GG, *PMOMARRAY;
  IDL_MEMINT dim[2];
  IDL_VPTR ivPmomArray;
  short getmom_error;

  getmom_error = 0;

  /* IDL inputs: Ensure we get scalars from IDL */
  IDL_ENSURE_SCALAR(argv[0]);
  IDL_ENSURE_SCALAR(argv[1]);
  IDL_ENSURE_SCALAR(argv[2]);

  if (argv[1]->type != IDL_TYP_FLOAT)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
	      "G must be FLOAT precision");

  IPHAS = IDL_LongScalar(argv[0]);
  GG = argv[1]->value.f;
 if (fabs(GG) >= 1.0){
     IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"abs GG must be less than 1.0");
  }
  NMOM=IDL_LongScalar(argv[2]);

  dim[0]=NMOM+1;
  dim[1]=1; /* single dimension */

  PMOMARRAY = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 1, dim,
					  IDL_ARR_INI_ZERO, &ivPmomArray);

  /* Call Fortran */
  getmom_error = getmom_(&IPHAS, &GG, &NMOM, PMOMARRAY);

  /* Copy the IDL_VPTR ivStatsArray to argv[3] */
  IDL_VarCopy(ivPmomArray, argv[3]);

  if (getmom_error == 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"GETMOM encountered an error and has aborted!");

}


/* =======================================================================
 * IDL Function:   plkavg
 * Description:    Calculates Plank emission across a wavenumber interval
 * =======================================================================
 */
IDL_VPTR IDL_CDECL plkavg(int argc, IDL_VPTR argv[], char *argk) {
  float wnlo, wnhi, temp;
  IDL_VPTR tmp0, tmp1, tmp2;
  IDL_VPTR ivplank;

  /* IDL inputs: Ensure we get scalars from IDL */
  IDL_ENSURE_SCALAR(argv[0]);
  IDL_ENSURE_SCALAR(argv[1]);
  IDL_ENSURE_SCALAR(argv[2]);

  /* Ensure that the input arguments are converted into floats */
  tmp0 = IDL_CvtFlt(1, &argv[0]);
  tmp1 = IDL_CvtFlt(1, &argv[1]);
  tmp2 = IDL_CvtFlt(1, &argv[2]);

  /* Extract data into local C variables */
  wnlo = tmp0->value.f;
  wnhi = tmp1->value.f;
  temp = tmp2->value.f;

  /* Create a temporary variable to pass back to IDL */
  ivplank = IDL_Gettmp();
  ivplank->type = IDL_TYP_FLOAT;

  /* Call Fortran */
  ivplank->value.f = plkavg_(&wnlo, &wnhi, &temp);

  /* Clean up local IDL variables created.. */
  IDL_DELTMP(tmp0);
  IDL_DELTMP(tmp1);
  IDL_DELTMP(tmp2);

  /* Return to IDL */
  return(ivplank);

}
