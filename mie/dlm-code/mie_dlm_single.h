/* Copyright (C) 1998-2017 University of Oxford
 *
 * This source code is licensed under the GNU General Public License (GPL),
 * Version 3.  See the file COPYING for more details.
 */


/* File:      mie_dlm_single.h
 * Purpose:   Header for mie_dlm_single.c
 * Author:    Gareth Thomas
 * Version:   1.1
 * Date:      June 2005
 */

#ifndef IDL_C_miesingle
#define IDL_C_miesingle

/* Message Numbers */
#define mie_single_ERROR          0
#define mie_single_NOSTRINGARRAY -1

/* Useful macro */
#define ARRLEN(arr) (sizeof(arr)/sizeof(arr[0]))

extern IDL_MSG_BLOCK msg_block;

/* Define the startup function that adds C functions to IDL
 * along with the exit handler */

/* IDL-Postres interface */
extern void mie_single_exit_handler(void);
extern int  mie_single_Startup(void);

/* Define the wrapper function which calls the fortran function */
extern void IDL_CDECL mie_dlm_single(int argc, IDL_VPTR argv[], char *argk);
extern void mieint_(const double *Dx, const IDL_DCOMPLEX *Cm, const IDL_LONG *inp, const double *Dqv, double *Qext, double *Qsca, double *Qbsc, double *g, double complex *S1, double complex *S2, double *F11, double *F33, double *F12, double *F34, IDL_LONG *mie_ok);
#endif
