/* Copyright (C) 1998-2017 University of Oxford
 *
 * This source code is licensed under the GNU General Public License (GPL),
 * Version 3.  See the file COPYING for more details.
 */


/* File:       mie_dlm_single.c
 * Purpose:    Wrapper for the FORTRAN routine mieint in the file mieint.f
 * Author:     Gareth Thomas
 * Version:    1.1
 * Date:       June 2005
 *
 * G McGarragh, 29, Jul 2015: Fixed returning of the phase function:
 *    outargv[32] must have enough elements to hold the maximum possible
 *    number of positional arguments.
 * G McGarragh, 29, Jul 2015: Add support to output the 11, 33, 12, and
 *    34 elements of the 4x4 single scattering phase matrix F.  F11 is
 *    identical to the old phase function output.
 */

/* ANSI */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

/* OpenMP */
#ifdef _OPENMP
#include <omp.h>
#endif

/* IDL */
#include "export.h"

/* Local */
#include "mie_dlm_single.h"

#ifdef __IDLPRE53__
    /* FOR IDL < 5.3 */
    /* Define the procedures */
    static IDL_SYSFUN_DEF mie_single_procedures[] = {
        {{(IDL_FUN_RET) mie_dlm_single},"MIE_DLM_SINGLE",2,12,IDL_SYSFUN_DEF_F_KEYWORDS},
    };
#else
    /* FOR IDL >= 5.3 */
    /* Define the procedures */
    static IDL_SYSFUN_DEF2 mie_single_procedures[] = {
        {{(IDL_FUN_RET) mie_dlm_single},"MIE_DLM_SINGLE",2,12,IDL_SYSFUN_DEF_F_KEYWORDS,0},
    };
#endif


/* Startup call when DLM is loaded */
int IDL_Load(void)
{
    /* IDL version 5.3 and greater use IDL_SYSFUN_DEF2 while earlier
     * versions use IDL_SYSFUN_DEF.  Note the addition of the final '0' in
     * each line for IDL_SYSFUN_DEF2. */

    /* If IDL is pre-5.3 then change IDL_SysRtnAdd to IDL_AddSystemRoutine,
     * (NB: the parameters stay the same) */

#ifdef __IDLPRE53__
    /* FOR IDL < 5.3 */
    /* Add procedures */
    if (!IDL_AddSystemRoutine(mie_single_procedures, FALSE,
                              ARRLEN(mie_single_procedures))) {
        return IDL_FALSE;
    }

#else
    /* FOR IDL >= 5.3 */
    /* Add procedures */
    if (!IDL_SysRtnAdd(mie_single_procedures, FALSE,
                       ARRLEN(mie_single_procedures))) {
        return IDL_FALSE;
    }

#endif

    /* Register the error handler */
    IDL_ExitRegister(mie_single_exit_handler);
    return(IDL_TRUE);
}

/* Called when IDL is shutdown */
void mie_single_exit_handler(void)
{
/* Nothing special to do in this case */
}

/* Wrapped Fortran routines below this point */

/* =======================================================================
 * IDL Procedure:  mie_dlm_single
 * Description:    Performs Mie computation for an array of particle sizes
 *                 for a given refractive index. Expected values are:
 *                 Dx:        Array of size parameters (MUST be an ARRAY of
 *                            type DOUBLE, even if it has one element!!!)
 *                 Cm:        Complex refractive index (MUST be be a scalar
 *                            of type DOUBLE COMPLEX)
 *                 Keywords:
 *                 Dqv = Dqv  An ARRAY of cos(scattering angles) at which
 *                            intensity functions and phase function are
 *                            evaluated. If not specified, a single value
 *                            of 1.0 (forward scatter) is used.
 *                 Returns:
 *                 Qext:      Array of extinction efficiencies
 *                 Qsca:      Array of scattering efficiencies
 *                 Qbsc:      Array of back-scatter efficiencies
 *                 g:         Array of asymetry parameters
 *                 s1 and s2: Arrays of scattering functions (in theta
 *                            and size)
 *                 Phase:     Array of Phase function (one for each size)
 *
 * G Thomas, November 2003
 * =======================================================================
 */
void IDL_CDECL mie_dlm_single(int argc, IDL_VPTR argv[], char *argk)
{
    /* local */
    long i, j;
    IDL_LONG64 dims[2];
    double dqvtmp[1];
#ifdef _OPENMP
    int nthreads;
#endif
    /* Variables to be passed to and from Mie routine itself */
    IDL_LONG mie_ok;
    IDL_LONG Npts;
    IDL_LONG Inp;
    IDL_DCOMPLEX Cm;
    double complex *S1, *S2;
    double Qext, Qsca, Qbsc, g;
/*
    double *F11, *F33, *F12, *F34;
*/
    /* Variables to hold values for different sizes */
    double *DxARR, *QextARR, *QscaARR, *QbscARR, *gARR;
    IDL_DCOMPLEX *S1ARR, *S2ARR;
    double *F11ARR, *F33ARR, *F12ARR, *F34ARR;

    /* Arrays for passing back to IDL */
    IDL_VPTR ivQextArr, ivQscaArr, ivQbscArr, ivgArr, ivS1Arr, ivS2Arr,
             ivF11Arr, ivF33Arr, ivF12Arr, ivF34Arr;

    /* Definition of keyword parameter Dqv - description stored in Dqvdesc */
    /*                                     - data stored in Dqvdata */
    IDL_VPTR outargv[32];

    static int dqvexists;
    static double Dqvdata[10000];
    static IDL_KW_ARR_DESC Dqvdesc = {(char*)Dqvdata,1,10000,0};

    static int mthread;
    static int mthreadexists;

    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        {"DQV",    IDL_TYP_DOUBLE,1,IDL_KW_ARRAY,&dqvexists,IDL_CHARA(Dqvdesc)},
        {"MTHREAD",IDL_TYP_LONG,  1,IDL_KW_ZERO, &mthreadexists, IDL_CHARA(mthread)},
        {NULL}
    };

    /* Ensure all temporary stuff resulting from keyword is cleaned up */
    IDL_KWCleanup(IDL_KW_MARK);
    /* Extraction of keyword and normal parameters */
    IDL_KWGetParams(argc,argv,argk,kw_pars,outargv,1);
    /* Ensure the refractive index is given as a scalar */
    IDL_ENSURE_SCALAR(argv[1]);

    if (argv[0]->type != IDL_TYP_DOUBLE)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Size parameters must be DOUBLE precision");
    if (argv[1]->type != IDL_TYP_DCOMPLEX)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Refractive index must be DOUBLE precision COMPLEX");
    IDL_ENSURE_ARRAY(argv[0]);

    /* Assign values to local variables */
    mie_ok  = 0;
    Npts    = argv[0]->value.arr->n_elts;
    DxARR   = (double *) argv[0]->value.arr->data;
    Cm      = argv[1]->value.dcmp;
    QextARR = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                             IDL_ARR_INI_ZERO, &ivQextArr);
    QscaARR = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                             IDL_ARR_INI_ZERO, &ivQscaArr);
    QbscARR = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                             IDL_ARR_INI_ZERO, &ivQbscArr);
    gARR    = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                             IDL_ARR_INI_ZERO, &ivgArr);

    /* Set up parallel processing of the "mieint" loop: How many threads
       to use. */
#ifdef _OPENMP
    if (mthreadexists) {
        if (mthread < 1)
            nthreads = omp_get_max_threads();
        else
            nthreads = mthread;
    }
    else
        nthreads = 1;
    omp_set_num_threads(nthreads);
#endif
    /* Copy the complex refractive index value to the C
       complex<double> variable Cm_c */

    if (dqvexists)
    {
        /* Dqv is defined, so use the user supplied scattering angles */
        Inp = Dqvdesc.n;
        dims[0] = Inp; dims[1] = Npts;

        S1ARR  = (IDL_DCOMPLEX *)   IDL_MakeTempArray(IDL_TYP_DCOMPLEX, 2,
                                        dims, IDL_ARR_INI_ZERO, &ivS1Arr);
        S2ARR  = (IDL_DCOMPLEX *)   IDL_MakeTempArray(IDL_TYP_DCOMPLEX, 2,
                                        dims, IDL_ARR_INI_ZERO, &ivS2Arr);
        F11ARR = (double *)         IDL_MakeTempArray(IDL_TYP_DOUBLE, 2,
                                        dims, IDL_ARR_INI_ZERO, &ivF11Arr);
        F33ARR = (double *)         IDL_MakeTempArray(IDL_TYP_DOUBLE, 2,
                                        dims, IDL_ARR_INI_ZERO, &ivF33Arr);
        F12ARR = (double *)         IDL_MakeTempArray(IDL_TYP_DOUBLE, 2,
                                        dims, IDL_ARR_INI_ZERO, &ivF12Arr);
        F34ARR = (double *)         IDL_MakeTempArray(IDL_TYP_DOUBLE, 2,
                                        dims, IDL_ARR_INI_ZERO, &ivF34Arr);

#pragma omp parallel num_threads(nthreads) \
    private(j, mie_ok, Qext, Qsca, Qbsc, g, S1, S2)
/*
#pragma omp parallel num_threads(nthreads) \
    private(j, mie_ok, Qext, Qsca, Qbsc, g, S1, S2, F11, F33, F12, F34)
*/
{
        S1  = (double complex *) malloc(Inp * sizeof(double complex));
        S2  = (double complex *) malloc(Inp * sizeof(double complex));
/*
        F11 = (double *)         malloc(Inp * sizeof(double));
        F33 = (double *)         malloc(Inp * sizeof(double));
        F12 = (double *)         malloc(Inp * sizeof(double));
        F34 = (double *)         malloc(Inp * sizeof(double));
*/
#pragma omp for schedule(guided)
        for (i = 0; i < Npts; i++)
        {
            /* Call the mieint procedure */
/*
            mieint_(&DxARR[i], &Cm, &Inp, Dqvdata, &Qext, &Qsca, &Qbsc, &g,
                    S1, S2, F11, F33, F12, F34, &mie_ok);
*/
            mieint_(&DxARR[i], &Cm, &Inp, Dqvdata, &Qext, &Qsca, &Qbsc, &g,
                    S1, S2, F11ARR+i*Inp, F33ARR+i*Inp, F12ARR+i*Inp, F34ARR+i*Inp, &mie_ok);
            if (mie_ok != 0)
                IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                            "Size parameter overflow!");
            QextARR[i] = Qext;
            QscaARR[i] = Qsca;
            QbscARR[i] = Qbsc;
            gARR[i] = g;

            for (j = 0; j < Inp; j++)
            {
                S1ARR[i*Inp+j].r = creal(S1[j]);
                S1ARR[i*Inp+j].i = cimag(S1[j]);
                S2ARR[i*Inp+j].r = creal(S2[j]);
                S2ARR[i*Inp+j].i = cimag(S2[j]);
/*
                F11ARR[i*Inp+j] = F11[j];
                F33ARR[i*Inp+j] = F33[j];
                F12ARR[i*Inp+j] = F12[j];
                F34ARR[i*Inp+j] = F34[j];
*/
            }
        }

        free(S1);
        free(S2);
/*
        free(F11);
        free(F33);
        free(F12);
        free(F34);
*/
}
        /* Check to see if each of the output parameters has been
           specified in the IDL call, if it has copy the value to the
           return argument, otherwise simply deallocate it.
           NOTE: The keyword is included as an argument in defining
           argc buy NOT argv - crazy!*/
        if (argc > 12)
        {
/*
            IDL_DELTMP(ivF34Arr);
            IDL_DELTMP(ivF12Arr);
            IDL_DELTMP(ivF33Arr);
            IDL_DELTMP(ivF11Arr);
*/
            IDL_VarCopy(ivF34Arr, argv[11]);
            IDL_VarCopy(ivF12Arr, argv[10]);
            IDL_VarCopy(ivF33Arr, argv[9]);
            IDL_VarCopy(ivF11Arr, argv[8]);

            IDL_VarCopy(ivS2Arr, argv[7]);
            IDL_VarCopy(ivS1Arr, argv[6]);
            IDL_VarCopy(ivgArr, argv[5]);
            IDL_VarCopy(ivQbscArr, argv[4]);
            IDL_VarCopy(ivQscaArr, argv[3]);
            IDL_VarCopy(ivQextArr, argv[2]);
        }
        else
        {
            IDL_DELTMP(ivF34Arr);
            IDL_DELTMP(ivF12Arr);
            IDL_DELTMP(ivF33Arr);
            IDL_DELTMP(ivF11Arr);
            if (argc > 8)
            {
                IDL_VarCopy(ivS2Arr, argv[7]);
                IDL_VarCopy(ivS1Arr, argv[6]);
                IDL_VarCopy(ivgArr, argv[5]);
                IDL_VarCopy(ivQbscArr, argv[4]);
                IDL_VarCopy(ivQscaArr, argv[3]);
                IDL_VarCopy(ivQextArr, argv[2]);
            }
            else
            {
                IDL_DELTMP(ivS2Arr);
                if (argc > 7)
                {
                   IDL_VarCopy(ivS1Arr, argv[6]);
                   IDL_VarCopy(ivgArr, argv[5]);
                   IDL_VarCopy(ivQbscArr, argv[4]);
                   IDL_VarCopy(ivQscaArr, argv[3]);
                   IDL_VarCopy(ivQextArr, argv[2]);
                }
                else
                {
                    IDL_DELTMP(ivS1Arr);
                    if (argc > 6)
                    {
                       IDL_VarCopy(ivgArr, argv[5]);
                       IDL_VarCopy(ivQbscArr, argv[4]);
                       IDL_VarCopy(ivQscaArr, argv[3]);
                       IDL_VarCopy(ivQextArr, argv[2]);
                    }
                    else
                    {
                        IDL_DELTMP(ivgArr);
                        if (argc > 5)
                        {
                           IDL_VarCopy(ivQbscArr, argv[4]);
                           IDL_VarCopy(ivQscaArr, argv[3]);
                           IDL_VarCopy(ivQextArr, argv[2]);
                        }
                        else
                        {
                            IDL_DELTMP(ivQbscArr);
                            if (argc > 4)
                            {
                               IDL_VarCopy(ivQscaArr, argv[3]);
                               IDL_VarCopy(ivQextArr, argv[2]);
                            }
                            else
                            {
                                IDL_DELTMP(ivQscaArr);
                                if (argc > 3)
                                   IDL_VarCopy(ivQextArr, argv[2]);
                                else
                                   IDL_DELTMP(ivQextArr)
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        /* Dqv has not been defined, so set Inp = 1 and do calculations
           at 0 degrees */
        Inp = 1;
        dqvtmp[0] = 1;

        S1ARR  = (IDL_DCOMPLEX *)   IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                        Npts, IDL_ARR_INI_ZERO, &ivS1Arr);
        S2ARR  = (IDL_DCOMPLEX *)   IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                        Npts, IDL_ARR_INI_ZERO, &ivS2Arr);
        F11ARR = (double *)         IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                        Npts, IDL_ARR_INI_ZERO, &ivF11Arr);
        F33ARR = (double *)         IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                        Npts, IDL_ARR_INI_ZERO, &ivF33Arr);
        F12ARR = (double *)         IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                        Npts, IDL_ARR_INI_ZERO, &ivF12Arr);
        F34ARR = (double *)         IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                        Npts, IDL_ARR_INI_ZERO, &ivF34Arr);

#pragma omp parallel num_threads(nthreads) \
    private(mie_ok, Qext, Qsca, Qbsc, g, S1, S2)
/*
#pragma omp parallel num_threads(nthreads) \
    private(mie_ok, Qext, Qsca, Qbsc, g, S1, S2, F11, F33, F12, F34)
*/
{
        S1  = (double complex *) malloc(1 * sizeof(double complex));
        S2  = (double complex *) malloc(1 * sizeof(double complex));
/*
        F11 = (double *)         malloc(1 * sizeof(double));
        F33 = (double *)         malloc(1 * sizeof(double));
        F12 = (double *)         malloc(1 * sizeof(double));
        F34 = (double *)         malloc(1 * sizeof(double));
*/
#pragma omp for schedule(guided)
        for (i = 0; i < Npts; i++)
        {
            /* Call the mieint procedure */
/*
            mieint_(&DxARR[i], &Cm, &Inp, dqvtmp, &Qext, &Qsca, &Qbsc, &g,
                    S1, S2, F11, F33, F12, F34, &mie_ok);
*/
            mieint_(&DxARR[i], &Cm, &Inp, dqvtmp, &Qext, &Qsca, &Qbsc, &g,
                    S1, S2, F11ARR+i*Inp, F33ARR+i*Inp, F12ARR+i*Inp, F34ARR+i*Inp, &mie_ok);
            if (mie_ok != 0)
                IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                            "Size parameter overflow!");

             QextARR[i] = Qext;
             QscaARR[i] = Qsca;
             QbscARR[i] = Qbsc;
             gARR[i] = g;
             S1ARR[i].r = creal(S1[0]);
             S1ARR[i].i = cimag(S1[0]);
             S2ARR[i].r = creal(S2[0]);
             S2ARR[i].i = cimag(S2[0]);
/*
             F11ARR[i] = F11[0];
             F33ARR[i] = F33[0];
             F12ARR[i] = F12[0];
             F34ARR[i] = F34[0];
*/
        }

        free(S1);
        free(S2);
/*
        free(F11);
        free(F33);
        free(F12);
        free(F34);
*/
}
        /* Check to see if each of the output parameters has been
           specified in the IDL call, if it has copy the value to the
           return argument, otherwise simply deallocate it. */

        if (argc > 8)
        {
/*
            IDL_DELTMP(ivF34Arr);
            IDL_DELTMP(ivF12Arr);
            IDL_DELTMP(ivF33Arr);
            IDL_DELTMP(ivF11Arr);
*/
            IDL_VarCopy(ivF34Arr, argv[11]);
            IDL_VarCopy(ivF12Arr, argv[10]);
            IDL_VarCopy(ivF33Arr, argv[9]);
            IDL_VarCopy(ivF11Arr, argv[8]);

            IDL_VarCopy(ivS2Arr, argv[7]);
            IDL_VarCopy(ivS1Arr, argv[6]);
            IDL_VarCopy(ivgArr, argv[5]);
            IDL_VarCopy(ivQbscArr, argv[4]);
            IDL_VarCopy(ivQscaArr, argv[3]);
            IDL_VarCopy(ivQextArr, argv[2]);
        }
        else
        {
            IDL_DELTMP(ivF34Arr);
            IDL_DELTMP(ivF12Arr);
            IDL_DELTMP(ivF33Arr);
            IDL_DELTMP(ivF11Arr);
            if (argc > 7)
            {
               IDL_VarCopy(ivS2Arr, argv[7]);
               IDL_VarCopy(ivS1Arr, argv[6]);
               IDL_VarCopy(ivgArr, argv[5]);
               IDL_VarCopy(ivQbscArr, argv[4]);
               IDL_VarCopy(ivQscaArr, argv[3]);
               IDL_VarCopy(ivQextArr, argv[2]);
            }
            else
            {
                IDL_DELTMP(ivS2Arr);
                if (argc > 6)
                {
                   IDL_VarCopy(ivS1Arr, argv[6]);
                   IDL_VarCopy(ivgArr, argv[5]);
                   IDL_VarCopy(ivQbscArr, argv[4]);
                   IDL_VarCopy(ivQscaArr, argv[3]);
                   IDL_VarCopy(ivQextArr, argv[2]);
                }
                else
                {
                    IDL_DELTMP(ivS1Arr);
                    if (argc > 5)
                    {
                       IDL_VarCopy(ivgArr, argv[5]);
                       IDL_VarCopy(ivQbscArr, argv[4]);
                       IDL_VarCopy(ivQscaArr, argv[3]);
                       IDL_VarCopy(ivQextArr, argv[2]);
                    }
                    else
                    {
                        IDL_DELTMP(ivgArr);
                        if (argc > 4)
                        {
                           IDL_VarCopy(ivQbscArr, argv[4]);
                           IDL_VarCopy(ivQscaArr, argv[3]);
                           IDL_VarCopy(ivQextArr, argv[2]);
                        }
                        else
                        {
                            IDL_DELTMP(ivQbscArr);
                            if (argc > 3)
                            {
                               IDL_VarCopy(ivQscaArr, argv[3]);
                               IDL_VarCopy(ivQextArr, argv[2]);
                            }
                            else
                            {
                               IDL_DELTMP(ivQscaArr);
                               if (argc > 2)
                                  IDL_VarCopy(ivQextArr, argv[2]);
                               else
                                  IDL_DELTMP(ivQextArr)
                            }
                        }
                    }
                }
            }
        }
    }

    /* Cleanup any temporaries due to the keyword and deallocate the
       remain temporary arrays.*/
    IDL_KWCleanup(IDL_KW_CLEAN);

    return;
}
