/* DISORTtoIDL_DP: diagnostic double-precision DISORT DLM for IDL.
 *
 * The Fortran kernel is the unmodified production DISORTfunctions.f compiled
 * with default REAL promoted to 8 bytes (ifort -real-size 64).  This wrapper
 * replaces DISORTwrapper.c for that kernel only.  It exposes
 *
 *    DISORT_DP, <the same 44 positional arguments as DISORT>
 *    info = DISORT_DP_INFO()
 *
 * Differences from the historical single-precision DISORT procedure are
 * deliberate and limited to precision handling:
 *   - every REAL scalar is read with IDL_DoubleScalar, every REAL array of
 *     any numeric IDL type is converted element by element to double (float
 *     to double is exact), so IDL float32 inputs reach the kernel unchanged;
 *   - all eight outputs are returned as IDL DOUBLE arrays;
 *   - array sizes and the PMOM leading dimension are checked, because the
 *     kernel indexes PMOM(0:MAXMOM, MAXCLY) with MAXMOM = NMOM;
 *   - LAMBER must be true (the BRDF path calls a float-only C routine and is
 *     disabled in this build);
 *   - TEMPER is copied into a zero-padded buffer of NLYR + 1 elements.
 * As in DISORTwrapper.c, MAXCLY, MAXULV, MAXUMU, MAXPHI and MAXMOM are set to
 * NLYR, NTAU, NUMU, NPHI and NMOM so that arrays are passed at exactly the
 * sizes the kernel indexes, and a non-zero kernel status is a fatal error.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "idl_export.h"

#ifndef DP_BUILD_INFO
#error "DP_BUILD_INFO must be defined at compile time"
#endif

extern int32_t disort_(int32_t *nlyr, double *dtauc, double *ssalb, int32_t *nmom, double *pmom,
                       double *temper, double *wvnmlo, double *wvnmhi, int32_t *usrtau, int32_t *ntau,
                       double *utau, int32_t *nstr, int32_t *usrang, int32_t *numu, double *umu,
                       int32_t *nphi, double *phi, int32_t *ibcnd, double *fbeam, double *umu0,
                       double *phi0, double *fisot, int32_t *lamber, double *albedo, double *btemp,
                       double *ttemp, double *temis, int32_t *plank, int32_t *onlyfl, double *accur,
                       int32_t *prnt, int32_t *maxcly, int32_t *maxulv, int32_t *maxumu, int32_t *maxphi,
                       int32_t *maxmom, double *rfldir, double *rfldn, double *flup, double *dfdt,
                       double *uavg, double *uu, double *albmed, double *trnmed);
extern int32_t dprealbytes_(void);

static void fail(const char *message)
{
   IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, message);
}

static int32_t get_int(IDL_VPTR v)
{
   IDL_ENSURE_SCALAR(v);
   return (int32_t) IDL_LongScalar(v);
}

static double get_real(IDL_VPTR v)
{
   IDL_ENSURE_SCALAR(v);
   return IDL_DoubleScalar(v);
}

static int32_t get_logical(IDL_VPTR v)
{
   return get_int(v) != 0 ? 1 : 0;
}

/* Element-by-element conversion of the first n elements of a numeric array
 * (or scalar when n == 1) to a newly allocated double buffer of length
 * `size` (>= n), zero-padded. */
static double *get_real_array(IDL_VPTR v, IDL_MEMINT n, IDL_MEMINT size, const char *name)
{
   char message[256];
   IDL_MEMINT available, i;
   char *data;
   double *out;

   if (v->type == IDL_TYP_STRING || v->type == IDL_TYP_STRUCT || v->type == IDL_TYP_COMPLEX ||
       v->type == IDL_TYP_DCOMPLEX || v->type == IDL_TYP_PTR || v->type == IDL_TYP_OBJREF ||
       v->type == IDL_TYP_UNDEF) {
      snprintf(message, sizeof message, "DISORT_DP: %s must be a real numeric array", name);
      fail(message);
   }
   if (v->flags & IDL_V_ARR) {
      available = v->value.arr->n_elts;
      data = (char *) v->value.arr->data;
   } else {
      available = 1;
      data = (char *) &v->value;
   }
   if (available < n) {
      snprintf(message, sizeof message, "DISORT_DP: %s has %lld elements, %lld required", name,
               (long long) available, (long long) n);
      fail(message);
   }
   out = (double *) calloc((size_t) (size > 0 ? size : 1), sizeof(double));
   if (out == NULL) fail("DISORT_DP: out of memory");
   for (i = 0; i < n; i++) {
      switch (v->type) {
         case IDL_TYP_BYTE:    out[i] = ((UCHAR *) data)[i]; break;
         case IDL_TYP_INT:     out[i] = ((IDL_INT *) data)[i]; break;
         case IDL_TYP_LONG:    out[i] = ((IDL_LONG *) data)[i]; break;
         case IDL_TYP_FLOAT:   out[i] = ((float *) data)[i]; break;
         case IDL_TYP_DOUBLE:  out[i] = ((double *) data)[i]; break;
         case IDL_TYP_UINT:    out[i] = ((IDL_UINT *) data)[i]; break;
         case IDL_TYP_ULONG:   out[i] = ((IDL_ULONG *) data)[i]; break;
         case IDL_TYP_LONG64:  out[i] = (double) ((IDL_LONG64 *) data)[i]; break;
         case IDL_TYP_ULONG64: out[i] = (double) ((IDL_ULONG64 *) data)[i]; break;
         default:
            free(out);
            snprintf(message, sizeof message, "DISORT_DP: unsupported type for %s", name);
            fail(message);
      }
   }
   return out;
}

static double *output_array(IDL_VPTR target, int n_dim, IDL_MEMINT *dims, IDL_VPTR *temp)
{
   IDL_EXCLUDE_EXPR(target);
   return (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, n_dim, dims, IDL_ARR_INI_ZERO, temp);
}

static void IDL_CDECL disort_dp(int argc, IDL_VPTR argv[], char *argk)
{
   int32_t nlyr, nmom, usrtau, ntau, nstr, usrang, numu, nphi, ibcnd, lamber, plank, onlyfl;
   int32_t maxcly, maxulv, maxumu, maxphi, maxmom, status, prnt[5];
   double wvnmlo, wvnmhi, fbeam, umu0, phi0, fisot, albedo, btemp, ttemp, temis, accur;
   double *dtauc, *ssalb, *pmom, *temper, *utau, *umu, *phi, *prnt_values;
   double *rfldir, *rfldn, *flup, *dfdt, *uavg, *uu, *albmed, *trnmed;
   IDL_VPTR out[8];
   IDL_MEMINT d1[1], d3[3], i, temper_available;
   char message[256];

   (void) argc;
   (void) argk;
   nlyr = get_int(argv[0]);
   nmom = get_int(argv[3]);
   wvnmlo = get_real(argv[6]);
   wvnmhi = get_real(argv[7]);
   usrtau = get_logical(argv[8]);
   ntau = get_int(argv[9]);
   nstr = get_int(argv[11]);
   usrang = get_logical(argv[12]);
   numu = get_int(argv[13]);
   nphi = get_int(argv[15]);
   ibcnd = get_int(argv[17]);
   fbeam = get_real(argv[18]);
   umu0 = get_real(argv[19]);
   phi0 = get_real(argv[20]);
   fisot = get_real(argv[21]);
   lamber = get_logical(argv[22]);
   albedo = get_real(argv[23]);
   btemp = get_real(argv[24]);
   ttemp = get_real(argv[25]);
   temis = get_real(argv[26]);
   plank = get_logical(argv[27]);
   onlyfl = get_logical(argv[28]);
   accur = get_real(argv[29]);

   if (nlyr < 1 || nmom < 0 || ntau < 1 || numu < 1 || nphi < 1 || nstr < 2)
      fail("DISORT_DP: invalid dimension argument");
   if (!lamber)
      fail("DISORT_DP: only LAMBER = 1 is supported (BRDF path disabled in this build)");
   if (!usrtau || !usrang)
      fail("DISORT_DP: USRTAU and USRANG must be set, as in all ORAC calls");

   /* PMOM must be laid out as (NMOM+1, >= NLYR) for PMOM(0:MAXMOM, MAXCLY). */
   if ((argv[4]->flags & IDL_V_ARR) && argv[4]->value.arr->n_dim >= 2 &&
       argv[4]->value.arr->dim[0] != (IDL_MEMINT) nmom + 1) {
      snprintf(message, sizeof message, "DISORT_DP: PMOM leading dimension %lld, NMOM+1 = %d required",
               (long long) argv[4]->value.arr->dim[0], nmom + 1);
      fail(message);
   }

   dtauc = get_real_array(argv[1], nlyr, nlyr, "DTAUC");
   ssalb = get_real_array(argv[2], nlyr, nlyr, "SSALB");
   for (i = 0; i < nlyr; i++)
      if (ssalb[i] > 1.0 || ssalb[i] < 0.0) fail("DISORT_DP: SSALB outside [0, 1]");
   pmom = get_real_array(argv[4], (IDL_MEMINT) (nmom + 1) * nlyr, (IDL_MEMINT) (nmom + 1) * nlyr, "PMOM");
   temper_available = (argv[5]->flags & IDL_V_ARR) ? argv[5]->value.arr->n_elts : 1;
   if (plank && temper_available < nlyr + 1)
      fail("DISORT_DP: TEMPER needs NLYR + 1 elements when PLANK is set");
   temper = get_real_array(argv[5], temper_available < nlyr + 1 ? temper_available : nlyr + 1, nlyr + 1,
                           "TEMPER");
   utau = get_real_array(argv[10], ntau, ntau, "UTAU");
   umu = get_real_array(argv[14], numu, numu, "UMU");
   phi = get_real_array(argv[16], nphi, nphi, "PHI");
   prnt_values = get_real_array(argv[30], 5, 5, "PRNT");
   for (i = 0; i < 5; i++) prnt[i] = prnt_values[i] != 0.0 ? 1 : 0;
   free(prnt_values);

   maxcly = nlyr;
   maxulv = ntau;
   maxumu = numu;
   maxphi = nphi;
   maxmom = nmom;

   d1[0] = ntau;
   rfldir = output_array(argv[36], 1, d1, &out[0]);
   rfldn = output_array(argv[37], 1, d1, &out[1]);
   flup = output_array(argv[38], 1, d1, &out[2]);
   dfdt = output_array(argv[39], 1, d1, &out[3]);
   uavg = output_array(argv[40], 1, d1, &out[4]);
   d3[0] = numu;
   d3[1] = ntau;
   d3[2] = nphi;
   uu = output_array(argv[41], 3, d3, &out[5]);
   d1[0] = numu;
   albmed = output_array(argv[42], 1, d1, &out[6]);
   trnmed = output_array(argv[43], 1, d1, &out[7]);

   status = disort_(&nlyr, dtauc, ssalb, &nmom, pmom, temper, &wvnmlo, &wvnmhi, &usrtau, &ntau, utau,
                    &nstr, &usrang, &numu, umu, &nphi, phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot,
                    &lamber, &albedo, &btemp, &ttemp, &temis, &plank, &onlyfl, &accur, prnt,
                    &maxcly, &maxulv, &maxumu, &maxphi, &maxmom, rfldir, rfldn, flup, dfdt, uavg, uu,
                    albmed, trnmed);

   free(dtauc);
   free(ssalb);
   free(pmom);
   free(temper);
   free(utau);
   free(umu);
   free(phi);
   if (status != 0) {
      for (i = 0; i < 8; i++) IDL_Deltmp(out[i]);
      fail("DISORT_DP: fatal error reported by the double-precision DISORT kernel");
   }
   for (i = 0; i < 8; i++) IDL_VarCopy(out[i], argv[36 + i]);
}

static IDL_VPTR IDL_CDECL disort_dp_info(int argc, IDL_VPTR argv[], char *argk)
{
   char text[1024];
   (void) argc;
   (void) argv;
   (void) argk;
   snprintf(text, sizeof text, "DISORTtoIDL_DP fortran_real_bytes=%d c_double_bytes=%d %s",
            (int) dprealbytes_(), (int) sizeof(double), DP_BUILD_INFO);
   return IDL_StrToSTRING(text);
}

int IDL_Load(void)
{
   static IDL_SYSFUN_DEF2 procedures[] = {
      {{(IDL_SYSRTN_GENERIC) disort_dp}, "DISORT_DP", 44, 44, 0, 0},
   };
   static IDL_SYSFUN_DEF2 functions[] = {
      {{(IDL_SYSRTN_GENERIC) disort_dp_info}, "DISORT_DP_INFO", 0, 0, 0, 0},
   };
   if (dprealbytes_() != 8) return IDL_FALSE;
   return IDL_SysRtnAdd(procedures, IDL_FALSE, IDL_CARRAY_ELTS(procedures)) &&
          IDL_SysRtnAdd(functions, IDL_TRUE, IDL_CARRAY_ELTS(functions));
}
