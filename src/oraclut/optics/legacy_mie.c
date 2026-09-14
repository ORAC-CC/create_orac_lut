/* Thin non-IDL adapter around the preserved ORAC Mie Fortran kernel.
 *
 * The historical IDL DLM calls mieint_ one particle size at a time.  This
 * adapter exposes the same call through a small C ABI so Python can use the
 * scientific kernel without loading IDL runtime symbols.
 */

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

extern void mieint_(const double *dx, const double _Complex *cm,
                    const int32_t *inp, const double *dqv,
                    double *qext, double *qsca, double *qbsc, double *g,
                    double _Complex *xs1, double _Complex *xs2,
                    double *f11, double *f33, double *f12, double *f34,
                    int32_t *error);

int oraclut_mie_batch(int32_t nsize, const double *dx,
                      double real_index, double imaginary_index,
                      int32_t nang, const double *dqv,
                      double *qext, double *qsca, double *qbsc, double *asym,
                      double *f11, double *f33, double *f12, double *f34) {
    const double _Complex cm = real_index + imaginary_index * I;
    int32_t error = 0;

    for (int32_t i = 0; i < nsize; ++i) {
        double _Complex xs1[nang];
        double _Complex xs2[nang];
        error = 0;
        mieint_(&dx[i], &cm, &nang, dqv,
                &qext[i], &qsca[i], &qbsc[i], &asym[i],
                xs1, xs2,
                f11 + ((size_t)i * (size_t)nang),
                f33 + ((size_t)i * (size_t)nang),
                f12 + ((size_t)i * (size_t)nang),
                f34 + ((size_t)i * (size_t)nang), &error);
        if (error != 0) {
            return (int)error;
        }
    }
    return 0;
}
