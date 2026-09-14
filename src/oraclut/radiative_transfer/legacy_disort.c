/* Non-IDL ABI adapter for the preserved production ORAC DISORT sources. */

#include <stdint.h>

/* An old-style declaration is intentional: the preserved Fortran entry point
 * has 44 pass-by-reference arguments, and the C adapter below keeps that ABI
 * explicit at the call site. */
extern int32_t disort_();
extern int32_t getmom_(int32_t *, float *, int32_t *, float *);
extern float plkavg_(float *, float *, float *);

int oraclut_disort(int32_t nlyr, float *dtauc, float *ssalb, int32_t nmom,
                   float *pmom, float *temper, float wvnlo, float wvnhi,
                   int32_t plank, int32_t ntau, float *utau, int32_t nstr,
                   int32_t numu, float *umu, int32_t nphi, float *phi,
                   float fbeam, float umu0, float fisot,
                   float *rfldir, float *rfldn, float *flup, float *dfdt,
                   float *uavg, float *uu, float *albmed, float *trnmed) {
    int32_t usrtau = 1;
    int32_t usrang = 1;
    int32_t ibcnd = 0;
    int32_t lamber = 1;
    float albedo = 0.0f;
    float btemp = 0.0f;
    float ttemp = 0.0f;
    float temis = 0.0f;
    int32_t onlyfl = 0;
    float accur = 1.0e-8f;
    int32_t prnt[5] = {0, 0, 0, 0, 0};
    int32_t maxcly = nlyr;
    int32_t maxulv = ntau;
    int32_t maxumu = numu;
    int32_t maxphi = nphi;
    int32_t maxmom = nmom;
    float phi0 = 0.0f;

    return (int)disort_(&nlyr, dtauc, ssalb, &nmom, pmom, temper,
                        &wvnlo, &wvnhi, &usrtau, &ntau, utau, &nstr,
                        &usrang, &numu, umu, &nphi, phi, &ibcnd, &fbeam,
                        &umu0, &phi0, &fisot, &lamber, &albedo, &btemp,
                        &ttemp, &temis, &plank, &onlyfl, &accur, prnt,
                        &maxcly, &maxulv, &maxumu, &maxphi, &maxmom,
                        rfldir, rfldn, flup, dfdt, uavg, uu, albmed,
                        trnmed);
}

int oraclut_getmom(int32_t iphas, float gg, int32_t nmom, float *pmom) {
    return (int)getmom_(&iphas, &gg, &nmom, pmom);
}

float oraclut_plkavg(float wvnlo, float wvnhi, float temperature) {
    return plkavg_(&wvnlo, &wvnhi, &temperature);
}
