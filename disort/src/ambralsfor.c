#include <math.h>
#include "ambralsfor.h"

/* All the following functions are taken from ambmodels.c.reciprocal.
 * The functions have been changed to ANSI C from K&R and trig functions
 * have been replaced by their float equivs
 */

void LiKernel(float, float, float, float, float, float, float *);
void GetPhaang(float, float, float, float, float, float *, float *,float *);
void GetpAngles(float, float, float *, float *, float *);
void GetDistance(float , float , float ,float *);
void GetOverlap(float, float, float, float, float, float, float, float *,
                float *);

float forward(geom_t geom, param_t params)
{
  /* this function runs the RossThickLiSparseReciprocal model in the
   * forward mode and returns the calculated reflectance.
   *
   * Parameters:
   * geom_t geom	view geometry
   * param_t params	BRDF kernel weights parameters
   */

  float phi, cosphi, sinphi, tv, costv, sintv, ti, costi, sinti;
  float cosphaang, phaang, sinphaang, rosskernel, tantv, tanti;
  float likernel, refl;

  ti = geom.szen;
  tv = geom.vzen;
  phi = geom.relaz;

  ti = D2R(ti);
  tv = D2R(tv);
  phi = D2R(phi);

  cosphi = cosf(phi);
  sinphi = sinf(phi);
  costv = cos(tv);
  costi = cos(ti);
  sintv = sinf(tv);
  sinti = sinf(ti);

  GetPhaang(costv, costi, sintv, sinti, cosphi, &cosphaang,
                  &phaang, &sinphaang);

  rosskernel = ((M_PI_2 - phaang) * cosphaang + sinphaang)/(costi + costv)
          -   M_PI_4;

  tantv = sintv/costv;
  tanti = sinti/costi;

  LiKernel(HB,BR,tantv,tanti,sinphi,cosphi, &likernel);

  refl = params.iso + params.vol * rosskernel + params.geo * likernel;

  return refl;
}

void LiKernel(float hbratio, float brratio, float tantv, float tanti,
              float sinphi, float cosphi, float *result)
{

  /* This func calculates the LiSparseModis kernel */

  float sintvp, costvp, tantvp, sintip, costip, tantip;
  float phaangp, cosphaangp, sinphaangp, distancep, overlap, temp;

  GetpAngles(brratio,tantv,&sintvp,&costvp,&tantvp);
  GetpAngles(brratio,tanti,&sintip,&costip,&tantip);
  GetPhaang(costvp,costip,sintvp,sintip,cosphi,&cosphaangp,&phaangp,
             &sinphaangp);
  GetDistance(tantvp,tantip,cosphi,&distancep);
  GetOverlap(hbratio,distancep,costvp,costip,tantvp,tantip,sinphi,
              &overlap,&temp);

  *result = overlap - temp + 0.5 * (1.+cosphaangp)/costvp/costip;
}

void GetPhaang(float cos1, float cos2, float sin1, float sin2, float cos3,
               float *cosres, float *res,float *sinres)
{
  /* This func calculates the phase angle */

  *cosres = cos1*cos2 + sin1*sin2*cos3;
  *res = acosf( MAX(-1., MIN(1.,*cosres)) );
  *sinres = sinf(*res);
}

void GetpAngles(float brratio, float tan1, float *sinp, float *cosp,
                float *tanp)
{
  /* This func calculates the 'prime' angles */

  float angp;

    *tanp = brratio*tan1;
  if(*tanp < 0) *tanp = 0.;
  angp = atanf(*tanp);
  *sinp = sinf(angp);
  *cosp = cosf(angp);
}

void GetDistance(float tan1, float tan2, float cos3,float *res)
{
  /* This func calculates the D distance */

  float temp;

  temp = tan1*tan1+tan2*tan2-2.*tan1*tan2*cos3;
  *res = sqrtf(MAX(0.,temp));
}

void GetOverlap(float hbratio, float distance, float cos1, float cos2,
                float tan1, float tan2, float sin3, float *overlap,
                float *temp)
{
  /* This func calculates the Overlap distance */

  float cost, sint, tvar;

  *temp = 1./cos1 + 1./cos2;
   cost = hbratio*sqrtf(distance*distance+tan1*tan1*tan2*tan2*sin3*sin3)/
          *temp;
   cost = MAX(-1., MIN(1.,cost));
   tvar = acosf(cost);
   sint = sinf(tvar);
   *overlap = 1./M_PI *(tvar-sint*cost)*(*temp);
   *overlap = MAX(0.,*overlap);
}

