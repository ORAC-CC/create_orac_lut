#ifndef __AMBRALSFOR_H__
#define __AMBRALSFOR_H__

#include <math.h>

#define D2R(p) M_PI*(p)/180.    /* degrees->radians */
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define BR 1.0                  /* LiSparse b/r */
#define HB 2.0                  /* LiSparse h/b */

typedef struct {
  float iso;
  float vol;
  float geo;
} param_t;

typedef struct {
  float szen;
  float vzen;
  float relaz;
} geom_t;

float forward(geom_t geom,  param_t params);

#endif /* __AMBRALSFOR_H__ */























