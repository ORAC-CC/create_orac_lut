#include "ambralsfor.h"
#include <stdio.h>
#include <stdlib.h>

void ambralsfortran_(float * pwvlo, float * pwvhi,
		     float * psz, float * pvz, float * pdphi,
		     float * pbdr)
{
  /* This function is a wrapper for the MODIS BRDF forward model for use
     with DISORT (it can be called from the DISORT function BDREF),
     which will use a sample of MODIS BRDF data, call the Ambrals
     Forward model and return a reflectance.

     G Thomas 02/02/2007 Original version
     G Thomas 05/01/2009 Code now looks for AMBRALS model parameters
                         in a file specified by the AMBRALS_BRDF_FILE
                         environment variable, and checks to see if
                         parameters should be reloaded using the
                         AMBRALS_BRDF_FILE_REFRESH environment variable
  */

  param_t params;
  geom_t  geom;
  char * FPath, * refresh;
  FILE * FFile;
  char cnt, rtest, dflt = 0;
  long i, wvl_i = 0;
  float wvldiff = 1e6;
  static long FNwvl = 0;
  static float Fwvl[1000], Fiso[1000], Fvol[1000], Fgeo[1000];

  /* Check to see if the AMBRALS_BRDF_FILE_REFRESH has been set to
     either "1" or has "t" or "T" as it's first character. If so, we
     reload the AMBRALS_BRDF_FILE...
     We then reset the AMBRALS_BRDF_FILE_REFRESH variable to false so
     that we don't keep re-reading the parameters (until re-instructed
     to do so by whichever routine set the environment variable in the
     first place).
   */
  refresh = getenv("AMBRALS_BRDF_FILE_REFRESH");
  if (refresh != NULL)
    {
      rtest = refresh[0];
      if (rtest == 'T' || rtest == 't' || (int) rtest == 1)
	{
	  FNwvl = 0;
	  setenv("AMBRALS_BRDF_FILE_REFRESH","false",1);
	}
    }
  if (FNwvl == 0)
    {

      /* Get the path to the file which contains the Ambrals BRDF model
	 parameters. This is specified by the AMBRALS_BRDF_FILE and is
	 expected to be a text file which contains the number of entries
	 on the first line, followed by space separated values:
	 nominal_wavelength (in microns), iso, vol, geo
      */
      printf ("DISORT BRDF surface: Atempting to load AMBRALS parameters...\n");
      FPath = getenv ("AMBRALS_BRDF_FILE");
      if (FPath == NULL)
	{
	  printf ("Warning: AMBRALS parameter file not found, using default parameters\n");
	  dflt = 1;
	  goto dfltpt;
	}
      /* Open the file and read it! */
      FFile = fopen(FPath,"r");
      if (FFile == NULL)
	{
	  printf ("Warning: Couldn't open AMBRALS parameter file, using default parameters\n");
	  dflt = 1;
	  goto dfltpt;
	}

      cnt = fscanf(FFile, "%ld", &FNwvl);
      if (cnt == 0)
	{
	  printf ("Warning: Error reading AMBRALS parameter file, using default parameters\n");
	  dflt = 1;
	  goto dfltpt;
	}
      if (FNwvl > 1000)
	{
	  printf ("Warning: AMBRALS parameter file contains more than 1000 entries\n         Only first 1000 will be read\n");
	  FNwvl = 1000;
	}
      for (i = 0; i < FNwvl; i++)
	{
	  cnt = fscanf(FFile, "%g %g %g %g", &Fwvl[i], &Fiso[i],
		       &Fvol[i], &Fgeo[i]);
	  if (cnt == 0)
	    {
	      printf ("Warning: Error reading AMBRALS parameter file, using default parameters\n");
	      dflt = 1;
	      goto dfltpt;
	    }
	  /* Figure out which AMBRALS wavelength is most appropriate for
	     our purposes */
	  if (wvldiff < ((1e4 / *pwvlo - Fwvl[i])*(1e4 / *pwvlo - Fwvl[i]) +
			 (1e4 / *pwvhi - Fwvl[i])*(1e4 / *pwvhi - Fwvl[i])))
	    wvl_i = i;
	}
      fclose(FFile);
    }
  /* Build our params structure using the values read from the file that
     best match our wavelength range */
  params.iso = Fiso[wvl_i];
  params.vol = Fvol[wvl_i];
  params.geo = Fgeo[wvl_i];

  /* printf ("Parameters are: %g, %g, %g\n",
             params.iso,params.vol,params.geo); */

  /* If we encountered a problem reading the file, substitute some
     default values */
 dfltpt:
  if (dflt == 1)
    {
      /* Dark forest (i). Amazon Basin. */
      params.iso = 0.040;
      params.vol = 0.036;
      params.geo = 0.0020;
    }
  /* Here we define a series our model parameters (which are usually
     extracted from the MODIS BRDF product) for different surface types
     (uncomment the one you would like to use:
  */

  /* Dark forest (i). Amazon Basin. */
  /*
  params.iso = 0.040;
  params.vol = 0.036;
  params.geo = 0.0020;
  */

  /* Dark forest (ii). Congo */
  /*
  params.iso = 0.040;
  params.vol = 0.037;
  params.geo = 0.0020;
  */

  /* Desert. Sahara */
  /*
  params.iso = 0.260;
  params.vol = 0.080;
  params.geo = 0.0010;
  */

  /* Mid-latitude cultivated summer (i). Central Europe */
  /*
  params.iso = 0.064;
  params.vol = 0.043;
  params.geo = 0.0030;
  */

  /* Mid-latitude cultivated summer (ii). S.W. USA */
  /*
  params.iso = 0.071;
  params.vol = 0.060;
  params.geo = 0.0030;
  */

  /* Set up the geometry variable */
  geom.szen  = *psz;
  geom.vzen  = *pvz;
  geom.relaz = *pdphi;

  /* Call the Ambrals forward model */
  *pbdr = forward(geom, params);

  /* That's it, return control to FORTRAN */
  return;
}
