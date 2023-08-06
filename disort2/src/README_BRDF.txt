This distribution of the DIScrete Ordinates Radiative Transfer code (DISORT)
includes a new version of the DISORT BDREF subroutine that allows the code to
use the MODIS AMBRALS BRDF surface reflectance model to define the bottom
boundary of the scattering medium. No changes have been made to the DISORT code
itself, which was written to allow a BRDF module to be swapped in for the
BDREF.f stub subroutine included with the original code (contained in the
printine disort/ directory).

The AMBRALS model itself is written in C and can be found in the ambralsfor.h
and .c files. These are essentially unchanged from the code which can be
obtained on request from the NASA MODIS surface properties team. This is called
by the ambrals-fortran.c code, which provides the interface to the FORTRAN code
and loads in the model parameters needed to run the AMBRALS model. Model
parameters are provided to the code via a plain text file, the (fully
qualified) path of which should be specified with the "AMBRALS_BRDF_FILE"
environment variable (which should be specified by the user prior to calling
the code). The parameter file must have the following format:

n
wvl fiso fvol fgeo
.
.
wvl fiso fvol fgeo

where "n" is the number of wavelengths supplied by the parameters file. Each of
the "n" subsequent lines gives the wavelength in microns, followed by the
isotropic, volumetric and geometric AMBRALS model parameters as decimal
numbers. Note: The parameters file should not contain more than 1000 sets of
parameters. If it does, only the first 1000 will be read. If no parameter file
is specified, the code cannot find the file, or there is a problem reading the
file, the code will issue a warning message and use default parameters in the
AMBRALS model.

By default 'AmbralsFortran' will only read the parameters file once per IDL
session (the AMBRALS parameters are then stored in static variables, where the
code can access them at will). To reload the parameters file (if, for instance,
you'd like to perform a forward model calculation over a different surface),
you can set the 'AMBRALS_BRDF_FILE_REFRESH' environment variable to "True" (or
"true", any string beginning with a t/T, or the number 1). This will cause
AmbralsFortran to read the parameters file specified by 'AMBRALS_BRDF_FILE'
again when DISORT is next called. Note: 'AmbralsFortran' then sets
'AMBRALS_BRDF_FILE_REFRESH' to "false", so that the file isn't then re-read
over and over.

Enjoy.

Gareth Thomas 05/01/2010
