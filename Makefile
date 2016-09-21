NAME    = create_orac_lut
VERSION = 0.01


include ${HOME}/macros.inc.new


INCDIRS += ${GRMLIB_INCDIR}
LIBDIRS += ${GRMLIB_LIBDIR}
LINKS   += ${GRMLIB_LINK} ${MATH_LINK}


SUBDIRS   =

OBJECTS   =

LIBRARIES =
BINARIES  =
PRODUCTS  =

TARGETS   = ${LIBRARIES} ${BINARIES} ${PRODUCTS}

EXTRA_CLEANS =


include ${INSTALL}/targets.inc


all:
	cd disort/src && $(MAKE) ifort

clean:
	cd disort/src && $(MAKE) clean


include dep.inc
