# DISORToIDL DLM Makefile (August 2001)
#
# HISTORY:
# 05/02/13, G Thomas: Steven Marsh.
# 04/12/06, G Thomas: Modifications for compiling with the Intel Fortran
#                     compiler v. 9 and IDL 6.2 (needed for 64 bit
#                     systems). At the moment 32 bit systems still use
#                     v. 7. To compile on a 64 bit system the ifort
#                     conditional should be set ALSO, gcc conditional
#                     now tries a fully gcc compile (i.e. using
#                     gfortran). This doesn't actually work at present!
# 02/01/07, G Thomas: Added in a condition for compiling the BRDF
#                     version of DISORT.
# 05/01/10, G Thomas: Tidied up and made code provide a statically
#                     linked library by default when called with ifort.
# 09/08/16, G McGarragh: A fully gcc compile (i.e. using gfortran) now
#                        works.


# Location specific configurations
include make.inc

# If the version number of IDL is prior to 5.2 then set the
# IDL flag below to "-D__IDLPRE53__".
#IDL = -D__IDLPRE53__


# VARIABLES ################################################

INC = -I$(IDL_DIR)/external

OBJS = DISORTfunctions.o DISORTtoIDL.o DISORTwrapper.o \
       BDREF.o PRTFIN.o ambralsfor.o ambrals-fortran.o \
       disort/ErrPack.o disort/LINPAK.o disort/RDI1MACH.o \

SO_EXT = so


# CONDITIONALS #############################################

dec: CC = cc
dec: FC = f77
dec: FLAGS = -float
dec: OFLAGS =
dec: CFLAGS = $(FLAGS) $(OFLAGS) $(INC) $(IDL)
dec: FFLAGS =
dec: LDFLAGS = -expect_unresolved '*' -shared -all

gcc: CC = gcc
gcc: FC = gfortran
gcc: CPL = -D__GFORTRAN__
gcc: FLAGS = -fPIC -Wall # -pedantic
gcc: OFLAGS = -O0
gcc: CFLAGS = $(FLAGS) $(OFLAGS) $(INC) $(IDL) $(CPL)
gcc: FFLAGS = $(FLAGS) $(OFLAGS) $(CPL)
gcc: LDFLAGS = -shared -lm -lgfortran

fujitsu: CC = gcc
fujitsu: FC = frt
fujitsu: CPL = -D__FUJITSU__
fujitsu: LIBS = -L/usr/local/PACK/FFC/lib
fujitsu: FLAGS = -Wall -Wimplicit -fpic -pedantic
fujitsu: OFLAGS = -O3
fujitsu: CFLAGS = $(FLAGS) $(OFLAGS) $(INC) $(IDL) $(CPL)
fujitsu: FFLAGS = -Fixed -K PIC $(OFLAGS) $(CPL)
fujitsu: LDFLAGS = -shared -lfj9f6 -lfj9i6

ifort: CC = gcc
ifort: FC = ifort
ifort: CPL = -D__GFORTRAN__
ifort: FLAGS = -O3 -fpic
ifort: OFLAGS = -O3
ifort: CFLAGS = -pedantic $(FLAGS) $(OFLAGS) $(INC) $(IDL) $(CPL)
ifort: FFLAGS = -fpic $(FLAGS) $(OFLAGS) $(CPL)
# Dynamic linking with Fortran libraries
#ifort: LDFLAGS = -shared -lm -lifcore
# Static linking with Fortran libraries
ifort: LDFLAGS = $(IFORT_DIR)/libifcore_pic.a $(IFORT_DIR)/libimf.a \
                 $(IFORT_DIR)/libirc.a $(IFORT_DIR)/libsvml.a -shared -lm

nag: CC = gcc
nag: FC = f95
nag: CPL = -D__NAG__
nag: LIBS = -L/usr/local/lib/NAGWare
nag: FLAGS = -Wall -Wimplicit -fpic -pedantic
nag: OFLAGS = -O3
nag: CFLAGS = $(FLAGS) $(OFLAGS) $(INC) $(IDL) $(CPL)
nag: FFLAGS = -C -fixed $(OFLAGS) $(CPL)
nag: LDFLAGS = -shared -lf96

sun: CC = cc
sun: FC = f77
sun: FLAGS = -K pic
sun: OFLAGS =
sun: CFLAGS = $(FLAGS) $(OFLAGS) $(INC) $(IDL)
sun: FFLAGS =
sun: LDFLAGS = -B symbolic -G


# TARGETS ##################################################

default:
	@echo
	@echo " Specify one of the following:"
	@echo "   dec      DEC/Compaq C Suite"
	@echo "   gcc      GNU Compiler Collection"
	@echo "   fujitsu  GNU GCC + Fujitsu Fortran 95 Compiler"
	@echo "   ifort    GNU GCC + Intel Fortran 95 Compiler"
	@echo "   nag      GNU GCC + NAGWare Fortran 95 Compiler"
	@echo "   sun      Sun Workshop Compiler Suite"
	@echo "   clean    Delete object/core/editor files"
	@echo
	@echo " (i.e. if you're doing this on a AOPP machine, you probably want make ifort')"
	@echo

dec:		pre debug DISORTtoIDL.$(SO_EXT) post
gcc:		pre debug DISORTtoIDL.$(SO_EXT) post
fujitsu:	pre debug DISORTtoIDL.$(SO_EXT) post
ifort:		pre debug DISORTtoIDL.$(SO_EXT) post
nag:		pre debug DISORTtoIDL.$(SO_EXT) post
sun:		pre debug DISORTtoIDL.$(SO_EXT) post

DISORTtoIDL.$(SO_EXT): $(OBJS)
	@echo '*** If you get an error saying export.h cannot be found, ***'
	@echo '    *** check the IDL directory set in the Makefile ***'
	$(CC) -o DISORTtoIDL.$(SO_EXT) $(OBJS) $(LDFLAGS) $(LIB) $(LIBS)

.PHONY: clean debug pre post

clean:
	rm -f *.o disort/*.o DISORTtoIDL.$(SO_EXT) core

debug:
	@echo 'Compiling with:'
	@echo CC = $(CC)
	@echo FC = $(FC)
	@echo CFLAGS = $(CFLAGS)
	@echo FFLAGS = $(FFLAGS)
	@echo OBJS = $(OBJS)
	@echo LDFLAGS = $(LDFLAGS)
	@echo

pre:
	@echo; date; echo "Compiling for: " `uname -srm`; echo

post:
	@echo "All targets built."; echo
