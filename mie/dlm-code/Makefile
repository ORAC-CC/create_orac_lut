# mie_dlm_single IDL DLM makefile

# Location specific configurations
include ../make.inc

# If the version number of IDL is prior to 5.2 then set the IDL flag below to
# "-D__IDLPRE53__".
#IDL = -D__IDLPRE53__

# Uncomment to compile a multithreaded version
GCC_OPENMP = -fopenmp

# C99 is required for long long, complex and variable length arrrays
GCC_FLAGS = -O3 -fPIC -pedantic -std=c99 -Wall -Wno-unknown-pragmas $(GCC_OPENMP)


ifeq ($(FCOMP),gfortran)
        # gfortran flags
	FC_FLAGS = -O3 -fPIC -ffixed-line-length-0
        #Linking flags
	LD_FLAGS = -shared -lm $(GCC_OPENMP)
else
        # ifort flags
	FC_FLAGS = -O3 -extend_source -fpic
        # Linking flags
	#LD_FLAGS = -shared -lifcore -limf -lirc -L$(IFORT_DIR)/lib/intel64 $(GCC_OPENMP)
        # Linking flags
	LD_FLAGS = -shared -lifcore -limf -lirc  $(GCC_OPENMP)
endif


OBJS = mie_dlm_single.o mieint.o mieintnocmplx.o

INC = -I$(IDL_DIR)/external

OUT_FILE = mie_dlm_single.so


$(OUT_FILE): $(OBJS)
	$(FCOMP) $(LD_FLAGS) -o $(OUT_FILE) $(OBJS)

mieint.o: mieint.f
	$(FCOMP) $(FC_FLAGS) -c -o mieint.o mieint.f
mieintnocmplx.o: mieintnocmplx.f
	$(FCOMP) $(FC_FLAGS) -c -o mieintnocmplx.o mieintnocmplx.f
mie_dlm_single.o: mie_dlm_single.c
	gcc $(GCC_FLAGS) $(INC) -c -o mie_dlm_single.o mie_dlm_single.c

clean:
	rm -f *.o *.so
