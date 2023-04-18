EXEC = int2
#HOST = gcc
HOST = gcc-openmp
#HOST = intel
#HOST = intel-ompenmp

#
# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# and libsolvers3d.so, for Macosx, these .so files also need to be
# copied over /usr/local/lib


FMMBIE_INSTALL_DIR = $(PREFIX)
FMM_INSTALL_DIR = $(PREFIX_FMM)
LFMMLINKLIB = -lfmm3d
LLINKLIB = -lfmm3dbie

ifneq ($(OS),Windows_NT) 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        ifeq ($(PREFIX_FMM),)
            FMM_INSTALL_DIR=/usr/local/lib
        endif
        ifeq ($(PREFIX),)
            FMMBIE_INSTALL_DIR=/usr/local/lib
        endif
    endif
    ifeq ($(UNAME_S),Linux)
        ifeq ($(PREFIX_FMM),)
            FMM_INSTALL_DIR=${HOME}/lib
        endif
        ifeq ($(PREFIX),)
            FMMBIE_INSTALL_DIR=${HOME}/lib
        endif
    endif
endif


ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
endif

ifeq ($(HOST),intel)
    FC=ifort 
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort 
    FFLAGS= -O3 -fPIC -march=native -qopenmp
endif


FC = gfortran
FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy 
CC = gcc


FEND = -L${FMMBIE_INSTALL_DIR} $(LLINKLIB) -L${FMM_INSTALL_DIR} $(LFMMLINKLIB)


# For your OS, override the above by placing make variables in make.inc
-include make.inc




.PHONY: usage clean list

OBJECTS =  example_em_efie_rwg_pec.o



default: usage

usage:
	@echo ""
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for Maxwell examplee, please specify what to make:"
	@echo "  make efie_rwg_open ----------- scattering from an open surface"
	@echo "  make clean ------------------- remove all .o files and executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "-------------------------------------------------------------------------"
	@echo ""



#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f  
	$(FC) -c $(FFLAGS) $< -o $@ $(FEND)

%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@ $(FEND)






efie_rwg_open: efie_rwg_open.o
	$(FC) $(FFLAGS) -o efie_rwg_open efie_rwg_open.o $(FEND)
	./efie_rwg_open

clean:
	rm -f $(OBJECTS)
	rm -f efie_rwg_open
