
#
# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# and libfmm3d.so, for Macosx, these .so files also need to be
# copied over /usr/local/lib


LFMMLINKLIB = -lfmm3d
LLINKLIB = -lfmm3dbie


FC = gfortran
FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy 
CC = gcc


FMMBIE_INSTALL_DIR = $(HOME)/lib
FMM_INSTALL_DIR = $(HOME)/lib
FEND = -L${FMMBIE_INSTALL_DIR} $(LLINKLIB) -L${FMM_INSTALL_DIR} $(LFMMLINKLIB)


# For your OS, override the above by placing make variables in make.inc
-include make.inc




.PHONY: usage clean list




default: usage

usage:
	@echo ""
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for Maxwell examplee, please specify what to make:"
	@echo ""
	@echo "  make efie_rwg_open ........... scattering from an open surface"
	@echo "  make clean ................... remove all .o files and executables"
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
	rm -f efie_rwg_open.o
	rm -f efie_rwg_open
