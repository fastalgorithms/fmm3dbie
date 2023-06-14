
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
	@echo "Makefile for boundary integral examplee, please specify what to make:"
	@echo ""
	@echo "  make maxwell/"
	@echo "         pec/"
	@echo "           aumfie .................. solve the mfie, then resolve for charge"
	@echo "           efie_rwg_open ........... scattering from an open surface"
	@echo "                                     using the electric field int eq"
	@echo "           nrccie................... scattering using the non-resonant"
	@echo "                                     charge-current int eq"
	@echo "         dielectric/"
	@echo "         hybrid/"
	@echo ""
	@echo "  make clean ....................... remove all .o files and executables"
	@echo ""
	@echo "All executables are built into the 'build/' directory, and should be"
	@echo "run from that directory."
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






maxwell/pec/aumfie: maxwell/pec/em_aumfie_example.o
	mkdir -p build
	$(FC) $(FFLAGS) -o build/em_aumfie_example maxwell/pec/em_aumfie_example.o $(FEND)
	cd ./build && ./em_aumfie_example

maxwell/pec/efie_rwg_open: maxwell/pec/efie_rwg_open.o
	mkdir -p build
	$(FC) $(FFLAGS) -o build/efie_rwg_open maxwell/pec/efie_rwg_open.o $(FEND)
	cd ./build && ./efie_rwg_open

maxwell/pec/nrccie: maxwell/pec/em_nrccie_example.o
	mkdir -p build
	rm -f build/*
	$(FC) $(FFLAGS) -o build/em_nrccie_example maxwell/pec/em_nrccie_example.o $(FEND)
	cd ./build && ./em_nrccie_example

clean:
	rm -f maxwell/pec/*.o
	rm -f build/*
