
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
CFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp

INCDIR = ../include

FMMBIE_INSTALL_DIR = $(HOME)/lib
FMM_INSTALL_DIR = $(HOME)/lib

FEND = -L${FMMBIE_INSTALL_DIR} $(LLINKLIB) -L${FMM_INSTALL_DIR} $(LFMMLINKLIB) -lm


# For your OS, override the above by placing make variables in make.inc
-include make.inc




.PHONY: usage clean list




default: usage

usage:
	@echo ""
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for boundary integral tests, please specify what to make:"
	@echo ""
	@echo "  make cisurf/"
	@echo "          load ..................... load and plot various base meshes"
	@echo "          skeltri .................. construct a skeleton tri mesh"
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

%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@ -I$(INCDIR) -I.

%.o : %.f  
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@




cisurf/load: cisurf/test_loadmesh.o ../src/cisurf/cisurf_load_mesh.o \
        ../src/cisurf/cisurf_plots.o ../src/common/cprini_long.o
	mkdir -p build
	$(CC) -o build/test_loadmesh cisurf/test_loadmesh.o \
            ../src/cisurf/cisurf_load_mesh.o \
            ../src/common/cprini_long.o \
            ../src/cisurf/cisurf_plots.o $(FEND)
	cd ./build && ./test_loadmesh

cisurf/skeltri: cisurf/test_skeltri.o ../src/cisurf/cisurf_load_mesh.o \
        ../src/cisurf/cisurf_plots.o \
		../src/cisurf/cisurf_skeltri.o ../src/common/cprini_long.o
	mkdir -p build
	$(CC) -o build/test_skeltri cisurf/test_skeltri.o \
            ../src/cisurf/cisurf_load_mesh.o \
            ../src/cisurf/cisurf_plots.o \
            ../src/cisurf/cisurf_skeltri.o \
			../src/common/cprini_long.o $(FEND)
	cd ./build && ./test_skeltri


clean:
	rm -f cisurf/*.o
	rm -f build/*
