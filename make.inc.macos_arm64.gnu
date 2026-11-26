# makefile overrides
# OS:       macOS
# Compiler: gfortran, clang
# OpenMP:   enabled
# BLAS:     framework acceralate
#
# NOTE for user:
#           Check gfortran version number 
#

CC=gcc
CXX=g++
FC=gfortran

FFLAGS= -fPIC -O3 -arch arm64 -std=legacy -w -mno-outline-atomics
FFLAGS2= -fPIC -arch arm64 -std=legacy -w -mno-outline-atomics
CFLAGS= -fPIC -O3 -arch arm64 -std=c99
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1 -fPIC -O3 -arch arm64 

ifeq ($(PREFIX),)
    FMMBIE_INSTALL_DIR=/usr/local/lib
endif

ifeq ($(PREFIX_FMM),)
    FMM_INSTALL_DIR=/usr/local/lib
endif

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

LBLAS=-framework accelerate

# Wl flags for loading objects from static library into dynamic library
WLDL = -Wl,-force_load
WLDLEND = 


#MATLAB interface:
FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
MFLAGS +=-L${FDIR}
MATLAB_ROOT := $(shell ls -d /Applications/MATLAB_R* | sort | tail -1)
MEX = $(MATLAB_ROOT)/bin/mex

# Link the mex against MATLAB's own OpenMP runtime (LLVM libomp) instead
# of libgomp. MATLAB's OpenBLAS is built against this libomp; with libgomp
# it cannot see our OpenMP parallel regions and multithreads from inside
# them, which gives wrong results. libomp provides the GOMP_* entry points
# that gfortran -fopenmp code calls.
MOMPLIBS = -L$(MATLAB_ROOT)/bin/maca64 -lomp

#MWRAP location
MWRAP=~/git/mwrap/mwrap
