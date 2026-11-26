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

# Link the mex against MATLAB's own OpenMP runtime instead of libgomp.
# MATLAB's BLAS/LAPACK (MKL) is threaded with Intel's OpenMP runtime
# (libiomp5); with libgomp it cannot see our OpenMP parallel regions and
# multithreads from inside them (oversubscription). libiomp5 provides the
# GOMP_* entry points that gfortran -fopenmp code calls.
MOMPLIBS = -L$(MATLAB_ROOT)/bin/maci64 -liomp5

#MWRAP location
MWRAP=~/git/mwrap/mwrap
