# makefile overrides
# OS:       macOS
# Compiler: gfortran 9.X
# OpenMP:   enabled
# BLAS:     framework acceralate
#
# NOTE for user:
#           Check gfortran version number 
#

CC=gcc-9
CXX=g++-9
FC=gfortran-9
FFLAGS= -fPIC -O3 -march=native -funroll-loops 

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



