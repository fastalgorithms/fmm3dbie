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
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy 

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


#MATLAB interface:
FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
MFLAGS +=-L${FDIR}
MEX = $(shell ls -d /Applications/MATLAB_R* | sort | tail -1)/bin/mex

#MWRAP location
MWRAP=~/git/mwrap/mwrap
