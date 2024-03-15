# makefile overrides
# OS:       macOS
# Compiler: gfortran X.X, clang
# OpenMP:   enabled
# BLAS:     framework acceralate
#
# NOTE for user:
#           Check gfortran version number 
#

CC=gcc
CXX=g++
FC=gfortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops 

ifeq ($(PREFIX),)
    FMMBIE_INSTALL_DIR=/usr/local/lib
endif

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

LBLAS=-framework accelerate



