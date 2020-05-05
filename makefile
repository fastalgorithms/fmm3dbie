# Makefile for fmm3dbie
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops 


CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS) 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1
CXXFLAGS+=$(FFLAGS)

CLIBS = -lgfortran -lm -ldl 

# Update blas
ifneq ($(OS),Windows_NT) 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LDF = /usr/local/lib
    endif
    ifeq ($(UNAME_S),Linux)
        LDF = ${LD_LIBRARY_PATH}
    endif
endif


# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap-0.33/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

ifeq ($(FAST_KER),ON)
  LIBS += -lstdc++
  CLIBS += -lstdc++
  OMP = ON
endif

LBLAS = -lblas -llapack

# For your OS, override the above by placing make variables in make.inc
-include make.inc

LIBS = -L${LDF} -lm -lfmm3d ${LBLAS} ${LDBLAS}

# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
CFLAGS += $(OMPFLAGS)
FFLAGS += $(OMPFLAGS)
MFLAGS += $(MOMPFLAGS)
LIBS += $(OMPLIBS)
MEXLIBS += $(OMPLIBS)
endif


LIBNAME=libfmm3dbie
DYNAMICLIB = $(LIBNAME).so
STATICLIB = lib-static/$(LIBNAME).a

# vectorized kernel directory
SRCDIR = ./vec-kernels/src
INCDIR = ./vec-kernels/include
LIBDIR = lib-static

# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/cumsum.o $(COM)/dlaran.o $(COM)/dotcross3d.o \
	$(COM)/hkrand.o $(COM)/lapack_slow.o $(COM)/lapack_wrap.o \
	$(COM)/legeexps.o $(COM)/orthom.o $(COM)/prini_new.o \
	$(COM)/rotmat_gmres.o $(COM)/setops.o \
	$(COM)/sort.o $(COM)/sparse_reps.o $(COM)/tree_lr_3d.o

# FMM wrappers
FMML = src/fmm_wrappers
FOBJS = $(FMML)/hfmm3d_ndiv.o $(FMML)/lfmm3d_ndiv.o 

# Helmholtz wrappers
HELM = src/helm_wrappers
HOBJS = $(HELM)/helm_comb_dir.o

# Laplace wrappers
LAP = src/lap_wrappers
LOBJS = $(LAP)/lap_comb_dir.o

# Kernels
KER = src/kernels
KOBJS = $(KER)/helm_kernels.o $(KER)/lap_kernels.o

# Quadrature wrappers
QUAD = src/quadratures
QOBJS = $(QUAD)/far_field_routs.o \
	$(QUAD)/ggq-selfquad-routs.o $(QUAD)/ggq-quads.o \
	$(QUAD)/ggq-selfquad.o \
	$(QUAD)/near_field_routs.o

# Surface wrappers
SURF = src/surface_routs
SOBJS = $(SURF)/in_go3.o $(SURF)/surf_routs.o $(SURF)/vtk_routs.o \
	$(SURF)/xtri_routs/xtri_parameterizations.o \
	$(SURF)/xtri_routs/xtri_plot.o

# Triangle adaptive integration routines
TRIA = src/tria_routs
TOBJS = $(TRIA)/ctriaints_main.o $(TRIA)/koornexps.o \
	$(TRIA)/triaintrouts.o $(TRIA)/dtriaints_main.o \
	$(TRIA)/triasymq.o $(TRIA)/triatreerouts.o $(TRIA)/dtriaintrouts.o


OBJS = $(COMOBJS) $(FOBJS) $(HOBJS) $(KOBJS) $(LOBJS) $(QOBJS) $(SOBJS) $(TOBJS)

.PHONY: usage lib python python3 

default: usage


usage:
	@echo "Makefile for FMM3D. Specify what to make:"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take around 30 secs)"
	@echo "  make python - compile and test python interfaces"
	@echo "  make python3 - compile and test python interfaces using python3"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif
$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(OMPFLAGS) $(OBJS) -o $(DYNAMICLIB) $(LIBS) 
	mv $(DYNAMICLIB) lib/


#
# testing routines
#

test: $(STATICLIB) test/com test/hwrap test/tria test/lwrap test/surf
	cd test/common; ./int2-com
	cd test/helm_wrappers; ./int2-helm
	cd test/lap_wrappers; ./int2-lap
	cd test/surface_routs; ./int2-surf
	cd test/tria_routs; ./int2-tria
	cat print_testres.txt
	rm print_testres.txt

test/com:
	$(FC) $(FFLAGS) test/common/test_common.f -o test/common/int2-com $(STATICLIB) $(LIBS) 

test/hwrap:
	$(FC) $(FFLAGS) test/helm_wrappers/test_helm_wrappers_qg_lp.f -o test/helm_wrappers/int2-helm $(STATICLIB) $(LIBS) 

test/lwrap:
	$(FC) $(FFLAGS) test/lap_wrappers/test_lap_wrappers_qg_lp.f -o test/lap_wrappers/int2-lap $(STATICLIB) $(LIBS) 

test/surf:
	$(FC) $(FFLAGS) test/surface_routs/test_surf_routs.f -o test/surface_routs/int2-surf $(STATICLIB) $(LIBS) 

TTOBJS = test/tria_routs/test_triaintrouts.o test/tria_routs/test_dtriaintrouts.o test/tria_routs/test_koornexps.o

test/tria: $(TTOBJS)
	$(FC) $(FFLAGS) test/tria_routs/test_triarouts.f -o test/tria_routs/int2-tria $(TTOBJS) $(STATICLIB) $(LIBS) 

#python
python: $(DYNAMICLIB)
	cd python && pip install -e . 

#python
python3: $(STATICLIB)
	cd python && export FAST_KER=$(FAST_KER) && export FLIBS='$(LIBS)' && export FFLAGS='$(FFLAGS)' && pip3 install -e . 

clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/common/int2-com
	rm -f test/helm_wrappers/int2-helm
	rm -f test/tria_routs/int2-tria
	rm -f python/*.so
	rm -rf python/build
	rm -rf python/fmm3dpy.egg-info

objclean: 
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/helm_wrappers/*.o test/common/*.o 
	rm -f test/tria_routs/*.o examples/helm_dir/*.o 
