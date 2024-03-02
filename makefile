# Makefile for fmm3dbie
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it


# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -w 

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran


FMMBIE_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FMMBIE_INSTALL_DIR = ${HOME}/lib
endif

FMM_INSTALL_DIR=$(PREFIX_FMM)
ifeq ($(PREFIX_FMM),)
	FMM_INSTALL_DIR=${HOME}/lib
endif

LBLAS = -lblas -llapack
MLBLAS = -lmwblas -lmwlapack

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

LIBNAME=$(PREFIX_LIBNAME)
ifeq ($(LIBNAME),)
	LIBNAME=libfmm3dbie
endif

DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

# dynamic libraries for matlab compilation
MLIBNAME = $(LIBNAME)_matlab
MSTATICLIB = $(MLIBNAME).a
MDYNAMICLIB = $(MLIBNAME).so
MLIMPLIB = $(MDYNAMICLIB)

LFMMLINKLIB = -lfmm3d
LFMMSTATICLIB = $(FMM_INSTALL_DIR)/libfmm3d.a
LLINKLIB = $(subst lib, -l, $(LIBNAME))
MLLINKLIB = $(subst lib, -l, $(MLIBNAME))

WLDL = -Wl,--whole-archive
WLDLEND = -Wl,--no-whole-archive

# For your OS, override the above by placing make variables in make.inc
-include make.inc

# update libs and dynamic libs to include appropriate versions of
# fmm3d
#
# Note: the static library is used for DYLIBS, so that fmm3d 
# does not get bundled in with the fmm3dbie dynamic library
#
#LIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB) 
#DYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)
F2PYDYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)

# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  FFLAGS += $(OMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  F2PYDYLIBS += $(OMPLIBS)
  MEXLIBS += $(OMPLIBS)
endif

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)
MEXLIBS += $(MLBLAS) 


# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/hkrand.o $(COM)/dotcross3d.o \
	$(COM)/dlaran.o \
	$(COM)/rotmat_gmres.o $(COM)/setops.o \
	$(COM)/sort.o $(COM)/sparse_reps.o $(COM)/get_fmm_thresh.o \
	$(COM)/common_Maxwell.o $(COM)/incoming_fields.o \
	$(COM)/rigidbodies.o $(COM)/polytens.o \
	$(COM)/chebexps.o

# Helmholtz wrappers
HELM = src/helm_wrappers
HOBJS = $(HELM)/helm_comb_dir.o $(HELM)/helm_rpcomb_neu.o \
	$(HELM)/helm_comb_trans.o $(HELM)/helm_rpcomb_imp.o \
	$(HELM)/helm_s_neu.o 

# Laplace wrappers
LAP = src/lap_wrappers
LOBJS = $(LAP)/lap_comb_dir.o

# Maxwell wrappers
EM = src/maxwell
EMOBJS = $(EM)/em_mfie_pec.o $(EM)/em_aumfie_pec.o \
	$(EM)/em_nrccie_pec.o $(EM)/em_auCKi_pec.o \
	$(EM)/em_dfie_trans.o $(EM)/em_adpie_pec.o \
	$(EM)/em_sdpie_pec.o $(EM)/em_cfie_rwg_pec.o \
	$(EM)/fix_tri.o $(EM)/analytic_sphere_pw_pec.o

# Stokes wrappers
STOK = src/stok_wrappers
STOKOBJS = $(STOK)/stok_comb_vel.o 

# Kernels
KER = src/kernels
KOBJS = $(KER)/helm_kernels.o $(KER)/lap_kernels.o $(KER)/DPIE_kernels.o \
	$(KER)/yuk_kernels.o $(KER)/stok_kernels.o

# Quadrature wrappers
QUAD = src/quadratures
QOBJS = $(QUAD)/far_field_routs.o \
	$(QUAD)/ggq-selfquad-routs.o $(QUAD)/ggq-quads.o \
	$(QUAD)/ggq-selfquad.o \
	$(QUAD)/near_field_routs.o $(QUAD)/near_quad_sub.o

# Surface wrappers
SURF = src/surface_routs
SOBJS = $(SURF)/in_go3.o $(SURF)/surf_routs.o $(SURF)/vtk_routs.o \
	$(SURF)/xtri_routs/xtri_parameterizations.o \
	$(SURF)/xtri_routs/xtri_plot.o $(SURF)/write_go3.o $(SURF)/in_gidmsh2.o \
	$(SURF)/in_gmsh2.o $(SURF)/patch_basis_routs.o \
	$(SURF)/xquad_routs/xquad_parametrizations.o

# Triangle adaptive integration routines
TRIA = src/tria_routs
TOBJS = $(TRIA)/ctriaints_main.o $(TRIA)/koornexps.o \
	$(TRIA)/triaintrouts.o $(TRIA)/dtriaints_main.o \
	$(TRIA)/triasymq.o $(TRIA)/triatreerouts.o $(TRIA)/dtriaintrouts.o

# Triangle adaptive integration routines
QUAD2 = src/quad_routs
QOBJS2 = $(QUAD2)/cquadints_main.o \
	$(QUAD2)/cquadintrouts.o $(QUAD2)/dquadints_main.o \
	$(QUAD2)/squarearbq.o $(QUAD2)/quadtreerouts.o $(QUAD2)/dquadintrouts.o

SURFSM = src/multiscale_mesher
SURFSMOBJS = $(SURFSM)/cisurf_loadmsh.o \
	$(SURFSM)/cisurf_skeleton.o $(SURFSM)/cisurf_plottools.o \
	$(SURFSM)/cisurf_tritools.o $(SURFSM)/tfmm_setsub.o $(SURFSM)/surface_smoother.o 

SURFSM_MOD_OBJS = $(SURFSM)/Mod_TreeLRD.o \
	$(SURFSM)/ModType_Smooth_Surface.o $(SURFSM)/Mod_Fast_Sigma.o \
	$(SURFSM)/Mod_Plot_Tools_sigma.o $(SURFSM)/Mod_Feval.o $(SURFSM)/Mod_Smooth_Surface.o

# Add to FFLAGS so that modules get compiled in the .mod folder
FFLAGS += -J .mod/

OBJS = $(COMOBJS) $(EMOBJS) $(HOBJS) $(KOBJS) $(LOBJS) $(QOBJS) $(SOBJS) $(TOBJS) $(STOKOBJS) $(QOBJS2)

OBJS_64 = $(COMOBJS) $(EMOBJS) $(HOBJS) $(KOBJS) $(LOBJS) $(QOBJS) $(SOBJS) $(TOBJS) $(STOKOBJS) $(QOBJS2)
OBJS_64 += $(COM)/lapack_wrap_64.o

ifeq ($(BLAS_64),ON)
OBJS += $(COM)/lapack_wrap_64.o
endif

ifneq ($(BLAS_64),ON)
OBJS += $(COM)/lapack_wrap.o
endif


.PHONY: usage lib install test test-dyn python mex mex-dyn matlab-dyn matlab surf-smooth-objs surf-smooth-test 

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for fmm3dbie. Specify what to make:"
	@echo "  make install      compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR)  "
	@echo "                    compile and install the main library at custom"
	@echo "                    location given by PREFIX"
	@echo "  make lib          compile the main library (in lib/ and lib-static/)"
	@echo "  make test         compile and run validation tests (will take around 30 secs)"
	@echo "  make test-dyn     test successful installation by validation tests linked "
	@echo "                    to dynamic library (will take a couple of mins)"
	@echo "  make matlab       compile matlab interfaces with static library linking"
	@echo "  make matlab-dyn   compile matlab interfaces with dynamic library linking"
	@echo "                    (static linking is preferred over dynamic linking)"
	@echo "  make python       compile and test python interfaces using python"
	@echo "  make objclean     removal all object files, preserving lib & MEX"
	@echo "  make clean        also remove lib, MEX, py, and demo executables"
	@echo "  make mex          generate matlab interfaces"
	@echo "                    (for expert users only, requires mwrap)"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=OFF' for single-threaded"
	@echo "-------------------------------------------------------------------------"




#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@


#
# build the library...
#
lib: $(STATICLIB) $(DYNAMICLIB)

lib-fmm3dbie-only: STATICLIBFMM3DBIE DYNAMICLIBFMM3DBIE	

ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif

STATICLIBFMM:
ifneq ($(wildcard ./FMM3D/.*),)
	[ ! -f make.inc ] || cp make.inc ./FMM3D; 
	cd FMM3D && make lib -j; 
	echo "$(LFMMSTATICLIB)"; 
	$(eval LFMMSTATICLIB := $(shell pwd)/FMM3D/lib-static/libfmm3d.a) 
	echo "$(LFMMSTATICLIB)"; 
endif

STATICLIBFMM3DBIE: $(OBJS)
	ar rcs $(STATICLIB) $(OBJS) 
	mv $(STATICLIB) lib-static/
	echo $(LFMMSTATICLIB)
	cd lib-static && ar -x $(STATICLIB)
	cd lib-static && ar -x $(LFMMSTATICLIB)
	cd lib-static && ar rcs $(STATICLIB) *.o
	cd lib-static && rm -rf *.o *__*

MSTATICLIBFMM3DBIE: $(OBJS_64)
	echo "$(OBJS_64)"
	ar rcs $(MSTATICLIB) $(OBJS_64) 
	mv $(MSTATICLIB) lib-static/
	cd lib-static && ar -x $(MSTATICLIB)
	cd lib-static && ar -x $(LFMMSTATICLIB)
	cd lib-static && ar rcs $(MSTATICLIB) *.o
	cd lib-static && rm -rf *.o *__*



DYNAMICLIBFMM3DBIE: STATICLIBFMM3DBIE 
	$(FC) -shared -fPIC $(FFLAGS) $(WLDL) lib-static/$(STATICLIB) $(WLDLEND) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/

$(STATICLIB): STATICLIBFMM STATICLIBFMM3DBIE

$(DYNAMICLIB): STATICLIBFMM DYNAMICLIBFMM3DBIE

$(MSTATICLIB): STATICLIBFMM MSTATICLIBFMM3DBIE

$(MDYNAMICLIB): $(MSTATICLIB) 
	$(FC) -shared -fPIC $(FFLAGS) $(WLDL) lib-static/$(MSTATICLIB) $(WLDLEND) -o $(MDYNAMICLIB) $(DYLIBS) 
	mv $(MDYNAMICLIB) lib/
	[ ! -f $(MLIMPLIB) ] || mv $(MLIMPLIB) lib/
	mkdir -p $(FMMBIE_INSTALL_DIR)
	cp -f lib/$(MDYNAMICLIB) $(FMMBIE_INSTALL_DIR)/

install: $(STATICLIB) $(DYNAMICLIB)
	rm -rf tmp
	echo $(FMMBIE_INSTALL_DIR)
	mkdir -p $(FMMBIE_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FMMBIE_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FMMBIE_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FMMBIE_INSTALL_DIR)/
	@echo "Make sure to include " $(FMMBIE_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FMMBIE_INSTALL_DIR)  " "$(LLINKLIB) 

surf-smooth-mod-objs: $(SURFSM_MOD_OBJS) 

surf-smooth-objs: surf-smooth-mod-objs  $(SURFSMOBJS)


# matlab..
#
MWDIR = matlab
MWF = fmm3dbie_routs
GW = $(MWF)



matlab:	$(MSTATICLIB) $(MWDIR)/$(GW).c
	$(MEX) $(MWDIR)/$(GW).c lib-static/$(MSTATICLIB) $(MFLAGS) \
	-output $(MWDIR)/$(GW) $(MEXLIBS) 


matlab-dyn:	$(MDYNAMICLIB) $(MWDIR)/$(GW).c
	$(MEX) $(MWDIR)/$(GW).c $(MFLAGS) \
	-output $(MWDIR)/$(GW) $(MEXLIBS) -L$(FMMBIE_INSTALL_DIR) $(MLLINKLIB) 

mex: $(MSTATICLIB)
	cd $(MWDIR); $(MWRAP) $(MWFLAGS) -list -mex $(GW) -mb $(MWF).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GW) -c $(GW).c $(MWF).mw;\
	$(MEX) $(GW).c ../lib-static/$(MSTATICLIB) $(MFLAGS) \
	-output $(GW) $(MEXLIBS) 

mex-dyn: $(MDYNAMICLIB)
	cd $(MWDIR); $(MWRAP) $(MWFLAGS) -list -mex $(GW) -mb $(MWF).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GW) -c $(GW).c $(MWF).mw;\
	$(MEX) $(GW).c $(MFLAGS) \
	-output $(GW) $(MEXLIBS) -L$(FMMBIE_INSTALL_DIR) $(MLLINKLIB); 

#
# testing routines
#
test: $(STATICLIB) test/com test/hwrap test/tria test/lwrap test/surf test/quadrature test/quad 
	cd test/common; ./int2-com
	cd test/helm_wrappers; ./int2-helm
	cd test/lap_wrappers; ./int2-lap
	cd test/surface_routs; ./int2-surf
	cd test/tria_routs; ./int2-tria
	cd test/quad_routs; ./int2-quad
	cd test/quadratures; ./int2-quad
	cat print_testres.txt
	rm print_testres.txt

test-dyn: $(DYNAMICLIB) test/com-dyn test/hwrap-dyn test/tria-dyn test/lwrap-dyn test/surf-dyn test/quadrature-dyn test/quad-dyn
	cd test/common; ./int2-com
	cd test/helm_wrappers; ./int2-helm
	cd test/lap_wrappers; ./int2-lap
	cd test/surface_routs; ./int2-surf
	cd test/tria_routs; ./int2-tria
	cd test/quadratures; ./int2-quad
	cat print_testres.txt
	rm print_testres.txt

test-fmm3dbie-only: STATICLIBFMM3DBIE test/com test/hwrap test/tria test/lwrap test/surf test/quadrature test/quad 
	cd test/common; ./int2-com
	cd test/helm_wrappers; ./int2-helm
	cd test/lap_wrappers; ./int2-lap
	cd test/surface_routs; ./int2-surf
	cd test/tria_routs; ./int2-tria
	cd test/quad_routs; ./int2-quad
	cd test/quadratures; ./int2-quad
	cat print_testres.txt
	rm print_testres.txt

test-dyn-fmm3dbie-only: DYNAMICLIBFMM3DBIE test/com-dyn test/hwrap-dyn test/tria-dyn test/lwrap-dyn test/surf-dyn test/quadrature-dyn test/quad-dyn
	cd test/common; ./int2-com
	cd test/helm_wrappers; ./int2-helm
	cd test/lap_wrappers; ./int2-lap
	cd test/surface_routs; ./int2-surf
	cd test/tria_routs; ./int2-tria
	cd test/quadratures; ./int2-quad
	cat print_testres.txt
	rm print_testres.txt

test/com: 
	$(FC) $(FFLAGS) test/common/test_common.f -o test/common/int2-com lib-static/$(STATICLIB) $(LIBS) 

test/hwrap:
	$(FC) $(FFLAGS) test/helm_wrappers/test_helm_wrappers_qg_lp.f -o test/helm_wrappers/int2-helm lib-static/$(STATICLIB) $(LIBS) 

test/lwrap:
	$(FC) $(FFLAGS) test/lap_wrappers/test_lap_wrappers_qg_lp.f -o test/lap_wrappers/int2-lap lib-static/$(STATICLIB) $(LIBS) 

test/surf:
	$(FC) $(FFLAGS) test/surface_routs/test_surf_routs.f -o test/surface_routs/int2-surf lib-static/$(STATICLIB) $(LIBS) 

TTOBJS = test/tria_routs/test_triaintrouts.o test/tria_routs/test_dtriaintrouts.o test/tria_routs/test_koornexps.o

test/tria: $(TTOBJS)
	$(FC) $(FFLAGS) test/tria_routs/test_triarouts.f -o test/tria_routs/int2-tria $(TTOBJS) lib-static/$(STATICLIB) $(LIBS) 

QTOBJS = test/quad_routs/test_quadintrouts.o test/quad_routs/test_dquadintrouts.o 

test/quad: $(QTOBJS)
	$(FC) $(FFLAGS) test/quad_routs/test_quadrouts.f -o test/quad_routs/int2-quad $(QTOBJS) lib-static/$(STATICLIB) $(LIBS) 

test/quadrature:
	$(FC) $(FFLAGS) test/quadratures/test_find_near.f -o test/quadratures/int2-quad lib-static/$(STATICLIB) $(LIBS) 


#
# Linking test files to dynamic libraries
#


test/com-dyn:
	$(FC) $(FFLAGS) test/common/test_common.f -o test/common/int2-com -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS) $(LFMMLINKLIB) 

test/hwrap-dyn:
	$(FC) $(FFLAGS) test/helm_wrappers/test_helm_wrappers_qg_lp.f -o test/helm_wrappers/int2-helm -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS)

test/lwrap-dyn:
	$(FC) $(FFLAGS) test/lap_wrappers/test_lap_wrappers_qg_lp.f -o test/lap_wrappers/int2-lap -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS)

test/surf-dyn:
	$(FC) $(FFLAGS) test/surface_routs/test_surf_routs.f -o test/surface_routs/int2-surf -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS) 

test/tria-dyn: $(TTOBJS)
	$(FC) $(FFLAGS) test/tria_routs/test_triarouts.f -o test/tria_routs/int2-tria $(TTOBJS) -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS) 

test/quad-dyn: $(QTOBJS)
	$(FC) $(FFLAGS) test/quad_routs/test_quadrouts.f -o test/quad_routs/int2-quad $(QTOBJS) -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS)


test/quadrature-dyn:
	$(FC) $(FFLAGS) test/quadratures/test_find_near.f -o test/quadratures/int2-quad -L$(FMMBIE_INSTALL_DIR) $(LLINKLIB) $(LIBS) 


# 
#  Surface smoother test
#
surf-smooth-test: surf-smooth-objs test/surf-smooth
	cd test/multiscale_mesher; ./int2-surfsmooth 

test/surf-smooth:
	$(FC) $(FFLAGS) test/multiscale_mesher/test_surfsmooth.f90 -o test/multiscale_mesher/int2-surfsmooth $(SURFSMOBJS) $(SURFSM_MOD_OBJS) lib-static/$(STATICLIB) $(LIBS) 

#
# build the python bindings/interface
#
python: $(STATICLIB)
	cd python && export FMMBIE_LIBS='$(LIBS)' && pip install -e . 


#
# build python gmsh to go3
#
python-gmsh: $(DYNAMICLIB)
	f2py $(F2PYDYLIBS) -lfmm3dbie -c src/tria_routs/koornexps.f90 -m kexp
	f2py $(F2PYDYLIBS) -lfmm3dbie -c src/surface_routs/surf_routs.f90 -m srout
	mv kexp*.so python/
	mv srout*.so python/
	rm -rf *.dSYM

#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f .mod/*
	rm -f test/common/int2-com
	rm -f test/helm_wrappers/int2-helm
	rm -f test/tria_routs/int2-tria
	rm -f python/*.so
	rm -rf python/build
	rm -rf python/fmm3dpy.egg-info
	rm -rf python/kexp*.so
	rm -rf python/srout*.so

objclean: 
	rm -f $(OBJS) $(TTOBJS) $(QTOBJS) $(OBJS_64) $(SURFSMOBJS) $(SURFSM_MOD_OBJS) 
	rm -f test/helm_wrappers/*.o test/common/*.o 
	rm -f test/tria_routs/*.o 
	rm -f test/lap_wrappers/*.o
	rm -f test/maxwell_wrappers/*.o
	rm -f test/quadratures/*.o
	rm -f test/quad_routs/*.o
	rm -f test/stok_wrappers/*.o
	rm -f test/surface_routs/*.o
