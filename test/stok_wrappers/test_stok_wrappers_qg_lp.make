
EXEC = stok_wrappers

HOST = osx-gfortran
HOST=linux-gfortran
#HOST=linux-gfortran-debug
HOST=linux-gfortran-openmp
#HOST=linux-ifort

ifeq ($(HOST),osx-gfortran)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -c -w --openmp
FLINK = gfortran -w --openmp -o $(EXEC)
FEND = -L/usr/local/lib -lfmm3d -framework accelerate
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -L../../lib -lblas -llapack -lfmm3d -L../../lib
endif

ifeq ($(HOST),linux-gfortran-debug)
FC = gfortran
FFLAGS = -g -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -L../../lib -lblas -llapack -lfmm3d -L../../lib
endif


ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math --openmp -c -w  
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -L../../lib -lfmm3d -lblas -llapack
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O1 -g -xHost -c -w -xW -qopenmp
FLINK = ifort -qopenmp -w -mkl -o $(EXEC)
FEND = -L/usr/local/lib -lfmm3d
endif


COM = ../../src/common
TRIA = ../../src/tria_routs
QUAD = ../../src/quadratures
HELM = ../../src/helm_wrappers
STOK = ../../src/stok_wrappers
KER = ../../src/kernels
SURF = ../../src/surface_routs
FMM = ../../src/fmm_wrappers


.PHONY: all clean list

SOURCES =  test_stok_wrappers_qg_lp.f \
  $(COM)/prini_new.f \
  $(COM)/hkrand.f \
  $(COM)/dlaran.f \
  $(COM)/rigidbodies.f \
  $(SURF)/surf_routs.f90 \
  $(SURF)/xtri_routs/xtri_parameterizations.f90 \
  $(SURF)/xtri_routs/xtri_plot.f90 \
  $(QUAD)/far_field_routs.f90 \
  $(QUAD)/near_field_routs.f \
  $(QUAD)/ggq-quads.f \
  $(QUAD)/ggq-selfquad.f \
  $(QUAD)/ggq-selfquad-routs.f \
  $(KER)/stok_kernels.f90 \
  $(KER)/helm_kernels.f90 \
  $(STOK)/stok_comb_vel.f \
  $(STOK)/stok_s_trac.f \
  $(TRIA)/koornexps.f90 \
  $(COM)/dotcross3d.f90 \
  $(COM)/lapack_wrap.f90 \
  $(COM)/legeexps.f \
  $(COM)/sparse_reps.f \
  $(COM)/sort.f \
  $(COM)/cumsum.f \
  $(COM)/tree_lr_3d.f \
  $(COM)/rotmat_gmres.f \
  $(COM)/setops.f \
  $(TRIA)/triasymq.f \
  $(TRIA)/ctriaints_main.f \
  $(TRIA)/dtriaints_main.f \
  $(TRIA)/dtriaintrouts.f \
  $(TRIA)/triatreerouts.f \
  $(TRIA)/triaintrouts.f  


OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(EXEC)  

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



