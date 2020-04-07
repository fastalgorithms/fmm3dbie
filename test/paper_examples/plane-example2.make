
EXEC = int2-plane

HOST = osx-gfortran
#HOST=linux-gfortran
#HOST=linux-gfortran-lblas
#HOST=linux-gfortran-openmp
#HOST=linux-ifort

ifeq ($(HOST),osx-gfortran)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -c -w --openmp
FLINK = gfortran -w --openmp -o $(EXEC)
FEND = -L../../lib -lfmm3d -framework accelerate
endif

ifeq ($(HOST),linux-gfortran-lblas)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -L../../lib -lblas -llapack -lfmm3d
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -L../../lib -lopenblas -lfmm3d
endif


ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math --openmp -c -w  
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -L../../lib -lfmm3d -lopenblas
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O1 -g -xHost -c -w -xW -qopenmp
FLINK = ifort -qopenmp -w -mkl -o $(EXEC)
FEND = -L../lib -lfmm3d
endif


COM = ../../src/common
TRIA = ../../src/tria_routs
QUAD = ../../src/quadratures
HELM = ../../src/helm_wrappers
KER = ../../src/kernels
SURF = ../../src/surface_routs
FMM = ../../src/fmm_wrappers


.PHONY: all clean list

SOURCES =  plane-example2.f \
  $(COM)/prini_new.f \
  $(COM)/hkrand.f \
  $(COM)/dlaran.f \
  $(SURF)/surf_routs.f90 \
  $(SURF)/vtk_routs.f90 \
  $(SURF)/xtri_routs/xtri_parameterizations.f90 \
  $(SURF)/xtri_routs/xtri_plot.f90 \
  $(SURF)/in_go3.f90 \
  $(QUAD)/far_field_routs.f90 \
  $(QUAD)/near_field_routs.f \
  $(QUAD)/ggq-quads.f \
  $(QUAD)/ggq-selfquad.f \
  $(QUAD)/ggq-pvselfquad.f \
  $(QUAD)/ggq-radial.f \
  $(QUAD)/ggq-pvradial.f \
  $(KER)/helm_kernels.f90 \
  $(TRIA)/koornexps.f90 \
  $(TRIA)/ortho2eva.f90 \
  $(TRIA)/ortho2exps.f90 \
  $(COM)/dotcross3d.f90 \
  $(COM)/lapack_wrap.f90 \
  $(COM)/orthom.f \
  $(COM)/legeexps.f \
  $(COM)/lapack_slow.f \
  $(COM)/sort.f \
  $(COM)/sparse_reps.f \
  $(COM)/cumsum.f \
  $(COM)/tree_lr_3d.f \
  $(COM)/rotmat_gmres.f \
  $(HELM)/helm_comb_dir.f \
  $(FMM)/hfmm3d_ndiv.f \
  $(FMM)/lfmm3d_ndiv.f \
  $(COM)/setops.f \
  $(TRIA)/triasymq.f \
  $(TRIA)/ctriaints_main.f \
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



