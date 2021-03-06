
EXEC = int2-triaintrouts

#HOST = osx
HOST=linux-gfortran
HOST=linux-gfortran-openmp

ifeq ($(HOST),osx)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -c -w
FLINK = gfortran -w -o $(EXEC)
FEND = -framework accelerate
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math --openmp -c -w  
FLINK = gfortran --openmp -w -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O3 -c -w -xW 
FLINK = ifort -w -mkl -o $(EXEC)
WITH_SECOND = 1
endif


TRIA = ../../src/tria_routs
COM = ../../src/common


.PHONY: all clean list

SOURCES =  test_triaintrouts.f \
  $(COM)/prini_new.f \
  $(TRIA)/koornexps.f90 \
  $(TRIA)/ortho2eva.f \
  $(COM)/dotcross3d.f \
  $(TRIA)/ortho2exps.f \
  $(COM)/lapack_wrap.f90 \
  $(TRIA)/triasymq.f \
  $(TRIA)/ctriaints_main.f \
  $(TRIA)/triatreerouts.f \
  $(TRIA)/triaintrouts.f \

ifeq ($(WITH_SECOND),1)
SOURCES += $(HELLSKITCHEN)/Common/second-r8.f
endif

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
	./$(EXEC) 2 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



