EXEC = int2-iter
#HOST = gcc
#HOST = gcc-openmp
#HOST = intel
HOST = intel-openmp

#
# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# and libsolvers3d.so, for Macosx, these .so files also need to be
# copied over /usr/local/lib


ifneq ($(OS),Windows_NT) 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LDF = /usr/local/lib
    endif
    ifeq ($(UNAME_S),Linux)
        LDF = ${HOME}/lib
    endif
endif


LIBS = -lfmm3d -lfmm3dbie

ifeq ($(HOST),gcc)
    FC=gfortran -L${LDF} 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native  
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp 
endif

ifeq ($(HOST),intel)
    FC=ifort -L${LDF} 
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort -L${LDF} 
    FFLAGS= -O3 -fPIC -march=native -qopenmp
    FEND = -mkl
endif

SURF=../../src/surface_routs

.PHONY: all clean 

OBJECTS =  helm_dir_iter_example.o \


#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f  
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) -L${LDF} $(LIBS) $(FEND)
	./$(EXEC)  

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



