#location of mex, can be set to ``mex'' on linux systems

MEX=/Applications/MATLAB_R2021a.app/bin/mex

FMM3DINSTALLDIR=/usr/local/lib
#FMM3DINSTALLDIR=/mnt/home/skailasa/lib
FMM3DBIEINSTALLDIR=/usr/local/lib
#FMM3DBIEINSTALLDIR=/mnt/home/skailasa/lib

# Optional on linux systems
GCCPATH=/usr/local/lib/gcc/13

OBJECTS = helmquadcorr.o read_plane_geom.o ellipsoid_routs.o

FC = gfortran
FFLAGS = -fPIC -march=native -O3 -fopenmp
FEND = -L$(FMM3DINSTALLDIR) -lfmm3d -L$(FMM3DBIEINSTALLDIR) -lfmm3dbie_matlab

.PHONY: all clean

%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

all: matlab

matlab: $(OBJECTS) 
	$(MEX) -v fmm3dbierouts.c $(OBJECTS) -largeArrayDims -DMWF77_UNDERSCORE1 -D_OPENMP -L$(GCCPATH) -output fmm3dbierouts -L$(FMM3DINSTALLDIR) -lfmm3d -L$(FMM3DBIEINSTALLDIR) -lfmm3dbie_matlab -lgomp -lstdc++ -lm -ldl -lgfortran

clean:
	rm -f $(OBJECTS)
	rm -f *.mex*
    
