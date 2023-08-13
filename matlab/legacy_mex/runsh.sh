rm -rf helmquadcorr.o
gfortran -O3 -march=native -fPIC -fopenmp -c helmquadcorr.f -o helmquadcorr.o -L/usr/local/lib -lfmm3d -lfmm3dbie_matlab
gfortran -fPIC -march=native -O3 -c read_plane_geom.f -o read_plane_geom.o -L/mnt/home/mrachh/lib -lfmm3d -lfmm3dbie_matlab
gfortran -fPIC -march=native -O3 -c ellipsoid_routs.f90 -o ellipsoid_routs.o -L/mnt/home/mrachh/lib -lfmm3d -lfmm3dbie_matlab
../../mwrap/mwrap -c99complex -list -mex fmm3dbierouts -mb fmm3dbierouts.mw
../../mwrap/mwrap -c99complex -mex -fmm3dbierouts -c fmm3dbierouts.c fmm3dbierouts.mw
#/Applications/MATLAB_R2021a.app/bin/mex -v fmm3dbierouts.c helmquadcorr.o read_plane_geom.o -largeArrayDims -DMWF77_UNDERSCORE1 -D_OPENMP -L/usr/local/lib/gcc/13 -output fmm3dbierouts -L/usr/local/lib -lfmm3d -lfmm3dbie_matlab -lgomp -lstdc++ -lm -ldl -lgfortran
