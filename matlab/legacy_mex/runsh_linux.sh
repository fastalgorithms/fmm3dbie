rm -rf helmquadcorr.o
MWPATH=/cm/shared/sw/pkg/vendor/matlab/R2020a/bin/glnxa64
gfortran -fPIC -march=native -O3 -c helmquadcorr.f -o helmquadcorr.o -L/mnt/home/mrachh/lib -lfmm3d -lfmm3dbie_matlab
gfortran -fPIC -march=native -O3 -c read_plane_geom.f -o read_plane_geom.o -L/mnt/home/mrachh/lib -lfmm3d -lfmm3dbie_matlab
../../mwrap/mwrap -c99complex -list -mex fmm3dbierouts -mb fmm3dbierouts.mw
../../mwrap/mwrap -c99complex -mex -fmm3dbierouts -c fmm3dbierouts.c fmm3dbierouts.mw
mex -v fmm3dbierouts.c helmquadcorr.o read_plane_geom.o -largeArrayDims -DMWF77_UNDERSCORE1 -output fmm3dbierouts -L/mnt/home/mrachh/lib -lfmm3d -lfmm3dbie_matlab -lstdc++ -lm -ldl -lgfortran 
