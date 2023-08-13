# Uncomment the next 3 lines if not using modules
module load modules/0-traditional
module load gcc/9.3.0
module load matlab
# Set FMM3D download directory here
FMM_PATH=/mnt/home/mrachh/git/FMM3D
# Set fmm3dbie download directory here
FMMBIE_PATH=/mnt/home/mrachh/git/fmm3dbie
# Set install directory here
INSTALL_DIR=/mnt/home/mrachh/libtmp
CDIR=${PWD}
rm -rf helmquadcorr.o
rm -rf read_plane_geom.o
cd ${FMM_PATH}
#make clean
make install PREFIX=${INSTALL_DIR}
make test
cd ${FMMBIE_PATH}
#make clean
make install PREFIX=${INSTALL_DIR} PREFIX_FMM=${INSTALL_DIR}
make test PREFIX_FMM=${INSTALL_DIR}
make install PREFIX=${INSTALL_DIR} PREFIX_FMM=${INSTALL_DIR} PREFIX_LIBNAME=libfmm3dbie_matlab BLAS_64=ON
cd ${CDIR}
gfortran -fPIC -march=native -O3 -c helmquadcorr.f -o helmquadcorr.o -L${INSTALL_DIR} -lfmm3d -lfmm3dbie_matlab
gfortran -fPIC -march=native -O3 -c read_plane_geom.f -o read_plane_geom.o -L${INSTALL_DIR} -lfmm3d -lfmm3dbie_matlab
mex -v fmm3dbierouts.c helmquadcorr.o read_plane_geom.o -largeArrayDims -DMWF77_UNDERSCORE1 -output fmm3dbierouts -L${INSTALL_DIR} -lfmm3d -lfmm3dbie_matlab -lstdc++ -lm -ldl -lgfortran 
