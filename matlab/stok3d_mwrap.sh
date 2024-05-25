# #Compile the .f file:
# gfortran -c -O3 -march=native -std=legacy ../fortran/funcs1.f -o ../fortran/funcs1.o
# gfortran -c -O3 -march=native -std=legacy ../../Jeremy/lap_bel/test8.f -L/usr/local/lib -lfmm3dbie -L/usr/local/lib -lfmm3d  -o test8.o

MWF=stok3d_routs

#Generate the .m file from the .mw file
../../../FastAlgs/mwrap/mwrap -c99complex -list -mex $MWF -mb $MWF.mw
#Generate a .c gateway
../../../FastAlgs/mwrap/mwrap -c99complex -mex $MWF -c $MWF.c $MWF.mw
#Generate mexmaci file to run on Mac
#/Applications/MATLAB_R2023b.app/bin/mex -v $MWF.c  -largeArrayDims -DMWF77_UNDERSCORE1 -l/opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/current/libgfortran.dylib -output $MWF -L/usr/local/lib -lfmm3dbie -L/usr/local/lib -lfmm3d 
FDIR=$(dirname `gfortran --print-file-name libgfortran.dylib`)
/Applications/MATLAB_R2023b.app/bin/mex -v $MWF.c  -largeArrayDims -DMWF77_UNDERSCORE1 -L${FDIR} -output $MWF -L/usr/local/lib -lfmm3dbie -L/usr/local/lib -lfmm3d -lgfortran
