#FFLAGS = -O3 -fallow-invalid-boz
#FFLAGS = -fdec
#FFLAGS=-fopenmp -fno-range-check
## -fallow-invalid-boz is for newer versions of gcc
#FFLAGS=-fopenmp -fdefault-integer-8 -fallow-invalid-boz
FFLAGS=-fopenmp -fdefault-integer-8
GF = gfortran

all:
	${GF} -c $(FFLAGS) sort.f -o sort.o
	${GF} $(FFLAGS) bvh.f sort.o -o bvh
	./bvh

clean:
	rm -rf bvh
