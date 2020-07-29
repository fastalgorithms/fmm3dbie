make install OMP=OFF
cd test
cd helm_wrappers
make -f test_helm_neu.make clean
make -f test_helm_neu.make
