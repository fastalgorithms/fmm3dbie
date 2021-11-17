#include <sctl.hpp>
#include <kernels.h>

int main() {
  omp_set_num_threads(1);
  std::cout<<"num of max omp threads: "<<omp_get_max_threads()<<"\n";

  int64_t N = 30000007;

  sctl::Vector<double> X(N), V0(N), V1(N), V2(N), V3(N);
  for (auto& x : X) x = (double)rand()/RAND_MAX * 5.0 + 0.1;
  V0 = 0;
  V1 = 0;
  V2 = 0;
  V3 = 0;

  { // Compute Vs
    sctl::Profile::Enable(true);

    sctl::Profile::Tic("Unvectorized");
    for (long i = 0; i < N; i++)
      V0[i] = erf(X[i]);
    sctl::Profile::Toc();

    sctl::Profile::Tic("Vectorized");
    erf_svml_cpp_(&X[0], &N, &V1[0]);
    sctl::Profile::Toc();

    sctl::Profile::Tic("Vectorized [0,1)");
    erf0_cpp_(&X[0], &N, &V2[0]);
    sctl::Profile::Toc();

    sctl::Profile::Tic("Vectorized [1,+inf)");
    erf1_cpp_(&X[0], &N, &V3[0]);
    sctl::Profile::Toc();

    sctl::Profile::Tic("Vectorized [0,+inf)");
    erf2_cpp_(&X[0], &N, &V3[0]);
    sctl::Profile::Toc();

    sctl::Profile::print();
  }

  double max_err = 0, max_val = 0;
  for (long i = 0; i < V0.Dim(); i++) {
    max_err = std::max(max_err, fabs(V0[i]-V1[i]));
    max_val = std::max(max_val, fabs(V0[i]));
  }
  std::cout<<"Relative Error = "<<max_err<<"/"<<max_val<<" = "<<max_err/max_val<<'\n';

  return 0;
}

