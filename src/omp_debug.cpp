#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
int ebrel_omp_thread_limit() {
#ifdef _OPENMP
  // OpenMP 4.0+ function: max number of threads the runtime will ever use
  return omp_get_thread_limit();
#else
  return -1; // signal "no OpenMP"
#endif
}

// [[Rcpp::export]]
int ebrel_omp_for_threads() {
  int nthreads = 1;

#ifdef _OPENMP
#pragma omp parallel
{
  // mimic your usual loop structure
#pragma omp for schedule(dynamic)
  for (int i = 0; i < 1000; ++i) {
    // some dummy work
  }

#pragma omp single
{
  nthreads = omp_get_num_threads();
}
}
#endif

return nthreads;
}
