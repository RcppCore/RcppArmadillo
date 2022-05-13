
#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::export]]
void rcpparma375(int nthreads) {
    omp_set_num_threads(nthreads);
    #pragma omp parallel for
    for (int i = 0; i < 5; i++){
        int tid = omp_get_thread_num();
        std::cout << tid << "\t tid" << std::endl;
        int nThreads = omp_get_num_threads();
        std::cout << nThreads << "\t nThreads" << std::endl;
    }
}
