// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Blas][Blas:1]]
#pragma once
namespace atrip {

  using Complex = std::complex<double>;

  extern "C" {
    void dgemm_(
      const char *transa,
      const char *transb,
      const int *m,
      const int *n,
      const int *k,
      double *alpha,
      const double *a,
      const int *lda,
      const double *b,
      const int *ldb,
      double *beta,
      double *c,
      const int *ldc
    );

    void zgemm_(
      const char *transa,
      const char *transb,
      const int *m,
      const int *n,
      const int *k,
      Complex *alpha,
      const Complex *A,
      const int *lda,
      const Complex *B,
      const int *ldb,
      Complex *beta,
      Complex *C,
      const int *ldc
    );
  }


  template <typename F=double>
  void xgemm(const char *transa,
             const char *transb,
             const int *m,
             const int *n,
             const int *k,
             F *alpha,
             const F *A,
             const int *lda,
             const F *B,
             const int *ldb,
             F *beta,
             F *C,
             const int *ldc) {
    dgemm_(transa, transb,
           m, n, k,
           alpha, A, lda,
           B, ldb, beta,
           C, ldc);
  }

  template <>
  void xgemm(const char *transa,
             const char *transb,
             const int *m,
             const int *n,
             const int *k,
             Complex *alpha,
             const Complex *A,
             const int *lda,
             const Complex *B,
             const int *ldb,
             Complex *beta,
             Complex *C,
             const int *ldc) {
    zgemm_(transa, transb,
           m, n, k,
           alpha, A, lda,
           B, ldb, beta,
           C, ldc);
  }
}
// Blas:1 ends here
