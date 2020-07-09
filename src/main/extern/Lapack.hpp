/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_DEFINED
#define LAPACK_DEFINED

#include <math/Real.hpp>
#include <math/Complex.hpp>

// TODO: use define for name mangling: underscore or not

extern "C" {
  void dgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    double *a,
    const int *lda,
    double *wReal,
    double *wImag,
    double *vLeft,
    const int *ldvLeft,
    double *vRight,
    const int *ldvRight,
    double *work,
    const int *workSize,
    int *info
  );
  void zgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    const complex<double> *a,
    const int *lda,
    complex<double> *w,
    complex<double> *vLeft,
    const int *ldvLeft,
    complex<double> *vRight,
    const int *ldvRight,
    complex<double> *work,
    const int *workSize,
    complex<double> *rwork,
    int *info
  );
  void dsysv_(
    const char *uplo,
    const int *n,
    const int *m,
    double *a,
    const int *lda,
    const int *ipiv,
    double *b,
    const int *ldb,
    double *work,
    const int *lwork,
    const int *info
  );
  void dgetrf_(
    const int *m,
    const int *n,
    double *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void zgetrf_(
    const int *m,
    const int *n,
    complex<double> *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void dgetri_(
    const int *n,
    double *a,
    const int *lda,
    const int *rowPermutation,
    double *work,
    const int *workSize,
    int *info
  );
  void zgetri_(
    const int *n,
    complex<double> *a,
    const int *lda,
    const int *rowPermutation,
    complex<double> *work,
    const int *workSize,
    int *info
  );
};

#endif

