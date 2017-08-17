/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_DEFINED
#define LAPACK_DEFINED

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
    const int *lwork,
    int *info
  );
  void zgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    const cc4s::complex *a,
    const int *lda,
    cc4s::complex *w,
    cc4s::complex *vLeft,
    const int *ldvLeft,
    cc4s::complex *vRight,
    const int *ldvRight,
    cc4s::complex *work,
    const int *lwork,
    double *rwork,
    int *info
  );
};

#endif

