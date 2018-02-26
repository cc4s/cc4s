/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_DEFINED
#define LAPACK_DEFINED

#include <math/Complex.hpp>

// TODO: use define for name mangling: underscore or not

extern "C" {
  void dgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    cc4s::FloatTypes<8>::real *a,
    const int *lda,
    cc4s::FloatTypes<8>::real *wReal,
    cc4s::FloatTypes<8>::real *wImag,
    cc4s::FloatTypes<8>::real *vLeft,
    const int *ldvLeft,
    cc4s::FloatTypes<8>::real *vRight,
    const int *ldvRight,
    cc4s::FloatTypes<8>::real *work,
    const int *workSize,
    int *info
  );
  void zgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    const cc4s::FloatTypes<8>::complex *a,
    const int *lda,
    cc4s::FloatTypes<8>::complex *w,
    cc4s::FloatTypes<8>::complex *vLeft,
    const int *ldvLeft,
    cc4s::FloatTypes<8>::complex *vRight,
    const int *ldvRight,
    cc4s::FloatTypes<8>::complex *work,
    const int *workSize,
    cc4s::FloatTypes<8>::real *rwork,
    int *info
  );

  void dgetrf_(
    const int *m,
    const int *n,
    cc4s::FloatTypes<8>::real *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void zgetrf_(
    const int *m,
    const int *n,
    cc4s::FloatTypes<8>::complex *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void dgetri_(
    const int *n,
    cc4s::FloatTypes<8>::real *a,
    const int *lda,
    const int *rowPermutation,
    cc4s::FloatTypes<8>::real *work,
    const int *workSize,
    int *info
  );
  void zgetri_(
    const int *n,
    cc4s::FloatTypes<8>::complex *a,
    const int *lda,
    const int *rowPermutation,
    cc4s::FloatTypes<8>::complex *work,
    const int *workSize,
    int *info
  );
};

#endif

