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
    cc4s::Real<64> *a,
    const int *lda,
    cc4s::Real<64> *wReal,
    cc4s::Real<64> *wImag,
    cc4s::Real<64> *vLeft,
    const int *ldvLeft,
    cc4s::Real<64> *vRight,
    const int *ldvRight,
    cc4s::Real<64> *work,
    const int *workSize,
    int *info
  );
  void zgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    const cc4s::Complex<64> *a,
    const int *lda,
    cc4s::Complex<64> *w,
    cc4s::Complex<64> *vLeft,
    const int *ldvLeft,
    cc4s::Complex<64> *vRight,
    const int *ldvRight,
    cc4s::Complex<64> *work,
    const int *workSize,
    cc4s::Complex<64> *rwork,
    int *info
  );
  void dsysv_(
    const char *uplo,
    const int *n,
    const int *m,
    cc4s::Real<64> *a,
    const int *lda,
    const int *ipiv,
    cc4s::Real<64> *b,
    const int *ldb,
    cc4s::Real<64> *work,
    const int *lwork,
    const int *info
  );
  void zsysv_(
    const char *uplo,
    const int *n,
    const int *m,
    cc4s::Complex<64> *a,
    const int *lda,
    const int *ipiv,
    cc4s::Complex<64> *b,
    const int *ldb,
    cc4s::Complex<64> *work,
    const int *lwork,
    const int *info
  );
  void dgetrf_(
    const int *m,
    const int *n,
    cc4s::Real<64> *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void zgetrf_(
    const int *m,
    const int *n,
    cc4s::Complex<64> *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void dgetri_(
    const int *n,
    cc4s::Real<64> *a,
    const int *lda,
    const int *rowPermutation,
    cc4s::Real<64> *work,
    const int *workSize,
    int *info
  );
  void zgetri_(
    const int *n,
    cc4s::Complex<64> *a,
    const int *lda,
    const int *rowPermutation,
    cc4s::Complex<64> *work,
    const int *workSize,
    int *info
  );
}

#endif

