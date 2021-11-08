/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

