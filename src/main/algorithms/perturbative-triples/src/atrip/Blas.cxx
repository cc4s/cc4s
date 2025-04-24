// Copyright 2022 Alejandro Gallo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// [[file:~/cuda/atrip/atrip.org::*Blas][Blas:2]]
#include <atrip/Blas.hpp>
#include <atrip/Acc.hpp>
#include <atrip/Atrip.hpp>

#if defined(HAVE_ACC)
#  include <cstring>

static size_t dgem_call = 0;

#  if defined(HAVE_CUDA)

static inline cublasOperation_t char_to_accblasOperation(const char *trans) {
  if (strncmp("N", trans, 1) == 0) return CUBLAS_OP_N;
  else if (strncmp("T", trans, 1) == 0) return CUBLAS_OP_T;
  else return CUBLAS_OP_C;
}

#  elif defined(HAVE_HIP)

static inline hipblasOperation_t char_to_accblasOperation(const char *trans) {
  if (strncmp("N", trans, 1) == 0) return HIPBLAS_OP_N;
  else if (strncmp("T", trans, 1) == 0) return HIPBLAS_OP_T;
  else return HIPBLAS_OP_C;
}

#  endif

#endif

namespace atrip {

template <>
void xgemm<float>(const char *transa,
                  const char *transb,
                  const int *m,
                  const int *n,
                  const int *k,
                  float *alpha,
                  const typename DataField<float>::type *A,
                  const int *lda,
                  const typename DataField<float>::type *B,
                  const int *ldb,
                  float *beta,
                  typename DataField<float>::type *C,
                  const int *ldc) {
#if defined(HAVE_ACC)
  // TODO: remove this verbose checking
  const ACC_BLAS_STATUS error = ACC_BLAS_SGEMM(Atrip::cuda.handle,
                                               char_to_accblasOperation(transa),
                                               char_to_accblasOperation(transb),
                                               *m,
                                               *n,
                                               *k,
                                               alpha,
                                               A,
                                               *lda,
                                               B,
                                               *ldb,
                                               beta,
                                               C,
                                               *ldc);
  if (error != 0)
    printf(
        ":%-3ld (%4ldth) ERR<%4d> blasDgemm: "
        "A = %20ld "
        "B = %20ld "
        "C = %20ld "
        "\n",
        Atrip::rank,
        dgem_call++,
        error,
        A,
        B,
        C);
#else
  sgemm_(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}

template <>
void xgemm<double>(const char *transa,
                   const char *transb,
                   const int *m,
                   const int *n,
                   const int *k,
                   double *alpha,
                   const typename DataField<double>::type *A,
                   const int *lda,
                   const typename DataField<double>::type *B,
                   const int *ldb,
                   double *beta,
                   typename DataField<double>::type *C,
                   const int *ldc) {
#if defined(HAVE_ACC)
  // TODO: remove this verbose checking
  const ACC_BLAS_STATUS error = ACC_BLAS_DGEMM(Atrip::cuda.handle,
                                               char_to_accblasOperation(transa),
                                               char_to_accblasOperation(transb),
                                               *m,
                                               *n,
                                               *k,
                                               alpha,
                                               A,
                                               *lda,
                                               B,
                                               *ldb,
                                               beta,
                                               C,
                                               *ldc);
  if (error != 0)
    printf(
        ":%-3ld (%4ldth) ERR<%4d> blasDgemm: "
        "A = %20ld "
        "B = %20ld "
        "C = %20ld "
        "\n",
        Atrip::rank,
        dgem_call++,
        error,
        A,
        B,
        C);
#else
  dgemm_(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}

template <>
void xgemm<Complex>(const char *transa,
                    const char *transb,
                    const int *m,
                    const int *n,
                    const int *k,
                    Complex *alpha,
                    const typename DataField<Complex>::type *A,
                    const int *lda,
                    const typename DataField<Complex>::type *B,
                    const int *ldb,
                    Complex *beta,
                    typename DataField<Complex>::type *C,
                    const int *ldc) {
#if defined(HAVE_ACC)
  ACC_DOUBLE_COMPLEX cu_alpha = {std::real(*alpha), std::imag(*alpha)},
                     cu_beta = {std::real(*beta), std::imag(*beta)};

  ACC_CHECK_BLAS("cublasZgemm",
                 ACC_BLAS_ZGEMM(Atrip::cuda.handle,
                                char_to_accblasOperation(transa),
                                char_to_accblasOperation(transb),
                                *m,
                                *n,
                                *k,
                                &cu_alpha,
                                A,
                                *lda,
                                B,
                                *ldb,
                                &cu_beta,
                                C,
                                *ldc));
#else
  zgemm_(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}

template <>
void xcopy<float>(int *n,
                  const DataFieldType<float> *x,
                  int *incx,
                  DataFieldType<float> *y,
                  int *incy) {
#if defined(HAVE_ACC)
  ACC_BLAS_SCOPY(Atrip::cuda.handle, *n, x, *incx, y, *incy);
#else
  scopy_(n, x, incx, y, incy);
#endif
}

template <>
void xcopy<double>(int *n,
                   const DataFieldType<double> *x,
                   int *incx,
                   DataFieldType<double> *y,
                   int *incy) {
#if defined(HAVE_ACC)
  ACC_BLAS_DCOPY(Atrip::cuda.handle, *n, x, *incx, y, *incy);
#else
  dcopy_(n, x, incx, y, incy);
#endif
}

template <>
void xcopy<Complex>(int *n,
                    const DataFieldType<Complex> *x,
                    int *incx,
                    DataFieldType<Complex> *y,
                    int *incy) {
#if defined(HAVE_ACC)
  ACC_BLAS_ZCOPY(Atrip::cuda.handle, *n, x, *incx, y, *incy);
#else
  zcopy_(n, x, incx, y, incy);
#endif
}

} // namespace atrip
// Blas:2 ends here
