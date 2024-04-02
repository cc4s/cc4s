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

// [[file:~/cuda/atrip/atrip.org::*Blas][Blas:1]]
#pragma once

#include <atrip/Complex.hpp>
#include <atrip/Types.hpp>
#include "config.h"

namespace atrip {

#if !defined(HAVE_CUDA)
extern "C" {
void dgemm_(const char *transa,
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
            const int *ldc);

void zgemm_(const char *transa,
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
            const int *ldc);

void dcopy_(int *n, const double *x, int *incx, double *y, int *incy);

void zcopy_(int *n, const void *x, int *incx, void *y, int *incy);
}
#endif

template <typename F>
void xcopy(int *n,
           const DataFieldType<F> *x,
           int *incx,
           DataFieldType<F> *y,
           int *incy);

template <typename F>
void xgemm(const char *transa,
           const char *transb,
           const int *m,
           const int *n,
           const int *k,
           F *alpha,
           const DataFieldType<F> *A,
           const int *lda,
           const DataFieldType<F> *B,
           const int *ldb,
           F *beta,
           DataFieldType<F> *C,
           const int *ldc);
} // namespace atrip
// Blas:1 ends here
