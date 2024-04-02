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

#ifndef OPERATIONS_HPP_
#define OPERATIONS_HPP_

#include <atrip/CUDA.hpp>
#include <atrip/Types.hpp>
#include <atrip/Complex.hpp>

namespace atrip {
namespace acc {

// cuda kernels

////
template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F
maybe_conjugate_scalar(const F &a) {
  return a;
}

template <>
__MAYBE_HOST__ __INLINE__ Complex maybe_conjugate_scalar(const Complex &a) {
  return std::conj(a);
}

// TODO: instantiate for std::complex<double>

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ cuDoubleComplex
maybe_conjugate_scalar(const cuDoubleComplex &a) {
  return {a.x, -a.y};
}
#endif /*  defined(HAVE_CUDA) */

template <typename F>
__MAYBE_GLOBAL__ void maybe_conjugate(F *to, F *from, size_t n) {
  for (size_t i = 0; i < n; ++i) { to[i] = maybe_conjugate_scalar<F>(from[i]); }
}

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ void
reorder(F *to, F *from, size_t size, size_t I, size_t J, size_t K) {
  size_t idx = 0;
  const size_t IDX = I + J * size + K * size * size;
  for (size_t k = 0; k < size; k++)
    for (size_t j = 0; j < size; j++)
      for (size_t i = 0; i < size; i++, idx++) to[idx] += from[IDX];
}

// Multiplication operation
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F prod(const F &a, const F &b) {
  return a * b;
}

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ cuDoubleComplex
prod(const cuDoubleComplex &a, const cuDoubleComplex &b) {
  return cuCmul(a, b);
}
#endif /*  defined(HAVE_CUDA) */

// Division operation
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F div(const F &a, const F &b) {
  return a / b;
}

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ cuDoubleComplex
div(const cuDoubleComplex &a, const cuDoubleComplex &b) {
  return cuCdiv(a, b);
}
#endif /*  defined(HAVE_CUDA) */

// Real part
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_HOST__ __INLINE__ double real(F &a) {
  return std::real(a);
}

template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ double real(double &a) {
  return a;
}

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ double real(cuDoubleComplex &a) {
  return cuCreal(a);
}
#endif /*  defined(HAVE_CUDA) */

// Substraction operator
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F sub(const F &a, const F &b) {
  return a - b;
}

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ cuDoubleComplex
sub(const cuDoubleComplex &a, const cuDoubleComplex &b) {
  return cuCsub(a, b);
}
#endif /*  defined(HAVE_CUDA) */

// Addition operator
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F add(const F &a, const F &b) {
  return a + b;
}

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ cuDoubleComplex
add(const cuDoubleComplex &a, const cuDoubleComplex &b) {
  return cuCadd(a, b);
}
#endif /*  defined(HAVE_CUDA) */

// Sum in place operator
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ void sum_in_place(F *to, const F *from) {
  *to += *from;
}

#if defined(HAVE_CUDA)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ void sum_in_place(cuDoubleComplex *to,
                                                  const cuDoubleComplex *from) {
  to->x += from->x;
  to->y += from->y;
}
#endif /*  defined(HAVE_CUDA) */

} // namespace acc
} // namespace atrip

#endif
