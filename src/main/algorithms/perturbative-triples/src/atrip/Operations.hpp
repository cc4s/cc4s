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

#include <atrip/Acc.hpp>
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

// template <>
// #if defined(HAVE_HIP)
// // HIP Complains if the attributes are different in the template
// // declaration
// __MAYBE_DEVICE__
// #endif
//     __MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ Complex
//     maybe_conjugate_scalar(const Complex &a) {
//   return std::conj(a);
// }

// TODO: instantiate for std::complex<double>

// template <>
// __MAYBE_DEVICE__ __INLINE__ Complex maybe_conjugate_scalar(const Complex &a)
// {
//   return std::conj(a);
// }

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ ACC_DOUBLE_COMPLEX
maybe_conjugate_scalar(const ACC_DOUBLE_COMPLEX &a) {
  return {a.x, -a.y};
}
#endif /*  defined(HAVE_ACC) */

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

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ ACC_DOUBLE_COMPLEX
prod(const ACC_DOUBLE_COMPLEX &a, const ACC_DOUBLE_COMPLEX &b) {
#  if defined(HAVE_HIP) && defined(BLAHBLAH)
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
#  else
  return ACC_COMPLEX_MUL(a, b);
#  endif /* defined(HAVE_HIP) */
}
#endif /*  defined(HAVE_ACC) */

// Division operation
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F div(const F &a, const F &b) {
  return a / b;
}

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ ACC_DOUBLE_COMPLEX
div(const ACC_DOUBLE_COMPLEX &a, const ACC_DOUBLE_COMPLEX &b) {
  return ACC_COMPLEX_DIV(a, b);
}
#endif /*  defined(HAVE_ACC) */

// Real part
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ double real(F &a) {
  return std::real(a);
}

template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ double real(double &a) {
  return a;
}

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ double real(ACC_DOUBLE_COMPLEX &a) {
  return ACC_COMPLEX_REAL(a);
}
#endif /*  defined(HAVE_ACC) */

// Substraction operator
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F sub(const F &a, const F &b) {
  return a - b;
}

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ ACC_DOUBLE_COMPLEX
sub(const ACC_DOUBLE_COMPLEX &a, const ACC_DOUBLE_COMPLEX &b) {
  return ACC_COMPLEX_SUB(a, b);
}
#endif /*  defined(HAVE_ACC) */

// Addition operator
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ F add(const F &a, const F &b) {
  return a + b;
}

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ __INLINE__ ACC_DOUBLE_COMPLEX
add(const ACC_DOUBLE_COMPLEX &a, const ACC_DOUBLE_COMPLEX &b) {
  return ACC_COMPLEX_ADD(a, b);
}
#endif /*  defined(HAVE_ACC) */

// Sum in place operator
//////////////////////////////////////////////////////////////////////////////

template <typename F>
__MAYBE_DEVICE__ __MAYBE_HOST__ void sum_in_place(F *to, const F *from) {
  *to += *from;
}

#if defined(HAVE_ACC)
template <>
__MAYBE_DEVICE__ __MAYBE_HOST__ void
sum_in_place(ACC_DOUBLE_COMPLEX *to, const ACC_DOUBLE_COMPLEX *from) {
  to->x += from->x;
  to->y += from->y;
}
#endif /*  defined(HAVE_ACC) */

} // namespace acc
} // namespace atrip

#endif
