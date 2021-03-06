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

#ifndef MATH_FUNCTIONS_DEFINED
#define MATH_FUNCTIONS_DEFINED

#include <Real.hpp>
#include <Complex.hpp>

#include <cmath>
#include <string>
#include <ctf.hpp>

namespace cc4s {
  // Provides template typed wrappers for std:: math constants and functions,
  // allowing unified usage as strongly typed function objects, also
  // inside class or function templates.
  // non-standard 128 bit reals are also wrapped where available.

  inline __complex128 toComplex128(const Complex<128> &c) {
    return *reinterpret_cast<const __complex128 *>(&c);
  }

  inline Complex<128> fromComplex128(const __complex128 &c) {
    return *reinterpret_cast<const Complex<128> *>(&c);
  }

  // constants
  template <typename F=Real<>>
  constexpr F Pi() {
    return acos(F(-1));
  }
  template <typename F=Real<>>
  constexpr F Tau() {
    return 2 * Pi<F>();
  }

  // functions
  template <typename F>
  inline typename ComplexTraits<F>::RealType abs(const F x) {
    return std::abs(x);
  }
  template <>
  inline Real<128> abs(const Real<128> x) {
    return fabsq(x);
  }
  template <>
  inline Real<128> abs(const Complex<128> c) {
    return cabsq(toComplex128(c));
  }

  // real part, same if real
  template <typename F>
  inline typename ComplexTraits<F>::RealType real(const F c) {
    return std::real(c);
  }
  template <>
  inline Real<128> real(const Real<128> c) {
    return c;
  }
  template <>
  inline Real<128> real(const Complex<128> c) {
    return crealq(toComplex128(c));
  }

  // imaginary part, 0 if real
  template <typename F>
  inline typename ComplexTraits<F>::RealType imag(const F c) {
    return std::imag(c);
  }
  template <>
  inline Real<128> imag(const Real<128> c) {
    return 0;
  }
  template <>
  inline Real<128> imag(const Complex<128> c) {
    return cimagq(toComplex128(c));
  }

  // complex argument, 0 if real
  template <typename F>
  inline typename ComplexTraits<F>::RealType arg(const F c) {
    return std::arg(c);
  }
  template <>
  inline Real<128> arg(const Real<128> c) {
    return 0;
  }
  template <>
  inline Real<128> arg(const Complex<128> c) {
    return cargq(toComplex128(c));
  }

  // complex conjugate, identity function if real type
  template <typename F>
  inline F conj(const F x) {
    return std::conj(x);
  }
  template <>
  inline Real<64> conj(const Real<64> x) {
    return x;
  }
  template <>
  inline Real<128> conj(const Real<128> x) {
    return x;
  }
  template <>
  inline Complex<128> conj(const Complex<128> x) {
    return fromComplex128(conjq(toComplex128(x)));
  }

  // functions defined on real and complex domain
  template <typename F>
  inline F pow(const F x, const F e) {
    return std::pow(x,e);
  }
  template <>
  inline Real<128> pow(const Real<128> x, const Real<128> e) {
    return powq(x,e);
  }
  template <>
  inline Complex<128> pow(const Complex<128> x, const Complex<128> e) {
    return fromComplex128(cpowq(toComplex128(x),toComplex128(e)));
  }

#define COMPLEX_AND_REAL_FUNCTION_DEFINITION(NAME) \
  template <typename F> \
  inline F NAME(const F x) { \
    return std::NAME(x); \
  } \
  template <> \
  inline Real<128> NAME(const Real<128> x) { \
    return NAME##q(x); \
  } \
  template <> \
  inline Complex<128> NAME(const Complex<128> x) { \
    return fromComplex128(c##NAME##q(toComplex128(x))); \
  }
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(sqrt)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(exp)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(log)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(cos)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(sin)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(tan)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(acos)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(asin)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(atan)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(cosh)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(sinh)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(tanh)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(acosh)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(asinh)
  COMPLEX_AND_REAL_FUNCTION_DEFINITION(atanh)

  // functions on the real domain only
  template <typename F>
  inline F frexp(const F x, int &e) {
    return std::frexp(x,&e);
  }
  template <>
  inline Real<128> frexp(const Real<128> x, int &e) {
    return frexpq(x,&e);
  }

  template <typename F>
  inline F ldexp(const F x, const int e) {
    return std::ldexp(x,e);
  }
  template <>
  inline Real<128> ldexp(const Real<128> x, const int e) {
    return ldexpq(x,e);
  }

  template <typename F>
  inline F fmod(const F x, const F y) {
    return std::fmod(x,y);
  }
  template <>
  inline Real<128> fmod(const Real<128> x, const Real<128> y) {
    return fmodq(x,y);
  }

#define REAL_FUNCTION_DEFINITION(NAME) \
  template <typename F> \
  inline F NAME(const F x) { \
    return std::NAME(x); \
  } \
  template <> \
  inline Real<128> NAME(const Real<128> x) { \
    return NAME##q(x); \
  }
  REAL_FUNCTION_DEFINITION(cbrt)
  REAL_FUNCTION_DEFINITION(erf)
  REAL_FUNCTION_DEFINITION(erfc)
  REAL_FUNCTION_DEFINITION(lgamma)
  REAL_FUNCTION_DEFINITION(tgamma)
  REAL_FUNCTION_DEFINITION(floor)
  REAL_FUNCTION_DEFINITION(ceil)
  REAL_FUNCTION_DEFINITION(trunc)
  REAL_FUNCTION_DEFINITION(round)
}

#endif

