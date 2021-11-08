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

#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <math/Real.hpp>
#include <complex>

namespace cc4s {
  // use standard library complex number support
  template <int FloatSize=DEFAULT_FLOAT_BIT_SIZE>
  using Complex = std::complex<Real<FloatSize>>;

  template <typename F>
  inline F absSqr(const F x) {
    return x*x;
  }

  template <int FloatSize>
  inline Real<FloatSize> absSqr(const Complex<FloatSize> z) {
    return absSqr(z.real()) + absSqr(z.imag());
  }

  // type info allowing inference of
  // real type from given optinoally complex type. e.g.:
  // ComplexTraits<Complex<32>>::RealType = Real<32>
  // ComplexTraits<Real<64>>::RealType = Real<64>
  template <typename F>
  class ComplexTraits;

  template <>
  class ComplexTraits<Real<32>> {
  public:
    typedef Real<32> RealType;
  };
  template <>
  class ComplexTraits<Real<64>> {
  public:
    typedef Real<64> RealType;
  };
  template <>
  class ComplexTraits<Real<128>> {
  public:
    typedef Real<128> RealType;
  };

  template <>
  class ComplexTraits<Complex<32>> {
  public:
    typedef Real<32> RealType;
  };
  template <>
  class ComplexTraits<Complex<64>> {
  public:
    typedef Real<64> RealType;
  };
  template <>
  class ComplexTraits<Complex<128>> {
  public:
    typedef Real<128> RealType;
  };

  // narrowing numeric conversions
  // TODO: implement usage properly
/*
  template <typename G, typename F>
  class NarrowingConversion {
  public:
    static G from(const G x) {
      return G(x);
    }
  };

  template <>
  class NarrowingConversion<Real<32>, Complex<32>> {
  public:
    static Real<32> from(const Complex<32> x) {
      return real(x);
    }
  };
  template <>
  class NarrowingConversion<Real<64>, Complex<64>> {
  public:
    static Real<64> from(const Complex<64> x) {
      return real(x);
    }
  };
  template <>
  class NarrowingConversion<Real<128>, Complex<128>> {
  public:
    static Real<128> from(const Complex<128> x) {
      return real(x);
    }
  };
*/

  inline std::basic_ostream<char> &operator <<(
    std::basic_ostream<char> &stream, const cc4s::Complex<128> x
  ) {
    return stream << '(' << x.real() << ',' << x.imag() << ')';
  }
}

#endif

