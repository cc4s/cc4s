#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <math/Real.hpp>
#include <complex>

namespace cc4s {
  // use standard library complex number support
  template <int FloatSize=64>
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
  // extended type from given optinoally complex type. e.g.:
  // ComplexTraits<Complex<>>::ExtendedType = Real<>
  // ComplexTraits<Real<>>::ExtendedType = Real<>
  template <typename T>
  class ComplexTraits {
  public:
    typedef T ExtendedType;
  };

  template <>
  class ComplexTraits<Complex<32>> {
  public:
    typedef Real<32> ExtendedType;
  };
  template <>
  class ComplexTraits<Complex<64>> {
  public:
    typedef Real<64> ExtendedType;
  };
  template <>
  class ComplexTraits<Complex<128>> {
  public:
    typedef Real<128> ExtendedType;
  };

  // narrowing numeric conversions
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
      return std::real(x);
    }
  };
  template <>
  class NarrowingConversion<Real<64>, Complex<64>> {
  public:
    static Real<64> from(const Complex<64> x) {
      return std::real(x);
    }
  };
  template <>
  class NarrowingConversion<Real<128>, Complex<128>> {
  public:
    static Real<128> from(const Complex<128> x) {
      return std::real(x);
    }
  };

#ifdef INTEL_COMPILER
    // TODO: implement for intel
#else
  inline std::ostream &operator <<(
    std::ostream &stream, const Complex<128> z
  ) {
    return stream << '(' << std::real(z) << ',' << std::imag(z) << ')';
  }
#endif
}

#endif

