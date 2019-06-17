#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <math/Float.hpp>
#include <complex>

namespace cc4s {
  // use standard library complex number support
  template <typename Real>
  using Complex = std::complex<Real>;

  // define explicit size complex types
  typedef Complex<Float32> Complex32;
  typedef Complex<Float64> Complex64;
  typedef Complex<Float128> Complex128;

  // define complex field over machine supported reals as default complex type
  typedef Complex<real> complex;


  template <typename Real>
  inline Real absSqr(const Real x) {
    return x*x;
  }

  template <typename Real>
  inline Real absSqr(const Complex<Real> z) {
    return absSqr(z.real()) + absSqr(z.imag());
  }

  // type info allowing inference of
  // extended type from given optinoally complex type. e.g.:
  // ComplexTraits<complex>::ExtendedType = real
  // ComplexTraits<real>::ExtendedType = real
  // ComplexTraits<Complex<int>>::ExtendedType = int
  template <typename T>
  class ComplexTraits {
  public:
    typedef T ExtendedType;
  };

  template <typename T>
  class ComplexTraits<Complex<T>> {
  public:
    typedef T ExtendedType;
  };

  // numeric conversions
  template <typename Target, typename Source>
  class Conversion;

  template <typename Target, typename Real>
  class Conversion<Target, Complex<Real>> {
  public:
    static Target from(const Complex<Real> x) {
      return Target(x);
    }
  };

  template <typename Real>
  class Conversion<Real, Complex<Real>> {
  public:
    static Real from(const Complex<Real> x) {
      return std::real(x);
    }
  };

#ifdef INTEL_COMPILER
    // TODO: implement for intel
#else
  inline std::ostream &operator <<(
    std::ostream &stream, const Complex128 z
  ) {
    return stream << '(' << std::real(z) << ',' << std::imag(z) << ')';
  }
#endif
}

#endif

