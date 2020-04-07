#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <math/Float.hpp>
#include <complex>

namespace cc4s {
  // use standard library complex number support
  template <typename F=Real<>>
  using Complex = std::complex<F>;

  // define explicit size complex types
  // DEPRECATED: use Complex<Real<32>> instead
  typedef Complex<Float32> Complex32;
  typedef Complex<Float64> Complex64;
  typedef Complex<Float128> Complex128;

  // DEPRECATED:
  typedef Complex<> complex;

  template <typename F>
  inline F absSqr(const F x) {
    return x*x;
  }

  template <typename F>
  inline F absSqr(const Complex<F> z) {
    return absSqr(z.real()) + absSqr(z.imag());
  }

  // type info allowing inference of
  // extended type from given optinoally complex type. e.g.:
  // ComplexTraits<Complex<>>::ExtendedType = Real<>
  // ComplexTraits<Real<>>::ExtendedType = Real<>
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
  template <typename G, typename F>
  class Conversion;

  template <typename G, typename F>
  class Conversion<G, Complex<F>> {
  public:
    static G from(const Complex<F> x) {
      return G(x);
    }
  };

  template <typename F>
  class Conversion<F, Complex<F>> {
  public:
    static F from(const Complex<F> x) {
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

