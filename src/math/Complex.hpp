#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <complex>
#ifndef INTEL_COMPILER
#include <quadmath.h>
#endif

// TODO: use configuration for setting default real type sizes
#define DEFAULT_REAL_SIZE 8

/*
#ifdef INTEL_COMPILER
namespace cc4s {
  class complex: public std::complex<real> {
  public:
    void real(real value) {
      this->real(value);
    }
    void imag(real value) {
      this->imag(value);
    }
  };
}
namespace std {
  real real(cc4s::complex c) {
    return std::real(std::complex<real>(c));
  }
  real imag(cc4s::complex c) {
    return std::imag(std::complex<real>(c));
  }
}
#else
#endif
*/
namespace cc4s {
  template <int RealSize=DEFAULT_REAL_SIZE>
  class RealTypes;

  template <>
  class RealTypes<4> {
  public:
    typedef float real;
    typedef std::complex<float> complex;
  };

  template <>
  class RealTypes<8> {
  public:
    typedef double real;
    typedef std::complex<double> complex;
  };

  template <>
  class RealTypes<16> {
  public:
#ifdef INTEL_COMPILER
    typedef _Quad real;
    typedef std::complex<_Quad> complex;
#else
    typedef __float128 real;
    typedef std::complex<__float128> complex;
#endif
  };

  // define types of default real size:
  typedef RealTypes<>::real real;
  typedef RealTypes<>::complex complex;

  inline real absSqr(const real x) {
    return x*x;
  }

  inline real absSqr(const complex z) {
    return absSqr(z.real()) + absSqr(z.imag());
  }

  // base template
  template <typename F>
  class ComplexTraits {
  public:
    static F convert(const complex x) {
      return F(x);
    }
  };

  template <>
  class ComplexTraits<real> {
  public:
    static real convert(const complex x) {
      return std::real(x);
    }
  };
}

#endif

