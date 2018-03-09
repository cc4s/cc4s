#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <complex>
#ifndef INTEL_COMPILER
#include <quadmath.h>
#endif

// TODO: move to own header file
// TODO: use configuration for setting default real type sizes in bits
#define DEFAULT_REAL_BIT_SIZE 64

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
  template <int FloatSize=DEFAULT_REAL_BIT_SIZE>
  class FloatTypes;

  template <>
  class FloatTypes<32> {
  public:
    typedef float real;
    typedef std::complex<float> complex;
  };

  template <>
  class FloatTypes<64> {
  public:
    typedef double real;
    typedef std::complex<double> complex;
  };

  template <>
  class FloatTypes<128> {
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
  typedef FloatTypes<>::real real;
  typedef FloatTypes<>::complex complex;

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

// define stream output for quadruple precision numbers
#ifdef INTEL_COMPILER
  // TODO: implement for intel
#else
inline std::ostream &operator <<(
  std::ostream &stream, const cc4s::FloatTypes<128>::real x
) {
  char buffer[1024];
  quadmath_snprintf(buffer, sizeof(buffer), "%*.36Qe", x);
  return stream << buffer;
}
#endif

#endif

