#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <complex>

/*
#ifdef INTEL_COMPILER
namespace cc4s {
  class complex: public std::complex<double> {
  public:
    void real(double value) {
      this->real(value);
    }
    void imag(double value) {
      this->imag(value);
    }
  };
}
namespace std {
  double real(cc4s::complex c) {
    return std::real(std::complex<double>(c));
  }
  double imag(cc4s::complex c) {
    return std::imag(std::complex<double>(c));
  }
}
#else
#endif
*/
namespace cc4s {
  typedef std::complex<double> complex;

  inline double absSqr(const double x) {
    return x*x;
  }

  inline double absSqr(const complex z) {
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
  class ComplexTraits<double> {
  public:
    static double convert(const complex x) {
      return std::real(x);
    }
  };
}

#endif

