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

typedef std::complex<double> complex;

#endif

