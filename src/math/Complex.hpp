#ifndef COMPLEX_DEFINED
#define COMPLEX_DEFINED

#include <complex>
#ifndef INTEL_COMPILER
#include <quadmath.h>
#endif

// TODO: use configuration for setting default float type sizes
/**
 * \brief Default size of IEEE floating pointer number in bytes.
 **/
#define DEFAULT_FLOAT_SIZE 8

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
  template <int FloatSize=DEFAULT_FLOAT_SIZE>
  class FloatTypes;

  template <>
  class FloatTypes<4> {
  public:
    typedef float real;
    typedef std::complex<float> complex;
  };

  template <>
  class FloatTypes<8> {
  public:
    typedef double real;
    typedef std::complex<double> complex;
  };

  template <>
  class FloatTypes<16> {
  public:
#ifdef INTEL_COMPILER
    typedef _Quad real;
    typedef std::complex<_Quad> complex;
#else
    typedef __float128 real;
    typedef std::complex<__float128> complex;
#endif
  };

  // define types of default Float size:
  typedef FloatTypes<>::real real;
  typedef FloatTypes<>::complex complex;

  // default cast template
  template <typename TargetType>
  class Cast {
  public:
    template <typename SourceType>
    static TargetType from(const SourceType x) { return TargetType(x); }
  };

  #define CAST_DEFINITION(SIZE) \
  template <> \
  class Cast<typename FloatTypes<SIZE>::real> { \
  public: \
    static FloatTypes<SIZE>::real from( \
      const real x \
    ) { \
      return x; \
    } \
    static FloatTypes<SIZE>::real from( \
      const FloatTypes<SIZE>::complex x \
    ) { \
      return std::real(x); \
    } \
  };

  CAST_DEFINITION(4)
  CAST_DEFINITION(8)
  CAST_DEFINITION(16)

  inline real absSqr(const real x) {
    return x*x;
  }

  inline real absSqr(const complex z) {
    return absSqr(z.real()) + absSqr(z.imag());
  }
}

#endif

