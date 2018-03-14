#ifndef FLOAT_DEFINED
#define FLOAT_DEFINED

#ifndef INTEL_COMPILER
#include <quadmath.h>
#include <ostream>
#endif


// TODO: use configuration for setting default float type sizes in bits
#define DEFAULT_FLOAT_BIT_SIZE 64
#define MACHINE_FLOAT_BIT_SIZE 64

namespace cc4s {
  template <int FloatSize>
  class FloatTypes;

  template <>
  class FloatTypes<32> {
  public:
    typedef float type;
  };

  template <>
  class FloatTypes<64> {
  public:
    typedef double type;
  };

  template <>
  class FloatTypes<128> {
  public:
#ifdef INTEL_COMPILER
    typedef _Quad type;
#else
    typedef __float128 type;
#endif
  };

  // define explicit size float types
  typedef FloatTypes<32>::type Float32;
  typedef FloatTypes<64>::type Float64;
  typedef FloatTypes<128>::type Float128;

  // define machine supported float as real type
  typedef FloatTypes<
    MACHINE_FLOAT_BIT_SIZE < DEFAULT_FLOAT_BIT_SIZE ?
      MACHINE_FLOAT_BIT_SIZE : DEFAULT_FLOAT_BIT_SIZE
  >::type real;

// define stream output for quadruple precision numbers
#ifdef INTEL_COMPILER
    // TODO: implement for intel
#else
  inline std::ostream &operator <<(
    std::ostream &stream, const Float128 x
  ) {
    char buffer[1024];
    quadmath_snprintf(buffer, sizeof(buffer), "%*.36Qe", x);
    return stream << buffer;
  }
#endif
}

#endif

