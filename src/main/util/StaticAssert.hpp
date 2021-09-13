#ifndef STATIC_ASSERT_DEFINED
#define STATIC_ASSERT_DEFINED

#include <math/Real.hpp>
#include <math/Complex.hpp>

namespace cc4s {
  template <typename T>
  class StaticAssert {
  public:
    enum { FALSE = false };
  };

  template <typename A, typename B>
  class TypeRelations {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = false };
  };
  template <typename A>
  class TypeRelations<A,A> {
  public:
    enum { EQUALS = true, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <typename A>
  class TypeRelations<A *, A> {
  public:
    enum { EQUALS = false, POINTER_TO = true, CASTABLE_TO = false };
  };

  template <>
  class TypeRelations<int, Real<64>> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <>
  class TypeRelations<int, Complex<64>> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <>
  class TypeRelations<Real<64>, Complex<64>> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };

  // TODO: 128 bit real and complex
}

#endif

