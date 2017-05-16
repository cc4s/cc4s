/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STATIC_ASSERT_DEFINED
#define STATIC_ASSERT_DEFINED

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
  class TypeRelations<int, double> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <>
  class TypeRelations<int, complex> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <>
  class TypeRelations<double, complex> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
}

#endif

