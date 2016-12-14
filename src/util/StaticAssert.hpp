/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STATIC_ASSERT_DEFINED
#define STATIC_ASSERT_DEFINED

#include <math/Complex.hpp>

namespace cc4s {
  template <typename T>
  class StaticAssert {
  public:
    enum { False = false };
  };

  template <typename A, typename B>
  class TypeRelations {
  public:
    enum { Equals = false, PointerTo = false, CompatibleTo = false };
  };
  template <typename A>
  class TypeRelations<A,A> {
  public:
    enum { Equals = true, PointerTo = false, CompatibleTo = true };
  };
  template <typename A>
  class TypeRelations<A *, A> {
  public:
    enum { Equals = false, PointerTo = true, CompatibleTo = false };
  };

  template <>
  class TypeRelations<int, double> {
  public:
    enum { Equals = false, PointerTo = false, CompatibleTo = true };
  };
  template <>
  class TypeRelations<int, complex> {
  public:
    enum { Equals = false, PointerTo = false, CompatibleTo = true };
  };
  template <>
  class TypeRelations<double, complex> {
  public:
    enum { Equals = false, PointerTo = false, CompatibleTo = true };
  };
}

#endif

