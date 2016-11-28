/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STATIC_ASSERT_DEFINED
#define STATIC_ASSERT_DEFINED

namespace cc4s {
  template <typename T>
  class StaticAssert {
  public:
    enum { False = false };
  };

  template <typename A, typename B>
  class TypeRelations {
  public:
    enum { Equals = false, PointerTo = false };
  };
  template <typename A>
  class TypeRelations<A,A> {
  public:
    enum { Equals = true, PointerTo = false };
  };
  template <typename A>
  class TypeRelations<A *, A> {
  public:
    enum { Equals = false, PointerTo = true };
  };
}

#endif

