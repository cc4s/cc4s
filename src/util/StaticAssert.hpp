/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STATIC_ASSERT_DEFINED
#define STATIC_ASSERT_DEFINED

namespace cc4s {
  template <typename T>
  class StaticAssert {
  public:
    enum { False = false };
  };
}

#endif

