/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CTF_TEST_DEFINED
#define CTF_TEST_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class CtfTest: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CtfTest);
    CtfTest(
      std::vector<Argument> const &argumentList
    );
    virtual ~CtfTest();
    virtual void run();
  };
}

#endif

