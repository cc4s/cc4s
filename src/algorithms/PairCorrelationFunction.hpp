/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PAIR_CORRELATION_FUNCTION_DEFINED
#define PAIR_CORRELATION_FUNCTION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <vector>

namespace cc4s {
  class PairCorrelationFunction: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(PairCorrelationFunction);
    PairCorrelationFunction(
      std::vector<Argument> const &argumentList
    );
    virtual ~PairCorrelationFunction();
    virtual void run();

  protected:
    int NG;
    void GtoRFourier();
  };
}


#endif
