/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_NORM_DEFINED 
#define TENSOR_NORM_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class TensorNorm: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorNorm);
    TensorNorm(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorNorm();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new TensorNorm(argumentList);
    }
  };
}

#endif

