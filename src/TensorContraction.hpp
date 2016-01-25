/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_CONTRACTION_DEFINED 
#define TENSOR_CONTRACTION_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class TensorContraction: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorContraction);
    TensorContraction(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorContraction();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new TensorContraction(argumentList);
    }
  };
}

#endif

