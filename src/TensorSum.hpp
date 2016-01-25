/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_SUM_DEFINED 
#define TENSOR_SUM_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class TensorSum: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorSum);
    TensorSum(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorSum();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new TensorSum(argumentList);
    }
  };
}

#endif

