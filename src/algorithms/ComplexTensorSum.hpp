/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COMPLEX_TENSOR_SUM_DEFINED 
#define COMPLEX_TENSOR_SUM_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ComplexTensorSum: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ComplexTensorSum);
    ComplexTensorSum(
      std::vector<Argument> const &argumentList
    );
    virtual ~ComplexTensorSum();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new ComplexTensorSum(argumentList);
    }
  };
}

#endif

