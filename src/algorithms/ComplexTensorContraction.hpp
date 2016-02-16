/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COMPLEX_TENSOR_CONTRACTION_DEFINED 
#define COMPLEX_TENSOR_CONTRACTION_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ComplexTensorContraction: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ComplexTensorContraction);
    ComplexTensorContraction(
      std::vector<Argument> const &argumentList
    );
    virtual ~ComplexTensorContraction();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new ComplexTensorContraction(argumentList);
    }
  };
}

#endif

