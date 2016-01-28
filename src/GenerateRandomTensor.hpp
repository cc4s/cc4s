/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef GENERATE_RANDOM_TENSOR_DEFINED 
#define GENERATE_RANDOM_TENSOR_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class GenerateRandomTensor: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(GenerateRandomTensor);
    GenerateRandomTensor(
      std::vector<Argument> const &argumentList
    );
    virtual ~GenerateRandomTensor();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new GenerateRandomTensor(argumentList);
    }
  };
}

#endif

