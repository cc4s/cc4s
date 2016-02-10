/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef GENERATE_RANDOM_MATRIX_DEFINED 
#define GENERATE_RANDOM_MATRIX_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class GenerateRandomMatrix: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(GenerateRandomMatrix);
    GenerateRandomMatrix(
      std::vector<Argument> const &argumentList
    );
    virtual ~GenerateRandomMatrix();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new GenerateRandomMatrix(argumentList);
    }
  };
}

#endif

