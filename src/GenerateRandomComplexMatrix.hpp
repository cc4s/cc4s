/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef GENERATE_RANDOM_COMPLEX_MATRIX_DEFINED 
#define GENERATE_RANDOM_COMPLEX_MATRIX_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class GenerateRandomComplexMatrix: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(GenerateRandomComplexMatrix);
    GenerateRandomComplexMatrix(
      std::vector<Argument> const &argumentList
    );
    virtual ~GenerateRandomComplexMatrix();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new GenerateRandomComplexMatrix(argumentList);
    }
  };
}

#endif

