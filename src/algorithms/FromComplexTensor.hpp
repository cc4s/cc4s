/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FROM_COMPLEX_TENSOR_DEFINED 
#define FROM_COMPLEX_TENSOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class FromComplexTensor: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(FromComplexTensor);
    FromComplexTensor(
      std::vector<Argument> const &argumentList
    );
    virtual ~FromComplexTensor();
    virtual void run();
  };
}

#endif

