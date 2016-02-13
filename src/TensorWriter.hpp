/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_WRITER_DEFINED 
#define TENSOR_WRITER_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class TensorWriter: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorWriter);
    TensorWriter(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorWriter();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new TensorWriter(argumentList);
    }
  };
}

#endif

