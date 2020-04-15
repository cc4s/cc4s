/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_READER_DEFINED 
#define TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class TensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorReader);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  };
}

#endif

