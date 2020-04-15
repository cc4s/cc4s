/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_NETWORK_DEFINED 
#define TENSOR_NETWORK_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class TensorNetwork: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorNetwork);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    template <typename TE>
    Real<> getTrace(const Ptr<MapNode> &tensor);
  };
}

#endif

