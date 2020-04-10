/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_NETWORK_DEFINED 
#define TENSOR_NETWORK_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class TensorNetwork: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorNetwork);
    TensorNetwork(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorNetwork();
    virtual void run();
    virtual void dryRun();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new TensorNetwork(argumentList);
    }
  };
}

#endif

