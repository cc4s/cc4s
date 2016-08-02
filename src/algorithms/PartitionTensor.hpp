/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTITION_TENSOR_DEFINED 
#define PARTITION_TENSOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class PartitionTensor: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(PartitionTensor);
    PartitionTensor(
      std::vector<Argument> const &argumentList
    );
    virtual ~PartitionTensor();
    virtual void run();
  protected:
    void partitionLastDimension();
  };
}

#endif

