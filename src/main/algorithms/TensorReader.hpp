/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_READER_DEFINED 
#define TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class TensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorReader);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    Ptr<Node> readText(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const std::string &scalarType
    );
    template <typename F, typename TE>
    Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readText(
      const std::string &fileName,
      const std::vector<size_t> &lens
    );
  };
}

#endif

