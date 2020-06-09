/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_READER_DEFINED 
#define TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>

namespace cc4s {
  class TensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorReader);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    Ptr<Node> readData(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const std::string &scalarType,
      const bool binary,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readText(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readBinary(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const SourceLocation &sourceLocation
    );
  };
}

#endif

