#ifndef TENSOR_WRITER_DEFINED 
#define TENSOR_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>

namespace cc4s {
  class TensorWriter: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorWriter)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    void writeData(
      const Ptr<MapNode> &tensor,
      const std::string &fileName,
      const std::string &scalarType,
      const bool binary,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    void writeText(
      const Ptr<MapNode> &tensor,
      const Ptr<Tensor<F,TE>> &tensorData,
      const std::string &fileName,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    void writeBinary(
      const Ptr<MapNode> &tensor,
      const Ptr<Tensor<F,TE>> &tensorData,
      const std::string &fileName,
      const SourceLocation &sourceLocation
    );
  };
}

#endif

