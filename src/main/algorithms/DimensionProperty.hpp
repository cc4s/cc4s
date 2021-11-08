#ifndef DIMENSION_PROPERTY_DEFINED
#define DIMENSION_PROPERTY_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>
#include <tcc/Tensor.hpp>

namespace cc4s {
  class DimensionProperty: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DimensionProperty)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    template <typename F, typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);
  };
}

#endif

