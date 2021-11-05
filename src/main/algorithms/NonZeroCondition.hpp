#ifndef NON_ZERO_CONDITION_DEFINED 
#define NON_ZERO_CONDITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>
#include <tcc/Tensor.hpp>

namespace cc4s {
  class NonZeroCondition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(NonZeroCondition)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    template <typename F, typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);
  };
}

#endif

