#ifndef WRITE_DEFINED 
#define WRITE_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>

namespace cc4s {
  class Write: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Write)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  };
}

#endif

