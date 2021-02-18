/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef READ_DEFINED 
#define READ_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>

namespace cc4s {
  class Read: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Read)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  };
}

#endif

