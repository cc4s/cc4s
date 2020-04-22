/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SLICE_COULOMB_VERTEX_DEFINED
#define SLICE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class SliceCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(SliceCoulombVertex);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);
  };
}

#endif

