/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DEFINE_HOLES_AND_PARTICLES_DEFINED
#define DEFINE_HOLES_AND_PARTICLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class DefineHolesAndParticles: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DefineHolesAndParticles);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);
  };
}

#endif

