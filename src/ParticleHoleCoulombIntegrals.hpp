/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED
#define PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class ParticleHoleCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ParticleHoleCoulombIntegrals);
    ParticleHoleCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~ParticleHoleCoulombIntegrals();
    virtual void run();
  };
}

#endif

