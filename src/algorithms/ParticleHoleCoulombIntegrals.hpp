/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED
#define PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ParticleHoleCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ParticleHoleCoulombIntegrals);
    ParticleHoleCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~ParticleHoleCoulombIntegrals();
    /**
     * \brief Calculates Coulomb integrals from ParticleHole Coulomb Vertex
     * GammaGai.
     */
    virtual void run();
  };
}

#endif

