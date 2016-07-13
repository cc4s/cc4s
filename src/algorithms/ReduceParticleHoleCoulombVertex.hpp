/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REDUCE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED
#define REDUCE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  class ReduceParticleHoleCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ReduceParticleHoleCoulombVertex);
    ReduceParticleHoleCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~ReduceParticleHoleCoulombVertex();
    /**
     * \brief Reduces the particle hole Coulomb vertex by applying the
     * energy matrix reduction transform:
     * \f$\Gamma^{ag}_i = \Gamma^{aG}_i U_G^g\f$.
     */
    virtual void run();
    /**
     * \brief Dry run of reducing the Coulomb vertex.
     */
    virtual void dryRun();
  };
}

#endif

