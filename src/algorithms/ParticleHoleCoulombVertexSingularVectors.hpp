/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED
#define PARTICLE_HOLE_COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ParticleHoleCoulombVertexSingularVectors: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ParticleHoleCoulombVertexSingularVectors);
    ParticleHoleCoulombVertexSingularVectors(std::vector<Argument> const &argumentList);
    virtual ~ParticleHoleCoulombVertexSingularVectors();
    /**
     * \brief Calculates left singular vectors of the particle-hole Coulomb vertex
     * \f$\tilde\Gamma^a_{iG}\f$.
     */
    virtual void run();
    /**
     * \brief Dry run for calculating the left singular vectors of
     * \f$\Gamma^a_{iG}\f$.
     */
    virtual void dryRun();

    static double constexpr DEFAULT_REDUCTION = 0.5;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;
  };
}

#endif

