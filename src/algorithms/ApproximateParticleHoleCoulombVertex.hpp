/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef APPROXIMATE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED
#define APPROXIMATE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  /* @WIKI
   * \brief Approximates the particle-hole Coulomb vertex $\tilde\Gamma^a_{iG}$
   * using the given set of left singular vectors $U^F_G$ associated to the
   * $\Gamma^a_{iF} = {U^\ast}^G_F \tilde\Gamma^a_{iG}$.
   * largest singular values, The approximated particle-hole Coulomb vertex is given by
   * \param[in] ParticleHoleCoulombVertex (complex tensor) (none)
   *   The particle-hole Coulomb vertex $\tilde\Gamma$ to approximate.
   * \param[in] ParticleHoleCoulomgVertexSingularVectors tensor none
   *   The left singular vectors $U$ to be used in the transformation.
   * \param[out] ApproximatedParticleHoleCoulombVertex tensor
   *   The approximated particle-hole Coulomb vertex $\Gamma$.
   */
  class ApproximateParticleHoleCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ApproximateParticleHoleCoulombVertex);
    ApproximateParticleHoleCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~ApproximateParticleHoleCoulombVertex();
    /**
     * \brief Approximates the particle-hole Coulomb vertex using the left
     * singular vectors:
     * \f$\Gamma^a_{iF} = {U^\ast}^G_F \tile\Gamma^a_{iG}\f$.
     */
    virtual void run();
    /**
     * \brief Dry run of approximating the particle-hole Coulomb vertex.
     */
    virtual void dryRun();
  };
}

#endif

