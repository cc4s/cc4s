/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED
#define PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates the Coulomb Integrals \f$V_{ij}^{ab} from the Particle Hole
   * Coulomb Vertex \f$\Gamma_{iG}^a\f$ and stores them in a CTF Tensor Vabij
   * The argument of the integrals is PPHHCoulombIntegrals.
   */
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

    /** \brief The Coulomb Vertex GammaGai  */
    CTF::Tensor<complex> *GammaGai;
    /** \brief The Coulomb integrals Vabij  */
    CTF::Tensor<> *Vabij;

    /**
     * \brief Dry run for calculating Coulomb integrals from ParticleHole Coulomb Vertex
     * GammaGai.
     */
    virtual void dryRun();
  };
}

#endif

