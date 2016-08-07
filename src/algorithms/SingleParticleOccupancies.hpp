/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SINGLE_PARTICLE_OCCUPANCIES_DEFINED
#define SINGLE_PARTICLE_OCCUPANCIES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Evaluates
   * \f$\langle\Psi|\hat N_p|\Psi\rangle>/\langle\Psi|\Psi\rangle\f$
   * given the DoublesAmplitudes from a linearized coupled cluster theory.
   */
  class SingleParticleOccupancies: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(SingleParticleOccupancies);
    SingleParticleOccupancies(
      std::vector<Argument> const &argumentList
    );
    virtual ~SingleParticleOccupancies();
    virtual void run();
    virtual void dryRun();
  };
}

#endif

