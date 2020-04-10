/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates MP2 energy from the Coulomb Integrals \f$V_{ij}^{ab}.
   */
  class Mp2EnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2EnergyFromCoulombIntegrals);
    Mp2EnergyFromCoulombIntegrals(std::vector<Argument> const &argumentList);
    virtual ~Mp2EnergyFromCoulombIntegrals();
    /**
     * \brief Calculates MP2 energy from Coulomb integrals Vabij
     */
    virtual void run();
    /**
     * \brief Dry run for the MP2 energy from Coulomb integrals Vabij
     */
    virtual void dryRun();
  protected:
    template <typename F>
    F calculateMp2Energy(CTF::Tensor<F> &Vabij);
  };
}

#endif

