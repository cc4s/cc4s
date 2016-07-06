/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MP2_ENERGY_MATRIX_FROM_COULOMB_INTEGRALS_DEFINED 
#define MP2_ENERGY_MATRIX_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class Mp2EnergyMatrixFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2EnergyMatrixFromCoulombIntegrals);
    Mp2EnergyMatrixFromCoulombIntegrals(std::vector<Argument> const &argumentList);
    virtual ~Mp2EnergyMatrixFromCoulombIntegrals();
    /**
     * \brief Calculates MP2 energy from Coulomb integrals Vabij
     */
    virtual void run();
    /**
     * \brief Dry run for the MP2 energy from Coulomb integrals Vabij
     */
    virtual void dryRun();
  };
}

#endif

