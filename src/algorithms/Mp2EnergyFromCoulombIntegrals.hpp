/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class Mp2EnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2EnergyFromCoulombIntegrals);
    Mp2EnergyFromCoulombIntegrals(std::vector<Argument> const &argumentList);
    virtual ~Mp2EnergyFromCoulombIntegrals();
    virtual void run();
  };
}

#endif

