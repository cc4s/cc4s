/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define THERMAL_MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ThermalMp2EnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ThermalMp2EnergyFromCoulombIntegrals);
    ThermalMp2EnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalMp2EnergyFromCoulombIntegrals();
    /**
     * \brief Calculates the thermal MP2 energy from Coulomb integrals Vabij
     */
    virtual void run();
    /**
     * \brief Dry run for the thermal MP2 energy from Coulomb integrals Vabij
     */
    virtual void dryRun();
  };
}

#endif

