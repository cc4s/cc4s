/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRCCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DRCCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <Algorithm.hpp>
#include <Options.hpp>

namespace cc4s {
  class DrccdEnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DrccdEnergyFromCoulombIntegrals);
    DrccdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~DrccdEnergyFromCoulombIntegrals();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new DrccdEnergyFromCoulombIntegrals(argumentList);
    }

  protected:
    void iterate();
  };
}

#endif

