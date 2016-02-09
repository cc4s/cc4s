/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <Algorithm.hpp>
#include <Options.hpp>

namespace cc4s {
  class DcdEnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DcdEnergyFromCoulombIntegrals);
    DcdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~DcdEnergyFromCoulombIntegrals();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new DcdEnergyFromCoulombIntegrals(argumentList);
    }

  protected:
    void iterateHirata();
    void iterateBartlett();
  };
}

#endif

