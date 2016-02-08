/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define CCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <Algorithm.hpp>
#include <Options.hpp>

namespace cc4s {
  class CcdEnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcdEnergyFromCoulombIntegrals);
    CcdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcdEnergyFromCoulombIntegrals();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new CcdEnergyFromCoulombIntegrals(argumentList);
    }

  protected:
    void iterateHirata();
    void iterateBartlett();
  };
}

#endif

