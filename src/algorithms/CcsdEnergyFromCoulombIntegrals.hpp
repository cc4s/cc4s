/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define CCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <Options.hpp>

namespace cc4s {
  class CcsdEnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEnergyFromCoulombIntegrals);
    CcsdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdEnergyFromCoulombIntegrals();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new CcsdEnergyFromCoulombIntegrals(argumentList);
    }

  protected:
    void iterate();
  };
}

#endif

