/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class DcsdEnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DcsdEnergyFromCoulombIntegrals);
    DcsdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~DcsdEnergyFromCoulombIntegrals();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new DcsdEnergyFromCoulombIntegrals(argumentList);
    }

  protected:
    void iterate();
  };
}

#endif

