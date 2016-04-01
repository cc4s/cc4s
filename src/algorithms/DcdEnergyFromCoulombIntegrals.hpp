/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

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

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 16;

  protected:
    void iterateHirata(int i);
    void iterateBartlett();
  };
}

#endif

