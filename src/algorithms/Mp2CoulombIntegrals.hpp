/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MP2_COULOMB_INTEGRALS_DEFINED
#define MP2_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class Mp2CoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2CoulombIntegrals);
    Mp2CoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~Mp2CoulombIntegrals();
    virtual void run();

  };
}

#endif

