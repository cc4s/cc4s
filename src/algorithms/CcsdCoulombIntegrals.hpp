/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_COULOMB_INTEGRALS_DEFINED
#define CCSD_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class CcsdCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdCoulombIntegrals);
    CcsdCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdCoulombIntegrals();
    virtual void run();

  };
}

#endif

