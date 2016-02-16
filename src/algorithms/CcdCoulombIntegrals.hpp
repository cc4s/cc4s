/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCD_COULOMB_INTEGRALS_DEFINED
#define CCD_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class CcdCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcdCoulombIntegrals);
    CcdCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcdCoulombIntegrals();
    virtual void run();

  };
}

#endif

