/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_FROM_VERTEX_DEFINED
#define COULOMB_INTEGRALS_FROM_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class CoulombIntegralsFromVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombIntegralsFromVertex);
    CoulombIntegralsFromVertex(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombIntegralsFromVertex();
    virtual void run();
    virtual void dryRun();

  };
}

#endif

