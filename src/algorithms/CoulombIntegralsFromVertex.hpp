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
    /**
     * \brief Calculates Coulomb integrals Vabcd,Vabij,Vaibj,Vabci,Vijka,Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices. Arguments can be
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    virtual void run();
    /**
     * \brief Dry run for calculating Coulomb integrals
     * Vabcd,Vabij,Vaibj,Vabci,Vijka Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices. Arguments can be
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    virtual void dryRun();

  };
}

#endif

