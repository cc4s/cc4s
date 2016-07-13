/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SLICE_COULOMB_VERTEX_DEFINED
#define SLICE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Slices the particle hole part \f$\Gamma^{aG}_i\f$ from
   * the full Coulomb vertex \f$\Gamma^{qG}_r\f$.
   */
  class SliceCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(SliceCoulombVertex);
    SliceCoulombVertex(
      std::vector<Argument> const &argumentList
    );
    virtual ~SliceCoulombVertex();
    /**
     * \brief Calculates Coulomb integrals Vabcd,Vabij,Vaibj,Vabci,Vijka,Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices. Arguments can be
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    virtual void run();


    /**
     * \brief Dry run for slicing the Coulomb vertex.
    */
    virtual void dryRun();
  };
}

#endif

