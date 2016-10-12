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
     * \brief Slices the particle hole part \f$\Gamma^{aG}_i\f$ from
     * the full Coulomb vertex \f$\Gamma^{qG}_r\f$.
     */
    virtual void run();


    /**
     * \brief Dry run for slicing the Coulomb vertex.
    */
    virtual void dryRun();
  };
}

#endif

