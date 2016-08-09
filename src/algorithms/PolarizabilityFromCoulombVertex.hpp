/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef POLARIZABILITY_FROM_COULOMB_VERTEX_DEFINED
#define POLARIZABILITY_FROM_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class PolarizabilityFromCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(PolarizabilityFromCoulombVertex);
    PolarizabilityFromCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~PolarizabilityFromCoulombVertex();
    /**
     * \brief Calculates the independent particle polarizability
     * \f${X_0}({\rm i}\nu)_G^H\f$ from the Coulomb vertex
     * \f$\Gamma^{qG}_r\f$.
     * At the moment, only the static contribution at \f${\rm i}\nu=0\f$
     * is calculated and parts of the Coulomb kernel \f$v^{1/2}_G\f$ and
     * \f$v^{1/2}^H\f$ are also still included.
     */
    virtual void run();
    /**
     * \brief Dry run for calculating the polarizability.
     */
    virtual void dryRun();
  };
}

#endif

