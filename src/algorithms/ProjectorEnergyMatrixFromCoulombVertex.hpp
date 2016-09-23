/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PROJECTOR_ENERGY_MATRIX_FROM_COULOMB_VERTEX_DEFINED
#define PROJECTOR_ENERGY_MATRIX_FROM_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ProjectorEnergyMatrixFromCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ProjectorEnergyMatrixFromCoulombVertex);
    ProjectorEnergyMatrixFromCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~ProjectorEnergyMatrixFromCoulombVertex();
    /**
     * \brief Calculates projector energy matrix from Coulomb vertex
     * \f$\Gamma_{pG}^q\f$
     */
    virtual void run();
    /**
     * \brief Dry run for projector energy matrix from Coulomb vertex
     * \f$\Gamma_{pG}^q\f$
     */
    virtual void dryRun();
  };
}

#endif

