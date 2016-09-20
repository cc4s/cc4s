/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef VIRTUAL_HARTREE_FOCK_ENERGY_MATRIX_FROM_COULOMB_VERTEX_DEFINED
#define VIRTUAL_HARTREE_FOCK_ENERGY_MATRIX_FROM_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class VirtualHartreeFockEnergyMatrixFromCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(VirtualHartreeFockEnergyMatrixFromCoulombVertex);
    VirtualHartreeFockEnergyMatrixFromCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~VirtualHartreeFockEnergyMatrixFromCoulombVertex();
    /**
     * \brief Calculates virtual Hartree-Fock energy from Coulomb vertex
     * \f$\Gamma_{pG}^q\f$
     */
    virtual void run();
    /**
     * \brief Dry run for virtual Hartree-Fock energy from Coulomb vertex
     * \f$\Gamma_{pG}^q\f$
     */
    virtual void dryRun();
  };
}

#endif

