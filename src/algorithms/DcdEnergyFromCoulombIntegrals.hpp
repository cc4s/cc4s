/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/LegacyClusterDoublesAlgorithm.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the Dcd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
   * V_{kl}^{ij}\f$, and \f$V_{cd}^{ab}\f$ (if given, else slicing and the Coulomb
   * Vertex \f$\Gamma_{pG}^q\f$  is used).
   */
  class DcdEnergyFromCoulombIntegrals: public LegacyClusterDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DcdEnergyFromCoulombIntegrals);
    DcdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~DcdEnergyFromCoulombIntegrals();
    /**
     * \brief Returns the abbreviation of the routine (DCD).
     * \return abbreviation of the routine
     */
    virtual std::string getAbbreviation() { return "Dcd"; }

  protected:
    /**
     * \brief Implements the iterate method with the DCD iteration. Iteration
     * routine taken from So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001).
     * \param[in] i Iteration number
     */
    virtual void iterate(int i);
    /**
     * \brief Implements the iterate method with the DCD iteration. Iteration
     * routine taken from Rev. Mod. Phys. 79, 291  Page 305, Figure 8.
     * \param[in] i Iteration number
     */
    virtual void iterateBartlett(int i);
    /**
     * \brief Implements the dry iterate method with the DCD iteration.
     */
    virtual void dryIterate();
  };
}

#endif

