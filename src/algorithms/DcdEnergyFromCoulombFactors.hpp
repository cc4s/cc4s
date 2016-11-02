/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCD_ENERGY_FROM_COULOMB_FACTORS_DEFINED 
#define DCD_ENERGY_FROM_COULOMB_FACTORS_DEFINED

#include <algorithms/ClusterDoublesAlgorithm.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the Dcd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
   * V_{kl}^{ij}\f$, and the factor orbitals \f$\Pi_{qR}\f$
   * and the Coulom factors \f$\Lambda_{GR}\f$.
   */
  class DcdEnergyFromCoulombFactors: public ClusterDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DcdEnergyFromCoulombFactors);
    DcdEnergyFromCoulombFactors(
      std::vector<Argument> const &argumentList
    );
    virtual ~DcdEnergyFromCoulombFactors();
    /**
     * \brief Returns the abbreviation of the routine (DCD).
     * \return abbreviation of the routine
     */
    virtual std::string getAbbreviation() { return "Dcd"; }

  protected:
    /**
     * \brief Implements the iterate method with the DCD iteration from So
     * Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001), modified
     * according to D. Kats, J. Chem. Phys. 139, 021102 (2013)
     */
    virtual void iterate(int i);
    /**
     * \brief Implements the dry iterate method with the DCD iteration
     */
    virtual void dryIterate();
  };
}

#endif

