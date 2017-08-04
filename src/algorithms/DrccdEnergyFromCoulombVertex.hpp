/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRCCD_ENERGY_FROM_COULOMB_VERTEX_DEFINED 
#define DRCCD_ENERGY_FROM_COULOMB_VERTEX_DEFINED

#include <algorithms/LegacyClusterDoublesAlgorithm.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the Drccd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb vertex \f$ \Gamma^{a}_{iG} \f$
   * in a \f$ \mathcal{O}(N^{5}) \f$ implementation.
   */
  class DrccdEnergyFromCoulombVertex: public LegacyClusterDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DrccdEnergyFromCoulombVertex);
    DrccdEnergyFromCoulombVertex(
      std::vector<Argument> const &argumentList
    );
    virtual ~DrccdEnergyFromCoulombVertex();
    /**
     * \brief Returns the abbreviation of the routine (DRCCD).
     * \return abbreviation of the routine
     */
    virtual std::string getAbbreviation() { return "Drccd"; }

  protected:
    /**
     * \brief Implements the iterate method with the DRCCD iteration.
     * in a \f$ \mathcal{O}(N^{5}) \f$ implementation, using only the Coulomb
     * vertex \f$ \Gamma^{a}_{iG} \f$.
     * \param[in] i Iteration number
     */
    virtual void iterate(int i);
    /**
     * \brief Implements the dry iterate method with the DRCCD iteration.
     */
    virtual void dryIterate();
  };
}

#endif

