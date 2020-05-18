/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_ENERGY_FROM_COULOMB_INTEGRALS_REFERENCE_DEFINED
#define CCSD_ENERGY_FROM_COULOMB_INTEGRALS_REFERENCE_DEFINED


#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterSinglesDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the Drccd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class CcsdEnergyFromCoulombIntegralsReference: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEnergyFromCoulombIntegralsReference);
    /**
     * \brief Returns the abbreviation of the routine (CCSD).
     * \return abbreviation of the routine
     */
    std::string getAbbreviation() override { return "Ccsd"; }

  protected:
    /**
     * \brief Implements the iterate method with the CCSD iteration.
     * \param[in] i Iteration number
     */
    Ptr<FockVector<Real<>, DryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Real<>, DryTensorEngine>> &amplitudes
    ) override;
    Ptr<FockVector<Real<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
    ) override;

    template <typename F, typename TE>
    Ptr<FockVector<F,TE>> getResiduum(
      const int iteration, const Ptr<const FockVector<F,TE>> &amplitudes
    );
  };
}

#endif

