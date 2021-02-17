/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRCCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DRCCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

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
  class DrccdEnergyFromCoulombIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DrccdEnergyFromCoulombIntegrals)
    /**
     * \brief Returns the abbreviation of the routine (DRCCD).
     * \return abbreviation of the routine
     */
    std::string getAbbreviation() override { return "Drccd"; }

  protected:
    /**
     * \brief Implements the iterate method with the DRCCD iteration.
     * \param[in] i Iteration number
     */
    Ptr<FockVector<Real<>, DefaultDryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Real<>, DefaultDryTensorEngine>> &amplitudes
    ) override;
    Ptr<FockVector<Complex<>, DefaultDryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Complex<>, DefaultDryTensorEngine>> &amplitudes
    ) override;
    Ptr<FockVector<Real<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
    ) override;
    Ptr<FockVector<Complex<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Complex<>, DefaultTensorEngine>> &amplitudes
    ) override;

    template <typename F, typename TE>
    Ptr<FockVector<F,TE>> getResiduum(
      const int iteration, const Ptr<const FockVector<F,TE>> &amplitudes
    );
  };
}

#endif

