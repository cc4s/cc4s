#ifndef CCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED
#define CCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED


#include <algorithms/CcsdEnergyFromCoulombIntegralsReference.hpp>

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
  class CcsdEnergyFromCoulombIntegrals:
    public CcsdEnergyFromCoulombIntegralsReference {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEnergyFromCoulombIntegrals)
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
    Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const TensorUnion<Real<>, DefaultDryTensorEngine>> &amplitudes
    ) override;
    Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const TensorUnion<Complex<>, DefaultDryTensorEngine>> &amplitudes
    ) override;
    Ptr<TensorUnion<Real<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const TensorUnion<Real<>, DefaultTensorEngine>> &amplitudes
    ) override;
    Ptr<TensorUnion<Complex<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const TensorUnion<Complex<>, DefaultTensorEngine>> &amplitudes
    ) override;

    template <typename TE>
    Ptr<TensorUnion<Real<>,TE>> getResiduum(
      const int iteration, const Ptr<const TensorUnion<Real<>,TE>> &amplitudes
    );
    template <typename TE>
    Ptr<TensorUnion<Complex<>,TE>> getResiduum(
      const int iteration, const Ptr<const TensorUnion<Complex<>,TE>> &amplitudes
    );
  };
}

#endif

