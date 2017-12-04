/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all rights reserved.*/
#ifndef MP2_EOM_DEFINED
#define MP2_EOM_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the Mp2 method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class UccsdAmplitudesFromCoulombIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UccsdAmplitudesFromCoulombIntegrals);
    UccsdAmplitudesFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~UccsdAmplitudesFromCoulombIntegrals();

    // Mask, one and two body parts
    PTR(CTF::Tensor<>) *Mai, *Mabij;

    virtual void run();
    void createMask();
    virtual std::string getAbbreviation() { return "Uccsd"; }

  protected:
    /**
     * \brief Implements the iterate method with the DRCCD iteration.
     * \param[in] i Iteration number
     */
    virtual PTR(FockVector<double>) getResiduum(
      const int iteration, const PTR(const FockVector<double>) &amplitudes
    );

    virtual PTR(FockVector<complex>) getResiduum(
      const int iteration, const PTR(const FockVector<complex>) &amplitudes
    );

  };
}

#endif

