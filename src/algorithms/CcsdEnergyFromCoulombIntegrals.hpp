/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED
#define CCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterSinglesDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the Ccsd method. Calculates the
   * amplitudes \f$T_{a}^{i}\f$ and \f$T_{ab}^{ij}\f$ from the Coulomb
   * integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
   * V_{kl}^{ij}, V_{ka}^{ij}, V_{ci}^{ab}\f$ and \f$V_{cd}^{ab}\f$ (if given, else slicing and the Coulomb
   * Vertex \f$\Gamma_{pG}^q\f$  is used).
   */
  class CcsdEnergyFromCoulombIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEnergyFromCoulombIntegrals);
    CcsdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdEnergyFromCoulombIntegrals();

    /**
     * \brief Returns the abbreviation of the routine (CCSD).
     * \return abbreviation of the routine
     */
    virtual std::string getAbbreviation() { return "Ccsd"; }

    static int64_t constexpr DEFAULT_SLICE_SIZE = -1;

  protected:
    /**
     * \brief Implements the iterate method with the CCSD iteration. Iteration
     * routine taken from So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001).
     * \param[in] i Iteration number
     */
    virtual void iterate(
      int i, Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
    );
    virtual void iterate(
      int i, Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
    );

    template <typename F>
    void iterate(
      int i, Mixer<F> *TaiMixer, Mixer<F> *TabijMixer
    );
  };
}

#endif

