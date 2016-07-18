/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_ENERGY_FROM_SLICED_INTEGRALS_DEFINED
#define CCSD_ENERGY_FROM_SLICED_INTEGRALS_DEFINED

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
   * V_{kl}^{ij}, V_{ka}^{ij} and  the Coulomb
   * Vertex \f$\Gamma_{pG}^q\f$.
   */
  class CcsdEnergyFromSlicedIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEnergyFromSlicedIntegrals);
    CcsdEnergyFromSlicedIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdEnergyFromSlicedIntegrals();

    /**
     * \brief Returns the abbreviation of the routine (CCSD).
     * \return abbreviation of the routine
     */
    virtual std::string getAbbreviation() { return "Ccsd"; }

    /** \brief The occupied orbital energies  */
    CTF::Tensor<> *epsi;
    /** \brief The virtual orbital energies  */
    CTF::Tensor<> *epsa;
    /** \brief The Coulomb integrals Vabij  */
    CTF::Tensor<> *Vabij;
    /** \brief The singles amplitudes Tai  */
    CTF::Tensor<> *Tai;
    /** \brief The doubles amplitudes Tabij  */
    CTF::Tensor<> *Tabij;

  protected:
    /**
     * \brief Implements the iterate method with the CCSD iteration. Iteration
     * routine taken from So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001).
     * \param[in] i Iteration number
     */
    virtual void iterate(int i);
    /**
     * \brief Implements the dry iterate method with the CCSD iteration.
     */
    virtual void dryIterate();
  };
}

#endif

