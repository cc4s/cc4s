/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FINITE_SIZE_DEFINED
#define FINITE_SIZE_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates the Coulomb Integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
   * V_{kl}^{ij}, V_{cd}^{ab}, V_{ka}^{ij}, V_{ci}^{ab}\f$ (if given) from the
   * Coulomb Vertex \f$\Gamma_{pG}^q\f$ and stores them in CTF Tensors Vabij,
   * Vaibj, Vijkl, Vabcd, Vijka, and Vabci respectively. The arguments of the
   * integrals are PPPP, PPHH, HHHH, PHPH, HHHP, and PPPHCoulombIntegrals.
   */
  class CoulombIntegralsFromVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombIntegralsFromVertex);
    CoulombIntegralsFromVertex(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombIntegralsFromVertex();
    /**
     * \brief Calculates Coulomb integrals Vabcd,Vabij,Vaibj,Vabci,Vijka,Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices. Arguments can be
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    virtual void run();

    /** \brief The occupied orbital energies  */
    CTF::Tensor<> *epsi;
    /** \brief The virtual orbital energies  */
    CTF::Tensor<> *epsa;
    /** \brief The Coulomb Vertex GammaGpq  */
    CTF::Tensor<complex> *GammaGpq;

    /**
     * \brief Dry run for calculating Coulomb integrals
     * Vabcd,Vabij,Vaibj,Vabci,Vijka Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices. Arguments can be
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    virtual void dryRun();

  };
}

#endif

