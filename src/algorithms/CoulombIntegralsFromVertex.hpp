/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_FROM_VERTEX_DEFINED
#define COULOMB_INTEGRALS_FROM_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <array>

namespace cc4s {
  /**
   * \brief Caclulates the Coulomb Integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
   * V_{kl}^{ij}, V_{cd}^{ab}, V_{ka}^{ij}, V_{ci}^{ab}\f$ (if given) from the
   * Coulomb Vertex \f$\Gamma_{rG}^q\f$ and stores them in CTF Tensors Vabij,
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
    /**
     * \brief Dry run for calculating Coulomb integrals
     * Vabcd,Vabij,Vaibj,Vabci,Vijka Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices. Arguments can be
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    virtual void dryRun();
  protected:
    void calculateRealIntegrals();
    void calculateComplexIntegrals();

    CTF::Tensor<complex> *GammaGai;
    CTF::Tensor<complex> *GammaGia;
    CTF::Tensor<complex> *GammaGab;
    CTF::Tensor<complex> *GammaGij;
    std::array<int,4> syms, vvvv, vovo, vvoo, voov, oovv, oooo, ooov, vvvo;
  };
}

#endif

