/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all
 * rights reserved.*/
#ifndef MP2_EOM_DAVIDSON_DEFINED
#define MP2_EOM_DAVIDSON_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>

namespace cc4s {

  template <typename F = complex>
  class CcsdSimilarityTransformedHamiltonian {
  public:
    CcsdSimilarityTransformedHamiltonian(
      CTF::Tensor<F> *Tai_,
      CTF::Tensor<F> *Tabij_,
      CTF::Tensor<F> *Fij_,
      CTF::Tensor<F> *Fab_,
      CTF::Tensor<F> *Vabcd_,
      CTF::Tensor<F> *Viajb_,
      CTF::Tensor<F> *Vijab_,
      CTF::Tensor<F> *Vijkl_,
      CTF::Tensor<F> *Vijka_,
      CTF::Tensor<F> *Viabc_,
      CTF::Tensor<F> *Viajk_,
      CTF::Tensor<F> *Vabic_,
      CTF::Tensor<F> *Vaibc_,
      CTF::Tensor<F> *Vaibj_,
      CTF::Tensor<F> *Viabj_,
      CTF::Tensor<F> *Vijak_,
      CTF::Tensor<F> *Vaijb_,
      CTF::Tensor<F> *Vabci_
    );

    FockVector<F> rightApplyIntermediates(FockVector<F> &v);
    FockVector<F> rightApplyHirata(FockVector<F> &v);
    FockVector<F> rightApply(FockVector<F> &v);
    FockVector<F> leftApply(FockVector<F> &v);

    /**
     * \brief This method should initialize the intermediates.
     *
     * \param[in] flag If true, then the rightApply method will be used
     * with intermediates, else without.
     */
    void buildIntermediates(bool flag);

  protected:
    CTF::Tensor<F> *Tai, *Tabij;
    CTF::Tensor<F> *Fij, *Fab;
    CTF::Tensor<F> *Vabcd;
    CTF::Tensor<F> *Viajb;
    CTF::Tensor<F> *Vijab;
    CTF::Tensor<F> *Vijkl;
    CTF::Tensor<F> *Vijka;
    CTF::Tensor<F> *Viabc;
    CTF::Tensor<F> *Viajk;
    CTF::Tensor<F> *Vabic;
    CTF::Tensor<F> *Vaibc;
    CTF::Tensor<F> *Vaibj;
    CTF::Tensor<F> *Viabj;
    CTF::Tensor<F> *Vijak;
    CTF::Tensor<F> *Vaijb;
    CTF::Tensor<F> *Vabci;
    bool withIntermediates;
    PTR(CTF::Tensor<F>) Wab, Wia, Wabcd, Wabci, Waibc,
                        Wiabj, Wiajk, Wij, Wijka, Wijkl;
  };

  /**
   * \brief Implements the diagonal preconditionar for the davidson method
   * \tparam F It is the field variable to be used, in general it will be
   * complex
   */
  template <typename F = complex>
  class CcsdPreConditioner {
  public:
    typedef FockVector<F> V;

    /**
     * \brief Constructor for the preconditioner.
     */
    CcsdPreConditioner (
      CTF::Tensor<F> &Tai,
      CTF::Tensor<F> &Tabij,
      CTF::Tensor<F> &Fij,
      CTF::Tensor<F> &Fab,
      CTF::Tensor<F> &Vabcd,
      CTF::Tensor<F> &Viajb,
      CTF::Tensor<F> &Vijab,
      CTF::Tensor<F> &Vijkl
    );

    std::vector<V> getInitialBasis(int eigenVectorsCount);

    V getCorrection(const complex eigenValue, V &residuum);

    V getDiagonalH() const {
      return diagonalH;
    }

  protected:
    V diagonalH;
  };

  class CcsdEquationOfMotionDavidsonIntermediates: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEquationOfMotionDavidsonIntermediates);
    CcsdEquationOfMotionDavidsonIntermediates(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdEquationOfMotionDavidsonIntermediates();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

    template <typename F = double>
    void getCanonicalPerturbationBasis(
        CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, size_t i
    );

  };
}

#endif

