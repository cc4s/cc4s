/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all
 * rights reserved.*/
#ifndef MP2_EOM_DAVIDSON_DEFINED
#define MP2_EOM_DAVIDSON_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>

namespace cc4s {

  template <typename F = complex>
  class Mp2SimilarityTransformedHamiltonian {
  public:
    Mp2SimilarityTransformedHamiltonian(
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
      CTF::Tensor<F> *Vabic_
    );

    FockVector<F> rightApply(FockVector<F> &v);

  protected:
    CTF::Tensor<F> *Tai, *Tabij;
    CTF::Tensor<F> *Fij, *Fab;
    CTF::Tensor<F> *Vabcd, *Viajb, *Vijab, *Vijkl;
    CTF::Tensor<F> *Vijka, *Viabc, *Viajk, *Vabic;
  };

  /**
   * \brief Implements the diagonal preconditionar for the davidson method
   * \tparam F It is the field variable to be used, in general it will be
   * complex
   */
  template <typename F = complex>
  class Mp2PreConditioner {
  public:
    typedef FockVector<F> V;

    /**
     * \brief Constructor for the preconditioner.
     */
    Mp2PreConditioner (
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

  protected:
    V diagonalH;
  };

  class Mp2EquationOfMotionDavidson: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2EquationOfMotionDavidson);
    Mp2EquationOfMotionDavidson(
      std::vector<Argument> const &argumentList
    );
    virtual ~Mp2EquationOfMotionDavidson();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

    template <typename F = double>
    void getCanonicalPerturbationBasis(
        CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, int64_t i
    );

  };
}

#endif

